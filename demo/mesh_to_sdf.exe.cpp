#include "smesh_buffer.hpp"
#include "smesh_context.hpp"
#include "smesh_env.hpp"
#include "smesh_grid.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"

#include "ssdf.hpp"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {

    static ptrdiff_t parse_extent(const char *const str) {
        char *end = nullptr;
        const long long value = std::strtoll(str, &end, 10);
        return (end != str && *end == '\0') ? static_cast<ptrdiff_t>(value) : ptrdiff_t(0);
    }

    // The box grown by the margins, then scaled about its centre.
    template <typename G>
    static void expand_box(G &xmin,
                           G &xmax,
                           G &ymin,
                           G &ymax,
                           G &zmin,
                           G &zmax,
                           const G scale,
                           const G margin_x,
                           const G margin_y,
                           const G margin_z) {
        const G cx = (xmin + xmax) * G(0.5);
        const G cy = (ymin + ymax) * G(0.5);
        const G cz = (zmin + zmax) * G(0.5);

        const G hx = ((xmax - xmin) * G(0.5) + margin_x) * scale;
        const G hy = ((ymax - ymin) * G(0.5) + margin_y) * scale;
        const G hz = ((zmax - zmin) * G(0.5) + margin_z) * scale;

        xmin = cx - hx;
        xmax = cx + hx;
        ymin = cy - hy;
        ymax = cy + hy;
        zmin = cz - hz;
        zmax = cz + hz;
    }

    static std::vector<std::string> split_names(const std::string &list) {
        std::vector<std::string> names;
        std::stringstream        ss(list);
        std::string              name;
        while (std::getline(ss, name, ',')) {
            if (!name.empty()) names.push_back(name);
        }
        return names;
    }

    // Points per axis so that the spacing is at most `hmax`.
    static ptrdiff_t points_for_spacing(const double range, const double hmax) {
        return static_cast<ptrdiff_t>(std::ceil(range / hmax)) + 1;
    }

}  // namespace

int main(int argc, char **argv) {
    auto ctx = smesh::initialize(argc, argv);
    {
        SMESH_TRACE_SCOPE("mesh_to_sdf.exe");

        if (argc != 6) {
            auto comm = ctx->communicator();
            if (!comm->rank()) {
                std::cerr << "Usage: " << argv[0] << " <mesh_folder> <nx> <ny> <nz> <output_folder>\n"
                          << "  The signed distance to the triangles of every block of <mesh_folder> on a grid.\n"
                          << "  SSDF_MARGIN[_X|_Y|_Z]  grow the box by this much (default 0)\n"
                          << "  SSDF_SCALE_BOX         then scale it about its centre (default 1)\n"
                          << "  SSDF_BOX_MESH          take the box from this mesh, not the surface\n"
                          << "  SSDF_BOX_BLOCKS        only these blocks of it (comma separated)\n"
                          << "  SSDF_HMAX              grid spacing at most this; then pass 0 0 0\n"
                          << "  SSDF_SCALE             multiply the distances (e.g. -1 flips the sign)" << std::endl;
            }
            return 1;
        }

        using G = smesh::geom_t;
        using T = smesh::geom_t;
        using I = smesh::idx_t;

        ptrdiff_t nx = parse_extent(argv[2]);
        ptrdiff_t ny = parse_extent(argv[3]);
        ptrdiff_t nz = parse_extent(argv[4]);
        const double hmax = smesh::Env::read("SSDF_HMAX", 0.0);

        if (hmax <= 0 && (nx < 2 || ny < 2 || nz < 2)) {
            std::cerr << "Error: nx, ny, nz must be integers >= 2 (or set SSDF_HMAX)" << std::endl;
            return 1;
        }

        auto mesh = smesh::Mesh::create_from_file(smesh::Communicator::self(), smesh::Path(argv[1]));
        auto points = mesh->points()->data();

        // The triangles of every block, as one list.
        ptrdiff_t n_triangles = 0;
        for (const auto &block : mesh->blocks()) {
            if (block->element_type() != smesh::TRI3 && block->element_type() != smesh::TRISHELL3) {
                std::cerr << "Error: block " << block->name() << " of " << argv[1] << " is not made of triangles"
                          << std::endl;
                return 1;
            }
            n_triangles += block->n_elements();
        }
        auto triangles = smesh::create_host_buffer<I>(3, n_triangles);
        auto elements  = triangles->data();
        {
            ptrdiff_t offset = 0;
            for (const auto &block : mesh->blocks()) {
                for (int v = 0; v < 3; ++v) {
                    const I *const src = block->elements()->data()[v];
                    for (ptrdiff_t e = 0; e < block->n_elements(); ++e) {
                        elements[v][offset + e] = src[e];
                    }
                }
                offset += block->n_elements();
            }
        }

        // The box: of the surface, or of another mesh (only some of its blocks).
        G xmin, xmax, ymin, ymax, zmin, zmax;
        const char *const box_mesh_path = std::getenv("SSDF_BOX_MESH");
        if (box_mesh_path && box_mesh_path[0]) {
            const char *const box_blocks = std::getenv("SSDF_BOX_BLOCKS");
            auto              box_mesh   = smesh::Mesh::create_from_file(smesh::Communicator::self(),
                                                                smesh::Path(box_mesh_path),
                                                                split_names(box_blocks ? box_blocks : ""));
            if (!box_mesh) {
                std::cerr << "Error: unable to read " << box_mesh_path << std::endl;
                return 1;
            }
            auto bp = box_mesh->points()->data();
            ssdf::compute_aabb(box_mesh->n_nodes(), bp[0], bp[1], bp[2], &xmin, &xmax, &ymin, &ymax, &zmin, &zmax);
        } else {
            ssdf::compute_aabb(mesh->n_nodes(), points[0], points[1], points[2], &xmin, &xmax, &ymin, &ymax, &zmin, &zmax);
        }

        const G margin = smesh::Env::read("SSDF_MARGIN", G(0));
        const G margin_x = smesh::Env::read("SSDF_MARGIN_X", margin);
        const G margin_y = smesh::Env::read("SSDF_MARGIN_Y", margin);
        const G margin_z = smesh::Env::read("SSDF_MARGIN_Z", margin);
        const G scale_box = smesh::Env::read("SSDF_SCALE_BOX", G(1));

        expand_box(xmin, xmax, ymin, ymax, zmin, zmax, scale_box, margin_x, margin_y, margin_z);

        if (hmax > 0) {
            nx = points_for_spacing(xmax - xmin, hmax);
            ny = points_for_spacing(ymax - ymin, hmax);
            nz = points_for_spacing(zmax - zmin, hmax);
        }
        std::cout << "mesh_to_sdf: " << n_triangles << " triangles, grid " << nx << " x " << ny << " x " << nz
                  << " on [" << xmin << ", " << xmax << "] x [" << ymin << ", " << ymax << "] x [" << zmin << ", "
                  << zmax << "]" << std::endl;

        auto sdf = smesh::Grid<T>::create(smesh::Communicator::self(), nx, ny, nz, xmin, ymin, zmin, xmax, ymax, zmax);
        T *const out = sdf->data();

        // Signed by the pseudonormal at the closest point: exact for a closed, consistently oriented
        // surface; the ray-crossing vote (sample_inside_bvh) misclassifies points whose rays pass
        // through edges and vertices.
        const int rc = ssdf::sample_sdf_bvh<G, T, I>(n_triangles,
                                                     elements[0],
                                                     elements[1],
                                                     elements[2],
                                                     mesh->n_nodes(),
                                                     points[0],
                                                     points[1],
                                                     points[2],
                                                     nx,
                                                     ny,
                                                     nz,
                                                     xmin,
                                                     ymin,
                                                     zmin,
                                                     xmax,
                                                     ymax,
                                                     zmax,
                                                     out);
        if (rc) {
            std::cerr << "Error: sample_sdf_bvh failed with code " << rc << std::endl;
            return rc;
        }

        const G scale = smesh::Env::read("SSDF_SCALE", G(1));
        sdf->scale(scale);

        const int rc_write = sdf->to_file(smesh::Path(argv[5]));
        if (rc_write) {
            std::cerr << "Error: failed to write grid to " << argv[5] << std::endl;
            return rc_write;
        }
    }

    return 0;
}
