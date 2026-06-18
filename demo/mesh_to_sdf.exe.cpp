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

namespace {

    static ptrdiff_t parse_extent(const char *const str) {
        char *end = nullptr;
        const long long value = std::strtoll(str, &end, 10);
        return (end != str && *end == '\0') ? static_cast<ptrdiff_t>(value) : ptrdiff_t(0);
    }

    template <typename G>
    static void expand_box(G &xmin, G &xmax, G &ymin, G &ymax, G &zmin, G &zmax, const G scale, const G margin) {
        const G cx = (xmin + xmax) * G(0.5);
        const G cy = (ymin + ymax) * G(0.5);
        const G cz = (zmin + zmax) * G(0.5);

        const G hx = (xmax - xmin) * G(0.5) * scale + margin;
        const G hy = (ymax - ymin) * G(0.5) * scale + margin;
        const G hz = (zmax - zmin) * G(0.5) * scale + margin;

        xmin = cx - hx;
        xmax = cx + hx;
        ymin = cy - hy;
        ymax = cy + hy;
        zmin = cz - hz;
        zmax = cz + hz;
    }

}  // namespace

int main(int argc, char **argv) {
    auto ctx = smesh::initialize(argc, argv);
    SMESH_TRACE_SCOPE("mesh_to_sdf.exe");

    if (argc != 6) {
        auto comm = ctx->communicator();
        if (!comm->rank()) {
            std::cerr << "Usage: " << argv[0] << " <mesh_folder> <nx> <ny> <nz> <output_folder>" << std::endl;
        }
        return 1;
    }

    using G = smesh::geom_t;
    using T = smesh::geom_t;
    using I = smesh::idx_t;

    const ptrdiff_t nx = parse_extent(argv[2]);
    const ptrdiff_t ny = parse_extent(argv[3]);
    const ptrdiff_t nz = parse_extent(argv[4]);

    if (nx < 2 || ny < 2 || nz < 2) {
        std::cerr << "Error: nx, ny, nz must be integers >= 2" << std::endl;
        return 1;
    }

    auto mesh = smesh::Mesh::create_from_file(smesh::Communicator::self(), smesh::Path(argv[1]));
    auto elements = mesh->elements(0)->data();
    auto points = mesh->points()->data();

    G xmin, xmax, ymin, ymax, zmin, zmax;
    ssdf::compute_aabb(mesh->n_nodes(), points[0], points[1], points[2], &xmin, &xmax, &ymin, &ymax, &zmin, &zmax);

    const G margin = smesh::Env::read("SSDF_MARGIN", G(0));
    const G scale = smesh::Env::read("SSDF_SCALE_FACTOR", G(1));
    expand_box(xmin, xmax, ymin, ymax, zmin, zmax, scale, margin);

    auto sdf = smesh::Grid<T>::create(smesh::Communicator::self(), nx, ny, nz, xmin, ymin, zmin, xmax, ymax, zmax);
    T *const out = sdf->data();
    const ptrdiff_t ngrid = nx * ny * nz;

#pragma omp parallel for
    for (ptrdiff_t i = 0; i < ngrid; ++i) {
        out[i] = static_cast<T>(std::numeric_limits<float>::max());
    }

    const int rc_dist = ssdf::sample_sqedf_bvh<G, T, I>(mesh->n_elements(0),
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
    if (rc_dist) {
        std::cerr << "Error: sample_sqedf_bvh failed with code " << rc_dist << std::endl;
        return rc_dist;
    }

    const bool closed = smesh::surface_is_closed(mesh);
    if (closed) {
        auto inside = smesh::create_host_buffer<uint8_t>(ngrid);
        uint8_t *const b_inside = inside->data();

        const int rc_sign = ssdf::sample_inside_bvh<G, I>(mesh->n_elements(0),
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
                                                          b_inside);
        if (rc_sign) {
            std::cerr << "Error: sample_inside_bvh failed with code " << rc_sign << std::endl;
            return rc_sign;
        }

#pragma omp parallel for
        for (ptrdiff_t i = 0; i < ngrid; ++i) {
            out[i] = std::sqrt(out[i]) * static_cast<T>(1 - 2 * int(b_inside[i]));
        }
    } else {
        auto qx = smesh::create_host_buffer<G>(ngrid);
        auto qy = smesh::create_host_buffer<G>(ngrid);
        auto qz = smesh::create_host_buffer<G>(ngrid);
        auto closest_dist = smesh::create_host_buffer<T>(ngrid);
        auto closest_x = smesh::create_host_buffer<G>(ngrid);
        auto closest_y = smesh::create_host_buffer<G>(ngrid);
        auto closest_z = smesh::create_host_buffer<G>(ngrid);
        auto closest_tri = smesh::create_host_buffer<I>(ngrid);

        G *const x = qx->data();
        G *const y = qy->data();
        G *const z = qz->data();
        T *const d = closest_dist->data();
        I *const tri = closest_tri->data();
        const G hx = (xmax - xmin) / G(nx - 1);
        const G hy = (ymax - ymin) / G(ny - 1);
        const G hz = (zmax - zmin) / G(nz - 1);

#pragma omp parallel for collapse(3)
        for (ptrdiff_t k = 0; k < nz; ++k) {
            for (ptrdiff_t j = 0; j < ny; ++j) {
                for (ptrdiff_t i = 0; i < nx; ++i) {
                    const ptrdiff_t idx = i + j * nx + k * nx * ny;
                    x[idx] = xmin + G(i) * hx;
                    y[idx] = ymin + G(j) * hy;
                    z[idx] = zmin + G(k) * hz;
                    d[idx] = static_cast<T>(std::numeric_limits<float>::max());
                    tri[idx] = I(-1);
                }
            }
        }

        const int rc_closest = ssdf::closest_point_bvh<G, T, I>(ngrid,
                                                                x,
                                                                y,
                                                                z,
                                                                mesh->n_elements(0),
                                                                elements[0],
                                                                elements[1],
                                                                elements[2],
                                                                mesh->n_nodes(),
                                                                points[0],
                                                                points[1],
                                                                points[2],
                                                                d,
                                                                closest_x->data(),
                                                                closest_y->data(),
                                                                closest_z->data(),
                                                                tri);
        if (rc_closest) {
            std::cerr << "Error: closest_point_bvh failed with code " << rc_closest << std::endl;
            return rc_closest;
        }

        G *const cx = closest_x->data();
        G *const cy = closest_y->data();
        G *const cz = closest_z->data();
        const I *const e0 = elements[0];
        const I *const e1 = elements[1];
        const I *const e2 = elements[2];
        const G *const px = points[0];
        const G *const py = points[1];
        const G *const pz = points[2];

#pragma omp parallel for
        for (ptrdiff_t idx = 0; idx < ngrid; ++idx) {
            const I t = tri[idx];
            const I i0 = e0[t];
            const I i1 = e1[t];
            const I i2 = e2[t];

            const G ux = px[i1] - px[i0];
            const G uy = py[i1] - py[i0];
            const G uz = pz[i1] - pz[i0];
            const G vx = px[i2] - px[i0];
            const G vy = py[i2] - py[i0];
            const G vz = pz[i2] - pz[i0];

            const G nxn = uy * vz - uz * vy;
            const G nyn = uz * vx - ux * vz;
            const G nzn = ux * vy - uy * vx;
            const G dot = (x[idx] - cx[idx]) * nxn + (y[idx] - cy[idx]) * nyn + (z[idx] - cz[idx]) * nzn;

            const T s = dot < G(0) ? T(-1) : T(1);
            out[idx] = std::sqrt(out[idx]) * s;
        }
    }

    const int rc_write = sdf->to_file(smesh::Path(argv[5]));
    if (rc_write) {
        std::cerr << "Error: failed to write grid to " << argv[5] << std::endl;
        return rc_write;
    }

    return 0;
}
