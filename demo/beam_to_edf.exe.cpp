#include "smesh_buffer.hpp"
#include "smesh_context.hpp"
#include "smesh_env.hpp"
#include "smesh_grid.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"

#include "ssdf.hpp"

#include <algorithm>
#include <cmath>
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
static void expand_box(G &xmin,
                       G &xmax,
                       G &ymin,
                       G &ymax,
                       G &zmin,
                       G &zmax,
                       const G margin_x,
                       const G margin_y,
                       const G margin_z) {
    xmin -= margin_x;
    xmax += margin_x;
    ymin -= margin_y;
    ymax += margin_y;
    zmin -= margin_z;
    zmax += margin_z;
}

}  // namespace

int main(int argc, char **argv) {
    auto ctx = smesh::initialize(argc, argv);
    {
        SMESH_TRACE_SCOPE("beam_to_edf.exe");

        if (argc != 6) {
            auto comm = ctx->communicator();
            if (!comm->rank()) {
                std::cerr << "Usage: " << argv[0] << " <beam_mesh_folder> <nx> <ny> <nz> <output_folder>" << std::endl;
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
        auto element_buffer = mesh->elements(0);
        if (mesh->spatial_dimension() != 3 || element_buffer->extent(0) != 2) {
            std::cerr << "Error: expected a 3D BEAM2 mesh" << std::endl;
            return 1;
        }

        auto elements = element_buffer->data();
        auto points = mesh->points()->data();
        const ptrdiff_t nbeams = mesh->n_elements(0);

        G xmin, xmax, ymin, ymax, zmin, zmax;
        ssdf::compute_aabb(mesh->n_nodes(), points[0], points[1], points[2], &xmin, &xmax, &ymin, &ymax, &zmin, &zmax);

        const G margin = smesh::Env::read("SSDF_MARGIN", G(0));
        expand_box(xmin,
                   xmax,
                   ymin,
                   ymax,
                   zmin,
                   zmax,
                   smesh::Env::read("SSDF_MARGIN_X", margin),
                   smesh::Env::read("SSDF_MARGIN_Y", margin),
                   smesh::Env::read("SSDF_MARGIN_Z", margin));

        xmin = smesh::Env::read("SSDF_XMIN", xmin);
        xmax = smesh::Env::read("SSDF_XMAX", xmax);
        ymin = smesh::Env::read("SSDF_YMIN", ymin);
        ymax = smesh::Env::read("SSDF_YMAX", ymax);
        zmin = smesh::Env::read("SSDF_ZMIN", zmin);
        zmax = smesh::Env::read("SSDF_ZMAX", zmax);

        auto edf = smesh::Grid<T>::create(smesh::Communicator::self(), nx, ny, nz, xmin, ymin, zmin, xmax, ymax, zmax);
        T *const out = edf->data();

        auto sax = smesh::create_host_buffer<G>(nbeams);
        auto say = smesh::create_host_buffer<G>(nbeams);
        auto saz = smesh::create_host_buffer<G>(nbeams);
        auto sabx = smesh::create_host_buffer<G>(nbeams);
        auto saby = smesh::create_host_buffer<G>(nbeams);
        auto szabz = smesh::create_host_buffer<G>(nbeams);
        auto sinv_ab2 = smesh::create_host_buffer<G>(nbeams);

        G *const ax = sax->data();
        G *const ay = say->data();
        G *const az = saz->data();
        G *const abx = sabx->data();
        G *const aby = saby->data();
        G *const abz = szabz->data();
        G *const inv_ab2 = sinv_ab2->data();

        const I *const e0 = elements[0];
        const I *const e1 = elements[1];
        const G *const px = points[0];
        const G *const py = points[1];
        const G *const pz = points[2];

#pragma omp parallel for
        for (ptrdiff_t b = 0; b < nbeams; ++b) {
            const I i0 = e0[b];
            const I i1 = e1[b];
            const G x0 = px[i0];
            const G y0 = py[i0];
            const G z0 = pz[i0];
            const G dx = px[i1] - x0;
            const G dy = py[i1] - y0;
            const G dz = pz[i1] - z0;
            const G len2 = dx * dx + dy * dy + dz * dz;

            ax[b] = x0;
            ay[b] = y0;
            az[b] = z0;
            abx[b] = dx;
            aby[b] = dy;
            abz[b] = dz;
            inv_ab2[b] = len2 > G(0) ? G(1) / len2 : G(0);
        }

        const G hx = (xmax - xmin) / G(nx - 1);
        const G hy = (ymax - ymin) / G(ny - 1);
        const G hz = (zmax - zmin) / G(nz - 1);

#pragma omp parallel for collapse(3)
        for (ptrdiff_t k = 0; k < nz; ++k) {
            for (ptrdiff_t j = 0; j < ny; ++j) {
                for (ptrdiff_t i = 0; i < nx; ++i) {
                    const G x = xmin + G(i) * hx;
                    const G y = ymin + G(j) * hy;
                    const G z = zmin + G(k) * hz;
                    G best = std::numeric_limits<G>::max();

                    for (ptrdiff_t b = 0; b < nbeams; ++b) {
                        const G aqx = x - ax[b];
                        const G aqy = y - ay[b];
                        const G aqz = z - az[b];
                        const G t_unclamped = (aqx * abx[b] + aqy * aby[b] + aqz * abz[b]) * inv_ab2[b];
                        const G t0 = t_unclamped < G(0) ? G(0) : t_unclamped;
                        const G t = t0 > G(1) ? G(1) : t0;
                        const G dx = aqx - t * abx[b];
                        const G dy = aqy - t * aby[b];
                        const G dz = aqz - t * abz[b];
                        const G d2 = dx * dx + dy * dy + dz * dz;
                        best = d2 < best ? d2 : best;
                    }

                    out[i + j * nx + k * nx * ny] = std::sqrt(best);
                }
            }
        }

        const G scale = smesh::Env::read("SSDF_SCALE", G(1));
        edf->scale(scale);

        const int rc_write = edf->to_file(smesh::Path(argv[5]));
        if (rc_write) {
            std::cerr << "Error: failed to write grid to " << argv[5] << std::endl;
            return rc_write;
        }
    }

    return 0;
}
