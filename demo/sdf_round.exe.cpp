#include "smesh_context.hpp"
#include "smesh_env.hpp"
#include "smesh_grid.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"

#include "bvh/sdf_rebuild.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char **argv) {
    auto ctx = smesh::initialize(argc, argv);
    {
        SMESH_TRACE_SCOPE("sdf_round.exe");

        if (argc != 4) {
            std::cerr << "Usage: " << argv[0] << " <grid_in> <radius> <grid_out>\n"
                      << "  Rounds the edges of the region where the grid's field is positive by <radius> (the\n"
                      << "  opening rounds its convex edges, the closing its concave ones) and writes the signed\n"
                      << "  distance to the result into <grid_out>, a new folder (round_sdf, sdf_rebuild.hpp).\n"
                      << "  Parts thinner than 2 radius are removed, gaps narrower than 2 radius filled: see the report.\n"
                      << "  SSDF_ROUND_OPEN, SSDF_ROUND_CLOSE  0 skips the opening or the closing (default 1)\n"
                      << "  SSDF_ROUND_BAND    exact within this of the surface, +-band past it (default 4 radius\n"
                      << "                     + 32 cells)\n"
                      << "  SSDF_ROUND_STRIDE  write every stride-th sample (default 1): the input on a finer grid\n"
                      << "                     than the output, which holds exact distances to the finer surface\n"
                      << "  SSDF_ROUND_REFINE  first refine the input by this factor (cubic, refine_field) and write\n"
                      << "                     on its own samples (default 1): the rounded surface's triangles, and\n"
                      << "                     the curvature their facets leave in the distance, shrink with it"
                      << std::endl;
            return 1;
        }

        using G = smesh::geom_t;
        using T = float;

        const std::string out_path = argv[3];
        if (std::filesystem::exists(out_path)) {
            std::cerr << "sdf_round: " << out_path << " exists; give a new folder" << std::endl;
            return 1;
        }

        auto in = smesh::Grid<T>::create_from_file(smesh::Communicator::self(), argv[1]);
        if (!in || in->block_size() != 1) {
            std::cerr << "sdf_round: unable to read a scalar grid from " << argv[1] << std::endl;
            return 1;
        }

        ssdf::GridGeometry<G> input;
        for (int d = 0; d < 3; ++d) {
            input.n[d] = in->extent(d);
            input.o[d] = in->origin(d);
            input.h[d] = in->delta(d);
        }
        const ptrdiff_t       refine = ptrdiff_t(smesh::Env::read("SSDF_ROUND_REFINE", int64_t(1)));
        const ssdf::GridGeometry<G> grid = refine > 1 ? input.refined(refine) : input;
        G                     hmax = 0;
        for (int d = 0; d < 3; ++d) hmax = std::max(hmax, grid.h[d]);
        const T         radius = T(std::atof(argv[2]));
        const bool      open   = smesh::Env::read("SSDF_ROUND_OPEN", 1) != 0;
        const bool      close  = smesh::Env::read("SSDF_ROUND_CLOSE", 1) != 0;
        const T         band   = smesh::Env::read("SSDF_ROUND_BAND", T(4 * radius + 32 * hmax));
        const ptrdiff_t stride =
                refine > 1 ? refine : ptrdiff_t(smesh::Env::read("SSDF_ROUND_STRIDE", int64_t(1)));
        if (!(radius > 0) || stride < 1 || refine < 1 || ((open || close) && !(band > 3 * radius))) {
            std::cerr << "sdf_round: needs radius > 0, stride >= 1, refine >= 1 and band > 3 radius" << std::endl;
            return 1;
        }

        const auto og = grid.strided(stride);
        std::cout << "sdf_round: " << argv[1] << ", " << grid.n[0] << " x " << grid.n[1] << " x " << grid.n[2]
                  << " samples at " << grid.h[0] << ", radius " << radius << (open ? ", opening" : "")
                  << (close ? ", closing" : "") << ", band " << band
                  << (refine > 1 ? ", refined by " + std::to_string(refine) : std::string()) << ", output " << og.n[0] << " x " << og.n[1]
                  << " x " << og.n[2] << " at " << og.h[0] << std::endl;

        std::vector<T> refined(refine > 1 ? grid.size() : 0);
        if (refine > 1) ssdf::refine_field<G, T>(input, in->data(), refine, refined.data());
        const T       *f = refine > 1 ? refined.data() : in->data();
        std::vector<T> f0(grid.size()), out(og.size());
        const auto     start = std::chrono::steady_clock::now();
        const int      rc    = ssdf::round_sdf<G, T>(grid, f, radius, band, open, close, stride, out.data(), f0.data());
        if (rc) {
            std::cerr << "sdf_round: round_sdf failed with code " << rc << std::endl;
            return rc;
        }
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

        // What changed, on the output's samples: the input against its distance f0 (a CSG of
        // distances is none inside), and the rounded surface against f0's within the band.
        const ptrdiff_t sy = grid.n[0], sz = grid.n[0] * grid.n[1];
        double          input_off = 0, moved = 0, moved_at[3] = {0, 0, 0};
        ptrdiff_t       removed = 0, filled = 0;
        for (ptrdiff_t k = 0; k < og.n[2]; ++k) {
            for (ptrdiff_t j = 0; j < og.n[1]; ++j) {
                for (ptrdiff_t i = 0; i < og.n[0]; ++i) {
                    const ptrdiff_t q = stride * (i + j * sy + k * sz), p = i + og.n[0] * (j + og.n[1] * k);
                    if (std::fabs(f0[q]) < band - radius) input_off = std::max(input_off, double(std::fabs(f[q] - f0[q])));
                    if (std::fabs(f0[q]) < band - radius && std::fabs(out[p] - f0[q]) > moved) {
                        moved       = std::fabs(out[p] - f0[q]);
                        moved_at[0] = double(og.o[0]) + double(i) * double(og.h[0]);
                        moved_at[1] = double(og.o[1]) + double(j) * double(og.h[1]);
                        moved_at[2] = double(og.o[2]) + double(k) * double(og.h[2]);
                    }
                    removed += f0[q] > 0 && out[p] < 0;
                    filled += f0[q] < 0 && out[p] > 0;
                }
            }
        }
        const double cell = double(og.h[0]) * double(og.h[1]) * double(og.h[2]);
        std::cout << "sdf_round: " << seconds << " s; the input differs from its distance by up to " << input_off
                  << "; the surface moved by up to " << moved << " (at " << moved_at[0] << ", " << moved_at[1]
                  << ", " << moved_at[2] << "; " << moved / radius << " radius: an edge is cut back by "
                  << "radius (1 / sin(angle / 2) - 1), 0.41 radius at 90 degrees, more where a part is "
                  << "thinner or a gap narrower than 2 radius); samples turned negative " << removed << " ("
                  << removed * cell << " volume), positive " << filled << " (" << filled * cell << ")" << std::endl;

        auto result = smesh::Grid<T>::create(smesh::Communicator::self(),
                                             og.n[0],
                                             og.n[1],
                                             og.n[2],
                                             og.o[0],
                                             og.o[1],
                                             og.o[2],
                                             G(og.o[0] + (og.n[0] - 1) * og.h[0]),
                                             G(og.o[1] + (og.n[1] - 1) * og.h[1]),
                                             G(og.o[2] + (og.n[2] - 1) * og.h[2]));
        std::copy(out.begin(), out.end(), result->data());
        const int rc_write = result->to_file(smesh::Path(out_path));
        if (rc_write) {
            std::cerr << "sdf_round: failed to write " << out_path << std::endl;
            return rc_write;
        }
    }
    return 0;
}
