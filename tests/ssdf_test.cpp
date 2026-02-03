#include "ssdf.hpp"
#include <cassert>
#include <cmath>
#include <iostream>

int main() {
    using G = float;
    using T = float;
    using I = int;

    // Test case: Single triangle surface
    const ptrdiff_t nspoints = 3;
    G sx[3] = {0.0f, 1.0f, 0.5f};
    G sy[3] = {0.0f, 0.0f, 1.0f};
    G sz[3] = {0.0f, 0.0f, 0.0f};

    // Triangle indices
    const ptrdiff_t nselements = 1;
    I s0[1] = {0};
    I s1[1] = {1};
    I s2[1] = {2};

    // Test points
    const ptrdiff_t npoints = 5;
    G x[5] = {0.5f, 0.0f, 1.0f, 0.5f, 0.25f};
    G y[5] = {0.5f, 0.0f, 0.0f, 0.5f, 0.25f};
    G z[5] = {0.0f, 0.0f, 0.0f, 1.0f, 0.0f};

    T out[5];
    for(int i = 0; i < 5; i++) {
        out[i] = 1111111;
    }

    // Call sdf function (use celllist variant for correctness check)
    int result = ssdf::edf_select(npoints, x, y, z, nselements, s0, s1, s2, nspoints, sx, sy, sz, out);

    assert(result == 0);

    auto brute = [&](ptrdiff_t p) {
        const G px = x[p], py = y[p], pz = z[p];
        const G dist_sq = ssdf::point_triangle_dist_sq(px,
                                                       py,
                                                       pz,
                                                       sx[s0[0]],
                                                       sy[s0[0]],
                                                       sz[s0[0]],
                                                       sx[s1[0]],
                                                       sy[s1[0]],
                                                       sz[s1[0]],
                                                       sx[s2[0]],
                                                       sy[s2[0]],
                                                       sz[s2[0]]);
        return std::sqrt(dist_sq);
    };

    // Check that distances are computed and match brute force
    for (ptrdiff_t i = 0; i < npoints; i++) {
        assert(std::isfinite(out[i]));
        assert(out[i] >= 0.0f);
        const G expected = brute(i);
        assert(std::abs(out[i] - expected) < 1e-5f);
        std::cout << "Point (" << x[i] << ", " << y[i] << ", " << z[i] 
                  << ") distance: " << out[i] << std::endl;
    }

    std::cout << "All tests passed!" << std::endl;
    return 0;
}

