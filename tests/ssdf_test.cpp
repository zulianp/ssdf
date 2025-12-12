#include "../ssdf.hpp"
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

    T out[5] = {1111111};

    // Call sdf function
    int result = ssdf::edf(npoints, x, y, z, nselements, s0, s1, s2, nspoints, sx, sy, sz, out);

    assert(result == 0);

    // Check that distances are computed (should be finite and non-negative)
    for (ptrdiff_t i = 0; i < npoints; i++) {
        assert(std::isfinite(out[i]));
        assert(out[i] >= 0.0f);
        std::cout << "Point (" << x[i] << ", " << y[i] << ", " << z[i] 
                  << ") distance: " << out[i] << std::endl;
    }

    // Point at (0.5, 0.5, 0.0) should be close to triangle center
    // Point at (0.5, 0.5, 1.0) should be approximately 1.0 away (vertical)
    assert(out[3] > 0.9f && out[3] < 1.1f);

    std::cout << "All tests passed!" << std::endl;
    return 0;
}

