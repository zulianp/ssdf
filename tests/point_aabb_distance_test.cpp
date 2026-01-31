// Basic correctness tests for point_aabb_distance_squared
#include "../edf/edf.hpp"
#include <cassert>
#include <cmath>
#include <iostream>

int main() {
    using T = double;

    const T minx = 0.0, maxx = 1.0;
    const T miny = 0.0, maxy = 1.0;
    const T minz = 0.0, maxz = 1.0;

    auto check = [&](T px, T py, T pz, T expected_sq) {
        const T d2 = ssdf::point_aabb_distance_squared(px, py, pz, minx, maxx, miny, maxy, minz, maxz);
        const T diff = std::abs(d2 - expected_sq);
        if (diff > 1e-12) {
            std::cerr << "Mismatch for (" << px << ", " << py << ", " << pz << ") expected "
                      << expected_sq << " got " << d2 << std::endl;
            assert(false);
        }
    };

    // Inside box → zero
    check(0.5, 0.5, 0.5, 0.0);
    check(0.0, 0.0, 0.0, 0.0);  // on a corner
    check(1.0, 0.2, 0.8, 0.0);  // on a face edge

    // Outside along one axis
    check(-1.0, 0.5, 0.5, 1.0);     // left of minx
    check(0.5, 2.0, 0.5, 1.0);      // above maxy
    check(0.5, 0.5, -0.25, 0.0625); // below minz

    // Outside along multiple axes
    check(-1.0, -2.0, 3.0, 9.0);  // dx=1, dy=2, dz=2 → 1+4+4

    std::cout << "point_aabb_distance_squared tests passed\n";
    return 0;
}
