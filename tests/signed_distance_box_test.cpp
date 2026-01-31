// Correctness tests for signed_distance_box
#include "../edf/edf.hpp"
#include <cassert>
#include <cmath>
#include <iostream>

int main() {
    using T = double;

    const T minx = -1.0, maxx = 1.0;
    const T miny = -2.0, maxy = 2.0;
    const T minz = 0.5, maxz = 1.5;

    auto check = [&](T px, T py, T pz, T expected, const char* label) {
        const T d = ssdf::signed_distance_box(px, py, pz, minx, maxx, miny, maxy, minz, maxz);
        const T diff = std::abs(d - expected);
        if (diff > 1e-12) {
            std::cerr << "Mismatch " << label << ": expected " << expected << " got " << d << std::endl;
            assert(false);
        }
    };

    // Inside: distance should be negative with magnitude equal to max penetration
    check(0.0, 0.0, 1.0, -0.5, "center depth to top face");
    check(0.0, 0.0, 1.4, -0.1, "near top face inside");
    check(0.9, 0.0, 1.0, -0.1, "near +X inside");
    check(-0.8, 1.9, 1.0, -0.1, "near corner inside");

    // On surface: should be exactly zero
    check(1.0, 0.0, 1.0, 0.0, "on +X face");
    check(0.0, -2.0, 1.0, 0.0, "on -Y face");
    check(0.0, 0.0, 1.5, 0.0, "on +Z face");

    // Outside: positive Euclidean distance to box
    check(2.0, 0.0, 1.0, 1.0, "outside +X");
    check(0.0, -3.0, 1.0, 1.0, "outside -Y");
    check(0.0, 0.0, 3.0, 1.5, "outside +Z");
    check(3.0, 4.0, -1.0, std::sqrt(2.0 * 2.0 + 2.0 * 2.0 + 1.5 * 1.5), "outside corner");

    std::cout << "signed_distance_box tests passed\n";
    return 0;
}
