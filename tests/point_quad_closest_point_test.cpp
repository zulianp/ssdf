#include "edf/edf.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

int main() {
    using T = double;

    T u = T(-1);
    const T segment_dist = ssdf::point_segment_parameter(T(0.25), T(1), T(0), T(0), T(0), T(0), T(1), T(0), T(0), &u);
    assert(std::abs(u - T(0.25)) < T(1e-12));
    assert(std::abs(segment_dist - T(1)) < T(1e-12));

    T wa = T(-1), wb = T(-1), wc = T(-1);
    const T triangle_dist = ssdf::point_triangle_barycentric(
        T(0.25), T(0.25), T(1), T(0), T(0), T(0), T(1), T(0), T(0), T(0), T(1), T(0), &wa, &wb, &wc);
    assert(std::abs(wa - T(0.5)) < T(1e-12));
    assert(std::abs(wb - T(0.25)) < T(1e-12));
    assert(std::abs(wc - T(0.25)) < T(1e-12));
    assert(std::abs(triangle_dist - T(1)) < T(1e-12));

    auto check = [](const char *label,
                    const T px,
                    const T py,
                    const T pz,
                    const T x0,
                    const T y0,
                    const T z0,
                    const T x1,
                    const T y1,
                    const T z1,
                    const T x2,
                    const T y2,
                    const T z2,
                    const T x3,
                    const T y3,
                    const T z3,
                    const T expected_s,
                    const T expected_t) {
        T s = T(-1);
        T t = T(-1);

        ssdf::point_to_quad_closest_point(px, py, pz, x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3, &s, &t);

        const T tol = T(1e-12);
        if (std::abs(s - expected_s) > tol || std::abs(t - expected_t) > tol) {
            std::cerr << label << ": expected (" << expected_s << ", " << expected_t << ") got (" << s << ", " << t
                      << ")\n";
            assert(false);
        }
    };

    check("basic",
          T(0.25),
          T(0.75),
          T(1),
          T(0),
          T(0),
          T(0),
          T(1),
          T(0),
          T(0),
          T(1),
          T(1),
          T(0),
          T(0),
          T(1),
          T(0),
          T(0.25),
          T(0.75));

    check("warped",
          T(0.7),
          T(0.2),
          T(0.225),
          T(0),
          T(0),
          T(0),
          T(1),
          T(0),
          T(0.25),
          T(1),
          T(1),
          T(0.5),
          T(0),
          T(1),
          T(-0.25),
          T(0.7),
          T(0.2));

    check("concave",
          T(0.3),
          T(1.0),
          T(0),
          T(0),
          T(0),
          T(0),
          T(2),
          T(0),
          T(0),
          T(0.75),
          T(0.5),
          T(0),
          T(0),
          T(2),
          T(0),
          T(0.4),
          T(0.8));

    check("degenerate",
          T(1.5),
          T(1),
          T(0),
          T(0),
          T(0),
          T(0),
          T(1),
          T(0),
          T(0),
          T(2),
          T(0),
          T(0),
          T(3),
          T(0),
          T(0),
          T(1),
          T(0.5));

    return 0;
}
