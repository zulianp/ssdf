#ifndef EDF_HPP
#define EDF_HPP

#include <stddef.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>

#ifndef SSDF_READ_ENV
#define SSDF_READ_ENV(name, conversion) \
    do {                                \
        char *var = getenv(#name);      \
        if (var) {                      \
            name = conversion(var);     \
        }                               \
    } while (0)
#endif

#ifndef SSDF_RESTRICT
#ifndef _WIN32
#define SSDF_RESTRICT __restrict__
#else
#define SSDF_RESTRICT __restrict
#endif
#endif

#define VECTOR_SIZE 32

#define MIN(a, b) ((a) < (b) ? (a) : (b))

namespace ssdf {

    template <typename T>
    inline static T point_triangle_dist_sq_approx(const T px,
                                                  const T py,
                                                  const T pz,
                                                  const T ax,
                                                  const T ay,
                                                  const T az,
                                                  const T bx,
                                                  const T by,
                                                  const T bz,
                                                  const T cx,
                                                  const T cy,
                                                  const T cz) {
        const T dx = std::min({std::abs(px - ax), std::abs(px - bx), std::abs(px - cx)});
        const T dy = std::min({std::abs(py - ay), std::abs(py - by), std::abs(py - cy)});
        const T dz = std::min({std::abs(pz - az), std::abs(pz - bz), std::abs(pz - cz)});
        return dx * dx + dy * dy + dz * dz;
    }

    template <typename T>
    inline static T point_aabb_distance_squared(const T px,
                                                const T py,
                                                const T pz,
                                                const T minx,
                                                const T maxx,
                                                const T miny,
                                                const T maxy,
                                                const T minz,
                                                const T maxz) {
        const T dx = (px < minx) ? (minx - px) : (px > maxx ? px - maxx : T(0));
        const T dy = (py < miny) ? (miny - py) : (py > maxy ? py - maxy : T(0));
        const T dz = (pz < minz) ? (minz - pz) : (pz > maxz ? pz - maxz : T(0));
        return dx * dx + dy * dy + dz * dz;
    }

    // 5.1.9 Closest Points of Two Line Segments. Real-Time Collision Detection by Christer Ericson.
    template <typename T>
    inline static void edge_to_edge_closest_points(
        //    First edge
        const T ap1x,
        const T ap1y,
        const T ap1z,
        const T ap2x,
        const T ap2y,
        const T ap2z,
        //    Second edge
        const T bp1x,
        const T bp1y,
        const T bp1z,
        const T bp2x,
        const T bp2y,
        const T bp2z,
        // Output line parameters
        T *const SSDF_RESTRICT s0,
        T *const SSDF_RESTRICT s1) {
        const T ux = ap2x - ap1x;
        const T uy = ap2y - ap1y;
        const T uz = ap2z - ap1z;

        const T vx = bp2x - bp1x;
        const T vy = bp2y - bp1y;
        const T vz = bp2z - bp1z;

        const T wx = ap1x - bp1x;
        const T wy = ap1y - bp1y;
        const T wz = ap1z - bp1z;

        const T a = ux * ux + uy * uy + uz * uz;
        const T b = ux * vx + uy * vy + uz * vz;
        const T c = vx * vx + vy * vy + vz * vz;
        const T d = ux * wx + uy * wy + uz * wz;
        const T e = vx * wx + vy * wy + vz * wz;

        // Treat very short edges as degenerate to avoid unstable divisions.
        const T eps = std::numeric_limits<T>::epsilon() * T(16);

        if (a <= eps) {
            // First edge collapses to a point, so only project that point onto edge B.
            *s0 = T(0);
            *s1 = (c <= eps) ? T(0) : std::clamp(e / c, T(0), T(1));
            return;
        }

        if (c <= eps) {
            // Second edge collapses to a point, so only project that point onto edge A.
            *s1 = T(0);
            *s0 = std::clamp(-d / a, T(0), T(1));
            return;
        }

        // Solve the unconstrained closest-point problem on the two supporting lines.
        const T D = a * c - b * b;

        T sN;
        T sD = D;
        T tN;
        T tD = D;

        if (D <= eps) {
            // Nearly parallel edges: pin the first parameter and project onto edge B.
            sN = T(0);
            sD = T(1);
            tN = e;
            tD = c;
        } else {
            sN = b * e - c * d;
            tN = a * e - b * d;

            if (sN <= T(0)) {
                // Closest point on edge A falls before ap1, clamp to the first endpoint.
                sN = T(0);
                tN = e;
                tD = c;
            } else if (sN >= sD) {
                // Closest point on edge A falls after ap2, clamp to the second endpoint.
                sN = sD;
                tN = e + b;
                tD = c;
            }
        }

        if (tN <= T(0)) {
            // Closest point on edge B falls before bp1, reproject the clamped endpoint onto edge A.
            tN = T(0);

            if (-d <= T(0)) {
                sN = T(0);
                sD = T(1);
            } else if (-d >= a) {
                sN = T(1);
                sD = T(1);
            } else {
                sN = -d;
                sD = a;
            }
        } else if (tN >= tD) {
            // Closest point on edge B falls after bp2, reproject the clamped endpoint onto edge A.
            tN = tD;

            const T bmd = b - d;
            if (bmd <= T(0)) {
                sN = T(0);
                sD = T(1);
            } else if (bmd >= a) {
                sN = T(1);
                sD = T(1);
            } else {
                sN = bmd;
                sD = a;
            }
        }

        // Convert the numerator/denominator form to segment-local coordinates in [0, 1].
        *s0 = (std::abs(sN) <= eps) ? T(0) : (sN / sD);
        *s1 = (std::abs(tN) <= eps) ? T(0) : (tN / tD);
    }

    template <typename T>
    inline static T point_aabb_signed_distance(const T px,
                                               const T py,
                                               const T pz,
                                               const T minx,
                                               const T maxx,
                                               const T miny,
                                               const T maxy,
                                               const T minz,
                                               const T maxz) {
        const T hx = (maxx - minx) * T(0.5);
        const T hy = (maxy - miny) * T(0.5);
        const T hz = (maxz - minz) * T(0.5);
        const T cx = (maxx + minx) * T(0.5);
        const T cy = (maxy + miny) * T(0.5);
        const T cz = (maxz + minz) * T(0.5);

        const T qx = std::abs(px - cx) - hx;
        const T qy = std::abs(py - cy) - hy;
        const T qz = std::abs(pz - cz) - hz;

        const T ox = qx > T(0) ? qx : T(0);
        const T oy = qy > T(0) ? qy : T(0);
        const T oz = qz > T(0) ? qz : T(0);

        const T outside = std::sqrt(ox * ox + oy * oy + oz * oz);
        const T inside = std::min<T>(T(0), std::max(qx, std::max(qy, qz)));
        return outside + inside;
    }

    template <typename T>
    inline static T point_segment_parameter(const T qx,
                                            const T qy,
                                            const T qz,
                                            const T ax,
                                            const T ay,
                                            const T az,
                                            const T bx,
                                            const T by,
                                            const T bz,
                                            T *const SSDF_RESTRICT out) {
        const T abx = bx - ax;
        const T aby = by - ay;
        const T abz = bz - az;
        const T aqx = qx - ax;
        const T aqy = qy - ay;
        const T aqz = qz - az;
        const T ab2 = abx * abx + aby * aby + abz * abz;
        const T u = (ab2 > T(0)) ? std::clamp((aqx * abx + aqy * aby + aqz * abz) / ab2, T(0), T(1)) : T(0);
        const T cx = ax + u * abx;
        const T cy = ay + u * aby;
        const T cz = az + u * abz;
        const T dx = qx - cx;
        const T dy = qy - cy;
        const T dz = qz - cz;
        *out = u;
        return dx * dx + dy * dy + dz * dz;
    }

    template <typename T>
    inline static T point_triangle_barycentric(const T qx,
                                               const T qy,
                                               const T qz,
                                               const T ax,
                                               const T ay,
                                               const T az,
                                               const T bx,
                                               const T by,
                                               const T bz,
                                               const T cx,
                                               const T cy,
                                               const T cz,
                                               T *const SSDF_RESTRICT wa,
                                               T *const SSDF_RESTRICT wb,
                                               T *const SSDF_RESTRICT wc) {
        const T abx = bx - ax, aby = by - ay, abz = bz - az;
        const T acx = cx - ax, acy = cy - ay, acz = cz - az;
        const T bcx = cx - bx, bcy = cy - by, bcz = cz - bz;
        const T ab2 = abx * abx + aby * aby + abz * abz;
        const T ac2 = acx * acx + acy * acy + acz * acz;
        const T bc2 = bcx * bcx + bcy * bcy + bcz * bcz;
        const T nx = aby * acz - abz * acy;
        const T ny = abz * acx - abx * acz;
        const T nz = abx * acy - aby * acx;
        const T max_edge2 = std::max(ab2, std::max(ac2, bc2));
        const T n2 = nx * nx + ny * ny + nz * nz;

        if (max_edge2 <= T(0) || n2 <= std::numeric_limits<T>::epsilon() * T(64) * max_edge2 * max_edge2) {
            T u;
            T best = point_segment_parameter(qx, qy, qz, ax, ay, az, bx, by, bz, &u);
            *wa = T(1) - u;
            *wb = u;
            *wc = T(0);

            T d = point_segment_parameter(qx, qy, qz, bx, by, bz, cx, cy, cz, &u);
            if (d < best) {
                best = d;
                *wa = T(0);
                *wb = T(1) - u;
                *wc = u;
            }

            d = point_segment_parameter(qx, qy, qz, cx, cy, cz, ax, ay, az, &u);
            if (d < best) {
                best = d;
                *wa = u;
                *wb = T(0);
                *wc = T(1) - u;
            }

            return best;
        }

        const T aqx = qx - ax, aqy = qy - ay, aqz = qz - az;
        const T d1 = abx * aqx + aby * aqy + abz * aqz;
        const T d2 = acx * aqx + acy * aqy + acz * aqz;
        if (d1 <= T(0) && d2 <= T(0)) {
            *wa = T(1);
            *wb = T(0);
            *wc = T(0);
            return aqx * aqx + aqy * aqy + aqz * aqz;
        }

        const T bqx = qx - bx, bqy = qy - by, bqz = qz - bz;
        const T d3 = abx * bqx + aby * bqy + abz * bqz;
        const T d4 = acx * bqx + acy * bqy + acz * bqz;
        if (d3 >= T(0) && d4 <= d3) {
            *wa = T(0);
            *wb = T(1);
            *wc = T(0);
            return bqx * bqx + bqy * bqy + bqz * bqz;
        }

        const T vc = d1 * d4 - d3 * d2;
        if (vc <= T(0) && d1 >= T(0) && d3 <= T(0)) {
            const T v = d1 / (d1 - d3);
            const T x = ax + v * abx;
            const T y = ay + v * aby;
            const T z = az + v * abz;
            const T dx = qx - x, dy = qy - y, dz = qz - z;
            *wa = T(1) - v;
            *wb = v;
            *wc = T(0);
            return dx * dx + dy * dy + dz * dz;
        }

        const T cqx = qx - cx, cqy = qy - cy, cqz = qz - cz;
        const T d5 = abx * cqx + aby * cqy + abz * cqz;
        const T d6 = acx * cqx + acy * cqy + acz * cqz;
        if (d6 >= T(0) && d5 <= d6) {
            *wa = T(0);
            *wb = T(0);
            *wc = T(1);
            return cqx * cqx + cqy * cqy + cqz * cqz;
        }

        const T vb = d5 * d2 - d1 * d6;
        if (vb <= T(0) && d2 >= T(0) && d6 <= T(0)) {
            const T w = d2 / (d2 - d6);
            const T x = ax + w * acx;
            const T y = ay + w * acy;
            const T z = az + w * acz;
            const T dx = qx - x, dy = qy - y, dz = qz - z;
            *wa = T(1) - w;
            *wb = T(0);
            *wc = w;
            return dx * dx + dy * dy + dz * dz;
        }

        const T va = d3 * d6 - d5 * d4;
        if (va <= T(0) && (d4 - d3) >= T(0) && (d5 - d6) >= T(0)) {
            const T w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
            const T x = bx + w * bcx;
            const T y = by + w * bcy;
            const T z = bz + w * bcz;
            const T dx = qx - x, dy = qy - y, dz = qz - z;
            *wa = T(0);
            *wb = T(1) - w;
            *wc = w;
            return dx * dx + dy * dy + dz * dz;
        }

        const T denom = T(1) / (va + vb + vc);
        const T v = vb * denom;
        const T w = vc * denom;
        const T u = T(1) - v - w;
        const T x = ax + abx * v + acx * w;
        const T y = ay + aby * v + acy * w;
        const T z = az + abz * v + acz * w;
        const T dx = qx - x, dy = qy - y, dz = qz - z;
        *wa = u;
        *wb = v;
        *wc = w;
        return dx * dx + dy * dy + dz * dz;
    }

    template <typename T>
    inline static void point_to_quad_closest_point(const T px,
                                                   const T py,
                                                   const T pz,
                                                   // v0
                                                   const T x0,
                                                   const T y0,
                                                   const T z0,
                                                   // v1
                                                   const T x1,
                                                   const T y1,
                                                   const T z1,
                                                   // v2
                                                   const T x2,
                                                   const T y2,
                                                   const T z2,
                                                   // v3
                                                   const T x3,
                                                   const T y3,
                                                   const T z3,
                                                   T *const SSDF_RESTRICT s,
                                                   T *const SSDF_RESTRICT t) {
        const T nwx = (y0 - y1) * (z0 + z1) + (y1 - y2) * (z1 + z2) + (y2 - y3) * (z2 + z3) +
                      (y3 - y0) * (z3 + z0);
        const T nwy = (z0 - z1) * (x0 + x1) + (z1 - z2) * (x1 + x2) + (z2 - z3) * (x2 + x3) +
                      (z3 - z0) * (x3 + x0);
        const T nwz = (x0 - x1) * (y0 + y1) + (x1 - x2) * (y1 + y2) + (x2 - x3) * (y2 + y3) +
                      (x3 - x0) * (y3 + y0);
        const T anx = std::abs(nwx);
        const T any = std::abs(nwy);
        const T anz = std::abs(nwz);

        T u0, v0, u1, v1, u2, v2, u3, v3;
        if (anx >= any && anx >= anz) {
            u0 = y0;
            v0 = z0;
            u1 = y1;
            v1 = z1;
            u2 = y2;
            v2 = z2;
            u3 = y3;
            v3 = z3;
        } else if (any >= anz) {
            u0 = x0;
            v0 = z0;
            u1 = x1;
            v1 = z1;
            u2 = x2;
            v2 = z2;
            u3 = x3;
            v3 = z3;
        } else {
            u0 = x0;
            v0 = y0;
            u1 = x1;
            v1 = y1;
            u2 = x2;
            v2 = y2;
            u3 = x3;
            v3 = y3;
        }

        auto cross2 = [](const T ax, const T ay, const T bx, const T by, const T cx, const T cy) -> T {
            return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
        };

        const T area2 = u0 * v1 - v0 * u1 + u1 * v2 - v1 * u2 + u2 * v3 - v2 * u3 + u3 * v0 - v3 * u0;
        const T turn1 = cross2(u0, v0, u1, v1, u2, v2);
        const T turn3 = cross2(u2, v2, u3, v3, u0, v0);

        bool diagonal02 = true;
        if (std::abs(area2) > std::numeric_limits<T>::epsilon() * T(64)) {
            const T orient = area2 > T(0) ? T(1) : T(-1);
            diagonal02 = !(turn1 * orient < T(0) || turn3 * orient < T(0));
        } else {
            const T d02x = x2 - x0, d02y = y2 - y0, d02z = z2 - z0;
            const T d13x = x3 - x1, d13y = y3 - y1, d13z = z3 - z1;
            diagonal02 = (d02x * d02x + d02y * d02y + d02z * d02z) <= (d13x * d13x + d13y * d13y + d13z * d13z);
        }

        T a, b, c;
        T best_s, best_t;
        T best;
        T candidate_s, candidate_t;
        T candidate;
        if (diagonal02) {
            best = point_triangle_barycentric(px, py, pz, x0, y0, z0, x1, y1, z1, x2, y2, z2, &a, &b, &c);
            best_s = b + c;
            best_t = c;

            candidate = point_triangle_barycentric(px, py, pz, x0, y0, z0, x2, y2, z2, x3, y3, z3, &a, &b, &c);
            candidate_s = b;
            candidate_t = b + c;
        } else {
            best = point_triangle_barycentric(px, py, pz, x0, y0, z0, x1, y1, z1, x3, y3, z3, &a, &b, &c);
            best_s = b;
            best_t = c;

            candidate = point_triangle_barycentric(px, py, pz, x1, y1, z1, x2, y2, z2, x3, y3, z3, &a, &b, &c);
            candidate_s = a + b;
            candidate_t = b + c;
        }

        if (candidate < best) {
            best_s = candidate_s;
            best_t = candidate_t;
        }

        *s = std::clamp(best_s, T(0), T(1));
        *t = std::clamp(best_t, T(0), T(1));
    }
    template <typename G, typename T>
    static void compute_aabb(const ptrdiff_t npoints,
                             const G *const SSDF_RESTRICT x,
                             const G *const SSDF_RESTRICT y,
                             const G *const SSDF_RESTRICT z,
                             T *const SSDF_RESTRICT xmin,
                             T *const SSDF_RESTRICT xmax,
                             T *const SSDF_RESTRICT ymin,
                             T *const SSDF_RESTRICT ymax,
                             T *const SSDF_RESTRICT zmin,
                             T *const SSDF_RESTRICT zmax) {
        // Reduce on scalars (OpenMP reduction cannot operate on pointer variables).
        T xmin_v = T(x[0]);
        T xmax_v = T(x[0]);
        T ymin_v = T(y[0]);
        T ymax_v = T(y[0]);
        T zmin_v = T(z[0]);
        T zmax_v = T(z[0]);

#pragma omp parallel for reduction(min : xmin_v, ymin_v, zmin_v) reduction(max : xmax_v, ymax_v, zmax_v)
        for (ptrdiff_t i = 0; i < npoints; ++i) {
            const T xi = T(x[i]);
            const T yi = T(y[i]);
            const T zi = T(z[i]);

            xmin_v = std::min(xmin_v, xi);
            xmax_v = std::max(xmax_v, xi);
            ymin_v = std::min(ymin_v, yi);
            ymax_v = std::max(ymax_v, yi);
            zmin_v = std::min(zmin_v, zi);
            zmax_v = std::max(zmax_v, zi);
        }

        *xmin = xmin_v;
        *xmax = xmax_v;
        *ymin = ymin_v;
        *ymax = ymax_v;
        *zmin = zmin_v;
        *zmax = zmax_v;
    }

    template <typename G, typename T>
    static void all_points_aabb_signed_distance(const ptrdiff_t npoints,
                                                const G *const SSDF_RESTRICT x,
                                                const G *const SSDF_RESTRICT y,
                                                const G *const SSDF_RESTRICT z,
                                                const T minx,
                                                const T maxx,
                                                const T miny,
                                                const T maxy,
                                                const T minz,
                                                const T maxz,
                                                T *const SSDF_RESTRICT out) {
#pragma omp parallel for
        for (ptrdiff_t i = 0; i < npoints; ++i) {
            out[i] = point_aabb_signed_distance(T(x[i]), T(y[i]), T(z[i]), minx, maxx, miny, maxy, minz, maxz);
        }
    }

}  // namespace ssdf

#endif  // EDF_HPP
