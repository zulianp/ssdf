#ifndef EDF_HPP
#define EDF_HPP

#include <algorithm>
#include <cmath>
#include <cstdlib>

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
        *xmin = x[0];
        *xmax = x[0];
        *ymin = y[0];
        *ymax = y[0];
        *zmin = z[0];
        *zmax = z[0];
#pragma omp parallel for reduction(min : xmin, ymin, zmin) reduction(max : xmax, ymax, zmax)
        for (ptrdiff_t i = 0; i < npoints; ++i) {
            *xmin = std::min({x[i], *xmin});
            *xmax = std::max({x[i], *xmax});
            *ymin = std::min({y[i], *ymax});
            *ymax = std::max({y[i], *ymax});
            *zmin = std::min({z[i], *zmax});
            *zmax = std::max({z[i], *zmax});
        }
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
