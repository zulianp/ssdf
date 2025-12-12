#ifndef SSDF_HPP
#define SSDF_HPP

#include <stddef.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <string>

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

    // Return current time in milliseconds
    inline double time_ms() {
        auto now = std::chrono::high_resolution_clock::now();
        auto duration = now.time_since_epoch();
        return std::chrono::duration<double, std::milli>(duration).count();
    }

    struct Timer {
        std::string name;
        double start_time;
        double end_time;
        double duration;
        inline Timer(const std::string &name) : name(name) { start_time = time_ms(); }
        inline ~Timer() {
            end_time = time_ms();
            duration = end_time - start_time;
            print();
        }
        inline void print() const { std::cout << name << " took " << duration << " ms" << std::endl; }
    };

// #define SSDF_TIMER(name) Timer t_##name(#name);
    // #define SSDF_TIMER(...)

    // Compute squared distance from point p to triangle (a, b, c)
    template <typename T>
    inline static T point_triangle_dist_sq(const T px,
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
        // Based on Real-Time Collision Detection (Christer Ericson)
        const T abx = bx - ax, aby = by - ay, abz = bz - az;
        const T acx = cx - ax, acy = cy - ay, acz = cz - az;
        const T apx = px - ax, apy = py - ay, apz = pz - az;

        const T d1 = abx * apx + aby * apy + abz * apz;
        const T d2 = acx * apx + acy * apy + acz * apz;
        if (d1 <= T(0) && d2 <= T(0)) return apx * apx + apy * apy + apz * apz;  // barycentric (1,0,0)

        const T bpx = px - bx, bpy = py - by, bpz = pz - bz;
        const T d3 = abx * bpx + aby * bpy + abz * bpz;
        const T d4 = acx * bpx + acy * bpy + acz * bpz;
        if (d3 >= T(0) && d4 <= d3) return bpx * bpx + bpy * bpy + bpz * bpz;  // barycentric (0,1,0)

        const T vc = d1 * d4 - d3 * d2;
        if (vc <= T(0) && d1 >= T(0) && d3 <= T(0)) {
            const T v = d1 / (d1 - d3);
            const T projx = ax + v * abx;
            const T projy = ay + v * aby;
            const T projz = az + v * abz;
            const T dx = px - projx, dy = py - projy, dz = pz - projz;
            return dx * dx + dy * dy + dz * dz;  // edge AB
        }

        const T cpx = px - cx, cpy = py - cy, cpz = pz - cz;
        const T d5 = abx * cpx + aby * cpy + abz * cpz;
        const T d6 = acx * cpx + acy * cpy + acz * cpz;
        if (d6 >= T(0) && d5 <= d6) return cpx * cpx + cpy * cpy + cpz * cpz;  // barycentric (0,0,1)

        const T vb = d5 * d2 - d1 * d6;
        if (vb <= T(0) && d2 >= T(0) && d6 <= T(0)) {
            const T w = d2 / (d2 - d6);
            const T projx = ax + w * acx;
            const T projy = ay + w * acy;
            const T projz = az + w * acz;
            const T dx = px - projx, dy = py - projy, dz = pz - projz;
            return dx * dx + dy * dy + dz * dz;  // edge AC
        }

        const T va = d3 * d6 - d5 * d4;
        if (va <= T(0) && (d4 - d3) >= T(0) && (d5 - d6) >= T(0)) {
            const T w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
            const T projx = bx + w * (cx - bx);
            const T projy = by + w * (cy - by);
            const T projz = bz + w * (cz - bz);
            const T dx = px - projx, dy = py - projy, dz = pz - projz;
            return dx * dx + dy * dy + dz * dz;  // edge BC
        }

        // Inside face region
        const T denom = T(1) / (va + vb + vc);
        const T v = vb * denom;
        const T w = vc * denom;
        const T projx = ax + abx * v + acx * w;
        const T projy = ay + aby * v + acy * w;
        const T projz = az + abz * v + acz * w;
        const T dx = px - projx, dy = py - projy, dz = pz - projz;
        return dx * dx + dy * dy + dz * dz;
    }

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

    // Conservative AABB check for a single element
    template <typename T>
    inline static bool aabb_can_improve(const T px,
                                        const T py,
                                        const T pz,
                                        const T best_dist_sq,
                                        const T minx,
                                        const T maxx,
                                        const T miny,
                                        const T maxy,
                                        const T minz,
                                        const T maxz) {
        const T dxmin = px - minx;
        const T dxmax = px - maxx;
        const T dymin = py - miny;
        const T dymax = py - maxy;
        const T dzmin = pz - minz;
        const T dzmax = pz - maxz;
        const T dx = MIN(dxmin * dxmin, dxmax * dxmax);
        const T dy = MIN(dymin * dymin, dymax * dymax);
        const T dz = MIN(dzmin * dzmin, dzmax * dzmax);
        const T dist_sq = dx + dy + dz;
        return dist_sq < best_dist_sq ||
               (dxmin >= 0 && dxmax <= 0 && dymin >= 0 && dymax <= 0 && dzmin >= 0 && dzmax <= 0);
    }

    /**
     * @brief Compute unsigned point-to-surface distances.
     *
     * @tparam G Geometry type (float/double) for coordinates.
     * @tparam T Output distance type.
     * @tparam I Index type for triangle vertices.
     * @param npoints Number of query points.
     * @param x,y,z Arrays of size npoints with point coordinates (SoA).
     * @param nselements Number of surface triangles.
     * @param s0,s1,s2 Arrays of size nselements with vertex indices (triangle list).
     * @param nspoints Number of surface vertices.
     * @param sx,sy,sz Arrays of size nspoints with surface vertex coordinates (SoA).
     * @param out Output array of size npoints. On input, values are treated as
     *            current best distances (initialize to large values if unused).
     * @return int 0 on success.
     *
     * Notes:
     * - Distances are unsigned (closest-point Euclidean distance).
     * - OpenMP is used when enabled at build time for parallelism.
     */
    template <typename G, typename T, typename I>
    int edf(const ptrdiff_t npoints,
            const G *const SSDF_RESTRICT x,
            const G *const SSDF_RESTRICT y,
            const G *const SSDF_RESTRICT z,
            const ptrdiff_t nselements,
            const I *const SSDF_RESTRICT s0,
            const I *const SSDF_RESTRICT s1,
            const I *const SSDF_RESTRICT s2,
            const ptrdiff_t nspoints,
            const G *const SSDF_RESTRICT sx,
            const G *const SSDF_RESTRICT sy,
            const G *const SSDF_RESTRICT sz,
            T *const SSDF_RESTRICT out) {
        SSDF_TIMER(total);
        // Allocate AABB arrays
        T *surf_minx = new T[nselements];
        T *surf_miny = new T[nselements];
        T *surf_minz = new T[nselements];
        T *surf_maxx = new T[nselements];
        T *surf_maxy = new T[nselements];
        T *surf_maxz = new T[nselements];
        I *sort_idx = new I[nselements];
        T *scratch = new T[nselements];
        T *cum_max = new T[nselements];

        // 1) Compute surface AABBs

        {
            SSDF_TIMER(aabb);
#pragma omp parallel for
            for (ptrdiff_t i = 0; i < nselements; i++) {
                const I i0 = s0[i], i1 = s1[i], i2 = s2[i];
                const G x0 = sx[i0], x1 = sx[i1], x2 = sx[i2];
                const G y0 = sy[i0], y1 = sy[i1], y2 = sy[i2];
                const G z0 = sz[i0], z1 = sz[i1], z2 = sz[i2];

                surf_minx[i] = std::min({x0, x1, x2});
                surf_maxx[i] = std::max({x0, x1, x2});
                surf_miny[i] = std::min({y0, y1, y2});
                surf_maxy[i] = std::max({y0, y1, y2});
                surf_minz[i] = std::min({z0, z1, z2});
                surf_maxz[i] = std::max({z0, z1, z2});
            }
        }

        // Process each dimension
        for (int dim = 0; dim < 3; dim++) {
            SSDF_TIMER(dim);

            T *surf_min = (dim == 0) ? surf_minx : (dim == 1) ? surf_miny : surf_minz;
            T *surf_max = (dim == 0) ? surf_maxx : (dim == 1) ? surf_maxy : surf_maxz;
            const G *pnt_coord = (dim == 0) ? x : (dim == 1) ? y : z;

            {
                SSDF_TIMER(sort);
                // Initialize sort index
                for (ptrdiff_t i = 0; i < nselements; i++) {
                    sort_idx[i] = i;
                }

                // Sort by min coordinate
                std::sort(sort_idx, sort_idx + nselements, [surf_min](I a, I b) { return surf_min[a] < surf_min[b]; });

                // Permute arrays
                memcpy(scratch, surf_min, sizeof(T) * nselements);
                for (ptrdiff_t i = 0; i < nselements; i++) {
                    surf_min[i] = scratch[sort_idx[i]];
                }

                memcpy(scratch, surf_max, sizeof(T) * nselements);
                for (ptrdiff_t i = 0; i < nselements; i++) {
                    surf_max[i] = scratch[sort_idx[i]];
                }

                {
                    T acc = surf_max[0];
#pragma omp parallel for reduction(inscan, max : acc)
                    for (ptrdiff_t i = 0; i < nselements; i++) {
                        acc = std::max(acc, surf_max[i]);

#pragma omp scan inclusive(acc)
                        cum_max[i] = acc;
                    }
                }
            }

// For each point, find closest surface element
#pragma omp parallel for
            for (ptrdiff_t p = 0; p < npoints; p++) {
                const G px = x[p], py = y[p], pz = z[p];
                const T pcoord = pnt_coord[p];
                T best_dist_sq = out[p] * out[p];

                // Binary search for insertion point
                ptrdiff_t left = std::lower_bound(surf_min, surf_min + nselements, pcoord) - surf_min;

                // Check elements to the left
                {
                    ptrdiff_t start = (left > 0 ? left - 1 : 0);
                    ptrdiff_t i = start;

                    for (; i >= 0; i--) {
                        const T margin = pcoord - cum_max[i];
                        if (margin * margin > best_dist_sq) break;
                        const T span = pcoord - surf_max[i];
                        if (span * span > best_dist_sq) continue;

                        const I orig_idx = sort_idx[i];
                        const I i0 = s0[orig_idx], i1 = s1[orig_idx], i2 = s2[orig_idx];
                        const G tx0 = sx[i0], tx1 = sx[i1], tx2 = sx[i2];
                        const G ty0 = sy[i0], ty1 = sy[i1], ty2 = sy[i2];
                        const G tz0 = sz[i0], tz1 = sz[i1], tz2 = sz[i2];

                        if (!aabb_can_improve<T>(px,
                                                 py,
                                                 pz,
                                                 best_dist_sq,
                                                 std::min({tx0, tx1, tx2}),
                                                 std::max({tx0, tx1, tx2}),
                                                 std::min({ty0, ty1, ty2}),
                                                 std::max({ty0, ty1, ty2}),
                                                 std::min({tz0, tz1, tz2}),
                                                 std::max({tz0, tz1, tz2}))) {
                            continue;
                        }

                        const T dist_sq =
                            point_triangle_dist_sq(px, py, pz, tx0, ty0, tz0, tx1, ty1, tz1, tx2, ty2, tz2);
                        best_dist_sq = std::min(best_dist_sq, dist_sq);
                    }
                }

                // Check elements to the right
                {
                    ptrdiff_t i = left;
                    ptrdiff_t end = nselements;
                    // Handle remainder (scalar)
                    for (; i < end; i++) {
                        const T margin = surf_min[i] - pcoord;
                        if (margin * margin > best_dist_sq) break;

                        const I orig_idx = sort_idx[i];

                        const I i0 = s0[orig_idx], i1 = s1[orig_idx], i2 = s2[orig_idx];
                        const G tx0 = sx[i0], tx1 = sx[i1], tx2 = sx[i2];
                        const G ty0 = sy[i0], ty1 = sy[i1], ty2 = sy[i2];
                        const G tz0 = sz[i0], tz1 = sz[i1], tz2 = sz[i2];

                        if (!aabb_can_improve<T>(px,
                                                 py,
                                                 pz,
                                                 best_dist_sq,
                                                 std::min({tx0, tx1, tx2}),
                                                 std::max({tx0, tx1, tx2}),
                                                 std::min({ty0, ty1, ty2}),
                                                 std::max({ty0, ty1, ty2}),
                                                 std::min({tz0, tz1, tz2}),
                                                 std::max({tz0, tz1, tz2}))) {
                            continue;
                        }

                        const T dist_sq =
                            point_triangle_dist_sq(px, py, pz, tx0, ty0, tz0, tx1, ty1, tz1, tx2, ty2, tz2);
                        best_dist_sq = std::min(best_dist_sq, dist_sq);
                    }
                }

                out[p] = std::sqrt(best_dist_sq);
            }
        }

        delete[] surf_minx;
        delete[] surf_miny;
        delete[] surf_minz;
        delete[] surf_maxx;
        delete[] surf_maxy;
        delete[] surf_maxz;
        delete[] sort_idx;
        delete[] scratch;
        delete[] cum_max;

        return 0;
    }

}  // namespace ssdf

#endif  // SSDF_HPP
