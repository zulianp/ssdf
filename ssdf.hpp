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
#include <tuple>
#include <vector>

#include "edf/edf.hpp"

#ifdef SSDF_ENABLE_CUBIQL
#include "bvh.hpp"
#endif


    
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

#define SSDF_TIMER(name) Timer t_##name(#name);
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

    // template <typename T>
    // inline static void point_triangle_distance_vector(const T px,
    //                                                   const T py,
    //                                                   const T pz,
    //                                                   const T ax,
    //                                                   const T ay,
    //                                                   const T az,
    //                                                   const T bx,
    //                                                   const T by,
    //                                                   const T bz,
    //                                                   const T cx,
    //                                                   const T cy,
    //                                                   const T cz,
    //                                                   T *const SSDF_RESTRICT outx,
    //                                                   T *const SSDF_RESTRICT outy,
    //                                                   T *const SSDF_RESTRICT outz) {
    //     //  Instead of returning the distance, store the vector to the closest point on the triangle
    //     const T abx = bx - ax, aby = by - ay, abz = bz - az;
    //     const T acx = cx - ax, acy = cy - ay, acz = cz - az;
    //     const T apx = px - ax, apy = py - ay, apz = pz - az;

    //     const T d1 = abx * apx + aby * apy + abz * apz;
    //     const T d2 = acx * apx + acy * apy + acz * apz;
    //     if (d1 <= T(0) && d2 <= T(0)) {
    //         outx[0] = ax - px;
    //         outy[0] = ay - py;
    //         outz[0] = az - pz;
    //     }
    //     const T bpx = px - bx, bpy = py - by, bpz = pz - bz;
    //     const T d3 = abx * bpx + aby * bpy + abz * bpz;
    //     const T d4 = acx * bpx + acy * bpy + acz * bpz;
    //     if (d3 >= T(0) && d4 <= d3) {
    //         outx[0] = bx - px;
    //         outy[0] = by - py;
    //         outz[0] = bz - pz;
    //     }
    //     const T cpx = px - cx, cpy = py - cy, cpz = pz - cz;
    //     const T d5 = abx * cpx + aby * cpy + abz * cpz;
    //     const T d6 = acx * cpx + acy * cpy + acz * cpz;
    //     if (d6 >= T(0) && d5 <= d6) {
    //         outx[0] = cx - px;
    //         outy[0] = cy - py;
    //         outz[0] = cz - pz;
    //     }
    //     const T vb = d5 * d2 - d1 * d6;
    //     if (vb <= T(0) && d2 >= T(0) && d6 <= T(0)) {
    //         outx[0] = ax + d2 / (d2 - d6) * (cx - ax);
    //         outy[0] = ay + d2 / (d2 - d6) * (cy - ay);
    //         outz[0] = az + d2 / (d2 - d6) * (cz - az);
    //     }
    //     const T va = d3 * d6 - d5 * d4;
    //     if (va <= T(0) && (d4 - d3) >= T(0) && (d5 - d6) >= T(0)) {
    //         outx[0] = bx + (d4 - d3) / ((d4 - d3) + (d5 - d6)) * (cx - bx);
    //         outy[0] = by + (d4 - d3) / ((d4 - d3) + (d5 - d6)) * (cy - by);
    //         outz[0] = bz + (d4 - d3) / ((d4 - d3) + (d5 - d6)) * (cz - bz);
    //     }
    //     const T denom = T(1) / (va + vb + vc);
    //     const T v = vb * denom;
    //     const T w = vc * denom;
    //     outx[0] = ax + abx * v + acx * w;
    //     outy[0] = ay + aby * v + acy * w;
    //     outz[0] = az + abz * v + acz * w;
    // }

    // template <typename T>
    // inline static T point_triangle_dist_sq_approx(const T px,
    //                                               const T py,
    //                                               const T pz,
    //                                               const T ax,
    //                                               const T ay,
    //                                               const T az,
    //                                               const T bx,
    //                                               const T by,
    //                                               const T bz,
    //                                               const T cx,
    //                                               const T cy,
    //                                               const T cz) {
    //     const T dx = std::min({std::abs(px - ax), std::abs(px - bx), std::abs(px - cx)});
    //     const T dy = std::min({std::abs(py - ay), std::abs(py - by), std::abs(py - cy)});
    //     const T dz = std::min({std::abs(pz - az), std::abs(pz - bz), std::abs(pz - cz)});
    //     return dx * dx + dy * dy + dz * dz;
    // }

    // template <typename T>
    // inline static T point_aabb_distance_squared(const T px,
    //                                             const T py,
    //                                             const T pz,
    //                                             const T minx,
    //                                             const T maxx,
    //                                             const T miny,
    //                                             const T maxy,
    //                                             const T minz,
    //                                             const T maxz) {
    //     const T dx = (px < minx) ? (minx - px) : (px > maxx ? px - maxx : T(0));
    //     const T dy = (py < miny) ? (miny - py) : (py > maxy ? py - maxy : T(0));
    //     const T dz = (pz < minz) ? (minz - pz) : (pz > maxz ? pz - maxz : T(0));
    //     return dx * dx + dy * dy + dz * dz;
    // }

    // Conservative AABB check for a single element
    template <typename T>
    inline static bool aabb_can_improve(const T px,
                                        const T py,
                                        const T pz,
                                        const T best_sq,
                                        const T minx,
                                        const T maxx,
                                        const T miny,
                                        const T maxy,
                                        const T minz,
                                        const T maxz) {
        return point_aabb_distance_squared(px, py, pz, minx, maxx, miny, maxy, minz, maxz) < best_sq;
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
        if (nselements == 0 || nspoints == 0) return 0;
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
                    // #pragma omp parallel for reduction(inscan, max : acc)
                    for (ptrdiff_t i = 0; i < nselements; i++) {
                        acc = std::max(acc, surf_max[i]);

                        // #pragma omp scan inclusive(acc)
                        cum_max[i] = acc;
                    }
                }
            }

// For each point, find closest surface element
#pragma omp parallel for
            for (ptrdiff_t p = 0; p < npoints; p++) {
                const G px = x[p], py = y[p], pz = z[p];
                const T pcoord = pnt_coord[p];
                T best_sq = out[p] * out[p];

                // Binary search for insertion point
                ptrdiff_t left = std::lower_bound(surf_min, surf_min + nselements, pcoord) - surf_min;

                // Check elements to the left
                for (ptrdiff_t i = (left > 0 ? left - 1 : 0); i >= 0; i--) {
                    const T margin = pcoord - cum_max[i];
                    if (margin >= T(0) && margin * margin >= best_sq) break;
                    // const T span = pcoord - surf_max[i];
                    // if (span * span > best_sq) continue;

                    const I orig_idx = sort_idx[i];
                    const I i0 = s0[orig_idx], i1 = s1[orig_idx], i2 = s2[orig_idx];
                    const G tx0 = sx[i0], tx1 = sx[i1], tx2 = sx[i2];
                    const G ty0 = sy[i0], ty1 = sy[i1], ty2 = sy[i2];
                    const G tz0 = sz[i0], tz1 = sz[i1], tz2 = sz[i2];

                    if (!aabb_can_improve<T>(px,
                                             py,
                                             pz,
                                             best_sq,
                                             std::min({tx0, tx1, tx2}),
                                             std::max({tx0, tx1, tx2}),
                                             std::min({ty0, ty1, ty2}),
                                             std::max({ty0, ty1, ty2}),
                                             std::min({tz0, tz1, tz2}),
                                             std::max({tz0, tz1, tz2}))) {
                        continue;
                    }

                    const T dist_sq = point_triangle_dist_sq(px, py, pz, tx0, ty0, tz0, tx1, ty1, tz1, tx2, ty2, tz2);
                    best_sq = std::min(best_sq, dist_sq);
                }

                // Check elements to the right
                for (ptrdiff_t i = left; i < nselements; i++) {
                    const T margin = surf_min[i] - pcoord;
                    if (margin >= T(0) && margin * margin >= best_sq) break;

                    const I orig_idx = sort_idx[i];

                    const I i0 = s0[orig_idx], i1 = s1[orig_idx], i2 = s2[orig_idx];
                    const G tx0 = sx[i0], tx1 = sx[i1], tx2 = sx[i2];
                    const G ty0 = sy[i0], ty1 = sy[i1], ty2 = sy[i2];
                    const G tz0 = sz[i0], tz1 = sz[i1], tz2 = sz[i2];

                    if (!aabb_can_improve<T>(px,
                                             py,
                                             pz,
                                             best_sq,
                                             std::min({tx0, tx1, tx2}),
                                             std::max({tx0, tx1, tx2}),
                                             std::min({ty0, ty1, ty2}),
                                             std::max({ty0, ty1, ty2}),
                                             std::min({tz0, tz1, tz2}),
                                             std::max({tz0, tz1, tz2}))) {
                        continue;
                    }

                    const T dist_sq = point_triangle_dist_sq(px, py, pz, tx0, ty0, tz0, tx1, ty1, tz1, tx2, ty2, tz2);
                    best_sq = std::min(best_sq, dist_sq);
                }

                out[p] = std::sqrt(best_sq);
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
     * - High-performance EDF using a 2D cell grid + sorted dimension
     */
    template <typename G, typename T, typename I>
    int edf_celllist(const ptrdiff_t npoints,
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
        SSDF_TIMER(edf_celllist);
        if (nselements == 0 || nspoints == 0) return 0;

        // Triangle AABBs
        std::vector<G> tminx(nselements), tminy(nselements), tminz(nselements);
        std::vector<G> tmaxx(nselements), tmaxy(nselements), tmaxz(nselements);

        G gminx = std::numeric_limits<G>::max();
        G gminy = std::numeric_limits<G>::max();
        G gminz = std::numeric_limits<G>::max();

        G gmaxx = std::numeric_limits<G>::lowest();
        G gmaxy = std::numeric_limits<G>::lowest();
        G gmaxz = std::numeric_limits<G>::lowest();

        G max_extent = G(0);

#pragma omp parallel
        {
            G lminx = std::numeric_limits<G>::max();
            G lminy = std::numeric_limits<G>::max();
            G lminz = std::numeric_limits<G>::max();

            G lmaxx = std::numeric_limits<G>::lowest();
            G lmaxy = std::numeric_limits<G>::lowest();
            G lmaxz = std::numeric_limits<G>::lowest();

            G lextent = G(0);
#pragma omp for nowait
            for (ptrdiff_t i = 0; i < nselements; ++i) {
                const I i0 = s0[i], i1 = s1[i], i2 = s2[i];
                const G x0 = sx[i0], x1 = sx[i1], x2 = sx[i2];
                const G y0 = sy[i0], y1 = sy[i1], y2 = sy[i2];
                const G z0 = sz[i0], z1 = sz[i1], z2 = sz[i2];

                const G xmin = std::min({x0, x1, x2});
                const G xmax = std::max({x0, x1, x2});
                const G ymin = std::min({y0, y1, y2});
                const G ymax = std::max({y0, y1, y2});
                const G zmin = std::min({z0, z1, z2});
                const G zmax = std::max({z0, z1, z2});

                tminx[i] = xmin;
                tmaxx[i] = xmax;
                tminy[i] = ymin;
                tmaxy[i] = ymax;
                tminz[i] = zmin;
                tmaxz[i] = zmax;

                lminx = std::min(lminx, xmin);
                lmaxx = std::max(lmaxx, xmax);
                lminy = std::min(lminy, ymin);
                lmaxy = std::max(lmaxy, ymax);
                lminz = std::min(lminz, zmin);
                lmaxz = std::max(lmaxz, zmax);
                lextent = std::max(lextent, std::max({xmax - xmin, ymax - ymin, zmax - zmin}));
            }
#pragma omp critical
            {
                gminx = std::min(gminx, lminx);
                gmaxx = std::max(gmaxx, lmaxx);
                gminy = std::min(gminy, lminy);
                gmaxy = std::max(gmaxy, lmaxy);
                gminz = std::min(gminz, lminz);
                gmaxz = std::max(gmaxz, lmaxz);
                max_extent = std::max(max_extent, lextent);
            }
        }

        // Choose grid axes: two widest spans
        const G span_x = gmaxx - gminx;
        const G span_y = gmaxy - gminy;
        const G span_z = gmaxz - gminz;
        int axis0 = 0, axis1 = 1, sort_axis = 2;
        {
            if (span_x >= span_y && span_x >= span_z) {
                axis0 = 0;
                if (span_y >= span_z) {
                    axis1 = 1;
                    sort_axis = 2;
                } else {
                    axis1 = 2;
                    sort_axis = 1;
                }
            } else if (span_y >= span_x && span_y >= span_z) {
                axis0 = 1;
                if (span_x >= span_z) {
                    axis1 = 0;
                    sort_axis = 2;
                } else {
                    axis1 = 2;
                    sort_axis = 0;
                }
            } else {
                axis0 = 2;
                if (span_x >= span_y) {
                    axis1 = 0;
                    sort_axis = 1;
                } else {
                    axis1 = 1;
                    sort_axis = 0;
                }
            }
        }

        const G cell_size = std::max(max_extent, G(1e-8));
        const G axis_min0 = (axis0 == 0) ? gminx : (axis0 == 1) ? gminy : gminz;
        const G axis_min1 = (axis1 == 0) ? gminx : (axis1 == 1) ? gminy : gminz;
        const G axis_span0 = (axis0 == 0) ? span_x : (axis0 == 1) ? span_y : span_z;
        const G axis_span1 = (axis1 == 0) ? span_x : (axis1 == 1) ? span_y : span_z;

        const I ncell0 = std::max(1, I(std::ceil(axis_span0 / cell_size)));
        const I ncell1 = std::max(1, I(std::ceil(axis_span1 / cell_size)));
        const ptrdiff_t ncells = ptrdiff_t(ncell0) * ptrdiff_t(ncell1);

        const G *tmin_axis0 = (axis0 == 0) ? tminx.data() : (axis0 == 1) ? tminy.data() : tminz.data();
        const G *tmin_axis1 = (axis1 == 0) ? tminx.data() : (axis1 == 1) ? tminy.data() : tminz.data();

        std::vector<I> cell_counts(ncells, 0);
        std::vector<G> cell_minx(ncells, std::numeric_limits<G>::max());
        std::vector<G> cell_miny(ncells, std::numeric_limits<G>::max());
        std::vector<G> cell_minz(ncells, std::numeric_limits<G>::max());
        std::vector<G> cell_maxx(ncells, std::numeric_limits<G>::lowest());
        std::vector<G> cell_maxy(ncells, std::numeric_limits<G>::lowest());
        std::vector<G> cell_maxz(ncells, std::numeric_limits<G>::lowest());

        auto cell_id = [&](G ax0, G ax1) -> ptrdiff_t {
            I ix = I((ax0 - axis_min0) / cell_size);
            I iy = I((ax1 - axis_min1) / cell_size);
            ix = std::max(0, std::min(ix, ncell0 - 1));
            iy = std::max(0, std::min(iy, ncell1 - 1));
            return ptrdiff_t(iy) * ncell0 + ix;
        };

        // Count and aggregate cell AABBs
        for (ptrdiff_t i = 0; i < nselements; ++i) {
            const ptrdiff_t cid = cell_id(tmin_axis0[i], tmin_axis1[i]);
            cell_counts[cid] += 1;
            cell_minx[cid] = std::min(cell_minx[cid], tminx[i]);
            cell_maxx[cid] = std::max(cell_maxx[cid], tmaxx[i]);
            cell_miny[cid] = std::min(cell_miny[cid], tminy[i]);
            cell_maxy[cid] = std::max(cell_maxy[cid], tmaxy[i]);
            cell_minz[cid] = std::min(cell_minz[cid], tminz[i]);
            cell_maxz[cid] = std::max(cell_maxz[cid], tmaxz[i]);
        }

        // Prefix sums
        std::vector<ptrdiff_t> cell_ptr(ncells + 1, 0);
        for (ptrdiff_t c = 0; c < ncells; ++c) cell_ptr[c + 1] = cell_ptr[c] + cell_counts[c];

        std::vector<I> cell_idx(nselements);
        std::vector<ptrdiff_t> fill_ptr = cell_ptr;
        for (ptrdiff_t i = 0; i < nselements; ++i) {
            const ptrdiff_t cid = cell_id(tmin_axis0[i], tmin_axis1[i]);
            cell_idx[fill_ptr[cid]++] = I(i);
        }

        const G *sorted_min = (sort_axis == 0) ? tminx.data() : (sort_axis == 1) ? tminy.data() : tminz.data();
        const G *sorted_max = (sort_axis == 0) ? tmaxx.data() : (sort_axis == 1) ? tmaxy.data() : tmaxz.data();

// Sort per cell on sort_axis
#pragma omp parallel for
        for (ptrdiff_t c = 0; c < ncells; ++c) {
            const ptrdiff_t begin = cell_ptr[c];
            const ptrdiff_t end = cell_ptr[c + 1];
            if (begin == end) continue;

            std::sort(cell_idx.begin() + begin, cell_idx.begin() + end, [&](I a, I b) {
                return sorted_min[a] < sorted_min[b];
            });
        }

        {  // Permute boxes to avoid aabb indirections during queries
            std::vector<G> temp(nselements);
            G *boxes[6] = {tminx.data(), tminy.data(), tminz.data(), tmaxx.data(), tmaxy.data(), tmaxz.data()};

            for (int d = 0; d < 6; d++) {
                memcpy(temp.data(), boxes[d], sizeof(G) * nselements);
                for (ptrdiff_t i = 0; i < nselements; i++) {
                    boxes[d][i] = temp[cell_idx[i]];
                }
            }
        }

        // build cumulative maxima per cell
        std::vector<G> cum_max(nselements, G(0));
        for (ptrdiff_t c = 0; c < ncells; ++c) {
            const ptrdiff_t begin = cell_ptr[c];
            const ptrdiff_t end = cell_ptr[c + 1];
            if (begin == end) continue;

            G acc = sorted_max[begin];
            for (ptrdiff_t i = begin; i < end; ++i) {
                acc = std::max(acc, sorted_max[i]);
                cum_max[i] = acc;
            }
        }

        {
            SSDF_TIMER(queries);

            auto process_cell = [&](const ptrdiff_t cid, const G px, const G py, const G pz, T &best_sq) {
                const ptrdiff_t begin = cell_ptr[cid];
                const ptrdiff_t end = cell_ptr[cid + 1];
                if (begin == end) return;

                const G pcoord = (sort_axis == 0) ? px : (sort_axis == 1) ? py : pz;
                const ptrdiff_t left = std::lower_bound(sorted_min + begin, sorted_min + end, pcoord) - sorted_min;

                auto evaluate_triangle = [&](const ptrdiff_t idx) {
                    const I tid = cell_idx[idx];
                    const I i0 = s0[tid], i1 = s1[tid], i2 = s2[tid];
                    const G dist_sq = point_triangle_dist_sq(
                        px, py, pz, sx[i0], sy[i0], sz[i0], sx[i1], sy[i1], sz[i1], sx[i2], sy[i2], sz[i2]);
                    if (dist_sq < best_sq) best_sq = dist_sq;
                };

                // Scan left
                for (ptrdiff_t i = (left > begin) ? left - 1 : begin; i >= begin; --i) {
                    const G margin = pcoord - cum_max[i];
                    if (margin >= G(0) && margin * margin >= best_sq) break;

                    if (aabb_can_improve<T>(
                            px, py, pz, best_sq, tminx[i], tmaxx[i], tminy[i], tmaxy[i], tminz[i], tmaxz[i])) {
                        evaluate_triangle(i);
                    }
                }

                // Scan right
                for (ptrdiff_t i = left; i < end; ++i) {
                    const G margin = sorted_min[i] - pcoord;
                    if (margin >= G(0) && margin * margin >= best_sq) break;

                    if (!aabb_can_improve<T>(
                            px, py, pz, best_sq, tminx[i], tmaxx[i], tminy[i], tmaxy[i], tminz[i], tmaxz[i])) {
                        continue;
                    }

                    evaluate_triangle(i);
                }
            };

// // warmup distance
// #pragma omp parallel for
//             for (ptrdiff_t p = 0; p < npoints; p++) {
//                 T sqdist =
//                     std::sqrt(point_aabb_distance_squared(x[p], y[p], z[p], gminx, gmaxx, gminy, gmaxy, gminz, gmaxz));
//                 out[p] = MIN(out[p], sqdist);
//             }

            const G *const xyz[3] = {x, y, z};

            // Query points
#pragma omp parallel for
            for (ptrdiff_t p = 0; p < npoints; ++p) {
                const G px = x[p], py = y[p], pz = z[p];
                T best_sq = out[p] * out[p];

                // TODO: change the brute force by
                // 1) finding cell containing the point
                // 2) process content of the cell
                // 3) Check sourrounding cells outword ring by ring
                // Make sure that the loop ends when no improvements are possible

                const ptrdiff_t first_cid = cell_id(xyz[axis0][p], xyz[axis1][p]);
                if (cell_counts[first_cid] != 0 && aabb_can_improve<T>(px,
                                                                       py,
                                                                       pz,
                                                                       best_sq,
                                                                       cell_minx[first_cid],
                                                                       cell_maxx[first_cid],
                                                                       cell_miny[first_cid],
                                                                       cell_maxy[first_cid],
                                                                       cell_minz[first_cid],
                                                                       cell_maxz[first_cid])) {
                    process_cell(first_cid, px, py, pz, best_sq);
                }

                // Iterate cells with conservative culling
                for (ptrdiff_t cid = 0; cid < ncells; ++cid) {
                    if (cell_counts[cid] == 0 ||
                        !aabb_can_improve<T>(px,
                                             py,
                                             pz,
                                             best_sq,
                                             cell_minx[cid],
                                             cell_maxx[cid],
                                             cell_miny[cid],
                                             cell_maxy[cid],
                                             cell_minz[cid],
                                             cell_maxz[cid]) ||
                        first_cid == cid) {
                        continue;
                    }

                    process_cell(cid, px, py, pz, best_sq);
                }

                out[p] = std::sqrt(best_sq);
            }
        }

        return 0;
    }

    template <typename G, typename T, typename I>
    int edf_select(const ptrdiff_t npoints,
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
        int SSDF_USE_CELL_LIST = 1;
        int SSDF_CELL_VALIDATE = 0;
        SSDF_READ_ENV(SSDF_USE_CELL_LIST, atoi);
        SSDF_READ_ENV(SSDF_CELL_VALIDATE, atoi);
#ifdef SSDF_ENABLE_CUBIQL
        int SSDF_USE_BVH = 0;
        SSDF_READ_ENV(SSDF_USE_BVH, atoi);
        if (SSDF_USE_BVH) {
            return edf_bvh<G, T, I>(npoints, x, y, z, nselements, s0, s1, s2, nspoints, sx, sy, sz, out);
        }
#endif

        if (SSDF_CELL_VALIDATE) {
            // Just for testing consistency
            std::vector<T> ref(out, out + npoints);
            std::vector<T> test(out, out + npoints);
            edf<G, T, I>(npoints, x, y, z, nselements, s0, s1, s2, nspoints, sx, sy, sz, ref.data());
            edf_celllist<G, T, I>(npoints, x, y, z, nselements, s0, s1, s2, nspoints, sx, sy, sz, test.data());
            T max_diff = 0;
            T diff_norm = 0;
            for (ptrdiff_t i = 0; i < npoints; ++i) {
                T rel_diff = std::abs(ref[i] - test[i]) / (std::abs(ref[i]) + 1e-16);
                max_diff = std::max<T>(max_diff, rel_diff);
                diff_norm += rel_diff * rel_diff;
            }
            if (max_diff > static_cast<T>(1e-5)) {
                std::cerr << "EDF cell-list validation max rel diff: " << max_diff << std::endl;
                std::cerr << "EDF cell-list validation rel diff norm: " << std::sqrt(diff_norm) << std::endl;
            }
            std::copy(test.begin(), test.end(), out);
            return 0;
        }

        if (SSDF_USE_CELL_LIST) {
            return edf_celllist<G, T, I>(npoints, x, y, z, nselements, s0, s1, s2, nspoints, sx, sy, sz, out);
        }
        return edf<G, T, I>(npoints, x, y, z, nselements, s0, s1, s2, nspoints, sx, sy, sz, out);
    }
}  // namespace ssdf

#endif  // SSDF_HPP
