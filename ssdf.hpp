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

    // Compute squared distance from point p to triangle (a, b, c)
    template <typename T>
    inline T point_triangle_dist_sq(const T px,
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
    inline T point_triangle_dist_sq_approx(const T px,
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
        const T dax = px - ax, day = py - ay, daz = pz - az;
        const T dbx = px - bx, dby = py - by, dbz = pz - bz;
        const T dcx = px - cx, dcy = py - cy, dcz = pz - cz;
        const T da_sq = dax * dax + day * day + daz * daz;
        const T db_sq = dbx * dbx + dby * dby + dbz * dbz;
        const T dc_sq = dcx * dcx + dcy * dcy + dcz * dcz;
        return std::min(da_sq, std::min(db_sq, dc_sq));
    }

    // Conservative check: whether any triangle in the chunk could improve best_dist_sq using AABB
    template <typename T, typename I>
    inline bool chunk_aabb_can_improve(const T px,
                                       const T py,
                                       const T pz,
                                       const T best_dist_sq,
                                       const ptrdiff_t chunk_start,
                                       const I *const SSDF_RESTRICT sort_idx,
                                       const T *const SSDF_RESTRICT surf_minx,
                                       const T *const SSDF_RESTRICT surf_maxx,
                                       const T *const SSDF_RESTRICT surf_miny,
                                       const T *const SSDF_RESTRICT surf_maxy,
                                       const T *const SSDF_RESTRICT surf_minz,
                                       const T *const SSDF_RESTRICT surf_maxz) {
        bool any = false;
#pragma omp simd reduction(| : any)
        for (ptrdiff_t lane = 0; lane < VECTOR_SIZE; ++lane) {
            const ptrdiff_t i = chunk_start + lane;
            const I orig_idx = sort_idx[i];

            const T dx = (px < surf_minx[orig_idx]) ? (surf_minx[orig_idx] - px)
                                                    : ((px > surf_maxx[orig_idx]) ? (px - surf_maxx[orig_idx]) : T(0));
            const T dy = (py < surf_miny[orig_idx]) ? (surf_miny[orig_idx] - py)
                                                    : ((py > surf_maxy[orig_idx]) ? (py - surf_maxy[orig_idx]) : T(0));
            const T dz = (pz < surf_minz[orig_idx]) ? (surf_minz[orig_idx] - pz)
                                                    : ((pz > surf_maxz[orig_idx]) ? (pz - surf_maxz[orig_idx]) : T(0));
            const T dist_sq = dx * dx + dy * dy + dz * dz;
            any = any || (dist_sq < best_dist_sq);
        }
        return any;
    }

    // Helper function to compute point-to-triangle distance for a fixed-size chunk
    // Uses fixed-size loop to enable compiler vectorization
    template <typename G, typename T, typename I>
    inline T compute_triangle_dist_chunk(const G px,
                                         const G py,
                                         const G pz,
                                         const ptrdiff_t chunk_start,
                                         const I *const SSDF_RESTRICT sort_idx,
                                         const I *const SSDF_RESTRICT s0,
                                         const I *const SSDF_RESTRICT s1,
                                         const I *const SSDF_RESTRICT s2,
                                         const G *const SSDF_RESTRICT sx,
                                         const G *const SSDF_RESTRICT sy,
                                         const G *const SSDF_RESTRICT sz) {
        T dist[VECTOR_SIZE];

// Fixed-size loop (VECTOR_SIZE) for compiler vectorization
#pragma omp simd
        for (ptrdiff_t lane = 0; lane < VECTOR_SIZE; ++lane) {
            const ptrdiff_t i = chunk_start + lane;
            const I orig_idx = sort_idx[i];
            const I i0 = s0[orig_idx], i1 = s1[orig_idx], i2 = s2[orig_idx];

            const G tx0 = sx[i0], tx1 = sx[i1], tx2 = sx[i2];
            const G ty0 = sy[i0], ty1 = sy[i1], ty2 = sy[i2];
            const G tz0 = sz[i0], tz1 = sz[i1], tz2 = sz[i2];

            // Compute minimum distance to triangle vertices
            const T dx = std::min({std::abs(px - tx0), std::abs(px - tx1), std::abs(px - tx2)});
            const T dy = std::min({std::abs(py - ty0), std::abs(py - ty1), std::abs(py - ty2)});
            const T dz = std::min({std::abs(pz - tz0), std::abs(pz - tz1), std::abs(pz - tz2)});
            dist[lane] = std::sqrt(dx * dx + dy * dy + dz * dz);
        }

        T best_dist = std::numeric_limits<T>::max();

#pragma omp simd reduction(min : best_dist)
        for (ptrdiff_t lane = 0; lane < VECTOR_SIZE; ++lane) {
            best_dist = MIN(best_dist, dist[lane]);
        }

        return best_dist;
    }

    struct Timer {
        std::string name;
        double start_time;
        double end_time;
        double duration;
        Timer(const std::string &name) : name(name) { start_time = time_ms(); }
        ~Timer() {
            end_time = time_ms();
            duration = end_time - start_time;
            print();
        }
        void print() const { std::cout << name << " took " << duration << " ms" << std::endl; }
    };

    // out is of size npoints
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
            Timer t_aabb("Compute surface AABBs");
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
            Timer t_dim("Process dimension " + std::to_string(dim));

            T *surf_min = (dim == 0) ? surf_minx : (dim == 1) ? surf_miny : surf_minz;
            T *surf_max = (dim == 0) ? surf_maxx : (dim == 1) ? surf_maxy : surf_maxz;
            const G *pnt_coord = (dim == 0) ? x : (dim == 1) ? y : z;

            {
                Timer t_sort("Sort surface elements");
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
                T best_dist = out[p];

                // Binary search for insertion point
                ptrdiff_t left = std::lower_bound(surf_min, surf_min + nselements, pcoord) - surf_min;

                // Check elements to the left (vectorized)
                {
                    ptrdiff_t start = (left > 0 ? left - 1 : 0);
                    ptrdiff_t i = start;

                    // Process in chunks of VECTOR_SIZE for vectorization
                    while (i >= VECTOR_SIZE - 1) {
                        // Early termination check for the chunk
                        if (cum_max[i - (VECTOR_SIZE - 1)] < pcoord - best_dist) break;

                        // Process chunk of VECTOR_SIZE elements
                        ptrdiff_t chunk_start = i - (VECTOR_SIZE - 1);
                        ptrdiff_t chunk_end = i + 1;

                        // Filter elements that can be skipped
                        bool process_chunk = false;
                        for (ptrdiff_t j = chunk_start; j < chunk_end; ++j) {
                            if (surf_max[j] >= pcoord - best_dist) {
                                process_chunk = true;
                                break;
                            }
                        }

                        if (process_chunk) {
                            T chunk_best = compute_triangle_dist_chunk<G, T, I>(
                                px, py, pz, chunk_start, sort_idx, s0, s1, s2, sx, sy, sz);
                            best_dist = std::min(best_dist, chunk_best);
                        }

                        i -= VECTOR_SIZE;
                    }

                    // Handle remainder (scalar)
                    for (; i >= 0; i--) {
                        if (cum_max[i] < pcoord - best_dist) break;
                        if (surf_max[i] < pcoord - best_dist) continue;

                        const I orig_idx = sort_idx[i];
                        const I i0 = s0[orig_idx], i1 = s1[orig_idx], i2 = s2[orig_idx];
                        const G tx0 = sx[i0], tx1 = sx[i1], tx2 = sx[i2];
                        const G ty0 = sy[i0], ty1 = sy[i1], ty2 = sy[i2];
                        const G tz0 = sz[i0], tz1 = sz[i1], tz2 = sz[i2];

                        const T dx = std::min({std::abs(px - tx0), std::abs(px - tx1), std::abs(px - tx2)});
                        const T dy = std::min({std::abs(py - ty0), std::abs(py - ty1), std::abs(py - ty2)});
                        const T dz = std::min({std::abs(pz - tz0), std::abs(pz - tz1), std::abs(pz - tz2)});
                        const T dist = std::sqrt(dx * dx + dy * dy + dz * dz);
                        best_dist = std::min(best_dist, dist);
                    }
                }

                // Check elements to the right (vectorized)
                {
                    ptrdiff_t i = left;
                    ptrdiff_t end = nselements;

                    // Process in chunks of VECTOR_SIZE for vectorization
                    while (i + VECTOR_SIZE <= end) {
                        // Early termination check
                        if (surf_min[i + VECTOR_SIZE - 1] > pcoord + best_dist) break;

                        // Process chunk of VECTOR_SIZE elements
                        T chunk_best =
                            compute_triangle_dist_chunk<G, T, I>(px, py, pz, i, sort_idx, s0, s1, s2, sx, sy, sz);
                        best_dist = std::min(best_dist, chunk_best);

                        i += VECTOR_SIZE;
                    }

                    // Handle remainder (scalar)
                    for (; i < end; i++) {
                        if (surf_min[i] > pcoord + best_dist) break;

                        const I orig_idx = sort_idx[i];
                        const I i0 = s0[orig_idx], i1 = s1[orig_idx], i2 = s2[orig_idx];
                        const G tx0 = sx[i0], tx1 = sx[i1], tx2 = sx[i2];
                        const G ty0 = sy[i0], ty1 = sy[i1], ty2 = sy[i2];
                        const G tz0 = sz[i0], tz1 = sz[i1], tz2 = sz[i2];

                        const T dx = std::min({std::abs(px - tx0), std::abs(px - tx1), std::abs(px - tx2)});
                        const T dy = std::min({std::abs(py - ty0), std::abs(py - ty1), std::abs(py - ty2)});
                        const T dz = std::min({std::abs(pz - tz0), std::abs(pz - tz1), std::abs(pz - tz2)});
                        const T dist = std::sqrt(dx * dx + dy * dy + dz * dz);
                        best_dist = std::min(best_dist, dist);
                    }
                }

                out[p] = best_dist;
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
