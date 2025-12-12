#ifndef SSDF_HPP
#define SSDF_HPP

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <stddef.h>
#include <string>

#define SSDF_RESTRICT __restrict

#define VECTOR_SIZE 16

namespace ssdf {

// Return current time in milliseconds
inline double time_ms() {
  auto now = std::chrono::high_resolution_clock::now();
  auto duration = now.time_since_epoch();
  return std::chrono::duration<double, std::milli>(duration).count();
}

// Helper function to compute point-to-triangle distance for a fixed-size chunk
// Uses fixed-size loop to enable compiler vectorization
template <typename G, typename T, typename I>
inline T compute_triangle_dist_chunk(const G px, const G py, const G pz,
                                     const ptrdiff_t chunk_start,
                                     const I *const SSDF_RESTRICT sort_idx,
                                     const I *const SSDF_RESTRICT s0,
                                     const I *const SSDF_RESTRICT s1,
                                     const I *const SSDF_RESTRICT s2,
                                     const G *const SSDF_RESTRICT sx,
                                     const G *const SSDF_RESTRICT sy,
                                     const G *const SSDF_RESTRICT sz) {
  T best_dist = std::numeric_limits<T>::max();

  // Fixed-size loop (VECTOR_SIZE) for compiler vectorization
  #pragma omp simd reduction(min:best_dist)
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
    const T dist = std::sqrt(dx * dx + dy * dy + dz * dz);
    best_dist = std::min(best_dist, dist);
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
  void print() const {
    std::cout << name << " took " << duration << " ms" << std::endl;
  }
};

// out is of size npoints
template <typename G, typename T, typename I>
int sdf(const ptrdiff_t npoints, const G *const SSDF_RESTRICT x,
        const G *const SSDF_RESTRICT y, const G *const SSDF_RESTRICT z,
        const ptrdiff_t nselements, const I *const SSDF_RESTRICT s0,
        const I *const SSDF_RESTRICT s1, const I *const SSDF_RESTRICT s2,
        const ptrdiff_t nspoints, const G *const SSDF_RESTRICT sx,
        const G *const SSDF_RESTRICT sy, const G *const SSDF_RESTRICT sz,
        T *const SSDF_RESTRICT out) {
  // TODO: Implement this
  // 1) Compute surf aabbs
  // 2) Find closets distance on the x dimension
  //  - Sort surface w.r.t x
  //  - compute cum max
  //  - for each point binary search cum max to find closest element from the
  //  left use computed distance to conservatively discard candidates
  //  - Then check elements on the right  use computed distance to
  //  conservatively discard candidates
  // 3) Repeat for y and z

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
      std::sort(sort_idx, sort_idx + nselements,
                [surf_min](I a, I b) { return surf_min[a] < surf_min[b]; });

      // Permute arrays
      memcpy(scratch, surf_min, sizeof(T) * nselements);
      for (ptrdiff_t i = 0; i < nselements; i++) {
        surf_min[i] = scratch[sort_idx[i]];
      }
      memcpy(scratch, surf_max, sizeof(T) * nselements);
      for (ptrdiff_t i = 0; i < nselements; i++) {
        surf_max[i] = scratch[sort_idx[i]];
      }

      // Compute cumulative max
      cum_max[0] = surf_max[0];
      for (ptrdiff_t i = 1; i < nselements; i++) {
        cum_max[i] = std::max(cum_max[i - 1], surf_max[i]);
      }
    }

// For each point, find closest surface element
#pragma omp parallel for
    for (ptrdiff_t p = 0; p < npoints; p++) {
      const G px = x[p], py = y[p], pz = z[p];
      const T pcoord = pnt_coord[p];
      T best_dist = out[p];

      // Binary search for insertion point
      ptrdiff_t left =
          std::lower_bound(surf_min, surf_min + nselements, pcoord) - surf_min;

      // Check elements to the left (vectorized)
      {
        ptrdiff_t start = (left > 0 ? left - 1 : 0);
        ptrdiff_t i = start;
        
        // Process in chunks of VECTOR_SIZE for vectorization
        while (i >= VECTOR_SIZE - 1) {
          // Early termination check for the chunk
          if (cum_max[i - (VECTOR_SIZE - 1)] < pcoord - best_dist)
            break;
          
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
          if (cum_max[i] < pcoord - best_dist)
            break;
          if (surf_max[i] < pcoord - best_dist)
            continue;

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
          if (surf_min[i + VECTOR_SIZE - 1] > pcoord + best_dist)
            break;
          
          // Process chunk of VECTOR_SIZE elements
          T chunk_best = compute_triangle_dist_chunk<G, T, I>(
              px, py, pz, i, sort_idx, s0, s1, s2, sx, sy, sz);
          best_dist = std::min(best_dist, chunk_best);
          
          i += VECTOR_SIZE;
        }
        
        // Handle remainder (scalar)
        for (; i < end; i++) {
          if (surf_min[i] > pcoord + best_dist)
            break;

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

} // namespace ssdf

#endif // SSDF_HPP
