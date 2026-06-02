#ifndef SSDF_BVH_HPP
#define SSDF_BVH_HPP

#include <stddef.h>
#include <stdint.h>

#include "ssdf_config.hpp"

#ifndef SSDF_RESTRICT
#ifndef _WIN32
#define SSDF_RESTRICT __restrict__
#else
#define SSDF_RESTRICT __restrict
#endif
#endif

namespace ssdf {

    template <typename G, typename T, typename I>
    int edf_bvh(const ptrdiff_t npoints,
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
                T *const SSDF_RESTRICT out);

    template <typename G, typename T, typename I>
    int closest_point_bvh(const ptrdiff_t npoints,
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
                          T *const SSDF_RESTRICT outd,
                          G *const SSDF_RESTRICT outx,
                          G *const SSDF_RESTRICT outy,
                          G *const SSDF_RESTRICT outz,
                          I *const SSDF_RESTRICT outtri);

    template <typename G, typename I>
    int points_inside_bvh(const ptrdiff_t npoints,
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
                          uint8_t *const SSDF_RESTRICT out);

    template <typename G, typename T, typename I>
    int closest_within_radius_bvh(const ptrdiff_t npoints,
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
                                  const ptrdiff_t radius_stride,
                                  const T *const SSDF_RESTRICT radius_squared,
                                  I *const SSDF_RESTRICT outtri,
                                  T *const SSDF_RESTRICT out_sqr_dist,
                                  T *const SSDF_RESTRICT closest_x,
                                  T *const SSDF_RESTRICT closest_y,
                                  T *const SSDF_RESTRICT closest_z,
                                  const bool skip_self = false);

    template <typename G, typename T, typename I, typename F>
    int potential_contact_faces_bvh(const ptrdiff_t nselements,
                                    const I *const SSDF_RESTRICT s0,
                                    const I *const SSDF_RESTRICT s1,
                                    const I *const SSDF_RESTRICT s2,
                                    const ptrdiff_t nspoints,
                                    const G *const SSDF_RESTRICT sx,
                                    const G *const SSDF_RESTRICT sy,
                                    const G *const SSDF_RESTRICT sz,
                                    const T extrusion,
                                    ptrdiff_t *const SSDF_RESTRICT pc_ptr,
                                    F **const SSDF_RESTRICT out_pc_idx);
}  // namespace ssdf

#endif  // SSDF_BVH_HPP
