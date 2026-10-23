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

    template <typename G, typename T, typename I>
    int closest_within_radius_local_bvh(const ptrdiff_t npoints,
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
                                        T *const SSDF_RESTRICT closest_s,
                                        T *const SSDF_RESTRICT closest_t,
                                        const bool skip_self = false);

    template <typename G, typename T, typename I>
    int closest_within_radius_quads_bvh(const ptrdiff_t npoints,
                                        const G *const SSDF_RESTRICT x,
                                        const G *const SSDF_RESTRICT y,
                                        const G *const SSDF_RESTRICT z,
                                        const ptrdiff_t nselements,
                                        const I *const SSDF_RESTRICT s0,
                                        const I *const SSDF_RESTRICT s1,
                                        const I *const SSDF_RESTRICT s2,
                                        const I *const SSDF_RESTRICT s3,
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

    template <typename G, typename T, typename I>
    int closest_within_radius_quads_local_bvh(const ptrdiff_t npoints,
                                              const G *const SSDF_RESTRICT x,
                                              const G *const SSDF_RESTRICT y,
                                              const G *const SSDF_RESTRICT z,
                                              const ptrdiff_t nselements,
                                              const I *const SSDF_RESTRICT s0,
                                              const I *const SSDF_RESTRICT s1,
                                              const I *const SSDF_RESTRICT s2,
                                              const I *const SSDF_RESTRICT s3,
                                              const ptrdiff_t nspoints,
                                              const G *const SSDF_RESTRICT sx,
                                              const G *const SSDF_RESTRICT sy,
                                              const G *const SSDF_RESTRICT sz,
                                              const ptrdiff_t radius_stride,
                                              const T *const SSDF_RESTRICT radius_squared,
                                              I *const SSDF_RESTRICT outtri,
                                              T *const SSDF_RESTRICT out_sqr_dist,
                                              T *const SSDF_RESTRICT closest_s,
                                              T *const SSDF_RESTRICT closest_t,
                                              const bool skip_self = false);

    template <typename G, typename T, typename I, typename F>
    int potential_contact_triangles_bvh(const ptrdiff_t nselements,
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

    template <typename G, typename T, typename I, typename F>
    int potential_contact_quads_bvh(const ptrdiff_t nselements,
                                    const I *const SSDF_RESTRICT s0,
                                    const I *const SSDF_RESTRICT s1,
                                    const I *const SSDF_RESTRICT s2,
                                    const I *const SSDF_RESTRICT s3,
                                    const ptrdiff_t nspoints,
                                    const G *const SSDF_RESTRICT sx,
                                    const G *const SSDF_RESTRICT sy,
                                    const G *const SSDF_RESTRICT sz,
                                    const T extrusion,
                                    ptrdiff_t *const SSDF_RESTRICT pc_ptr,
                                    F **const SSDF_RESTRICT out_pc_idx);

    template <typename G, typename T, typename I>
    int sample_sqedf_bvh(const ptrdiff_t nselements,
                         const I *const SSDF_RESTRICT s0,
                         const I *const SSDF_RESTRICT s1,
                         const I *const SSDF_RESTRICT s2,
                         const ptrdiff_t nspoints,
                         const G *const SSDF_RESTRICT sx,
                         const G *const SSDF_RESTRICT sy,
                         const G *const SSDF_RESTRICT sz,
                         const ptrdiff_t nx,
                         const ptrdiff_t ny,
                         const ptrdiff_t nz,
                         const G xmin,
                         const G ymin,
                         const G zmin,
                         const G xmax,
                         const G ymax,
                         const G zmax,
                         T *const SSDF_RESTRICT out);

    /// The signed distance to a triangulated surface on the grid's nodes: positive on the side its
    /// triangles' normals (right-handed in s0, s1, s2) point to. The sign is that of (p - c) . n at the
    /// closest point c, n the angle-weighted pseudonormal there (Baerentzen and Aanaes, IEEE TVCG
    /// 11(3), 2005): the face's normal inside a face, the sum of the two faces' on an edge, the
    /// angle-weighted sum of the incident faces' at a vertex. Exact for a closed, consistently
    /// oriented surface, wherever p is; a ray-crossing count is not, through edges and vertices.
    template <typename G, typename T, typename I>
    int sample_sdf_bvh(const ptrdiff_t nselements,
                       const I *const SSDF_RESTRICT s0,
                       const I *const SSDF_RESTRICT s1,
                       const I *const SSDF_RESTRICT s2,
                       const ptrdiff_t nspoints,
                       const G *const SSDF_RESTRICT sx,
                       const G *const SSDF_RESTRICT sy,
                       const G *const SSDF_RESTRICT sz,
                       const ptrdiff_t nx,
                       const ptrdiff_t ny,
                       const ptrdiff_t nz,
                       const G xmin,
                       const G ymin,
                       const G zmin,
                       const G xmax,
                       const G ymax,
                       const G zmax,
                       T *const SSDF_RESTRICT out);

    template <typename G, typename I>
    int sample_inside_bvh(const ptrdiff_t nselements,
                          const I *const SSDF_RESTRICT s0,
                          const I *const SSDF_RESTRICT s1,
                          const I *const SSDF_RESTRICT s2,
                          const ptrdiff_t nspoints,
                          const G *const SSDF_RESTRICT sx,
                          const G *const SSDF_RESTRICT sy,
                          const G *const SSDF_RESTRICT sz,
                          const ptrdiff_t nx,
                          const ptrdiff_t ny,
                          const ptrdiff_t nz,
                          const G xmin,
                          const G ymin,
                          const G zmin,
                          const G xmax,
                          const G ymax,
                          const G zmax,
                          uint8_t *const SSDF_RESTRICT out);
}  // namespace ssdf

#endif  // SSDF_BVH_HPP
