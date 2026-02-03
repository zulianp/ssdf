#include "bvh.hpp"

#include "cuBQL/queries/triangleData/closestPointOnAnyTriangle.h"
#include "cuBQL/queries/triangleData/pointInsideOutside.h"

#define MIN(a, b) ((a) < (b) ? (a) : (b))

namespace ssdf {

    template <typename T>
    static inline void triangle_local_coordinates(const T px,
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
                                                  const T cz,
                                                  T *const SSDF_RESTRICT w0,
                                                  T *const SSDF_RESTRICT w1,
                                                  T *const SSDF_RESTRICT w2) {
        const T v0x = bx - ax, v0y = by - ay, v0z = bz - az;
        const T v1x = cx - ax, v1y = cy - ay, v1z = cz - az;
        const T v2x = px - ax, v2y = py - ay, v2z = pz - az;

        const T d00 = v0x * v0x + v0y * v0y + v0z * v0z;
        const T d01 = v0x * v1x + v0y * v1y + v0z * v1z;
        const T d11 = v1x * v1x + v1y * v1y + v1z * v1z;
        const T d20 = v2x * v0x + v2y * v0y + v2z * v0z;
        const T d21 = v2x * v1x + v2y * v1y + v2z * v1z;

        const T inv_denom = T(1) / (d00 * d11 - d01 * d01);
        const T v = (d11 * d20 - d01 * d21) * inv_denom;
        const T w = (d00 * d21 - d01 * d20) * inv_denom;
        const T u = T(1) - v - w;

        w0[0] = u;
        w1[0] = v;
        w2[0] = w;
    }

    template <typename G, typename T, typename I>
    int interpolant_cpu(const ptrdiff_t npoints,
                        const G *const SSDF_RESTRICT x,
                        const G *const SSDF_RESTRICT y,
                        const G *const SSDF_RESTRICT z,
                        const I *const SSDF_RESTRICT triangles,
                        const ptrdiff_t nselements,
                        const I *const SSDF_RESTRICT s0,
                        const I *const SSDF_RESTRICT s1,
                        const I *const SSDF_RESTRICT s2,
                        const ptrdiff_t nspoints,
                        const G *const SSDF_RESTRICT sx,
                        const G *const SSDF_RESTRICT sy,
                        const G *const SSDF_RESTRICT sz,
                        T *const SSDF_RESTRICT w0,
                        T *const SSDF_RESTRICT w1,
                        T *const SSDF_RESTRICT w2) {
#pragma omp parallel for
        for (ptrdiff_t i = 0; i < npoints; i++) {
            const I tri = triangles[i];
            const I i0 = s0[tri];
            const I i1 = s1[tri];
            const I i2 = s2[tri];

            const G ax = sx[i0], ay = sy[i0], az = sz[i0];
            const G bx = sx[i1], by = sy[i1], bz = sz[i1];
            const G cx = sx[i2], cy = sy[i2], cz = sz[i2];
            triangle_local_coordinates<T>(x[i], y[i], z[i], ax, ay, az, bx, by, bz, cx, cy, cz, &w0[i], &w1[i], &w2[i]);
        }

        return 0;
    }

    template <typename G, typename T, typename I>
    int closest_point_bvh_cpu(const ptrdiff_t npoints,
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
                              I *const SSDF_RESTRICT outtri) {
        using cuBQL::box3f;
        using cuBQL::bvh3f;
        using cuBQL::divRoundUp;
        using cuBQL::getCurrentTime;
        using cuBQL::prettyDouble;
        using cuBQL::prettyNumber;
        using cuBQL::Triangle;
        using cuBQL::vec3f;
        using bvh_t = bvh3f;

        std::vector<Triangle> triangles(nselements);
        std::vector<box3f> boxes(nselements);
        for (ptrdiff_t i = 0; i < nselements; i++) {
            triangles[i] = Triangle(vec3f(sx[s0[i]], sy[s0[i]], sz[s0[i]]),
                                    vec3f(sx[s1[i]], sy[s1[i]], sz[s1[i]]),
                                    vec3f(sx[s2[i]], sy[s2[i]], sz[s2[i]]));
            boxes[i] = triangles[i].bounds();
        }

        bvh_t trianglesBVH;
        cuBQL::cpuBuilder(trianglesBVH, boxes.data(), boxes.size(), cuBQL::BuildConfig());

        auto runQueries = [&](bvh_t trianglesBVH,
                              const Triangle *triangles,
                              int numQueries,
                              const G *const SSDF_RESTRICT x,
                              const G *const SSDF_RESTRICT y,
                              const G *const SSDF_RESTRICT z,
                              T *const SSDF_RESTRICT outd,
                              G *const SSDF_RESTRICT outx,
                              G *const SSDF_RESTRICT outy,
                              G *const SSDF_RESTRICT outz,
                              I *const SSDF_RESTRICT outtri) {
#pragma omp parallel for
            for (int tid = 0; tid < numQueries; tid++) {
                vec3f queryPoint(x[tid], y[tid], z[tid]);
                cuBQL::triangles::CPAT cpat;
                cpat.sqrDist = outd[tid] * outd[tid];
                cpat.runQuery(triangles, trianglesBVH, queryPoint);
                const T dist = sqrt(cpat.sqrDist);
                if (outd[tid] > dist) {
                    outd[tid] = dist;
                    outx[tid] = cpat.P.x;
                    outy[tid] = cpat.P.y;
                    outz[tid] = cpat.P.z;
                    outtri[tid] = cpat.triangleIdx;
                }
            }
        };

        runQueries(trianglesBVH, triangles.data(), npoints, x, y, z, outd, outx, outy, outz, outtri);
        return 0;
    }

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
                          I *const SSDF_RESTRICT outtri) {
        return closest_point_bvh_cpu<G, T, I>(
            npoints, x, y, z, nselements, s0, s1, s2, nspoints, sx, sy, sz, outd, outx, outy, outz, outtri);
    }

    template <typename G, typename T, typename I>
    int edf_bvh_cpu(const ptrdiff_t npoints,
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
        using cuBQL::box3f;
        using cuBQL::bvh3f;
        using cuBQL::divRoundUp;
        using cuBQL::getCurrentTime;
        using cuBQL::prettyDouble;
        using cuBQL::prettyNumber;
        using cuBQL::Triangle;
        using cuBQL::vec3f;
        using bvh_t = bvh3f;

        std::vector<Triangle> triangles(nselements);
        std::vector<box3f> boxes(nselements);
        for (ptrdiff_t i = 0; i < nselements; i++) {
            triangles[i] = Triangle(vec3f(sx[s0[i]], sy[s0[i]], sz[s0[i]]),
                                    vec3f(sx[s1[i]], sy[s1[i]], sz[s1[i]]),
                                    vec3f(sx[s2[i]], sy[s2[i]], sz[s2[i]]));
            boxes[i] = triangles[i].bounds();
        }

        bvh_t trianglesBVH;
        cuBQL::cpuBuilder(trianglesBVH, boxes.data(), boxes.size(), cuBQL::BuildConfig());

        auto runQueries = [&](bvh_t trianglesBVH,
                              const Triangle *triangles,
                              int numQueries,
                              const G *const SSDF_RESTRICT x,
                              const G *const SSDF_RESTRICT y,
                              const G *const SSDF_RESTRICT z,
                              T *const SSDF_RESTRICT out) {
#pragma omp parallel for
            for (int tid = 0; tid < numQueries; tid++) {
                vec3f queryPoint(x[tid], y[tid], z[tid]);
                cuBQL::triangles::CPAT cpat;
                cpat.sqrDist = out[tid] * out[tid];
                cpat.runQuery(triangles, trianglesBVH, queryPoint);
                out[tid] = MIN(out[tid], sqrt(cpat.sqrDist));
            }
        };

        runQueries(trianglesBVH, triangles.data(), npoints, x, y, z, out);
        return 0;
    }

#ifdef SSDF_ENABLE_CUDA

    template <typename G, typename T, typename I>
    int edf_bvh_cuda(const ptrdiff_t npoints,
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
#endif

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
                T *const SSDF_RESTRICT out) {
#ifdef SSDF_ENABLE_CUDA
        return edf_bvh_cuda<G, T, I>(npoints, x, y, z, nselements, s0, s1, s2, nspoints, sx, sy, sz, out);
#endif
        return edf_bvh_cpu<G, T, I>(npoints, x, y, z, nselements, s0, s1, s2, nspoints, sx, sy, sz, out);
    }

    template <typename G, typename I>
    int points_inside_bvh_cpu(const ptrdiff_t npoints,
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
                              uint8_t *const SSDF_RESTRICT out) {
        using cuBQL::box3f;
        using cuBQL::bvh3f;
        using cuBQL::Triangle;
        using cuBQL::vec3f;
        using bvh_t = bvh3f;

        std::vector<box3f> boxes(nselements);
        for (ptrdiff_t i = 0; i < nselements; i++) {
            auto tri = Triangle(vec3f(sx[s0[i]], sy[s0[i]], sz[s0[i]]),
                                vec3f(sx[s1[i]], sy[s1[i]], sz[s1[i]]),
                                vec3f(sx[s2[i]], sy[s2[i]], sz[s2[i]]));
            boxes[i] = tri.bounds();
        }

        bvh_t trianglesBVH;
        cuBQL::cpuBuilder(trianglesBVH, boxes.data(), boxes.size(), cuBQL::BuildConfig());

        auto getTriangle = [s0, s1, s2, sx, sy, sz](uint32_t primID) {
            return Triangle(vec3f(sx[s0[primID]], sy[s0[primID]], sz[s0[primID]]),
                            vec3f(sx[s1[primID]], sy[s1[primID]], sz[s1[primID]]),
                            vec3f(sx[s2[primID]], sy[s2[primID]], sz[s2[primID]]));
        };

#pragma omp parallel for
        for (int tid = 0; tid < npoints; tid++) {
            vec3f queryPoint(x[tid], y[tid], z[tid]);
            bool inside = cuBQL::triangles::pointIsInsideSurface(trianglesBVH, getTriangle, queryPoint);
            out[tid] = inside;
        }

        return 0;
    }

    template int closest_point_bvh<float, double, int>(const ptrdiff_t,
                                                       const float *const SSDF_RESTRICT,
                                                       const float *const SSDF_RESTRICT,
                                                       const float *const SSDF_RESTRICT,
                                                       const ptrdiff_t,
                                                       const int *const SSDF_RESTRICT,
                                                       const int *const SSDF_RESTRICT,
                                                       const int *const SSDF_RESTRICT,
                                                       const ptrdiff_t,
                                                       const float *const SSDF_RESTRICT,
                                                       const float *const SSDF_RESTRICT,
                                                       const float *const SSDF_RESTRICT,
                                                       double *const SSDF_RESTRICT,
                                                       float *const SSDF_RESTRICT,
                                                       float *const SSDF_RESTRICT,
                                                       float *const SSDF_RESTRICT,
                                                       int *const SSDF_RESTRICT);

    template int closest_point_bvh<double, double, int>(const ptrdiff_t,
                                                        const double *const SSDF_RESTRICT,
                                                        const double *const SSDF_RESTRICT,
                                                        const double *const SSDF_RESTRICT,
                                                        const ptrdiff_t,
                                                        const int *const SSDF_RESTRICT,
                                                        const int *const SSDF_RESTRICT,
                                                        const int *const SSDF_RESTRICT,
                                                        const ptrdiff_t,
                                                        const double *const SSDF_RESTRICT,
                                                        const double *const SSDF_RESTRICT,
                                                        const double *const SSDF_RESTRICT,
                                                        double *const SSDF_RESTRICT,
                                                        double *const SSDF_RESTRICT,
                                                        double *const SSDF_RESTRICT,
                                                        double *const SSDF_RESTRICT,
                                                        int *const SSDF_RESTRICT);

    template int closest_point_bvh<float, float, int>(const ptrdiff_t,
                                                      const float *const SSDF_RESTRICT,
                                                      const float *const SSDF_RESTRICT,
                                                      const float *const SSDF_RESTRICT,
                                                      const ptrdiff_t,
                                                      const int *const SSDF_RESTRICT,
                                                      const int *const SSDF_RESTRICT,
                                                      const int *const SSDF_RESTRICT,
                                                      const ptrdiff_t,
                                                      const float *const SSDF_RESTRICT,
                                                      const float *const SSDF_RESTRICT,
                                                      const float *const SSDF_RESTRICT,
                                                      float *const SSDF_RESTRICT,
                                                      float *const SSDF_RESTRICT,
                                                      float *const SSDF_RESTRICT,
                                                      float *const SSDF_RESTRICT,
                                                      int *const SSDF_RESTRICT);

    template int edf_bvh<float, double, int>(const ptrdiff_t,
                                             const float *const SSDF_RESTRICT,
                                             const float *const SSDF_RESTRICT,
                                             const float *const SSDF_RESTRICT,
                                             const ptrdiff_t,
                                             const int *const SSDF_RESTRICT,
                                             const int *const SSDF_RESTRICT,
                                             const int *const SSDF_RESTRICT,
                                             const ptrdiff_t,
                                             const float *const SSDF_RESTRICT,
                                             const float *const SSDF_RESTRICT,
                                             const float *const SSDF_RESTRICT,
                                             double *const SSDF_RESTRICT);

    template int edf_bvh<float, float, int>(const ptrdiff_t,
                                            const float *const SSDF_RESTRICT,
                                            const float *const SSDF_RESTRICT,
                                            const float *const SSDF_RESTRICT,
                                            const ptrdiff_t,
                                            const int *const SSDF_RESTRICT,
                                            const int *const SSDF_RESTRICT,
                                            const int *const SSDF_RESTRICT,
                                            const ptrdiff_t,
                                            const float *const SSDF_RESTRICT,
                                            const float *const SSDF_RESTRICT,
                                            const float *const SSDF_RESTRICT,
                                            float *const SSDF_RESTRICT);

    template int edf_bvh<double, double, int>(const ptrdiff_t,
                                              const double *const SSDF_RESTRICT,
                                              const double *const SSDF_RESTRICT,
                                              const double *const SSDF_RESTRICT,
                                              const ptrdiff_t,
                                              const int *const SSDF_RESTRICT,
                                              const int *const SSDF_RESTRICT,
                                              const int *const SSDF_RESTRICT,
                                              const ptrdiff_t,
                                              const double *const SSDF_RESTRICT,
                                              const double *const SSDF_RESTRICT,
                                              const double *const SSDF_RESTRICT,
                                              double *const SSDF_RESTRICT);

    template int points_inside_bvh_cpu<float, int>(const ptrdiff_t npoints,
                                                   const float *const SSDF_RESTRICT x,
                                                   const float *const SSDF_RESTRICT y,
                                                   const float *const SSDF_RESTRICT z,
                                                   const ptrdiff_t nselements,
                                                   const int *const SSDF_RESTRICT s0,
                                                   const int *const SSDF_RESTRICT s1,
                                                   const int *const SSDF_RESTRICT s2,
                                                   const ptrdiff_t nspoints,
                                                   const float *const SSDF_RESTRICT sx,
                                                   const float *const SSDF_RESTRICT sy,
                                                   const float *const SSDF_RESTRICT sz,
                                                   uint8_t *const SSDF_RESTRICT out);
}  // namespace ssdf
