#include "bvh.hpp"

#include "cuBQL/queries/triangleData/closestPointOnAnyTriangle.h"

#define MIN(a, b) ((a) < (b) ? (a) : (b))

namespace ssdf {

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
        using cuBQL::vec3i;

        std::vector<Triangle> triangles(nselements);
        std::vector<box3f> boxes(nselements);
        for (ptrdiff_t i = 0; i < nselements; i++) {
            triangles[i] = Triangle(vec3f(sx[s0[i]], sy[s0[i]], sz[s0[i]]),
                                    vec3f(sx[s1[i]], sy[s1[i]], sz[s1[i]]),
                                    vec3f(sx[s2[i]], sy[s2[i]], sz[s2[i]]));
            boxes[i] = triangles[i].bounds();
        }

        bvh3f trianglesBVH;
        cuBQL::cpuBuilder(trianglesBVH, boxes.data(), boxes.size(), cuBQL::BuildConfig());

        auto runQueries = [&](bvh3f trianglesBVH,
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

    template int edf_bvh<float, double, int>(const ptrdiff_t npoints,
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

    template int edf_bvh<float, float, int>(const ptrdiff_t npoints,
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

    template int edf_bvh<double, double, int>(const ptrdiff_t npoints,
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
}  // namespace ssdf
