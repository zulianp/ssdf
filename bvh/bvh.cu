#include "bvh.hpp"

#include "cuBQL/bvh.h"
#include "cuBQL/queries/triangleData/closestPointOnAnyTriangle.h"

#define MIN(a, b) ((a) < (b) ? (a) : (b))

using cuBQL::box3f;
using cuBQL::bvh3f;
using cuBQL::divRoundUp;
using cuBQL::Triangle;
using cuBQL::vec3f;

/*! generate boxes (required for bvh builder) from prim type 'index line triangles' */
__global__ void edf_bvh_cuda_generate_boxes_kernel(box3f *boxForBuilder, const Triangle *triangles, int numTriangles) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= numTriangles) return;

    auto triangle = triangles[tid];
    boxForBuilder[tid] = triangle.bounds();
}

template <typename G, typename T>
__global__ void edf_bvh_cuda_queries_kernel(bvh3f trianglesBVH,
                                            const Triangle *triangles,
                                            int numQueries,
                                            const G *const SSDF_RESTRICT x,
                                            const G *const SSDF_RESTRICT y,
                                            const G *const SSDF_RESTRICT z,
                                            T *const SSDF_RESTRICT out) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= numQueries) return;

    vec3f queryPoint(x[tid], y[tid], z[tid]);
    cuBQL::triangles::CPAT cpat;
    cpat.runQuery(triangles, trianglesBVH, queryPoint);
    out[tid] = MIN(out[tid], sqrt(cpat.sqrDist));
}

namespace ssdf {

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
                     T *const SSDF_RESTRICT out) {
        G *d_x = 0;
        CUBQL_CUDA_CALL(Malloc((void **)&d_x, npoints * sizeof(G)));
        CUBQL_CUDA_CALL(Memcpy(d_x, x, npoints * sizeof(G), cudaMemcpyHostToDevice));
        G *d_y = 0;
        CUBQL_CUDA_CALL(Malloc((void **)&d_y, npoints * sizeof(G)));
        CUBQL_CUDA_CALL(Memcpy(d_y, y, npoints * sizeof(G), cudaMemcpyHostToDevice));
        G *d_z = 0;
        CUBQL_CUDA_CALL(Malloc((void **)&d_z, npoints * sizeof(G)));
        CUBQL_CUDA_CALL(Memcpy(d_z, z, npoints * sizeof(G), cudaMemcpyHostToDevice));

        T *d_out = 0;
        CUBQL_CUDA_CALL(Malloc((void **)&d_out, npoints * sizeof(T)));

        std::vector<Triangle> h_triangles(nselements);
        for (ptrdiff_t i = 0; i < nselements; i++) {
            h_triangles[i] = Triangle(vec3f(sx[s0[i]], sy[s0[i]], sz[s0[i]]),
                                      vec3f(sx[s1[i]], sy[s1[i]], sz[s1[i]]),
                                      vec3f(sx[s2[i]], sy[s2[i]], sz[s2[i]]));
        }

        Triangle *d_triangles = 0;
        CUBQL_CUDA_CALL(Malloc((void **)&d_triangles, nselements * sizeof(Triangle)));
        CUBQL_CUDA_CALL(Memcpy(d_triangles, h_triangles.data(), nselements * sizeof(Triangle), cudaMemcpyHostToDevice));

        box3f *d_boxes = 0;
        CUBQL_CUDA_CALL(Malloc((void **)&d_boxes, nselements * sizeof(box3f)));
        edf_bvh_cuda_generate_boxes_kernel<<<divRoundUp(nselements, ptrdiff_t(256)), 256>>>(
            d_boxes, d_triangles, nselements);
        CUBQL_CUDA_CALL(Free(d_boxes));

        bvh3f d_bvh;
        cuBQL::gpuBuilder(d_bvh, d_boxes, nselements, cuBQL::BuildConfig());

        edf_bvh_cuda_queries_kernel<<<divRoundUp(npoints, ptrdiff_t(256)), 256>>>(
            d_bvh, d_triangles, npoints, d_x, d_y, d_z, d_out);

        CUBQL_CUDA_CALL(Memcpy(out, d_out, npoints * sizeof(T), cudaMemcpyDeviceToHost));

        CUBQL_CUDA_CALL(Free(d_triangles));
        CUBQL_CUDA_CALL(Free(d_boxes));
        CUBQL_CUDA_CALL(Free(d_x));
        CUBQL_CUDA_CALL(Free(d_y));
        CUBQL_CUDA_CALL(Free(d_z));
        CUBQL_CUDA_CALL(Free(d_out));
        cuBQL::cuda::free(d_bvh);
        return 0;
    }

    template int edf_bvh_cuda<float, double, int>(const ptrdiff_t npoints,
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

    template int edf_bvh_cuda<float, float, int>(const ptrdiff_t npoints,
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

    template int edf_bvh_cuda<double, double, int>(const ptrdiff_t npoints,
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
