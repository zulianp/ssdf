#include "bvh/bvh.hpp"

#if defined(__HIP_PLATFORM_SPIRV__)
// chipStar (HIP on SPIR-V, Intel GPUs) has no signed 64-bit atomicMin and atomicMax, which cuBQL's GPU
// builder names (builder/cuda/builder_common.h; a float BVH such as bvh3f uses the 32-bit ones only).
// CUDA's, as compare-and-swap loops on chipStar's 64-bit atomicCAS (relaxed, device scope, as CUDA's):
// these need 64-bit base atomics only, not the extended ones (cl_khr_int64_extended_atomics) a native
// 64-bit min/max lowers to. cuBQL calls them unqualified, so they stay global; as templates, a
// chipStar release that declares these overloads takes precedence instead of colliding with them.
template <int = 0>
__device__ inline long long atomicMin(long long *const address, const long long value) {
    unsigned long long *const bits = (unsigned long long *)address;
    unsigned long long old = *bits;
    unsigned long long assumed;
    do {
        assumed = old;
        if ((long long)assumed <= value) break;
        old = atomicCAS(bits, assumed, (unsigned long long)value);
    } while (old != assumed);
    return (long long)old;
}

template <int = 0>
__device__ inline long long atomicMax(long long *const address, const long long value) {
    unsigned long long *const bits = (unsigned long long *)address;
    unsigned long long old = *bits;
    unsigned long long assumed;
    do {
        assumed = old;
        if ((long long)assumed >= value) break;
        old = atomicCAS(bits, assumed, (unsigned long long)value);
    } while (old != assumed);
    return (long long)old;
}
#endif

#define CUBQL_GPU_BUILDER_IMPLEMENTATION 1
#define CUBQL_TRIANGLE_CPAT_IMPLEMENTATION 1

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
    cpat.sqrDist = out[tid]*out[tid];
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
        CUBQL_CUDA_CALL(Memcpy(d_out, out, npoints * sizeof(T), cudaMemcpyHostToDevice));

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

        CUBQL_CUDA_SYNC_CHECK();

        bvh3f d_bvh;
#if defined(__HIP_PLATFORM_SPIRV__)
        // chipStar (HIP on SPIR-V, Intel GPUs) implements no memory pools, which cuBQL's default
        // allocator sets up under HIP: the BVH is built and freed with plain device allocations.
        // Each allocation then waits for the stream first, and each free for the device. chipStar's
        // hipMallocAsync and hipFreeAsync are hipMalloc and that synchronous free too, so the stream
        // waits are all this costs against the pooled path.
        cuBQL::DeviceMemoryResource mem_resource;
        cuBQL::gpuBuilder(d_bvh, d_boxes, nselements, cuBQL::BuildConfig(), 0, mem_resource);
#else
        cuBQL::gpuBuilder(d_bvh, d_boxes, nselements, cuBQL::BuildConfig());
#endif

        edf_bvh_cuda_queries_kernel<<<divRoundUp(npoints, ptrdiff_t(256)), 256>>>(
            d_bvh, d_triangles, npoints, d_x, d_y, d_z, d_out);

        CUBQL_CUDA_SYNC_CHECK();

        CUBQL_CUDA_CALL(Memcpy(out, d_out, npoints * sizeof(T), cudaMemcpyDeviceToHost));

        CUBQL_CUDA_CALL(Free(d_triangles));
        CUBQL_CUDA_CALL(Free(d_boxes));
        CUBQL_CUDA_CALL(Free(d_x));
        CUBQL_CUDA_CALL(Free(d_y));
        CUBQL_CUDA_CALL(Free(d_z));
        CUBQL_CUDA_CALL(Free(d_out));
#if defined(__HIP_PLATFORM_SPIRV__)
        cuBQL::cuda::free(d_bvh, 0, mem_resource);
#else
        cuBQL::cuda::free(d_bvh);
#endif
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
