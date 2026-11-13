#ifndef SSDF_GPU_RUNTIME_CUH
#define SSDF_GPU_RUNTIME_CUH

// The GPU runtime the device BVH is built against: CUDA, or HIP for AMD GPUs.
//
// bvh.cu and the cuBQL GPU builder it instantiates are written against the CUDA runtime API, which
// HIP mirrors name for name. Under SSDF_ENABLE_HIP this header includes HIP and maps the cuda*
// names they call onto their hip* counterparts; cuBQL's CUBQL_CUDA_CALL(X) pastes cuda##X, and the
// pasted name is mapped too. It must come before cuBQL's headers. Otherwise it includes the CUDA
// runtime.
//
// smesh, SCCD and sresample have headers like this one. A macro may be redefined only with the
// same replacement, so a name mapped in more than one header must map to the same tokens in each.

#include "ssdf_config.hpp"

#if defined(SSDF_ENABLE_HIP)

#include <hip/hip_runtime.h>

// The HIP runtime follows the CUDA 12 API: cuBQL then allocates with hipMallocAsync, where for older
// toolkits it takes managed memory, which an AMD GPU without XNACK faults on.
#ifndef CUDART_VERSION
#define CUDART_VERSION 12000
#endif

#define cudaError_t hipError_t
#define cudaSuccess hipSuccess
#define cudaGetErrorString hipGetErrorString
#define cudaGetLastError hipGetLastError

#define cudaMalloc hipMalloc
#define cudaFree hipFree
#define cudaMallocAsync hipMallocAsync
#define cudaFreeAsync hipFreeAsync
#define cudaMallocHost hipHostMalloc
#define cudaFreeHost hipHostFree
#define cudaMallocManaged hipMallocManaged
#define cudaMemPool_t hipMemPool_t
#define cudaMemPoolAttrReleaseThreshold hipMemPoolAttrReleaseThreshold
#define cudaMemPoolSetAttribute hipMemPoolSetAttribute
#define cudaDeviceGetDefaultMemPool hipDeviceGetDefaultMemPool

#define cudaMemcpy hipMemcpy
#define cudaMemcpyAsync hipMemcpyAsync
#define cudaMemcpyHostToDevice hipMemcpyHostToDevice
#define cudaMemcpyDeviceToHost hipMemcpyDeviceToHost
#define cudaMemcpyDeviceToDevice hipMemcpyDeviceToDevice
#define cudaMemcpyDefault hipMemcpyDefault
#define cudaMemset hipMemset

#define cudaStream_t hipStream_t
#define cudaStreamSynchronize hipStreamSynchronize
#define cudaDeviceSynchronize hipDeviceSynchronize
#define cudaEvent_t hipEvent_t
#define cudaEventCreate hipEventCreate
#define cudaEventRecord hipEventRecord
#define cudaEventSynchronize hipEventSynchronize
#define cudaEventDestroy hipEventDestroy

#define cudaGetDevice hipGetDevice
#define cudaSetDevice hipSetDevice
#define cudaGetDeviceCount hipGetDeviceCount

#else

#include <cuda_runtime.h>

#endif

#endif  // SSDF_GPU_RUNTIME_CUH
