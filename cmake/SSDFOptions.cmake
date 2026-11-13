include_guard(GLOBAL)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Build options
option(BUILD_SHARED_LIBS "Build shared library" ON)
option(SSDF_ENABLE_MPI "Enable MPI support" OFF)
option(SSDF_ENABLE_OPENMP "Enable OpenMP support" OFF)
option(SSDF_ENABLE_CUDA "Enable CUDA support" OFF)
# SSDF_ENABLE_HIP=ON compiles the device sources with HIP instead of CUDA, for AMD GPUs. The device
# path is what SSDF_ENABLE_CUDA turns on, so HIP turns it on too. gfx942 is MI300A.
option(SSDF_ENABLE_HIP "Build the device path with HIP, for AMD GPUs" OFF)
set(SSDF_HIP_ARCHITECTURES "gfx942" CACHE STRING "HIP (AMD GPU) architectures to build for")
if(SSDF_ENABLE_HIP AND NOT SSDF_ENABLE_CUDA)
    message(STATUS "SSDF_ENABLE_HIP=ON turns on SSDF_ENABLE_CUDA (the device path)")
    set(SSDF_ENABLE_CUDA ON CACHE BOOL "Enable CUDA support" FORCE)
    # A parent project may hold it in a normal variable, which hides the cache entry.
    set(SSDF_ENABLE_CUDA ON)
endif()
option(SSDF_ENABLE_CUBIQL "Enable cuBQL support" ON)
option(SSDF_ENABLE_SMESH "Enable smesh support" OFF)
option(SSDF_ENABLE_ASAN "Enable AddressSanitizer instrumentation" OFF)
option(SSDF_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
