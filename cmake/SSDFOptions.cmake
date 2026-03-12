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
option(SSDF_ENABLE_HIP "Enable HIP support" OFF)
option(SSDF_ENABLE_CUBIQL "Enable cuBQL support" OFF)
option(SSDF_ENABLE_SMESH "Enable smesh support" ON)
option(SSDF_ENABLE_ASAN "Enable AddressSanitizer instrumentation" OFF)
option(SSDF_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
