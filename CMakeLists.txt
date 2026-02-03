cmake_minimum_required(VERSION 3.15)
project(ssdf VERSION 0.1.0 LANGUAGES CXX)

# Standard install dirs (bin, lib, include, ...)
include(GNUInstallDirs)

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

# Set output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Build-time generated config header
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/ssdf_config.hpp.in"
    "${CMAKE_CURRENT_BINARY_DIR}/ssdf_config.hpp"
    @ONLY
)

# Optional dependencies
if(SSDF_ENABLE_MPI)
    find_package(MPI REQUIRED)
    message(STATUS "MPI found: ${MPI_FOUND}")
    message(STATUS "MPI CXX compiler: ${MPI_CXX_COMPILER}")
endif()

if(SSDF_ENABLE_OPENMP)
    find_package(OpenMP COMPONENTS CXX REQUIRED)
    message(STATUS "OpenMP found: ${OpenMP_FOUND}")
    message(STATUS "OpenMP CXX: ${OpenMP_CXX_VERSION}")
    message(STATUS "OpenMP CXX flags: ${OpenMP_CXX_FLAGS}")
endif()

if(SSDF_ENABLE_CUDA)
    enable_language(CUDA)
    message(STATUS "CUDA enabled")
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)

    if(NOT CMAKE_CUDA_ARCHITECTURES)
        set(CMAKE_CUDA_ARCHITECTURES "native")
    endif()
    
endif()

if(SSDF_ENABLE_HIP)
    enable_language(HIP)
    message(STATUS "HIP enabled")
endif()

if(SSDF_ENABLE_CUBIQL)
    if(NOT cuBQL_DIR)
        message(FATAL_ERROR "cuBQL_DIR is not set")
    endif()
    add_subdirectory(${cuBQL_DIR} "${CMAKE_BINARY_DIR}/_deps/cuBQL-build")
endif()

enable_testing()

add_subdirectory(serial)
add_subdirectory(distributed)
add_subdirectory(drivers)
add_subdirectory(tests)

# Install targets (library only; executables installed in drivers/)
install(TARGETS ssdf
    EXPORT ssdfTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# Install generated config header (flat include)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/ssdf_config.hpp"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Install CMake package files
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/ssdfConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/ssdfConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/ssdfConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ssdf
)

install(EXPORT ssdfTargets
    FILE ssdfTargets.cmake
    NAMESPACE ssdf::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ssdf
    EXPORT_LINK_INTERFACE_LIBRARIES
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/ssdfConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/ssdfConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ssdf
)



