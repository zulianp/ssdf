cmake_minimum_required(VERSION 3.15)
project(ssdf VERSION 0.1.0 LANGUAGES CXX)

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

# Set output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Optional dependencies
if(SSDF_ENABLE_MPI)
    find_package(MPI REQUIRED)
    message(STATUS "MPI found: ${MPI_FOUND}")
    message(STATUS "MPI CXX compiler: ${MPI_CXX_COMPILER}")
endif()

if(SSDF_ENABLE_OPENMP)
    find_package(OpenMP REQUIRED)
    message(STATUS "OpenMP found: ${OpenMP_FOUND}")
    if(OpenMP_CXX_FOUND)
        message(STATUS "OpenMP CXX: ${OpenMP_CXX_VERSION}")
    endif()
endif()

if(SSDF_ENABLE_CUDA)
    enable_language(CUDA)
    message(STATUS "CUDA enabled")
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
endif()

if(SSDF_ENABLE_HIP)
    enable_language(HIP)
    message(STATUS "HIP enabled")
endif()

# Library sources
set(SSDF_SOURCES
    ssdf.cpp
)


# Collect header files (only include existing files)
set(SSDF_HEADERS "")
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/ssdf.hpp")
    list(APPEND SSDF_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/ssdf.hpp")
endif()
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/cell_list.hpp")
    list(APPEND SSDF_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/cell_list.hpp")
endif()

if(SSDF_ENABLE_MPI)
    list(APPEND SSDF_HEADERS psdf.hpp)
endif()


# Create library target (headers are for installation only, not compilation)
add_library(ssdf ${SSDF_SOURCES})

# Set target properties
set_target_properties(ssdf PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
)

# Link dependencies

# Add optional dependencies
if(SSDF_ENABLE_MPI)
    target_link_libraries(ssdf PUBLIC MPI::MPI_CXX)
    target_compile_definitions(ssdf PUBLIC SSDF_ENABLE_MPI)
endif()

if(SSDF_ENABLE_OPENMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(ssdf PUBLIC OpenMP::OpenMP_CXX)
        target_compile_definitions(ssdf PUBLIC SSDF_ENABLE_OPENMP)
    endif()
endif()

if(SSDF_ENABLE_CUDA)
    target_compile_definitions(ssdf PUBLIC SSDF_ENABLE_CUDA)
    set_property(TARGET ssdf PROPERTY CUDA_SEPARABLE_COMPILATION ON)
endif()

if(SSDF_ENABLE_HIP)
    target_compile_definitions(ssdf PUBLIC SSDF_ENABLE_HIP)
endif()

# Include directories for public interface
target_include_directories(ssdf
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
)

# Install targets
install(TARGETS ssdf
    EXPORT ssdfTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
)

# Install headers
install(FILES ${SSDF_HEADERS}
    DESTINATION include/ssdf
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
    INSTALL_DESTINATION lib/cmake/ssdf
)

install(EXPORT ssdfTargets
    FILE ssdfTargets.cmake
    NAMESPACE ssdf::
    DESTINATION lib/cmake/ssdf
    EXPORT_LINK_INTERFACE_LIBRARIES
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/ssdfConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/ssdfConfigVersion.cmake"
    DESTINATION lib/cmake/ssdf
)

# Testing
enable_testing()

# Test executable
add_executable(ssdf_test tests/ssdf_test.cpp)
target_link_libraries(ssdf_test PRIVATE ssdf)

# Add test to CTest
add_test(NAME ssdf_basic_test COMMAND ssdf_test)

# SDF executable
add_executable(sdf_exe sdf.exe.cpp)
target_link_libraries(sdf_exe PRIVATE ssdf)


if(SSDF_ENABLE_MPI)
    add_executable(psdf_exe psdf.exe.cpp)
    target_link_libraries(psdf_exe PRIVATE ssdf)
endif()