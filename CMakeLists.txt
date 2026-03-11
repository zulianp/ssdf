cmake_minimum_required(VERSION 3.15)
project(ssdf VERSION 0.1.0 LANGUAGES CXX)

# TODO: fix this cmake to be high-quality

# Standard install dirs (bin, lib, include, ...)
include(GNUInstallDirs)
include(cmake/SSDFDependencies.cmake)

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

# ------------------------------------------------------------------------------
# File collection helpers
# ------------------------------------------------------------------------------

# Folders are provided manually; files are collected via glob automatically.
function(ssdf_glob_into out_var)
    set(result "")
    foreach(dir IN LISTS ARGN)
        if(NOT IS_DIRECTORY "${dir}")
            message(FATAL_ERROR "ssdf_glob_into(): directory does not exist: ${dir}")
        endif()
        file(GLOB _cpp CONFIGURE_DEPENDS
            "${dir}/*.c"
            "${dir}/*.cc"
            "${dir}/*.cxx"
            "${dir}/*.cpp"
        )
        list(APPEND result ${_cpp})
        if(SSDF_ENABLE_CUDA)
            file(GLOB _cu CONFIGURE_DEPENDS "${dir}/*.cu")
            list(APPEND result ${_cu})
        endif()
    endforeach()
    list(REMOVE_DUPLICATES result)
    set(${out_var} ${result} PARENT_SCOPE)
endfunction()

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

# ------------------------------------------------------------------------------
# Library: ssdf
# ------------------------------------------------------------------------------

set(SSDF_SERIAL_DIRS
    "${PROJECT_SOURCE_DIR}/serial"
    "${PROJECT_SOURCE_DIR}/serial/kdtree"
    "${PROJECT_SOURCE_DIR}/serial/graph"
)

if(SSDF_ENABLE_CUBIQL)
    list(APPEND SSDF_SERIAL_DIRS "${PROJECT_SOURCE_DIR}/serial/bvh")
endif()

ssdf_glob_into(SSDF_SERIAL_SOURCES ${SSDF_SERIAL_DIRS})

add_library(ssdf ${SSDF_SERIAL_SOURCES})

set_target_properties(ssdf PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
)

# Feature dependencies
if(SSDF_ENABLE_MPI)
    target_link_libraries(ssdf PUBLIC MPI::MPI_CXX)
endif()

if(SSDF_ENABLE_OPENMP)
    target_link_libraries(ssdf PUBLIC OpenMP::OpenMP_CXX)
endif()

if(SSDF_ENABLE_CUDA)
    set_property(TARGET ssdf PROPERTY CUDA_SEPARABLE_COMPILATION ON)
endif()

if(SSDF_ENABLE_CUBIQL)
    target_link_libraries(ssdf PRIVATE cuBQL)
endif()

target_link_libraries(ssdf PUBLIC smesh::smesh)

# Public include dirs:
# - build interface: serial headers + generated ssdf_config.hpp
# - install interface: allow both <ssdf_config.hpp> and "ssdf.hpp" (from include/ssdf)
target_include_directories(ssdf
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/serial>
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/ssdf>
)

# ------------------------------------------------------------------------------
# Drivers (executables)
# ------------------------------------------------------------------------------

file(GLOB SSDF_DRIVER_SOURCES CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/drivers/*.exe.cpp")
file(GLOB SSDF_DEMO_SOURCES CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/demo/*.exe.cpp")

list(APPEND SSDF_DRIVER_SOURCES ${SSDF_DEMO_SOURCES})

# message(STATUS "SSDF_DRIVER_SOURCES: ${SSDF_DRIVER_SOURCES}")

set(SSDF_DRIVER_TARGETS "")

foreach(driver_src IN LISTS SSDF_DRIVER_SOURCES)
    get_filename_component(_name_we "${driver_src}" NAME_WE) # e.g. sdf.exe
    string(REGEX REPLACE "\\.exe$" "" driver_name "${_name_we}")


    if(driver_name STREQUAL "psdf" AND NOT SSDF_ENABLE_MPI)
        continue()
    endif()
    if(driver_name STREQUAL "points_inside" AND NOT SSDF_ENABLE_CUBIQL)
        continue()
    endif()

    add_executable(${driver_name} "${driver_src}")
    target_link_libraries(${driver_name} PRIVATE ssdf)

    if(driver_name STREQUAL "psdf")
        target_include_directories(psdf PRIVATE "${PROJECT_SOURCE_DIR}")
    endif()

    list(APPEND SSDF_DRIVER_TARGETS ${driver_name})
endforeach()

if(SSDF_DRIVER_TARGETS)
    install(TARGETS ${SSDF_DRIVER_TARGETS}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()

# ------------------------------------------------------------------------------
# Tests
# ------------------------------------------------------------------------------

file(GLOB SSDF_TEST_SOURCES CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/tests/*_test.cpp")
set(SSDF_TEST_TARGETS "")

foreach(test_src IN LISTS SSDF_TEST_SOURCES)
    get_filename_component(test_name "${test_src}" NAME_WE)
    add_executable(${test_name} "${test_src}")
    target_link_libraries(${test_name} PRIVATE ssdf)
    add_test(NAME ${test_name} COMMAND ${test_name})
    list(APPEND SSDF_TEST_TARGETS ${test_name})
endforeach()

if(SSDF_TEST_TARGETS)
    add_custom_target(ssdf_tests ALL DEPENDS ${SSDF_TEST_TARGETS})
endif()

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

# Install headers (core)
install(FILES
    "${PROJECT_SOURCE_DIR}/serial/ssdf.hpp"
    "${PROJECT_SOURCE_DIR}/serial/sqedf.hpp"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/ssdf"
)

install(FILES
    "${PROJECT_SOURCE_DIR}/serial/edf/edf.hpp"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/ssdf/edf"
)

install(FILES
    "${PROJECT_SOURCE_DIR}/serial/kdtree/kdtree.hpp"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/ssdf/kdtree"
)

install(FILES
    "${PROJECT_SOURCE_DIR}/serial/bvh/bvh.hpp"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/ssdf/bvh"
)

# Distributed headers are installed only when MPI is enabled.
if(SSDF_ENABLE_MPI)
    install(FILES
        "${PROJECT_SOURCE_DIR}/distributed/psdf.hpp"
        "${PROJECT_SOURCE_DIR}/distributed/mpi_io.hpp"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/ssdf/distributed"
    )
endif()

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
