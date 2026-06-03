cmake_minimum_required(VERSION 3.15)
project(ssdf VERSION 0.1.0 LANGUAGES CXX C)

if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
endif()

if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

function(ssdf_enable_strict_warnings target_name)
    if(MSVC)
        target_compile_options(${target_name} PRIVATE /W4 /permissive-)
        if(SSDF_WARNINGS_AS_ERRORS)
            target_compile_options(${target_name} PRIVATE /WX)
        endif()
    else()
        target_compile_options(${target_name} PRIVATE -Wall -Wextra -Wpedantic)
        if(SSDF_WARNINGS_AS_ERRORS)
            target_compile_options(${target_name} PRIVATE -Werror)
        endif()
    endif()
endfunction()

# Standard install dirs (bin, lib, include, ...)
include(GNUInstallDirs)
include(cmake/SSDFOptions.cmake)
include(cmake/SSDFDependencies.cmake)

set(SSDF_BUILD_INCLUDES "")
list(APPEND SSDF_BUILD_INCLUDES
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/serial>
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/ssdf>
)

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
function(ssdf_glob_files out_var)
    set(options)
    set(one_value_args)
    set(multi_value_args DIRECTORIES PATTERNS)
    cmake_parse_arguments(SSDF_GLOB "${options}" "${one_value_args}" "${multi_value_args}" ${ARGN})

    if(NOT SSDF_GLOB_DIRECTORIES)
        message(FATAL_ERROR "ssdf_glob_files(): DIRECTORIES must not be empty")
    endif()
    if(NOT SSDF_GLOB_PATTERNS)
        message(FATAL_ERROR "ssdf_glob_files(): PATTERNS must not be empty")
    endif()

    set(result "")
    foreach(dir IN LISTS SSDF_GLOB_DIRECTORIES)
        get_filename_component(abs_dir "${dir}" ABSOLUTE BASE_DIR "${PROJECT_SOURCE_DIR}")
        if(NOT IS_DIRECTORY "${abs_dir}")
            message(FATAL_ERROR "ssdf_glob_files(): directory does not exist: ${dir}")
        endif()
        foreach(pattern IN LISTS SSDF_GLOB_PATTERNS)
            file(GLOB _matches CONFIGURE_DEPENDS "${abs_dir}/${pattern}")
            list(APPEND result ${_matches})
        endforeach()
    endforeach()
    list(REMOVE_DUPLICATES result)
    list(SORT result)
    set(${out_var} ${result} PARENT_SCOPE)
endfunction()

function(ssdf_configure_binary target_name)
    target_link_libraries(${target_name} PRIVATE ssdf)
    ssdf_enable_strict_warnings(${target_name})
    if(SSDF_BINARY_PRIVATE_COMPILE_OPTIONS)
        target_compile_options(${target_name} PRIVATE ${SSDF_BINARY_PRIVATE_COMPILE_OPTIONS})
    endif()
    if(SSDF_BINARY_PRIVATE_LINK_OPTIONS)
        target_link_options(${target_name} PRIVATE ${SSDF_BINARY_PRIVATE_LINK_OPTIONS})
    endif()

    if(target_name STREQUAL "psdf")
        if(SSDF_PSDF_PRIVATE_LIBRARIES)
            target_link_libraries(${target_name} PRIVATE ${SSDF_PSDF_PRIVATE_LIBRARIES})
        endif()
        if(SSDF_PSDF_PRIVATE_INCLUDES)
            target_include_directories(${target_name} PRIVATE ${SSDF_PSDF_PRIVATE_INCLUDES})
        endif()
    endif()
endfunction()

enable_testing()

# ------------------------------------------------------------------------------
# Library: ssdf
# ------------------------------------------------------------------------------

set(SSDF_SERIAL_DIRS
    "serial"
    "serial/kdtree"
    "serial/graph"
)

if(SSDF_ENABLE_CUBIQL)
    list(APPEND SSDF_SERIAL_DIRS "serial/bvh")
endif()

set(SSDF_SERIAL_SOURCE_PATTERNS "*.c" "*.cc" "*.cxx" "*.cpp")
if(SSDF_ENABLE_CUDA)
    list(APPEND SSDF_SERIAL_SOURCE_PATTERNS "*.cu")
endif()

ssdf_glob_files(
    SSDF_SERIAL_SOURCES
    DIRECTORIES ${SSDF_SERIAL_DIRS}
    PATTERNS ${SSDF_SERIAL_SOURCE_PATTERNS}
)

add_library(ssdf ${SSDF_SERIAL_SOURCES})
target_compile_features(ssdf PUBLIC cxx_std_17)
ssdf_enable_strict_warnings(ssdf)

set_target_properties(ssdf PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
)

if(SSDF_ENABLE_CUDA)
    set_property(TARGET ssdf PROPERTY CUDA_SEPARABLE_COMPILATION ON)
endif()

if(SSDF_PUBLIC_DEP_LIBRARIES)
    target_link_libraries(ssdf PUBLIC ${SSDF_PUBLIC_DEP_LIBRARIES})
endif()
if(SSDF_PRIVATE_DEP_LIBRARIES)
    target_link_libraries(ssdf PRIVATE ${SSDF_PRIVATE_DEP_LIBRARIES})
endif()
if(SSDF_PRIVATE_COMPILE_OPTIONS)
    target_compile_options(ssdf PRIVATE ${SSDF_PRIVATE_COMPILE_OPTIONS})
endif()
if(SSDF_PRIVATE_LINK_OPTIONS)
    target_link_options(ssdf PRIVATE ${SSDF_PRIVATE_LINK_OPTIONS})
endif()

# Public include dirs:
# - build interface: serial headers + generated ssdf_config.hpp
# - install interface: allow both <ssdf_config.hpp> and "ssdf.hpp" (from include/ssdf)
target_include_directories(ssdf PUBLIC ${SSDF_BUILD_INCLUDES})

# ------------------------------------------------------------------------------
# Drivers and demos (executables)
# ------------------------------------------------------------------------------

set(SSDF_APP_DIRS
    "drivers"
)

if(SSDF_ENABLE_SMESH)
    list(APPEND SSDF_APP_DIRS "demo")
endif()

if(SSDF_ENABLE_MPI)
    list(APPEND SSDF_APP_DIRS
        "drivers/distributed"
    )
endif()

ssdf_glob_files(
    SSDF_APP_SOURCES
    DIRECTORIES ${SSDF_APP_DIRS}
    PATTERNS "*.exe.cpp"
)

set(SSDF_APP_TARGETS "")

foreach(app_src IN LISTS SSDF_APP_SOURCES)
    get_filename_component(_name_we "${app_src}" NAME_WE) # e.g. sdf.exe
    string(REGEX REPLACE "\\.exe$" "" app_name "${_name_we}")

    add_executable(${app_name} "${app_src}")
    ssdf_configure_binary(${app_name})
    list(APPEND SSDF_APP_TARGETS ${app_name})
endforeach()

if(SSDF_APP_TARGETS)
    add_custom_target(ssdf_apps ALL DEPENDS ${SSDF_APP_TARGETS})
    install(TARGETS ${SSDF_APP_TARGETS}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()

# ------------------------------------------------------------------------------
# Tests
# ------------------------------------------------------------------------------

ssdf_glob_files(
    SSDF_TEST_SOURCES
    DIRECTORIES "tests"
    PATTERNS "*_test.cpp"
)

if(NOT SSDF_ENABLE_CUBIQL)
    list(FILTER SSDF_TEST_SOURCES EXCLUDE REGEX "/closest_within_radius_test\\.cpp$")
    list(FILTER SSDF_TEST_SOURCES EXCLUDE REGEX "/potential_contact_triangles_bvh_test\\.cpp$")
    list(FILTER SSDF_TEST_SOURCES EXCLUDE REGEX "/potential_contact_quads_bvh_test\\.cpp$")
endif()

set(SSDF_TEST_TARGETS "")

foreach(test_src IN LISTS SSDF_TEST_SOURCES)
    get_filename_component(test_name "${test_src}" NAME_WE)
    add_executable(${test_name} "${test_src}")
    ssdf_configure_binary(${test_name})
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

set(SSDF_PUBLIC_HEADER_SUBDIRS "" "edf" "kdtree")
if(SSDF_ENABLE_CUBIQL)
    list(APPEND SSDF_PUBLIC_HEADER_SUBDIRS "bvh")
endif()

foreach(header_subdir IN LISTS SSDF_PUBLIC_HEADER_SUBDIRS)
    set(header_dir "serial")
    set(install_dir "${CMAKE_INSTALL_INCLUDEDIR}/ssdf")

    if(header_subdir)
        string(APPEND header_dir "/${header_subdir}")
        string(APPEND install_dir "/${header_subdir}")
    endif()

    ssdf_glob_files(
        SSDF_PUBLIC_HEADERS
        DIRECTORIES "${header_dir}"
        PATTERNS "*.h" "*.hh" "*.hpp" "*.hxx"
    )

    if(SSDF_PUBLIC_HEADERS)
        install(FILES ${SSDF_PUBLIC_HEADERS}
            DESTINATION "${install_dir}"
        )
    endif()
endforeach()

# Distributed headers are installed only when MPI is enabled.
if(SSDF_ENABLE_MPI)
    ssdf_glob_files(
        SSDF_DISTRIBUTED_HEADERS
        DIRECTORIES "distributed"
        PATTERNS "*.h" "*.hh" "*.hpp" "*.hxx"
    )
    if(SSDF_DISTRIBUTED_HEADERS)
        install(FILES ${SSDF_DISTRIBUTED_HEADERS}
            DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/ssdf/distributed"
        )
    endif()
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
