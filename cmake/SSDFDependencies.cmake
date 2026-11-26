include_guard(GLOBAL)

set(SSDF_PUBLIC_DEP_LIBRARIES "")
set(SSDF_PRIVATE_DEP_LIBRARIES "")
set(SSDF_PSDF_PRIVATE_LIBRARIES "")
set(SSDF_PSDF_PRIVATE_INCLUDES "")
set(SSDF_PRIVATE_COMPILE_OPTIONS "")
set(SSDF_PRIVATE_LINK_OPTIONS "")
set(SSDF_BINARY_PRIVATE_COMPILE_OPTIONS "")
set(SSDF_BINARY_PRIVATE_LINK_OPTIONS "")


if(SSDF_ENABLE_SMESH)
    if(NOT TARGET smesh::smesh)
        find_package(smesh CONFIG REQUIRED)
    endif()
    list(APPEND SSDF_PUBLIC_DEP_LIBRARIES smesh::smesh)

    function(ssdf_inherit_feature_from_smesh feature_name dependency_token)
        if(NOT TARGET smesh::smesh)
            return()
        endif()

        get_target_property(SSDF_SMESH_LINK_LIBRARIES smesh::smesh INTERFACE_LINK_LIBRARIES)
        if(NOT SSDF_SMESH_LINK_LIBRARIES MATCHES "(^|;)${dependency_token}(;|$)")
            return()
        endif()

        if(NOT ${feature_name})
            message(STATUS "${feature_name} enabled to match smesh")
        endif()
        set(${feature_name} ON CACHE BOOL "" FORCE)
        set(${feature_name} ON PARENT_SCOPE)
    endfunction()

    ssdf_inherit_feature_from_smesh(SSDF_ENABLE_MPI "MPI::MPI_CXX")
    ssdf_inherit_feature_from_smesh(SSDF_ENABLE_OPENMP "OpenMP::OpenMP_CXX")

    if(TARGET smesh::smesh)
        get_target_property(SSDF_SMESH_LINK_LIBRARIES smesh::smesh INTERFACE_LINK_LIBRARIES)
        if(SSDF_SMESH_LINK_LIBRARIES MATCHES "(^|;)OpenMP::OpenMP_CXX(;|$)")
            set(SSDF_ENABLE_OPENMP ON)
        endif()
    endif()
endif()



if(SSDF_ENABLE_MPI)
    find_package(MPI REQUIRED)
    message(STATUS "MPI found: ${MPI_FOUND}")
    message(STATUS "MPI CXX compiler: ${MPI_CXX_COMPILER}")
    list(APPEND SSDF_PUBLIC_DEP_LIBRARIES MPI::MPI_CXX)
    list(APPEND SSDF_PSDF_PRIVATE_LIBRARIES MPI::MPI_CXX)
endif()

if(SSDF_ENABLE_OPENMP)
    find_package(OpenMP COMPONENTS CXX REQUIRED)
    message(STATUS "OpenMP found: ${OpenMP_FOUND}")
    message(STATUS "OpenMP CXX: ${OpenMP_CXX_VERSION}")
    message(STATUS "OpenMP CXX flags: ${OpenMP_CXX_FLAGS}")
    list(APPEND SSDF_PUBLIC_DEP_LIBRARIES OpenMP::OpenMP_CXX)
endif()

if(SSDF_ENABLE_HIP)
    # The platform as given (enable_language may detect it later, below).
    set(_ssdf_hip_platform_given "${CMAKE_HIP_PLATFORM}")
    if(NOT CMAKE_HIP_PLATFORM)
        # Not given: enable_language takes it from `hipconfig --platform` (CMake's
        # CMakeDetermineHIPCompiler.cmake), which follows the HIP_PLATFORM environment variable and
        # on chipStar reports spirv. Checked before the AMD architectures below are set, on which
        # CMake's own compiler test would fail first.
        set(_ssdf_hip_platform "$ENV{HIP_PLATFORM}")
        if(NOT _ssdf_hip_platform)
            execute_process(COMMAND hipconfig --platform
                OUTPUT_VARIABLE _ssdf_hipconfig_platform OUTPUT_STRIP_TRAILING_WHITESPACE
                RESULT_VARIABLE _ssdf_hipconfig_result
                ERROR_QUIET)
            if(_ssdf_hipconfig_result EQUAL 0)
                set(_ssdf_hip_platform "${_ssdf_hipconfig_platform}")
            endif()
        endif()
        if(_ssdf_hip_platform STREQUAL "spirv")
            message(FATAL_ERROR "The HIP installation is chipStar's (HIP_PLATFORM or hipconfig --platform: spirv), but "
                                "the platform was not given: configure with -DCMAKE_HIP_PLATFORM=spirv")
        endif()
    endif()
    # The .cu sources are compiled as HIP (see CMakeLists.txt). The architectures must be set
    # before enable_language.
    if(CMAKE_HIP_PLATFORM STREQUAL "spirv")
        # chipStar (Intel GPUs, e.g. Aurora's PVC): HIP compiled to SPIR-V, with no architecture list.
        # CMake knows this platform from 4.3.
        if(CMAKE_VERSION VERSION_LESS 4.3)
            message(FATAL_ERROR "CMAKE_HIP_PLATFORM=spirv (chipStar) needs CMake 4.3 or later; this is ${CMAKE_VERSION}")
        endif()
    elseif(NOT DEFINED CMAKE_HIP_ARCHITECTURES)
        set(CMAKE_HIP_ARCHITECTURES "${SSDF_HIP_ARCHITECTURES}")
    endif()
    enable_language(HIP)
    if(CMAKE_HIP_PLATFORM STREQUAL "spirv" AND NOT _ssdf_hip_platform_given STREQUAL "spirv")
        # CMake detected spirv after SSDF_HIP_ARCHITECTURES (gfx...) was set above: the two do not
        # mix. The check above catches this first; this one stays for a CMake that detects the
        # platform by other means.
        message(FATAL_ERROR "The HIP compiler is chipStar's (CMAKE_HIP_PLATFORM detected as spirv), but the "
                            "platform was not given: configure a new build directory with -DCMAKE_HIP_PLATFORM=spirv")
    endif()
    message(STATUS "HIP enabled for ${CMAKE_HIP_ARCHITECTURES}")
    if(NOT DEFINED CMAKE_HIP_STANDARD)
        set(CMAKE_HIP_STANDARD 17)
        set(CMAKE_HIP_STANDARD_REQUIRED ON)
    endif()
    # Executables made of C++ sources are linked by the C++ compiler, as smesh's and SCCD's are.
    set(CMAKE_HIP_LINKER_PREFERENCE_PROPAGATES OFF)
    if(NOT TARGET hip::host)
        find_package(hip REQUIRED)
    endif()
    # hipCUB (cuBQL's GPU builder sorts with it) is header-only and used by the .cu files only; its
    # target is not linked, because through rocPRIM it brings hip::device, which compiles every C++
    # source as HIP.
    if(NOT TARGET hip::hipcub)
        find_package(hipcub REQUIRED)
    endif()
    list(APPEND SSDF_PUBLIC_DEP_LIBRARIES hip::host)
    # cuBQL enables CUDA whenever it finds nvcc; with HIP its headers are compiled by ssdf's HIP
    # sources instead.
    set(CUBQL_DISABLE_CUDA ON CACHE BOOL "" FORCE)
elseif(SSDF_ENABLE_CUDA)
    enable_language(CUDA)
    message(STATUS "CUDA enabled")
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)

    if(NOT CMAKE_CUDA_ARCHITECTURES)
        set(CMAKE_CUDA_ARCHITECTURES "native")
    endif()
endif()

if(SSDF_ENABLE_CUBIQL)
    if(NOT cuBQL_DIR)
        set(cuBQL_DIR "${PROJECT_SOURCE_DIR}/external/cubql" CACHE PATH "Path to cuBQL (git submodule: external/cubql)")
    endif()
    if(NOT EXISTS "${cuBQL_DIR}/CMakeLists.txt")
        message(FATAL_ERROR "cuBQL not found at ${cuBQL_DIR}. Clone submodules: git submodule update --init --recursive")
    endif()
    add_subdirectory(${cuBQL_DIR} "${CMAKE_BINARY_DIR}/_deps/cuBQL-build")
    # Third-party headers: include as system so our -Werror does not apply to them.
    get_target_property(_ssdf_cubql_includes cuBQL INTERFACE_INCLUDE_DIRECTORIES)
    set_target_properties(cuBQL PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${_ssdf_cubql_includes}")
    # Only link cuBQL in the build tree so install(EXPORT ssdfTargets) does not
    # reference a non-exported target. Requires a shared libssdf (or otherwise
    # self-contained link) so consumers are not missing cuBQL symbols.
    list(APPEND SSDF_PRIVATE_DEP_LIBRARIES $<BUILD_INTERFACE:cuBQL>)
endif()

if(SSDF_ENABLE_ASAN)
    if(NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
        message(FATAL_ERROR "SSDF_ENABLE_ASAN requires Clang or GCC")
    endif()
    if(SSDF_ENABLE_CUDA OR SSDF_ENABLE_HIP)
        message(FATAL_ERROR "SSDF_ENABLE_ASAN does not support CUDA or HIP builds")
    endif()

    list(APPEND SSDF_PRIVATE_COMPILE_OPTIONS -fsanitize=address -fno-omit-frame-pointer)
    list(APPEND SSDF_PRIVATE_LINK_OPTIONS -fsanitize=address)
    list(APPEND SSDF_BINARY_PRIVATE_COMPILE_OPTIONS -fsanitize=address -fno-omit-frame-pointer)
    list(APPEND SSDF_BINARY_PRIVATE_LINK_OPTIONS -fsanitize=address)
endif()

list(APPEND SSDF_PSDF_PRIVATE_INCLUDES "${PROJECT_SOURCE_DIR}")
