list(APPEND SSDF_BVH_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/bvh.cpp")
list(APPEND SSDF_BVH_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/bvh.hpp")

if(SSDF_ENABLE_CUDA)
    message(STATUS "Building BVH CUDA code")
    list(APPEND SSDF_BVH_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/bvh.cu")
endif()

target_sources(ssdf PRIVATE ${SSDF_BVH_SOURCES})
target_include_directories(ssdf     PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

