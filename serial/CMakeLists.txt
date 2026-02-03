set(SSDF_SERIAL_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/ssdf.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/kdtree/kdtree.cpp"
)

if(SSDF_ENABLE_CUBIQL)
    list(APPEND SSDF_SERIAL_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/bvh/bvh.cpp")
    if(SSDF_ENABLE_CUDA)
        list(APPEND SSDF_SERIAL_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/bvh/bvh.cu")
    endif()
endif()

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

if(SSDF_ENABLE_HIP)
endif()

if(SSDF_ENABLE_CUBIQL)
    target_link_libraries(ssdf PRIVATE cuBQL)
endif()

# Public include dirs:
# - build interface: serial headers + generated ssdf_config.hpp + (repo root for distributed includes in-tree)
# - install interface: allow both <ssdf_config.hpp> and "ssdf.hpp" (from include/ssdf)
target_include_directories(ssdf
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/serial>
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/ssdf>
)

# Install headers (core)
install(FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/ssdf.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/sqedf.hpp"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/ssdf"
)

install(FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/edf/edf.hpp"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/ssdf/edf"
)

install(FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/kdtree/kdtree.hpp"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/ssdf/kdtree"
)

install(FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/bvh/bvh.hpp"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/ssdf/bvh"
)

