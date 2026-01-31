add_library(ssdf_edf INTERFACE)

target_sources(ssdf_edf INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/edf.hpp>
    $<INSTALL_INTERFACE:include/ssdf/edf/edf.hpp>
)

target_include_directories(ssdf_edf INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:include/ssdf>
)

install(TARGETS ssdf_edf
    EXPORT ssdfTargets
)

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/edf.hpp
    DESTINATION include/ssdf/edf
)
