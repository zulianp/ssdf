file(GLOB SSDF_TEST_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*_test.cpp")
set(SSDF_TEST_TARGETS "")

foreach(test_src ${SSDF_TEST_SOURCES})
    get_filename_component(test_name "${test_src}" NAME_WE)
    add_executable(${test_name} "${test_src}")
    target_link_libraries(${test_name} PRIVATE ssdf)
    add_test(NAME ${test_name} COMMAND ${test_name})
    list(APPEND SSDF_TEST_TARGETS ${test_name})
endforeach()

if(SSDF_TEST_TARGETS)
    add_custom_target(ssdf_tests ALL DEPENDS ${SSDF_TEST_TARGETS})
endif()

