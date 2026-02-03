add_executable(sdf "${CMAKE_CURRENT_SOURCE_DIR}/sdf.exe.cpp")
target_link_libraries(sdf PRIVATE ssdf)

add_executable(sqedf "${CMAKE_CURRENT_SOURCE_DIR}/sqedf.exe.cpp")
target_link_libraries(sqedf PRIVATE ssdf)

if(SSDF_ENABLE_MPI)
    add_executable(psdf "${CMAKE_CURRENT_SOURCE_DIR}/psdf.exe.cpp")
    target_link_libraries(psdf PRIVATE ssdf)
    target_include_directories(psdf PRIVATE "${PROJECT_SOURCE_DIR}")
endif()

install(TARGETS sdf sqedf
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

if(SSDF_ENABLE_MPI)
    install(TARGETS psdf
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()

