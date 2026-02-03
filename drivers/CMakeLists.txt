add_executable(sdf "${CMAKE_CURRENT_SOURCE_DIR}/sdf.exe.cpp")
target_link_libraries(sdf PRIVATE ssdf)

add_executable(sqedf "${CMAKE_CURRENT_SOURCE_DIR}/sqedf.exe.cpp")
target_link_libraries(sqedf PRIVATE ssdf)

add_executable(sdf_sign "${CMAKE_CURRENT_SOURCE_DIR}/sdf_sign.exe.cpp")
target_link_libraries(sdf_sign PRIVATE ssdf)

if(SSDF_ENABLE_CUBIQL)
    add_executable(points_inside "${CMAKE_CURRENT_SOURCE_DIR}/points_inside.exe.cpp")
    target_link_libraries(points_inside PRIVATE ssdf)
endif()

if(SSDF_ENABLE_MPI)
    add_executable(psdf "${CMAKE_CURRENT_SOURCE_DIR}/psdf.exe.cpp")
    target_link_libraries(psdf PRIVATE ssdf)
    target_include_directories(psdf PRIVATE "${PROJECT_SOURCE_DIR}")
endif()

install(TARGETS sdf sqedf sdf_sign
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

if(SSDF_ENABLE_CUBIQL)
    install(TARGETS points_inside
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()

if(SSDF_ENABLE_MPI)
    install(TARGETS psdf
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()

