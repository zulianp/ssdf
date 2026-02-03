# MPI-related headers are installed only when MPI is enabled.
if(SSDF_ENABLE_MPI)
    install(FILES
        "${CMAKE_CURRENT_SOURCE_DIR}/psdf.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/mpi_io.hpp"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/ssdf/distributed"
    )
endif()

