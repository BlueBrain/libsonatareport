set(sonatareport_SOURCES
    "library/reports.cpp"
    "library/report.cpp"
    "library/sonatareport.cpp"
    "library/soma_report.cpp"
    "library/element_report.cpp"
    "data/node.cpp"
    "data/soma_node.cpp"
    "data/sonata_data.cpp"
    "io/hdf5_writer.cpp"
    "utils/logger.cpp"
    "utils/imeutil.cpp"
    )

# Shared library
add_library(sonata_report SHARED ${sonatareport_SOURCES})
target_include_directories(sonata_report
    SYSTEM PRIVATE
        ${PROJECT_SOURCE_DIR}/extlib/spdlog/include
    PUBLIC
        ${MPI_INCLUDE_PATH}
        ${HDF5_INCLUDE_DIRS}
        $<BUILD_INTERFACE:${SONATA_REPORT_INCLUDE_DIR}>
        $<INSTALL_INTERFACE:include>
)
target_compile_options(sonata_report
    PRIVATE ${SONATA_REPORT_COMPILE_OPTIONS}
)
if (ENABLE_COVERAGE)
    target_compile_options(sonata_report
        PRIVATE -g -O0 --coverage -fprofile-arcs -ftest-coverage
    )
    target_link_libraries(sonata_report
        PRIVATE gcov
    )
endif()
set_target_properties(sonata_report
    PROPERTIES
    OUTPUT_NAME "sonatareport"
)
target_link_libraries(sonata_report
    PRIVATE spdlog::spdlog_header_only
    PRIVATE ${MPI_C_LIBRARIES}
    PRIVATE ${HDF5_C_LIBRARIES}
)

install(TARGETS sonata_report
    EXPORT sonata-targets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
