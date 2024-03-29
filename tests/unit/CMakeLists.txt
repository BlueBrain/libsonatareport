set(TEST_SOURCES
    tests.cpp
    test_node.cpp
    test_report.cpp
    test_sonatadata.cpp
    test_sonatareport.cpp
    )

add_executable(reports_unit_tests ${TEST_SOURCES})
target_include_directories(reports_unit_tests
    PRIVATE
        $<BUILD_INTERFACE:${SONATA_REPORT_INCLUDE_DIR}>
        $<BUILD_INTERFACE:${SONATA_REPORT_SRC_DIR}>
    SYSTEM
        ${PROJECT_SOURCE_DIR}/extlib/spdlog/include
)
target_link_libraries(reports_unit_tests
    PRIVATE
        sonata_report
        Catch2::Catch2
        spdlog::spdlog_header_only
        ${MPI_C_LIBRARIES}
        ${HDF5_C_LIBRARIES}
)

catch_discover_tests(reports_unit_tests)
