add_subdirectory(binary_reader)
# Convert binary reports to SONATA format
add_executable(spikes_converter spikes_converter.cpp)
target_include_directories(spikes_converter
    PRIVATE
        $<BUILD_INTERFACE:${SONATA_REPORT_INCLUDE_DIR}>
        $<BUILD_INTERFACE:${SONATA_REPORT_SRC_DIR}>
    SYSTEM
        ${PROJECT_SOURCE_DIR}/extlib/spdlog/include
)
target_compile_options(spikes_converter
    PRIVATE ${SONATA_REPORT_COMPILE_OPTIONS}
)
target_link_libraries(spikes_converter
    PUBLIC sonata_report
    PRIVATE spdlog::spdlog_header_only
    PRIVATE ${MPI_C_LIBRARIES}
)

add_executable(reports_converter reports_converter.cpp)
target_include_directories(reports_converter
    PRIVATE
        $<BUILD_INTERFACE:${SONATA_REPORT_INCLUDE_DIR}>
        $<BUILD_INTERFACE:${SONATA_REPORT_SRC_DIR}>
    SYSTEM
        ${PROJECT_SOURCE_DIR}/extlib/spdlog/include
)
target_compile_options(reports_converter
    PRIVATE ${SONATA_REPORT_COMPILE_OPTIONS}
)
target_link_libraries(reports_converter
    PUBLIC sonata_report
    PUBLIC binary_reader
    PRIVATE spdlog::spdlog_header_only
    PRIVATE ${MPI_C_LIBRARIES}
)

install(TARGETS spikes_converter DESTINATION bin)
install(TARGETS reports_converter DESTINATION bin)