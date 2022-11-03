configure_file(converter_test.sh.in converter_test.sh @ONLY)

add_test(NAME converter_test
    COMMAND "/bin/sh" ${CMAKE_CURRENT_BINARY_DIR}/converter_test.sh
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)
