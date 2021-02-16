set(SONATA_REPORT_PYTHON_BUILD ${CMAKE_CURRENT_BINARY_DIR}/libsonata)

if (SONATA_REPORT_ENABLE_SUBMODULES)
    add_subdirectory(pybind11 EXCLUDE_FROM_ALL)
else()
    find_package(pybind11)
endif()

pybind11_add_module(sonatareport_python SYSTEM bindings.cpp)

set_target_properties(sonatareport_python PROPERTIES
    OUTPUT_NAME "_libsonatareport"
)

target_link_libraries(sonatareport_python
    PRIVATE sonata_report
    PRIVATE pybind11::module
)
