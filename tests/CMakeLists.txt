cmake_minimum_required(VERSION 3.16)

project(nonlocal_tests)

enable_testing()

FileEmbedSetup()
FileEmbedAdd(config/required_fields.json)
FileEmbedAdd(config/save_fields.json)

add_executable(unit_tests 
    nonlocal_tests.cpp
    config/required_fields.cpp
    config/save_fields.cpp
    finite_elements/quadratures.cpp
    finite_elements/elements.cpp
    parallel_utils/uniform_ranges.cpp
)
target_include_directories(unit_tests PUBLIC ".")

target_link_libraries(unit_tests 
    gtest_main
    gtest
    embedded_files
    metamath_lib
    config_lib
    parallel_utils_lib
)