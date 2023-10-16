cmake_minimum_required(VERSION 3.16)

project(nonlocal_tests)

FileEmbedSetup()

add_subdirectory(config)
add_subdirectory(finite_elements)
add_subdirectory(parallel_utils)

add_executable(unit_tests nonlocal_tests.cpp)
target_include_directories(unit_tests PUBLIC ".")
target_link_libraries(unit_tests
    finite_elements_test_lib
    parallel_utils_test_lib
    config_test_lib
)