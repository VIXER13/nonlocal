cmake_minimum_required(VERSION 3.16)

project(nonlocal_tests)

FileEmbedSetup()

add_subdirectory(config)
add_subdirectory(finite_elements)
add_subdirectory(math_expression_tests)
add_subdirectory(nonstationary_thermal_1d)
add_subdirectory(parallel)


add_executable(unit_tests nonlocal_tests.cpp)
target_include_directories(unit_tests PUBLIC ${LIBRARIES_DIR})
target_link_libraries(unit_tests
    finite_elements_test_lib
    parallel_test_lib
    config_test_lib
    nonstationary_thermal_1d_lib
    math_expression_test_lib
)

add_test(NAME unit_tests COMMAND unit_tests)

