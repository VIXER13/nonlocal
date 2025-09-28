cmake_minimum_required(VERSION 3.16)

project(nonlocal_tests)

set(PROJECT_ROOT ${CMAKE_SOURCE_DIR})
message(STATUS "PROJECT_ROOT: ${PROJECT_ROOT}")

FileEmbedSetup()

add_subdirectory(config)
add_subdirectory(finite_elements)
add_subdirectory(math_expression_tests)
add_subdirectory(parallel)
add_subdirectory(thermal_1d)
add_subdirectory(solvers_2d)
add_subdirectory(stationary_thermal_2d)

add_executable(unit_tests nonlocal_tests.cpp)
target_include_directories(unit_tests PUBLIC ${LIBRARIES_DIR})
target_link_libraries(unit_tests
    embedded_files
    finite_elements_test_lib
    math_expression_test_lib
    parallel_test_lib
    config_test_lib
	thermal_1d_test_lib
    thermal_2d_test_lib
    solvers_2d_test_lib
    stationary_thermal_2d_lib
)

add_test(NAME unit_tests COMMAND unit_tests)

