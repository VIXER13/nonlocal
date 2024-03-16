cmake_minimum_required(VERSION 3.16)

project(nonlocal_tests)

FileEmbedSetup()

add_subdirectory(config)
add_subdirectory(finite_elements)
add_subdirectory(parallel_utils)
add_subdirectory(test_template)

add_executable(unit_tests nonlocal_tests.cpp)
target_include_directories(unit_tests PUBLIC ".")
target_link_libraries(unit_tests
    finite_elements_test_lib
    parallel_utils_test_lib
    config_test_lib
)

add_executable(test_template nonlocal_tests.cpp)
target_include_directories(test_template PUBLIC ".")
target_link_libraries(test_template
    test_template_lib
)

add_test(NAME unit_tests COMMAND unit_tests)
add_test(NAME test_template COMMAND test_template)





