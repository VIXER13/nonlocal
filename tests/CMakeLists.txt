cmake_minimum_required(VERSION 3.16)

project(nonlocal_tests)

enable_testing()

add_executable(unit_tests 
    nonlocal_tests.cpp
    finite_elements/quadratures.cpp
    finite_elements/elements.cpp
)
target_include_directories(unit_tests PUBLIC ".")
target_link_libraries(unit_tests gtest_main gtest metamath_lib)