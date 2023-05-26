cmake_minimum_required(VERSION 3.16)

project(finite_elements_tests)

add_library(finite_elements_test_lib OBJECT
    init_element.cpp
    element_1d.cpp
    quadratures.cpp
)

target_include_directories(finite_elements_test_lib PUBLIC ".")
target_link_libraries(finite_elements_test_lib
    Boost::ut
    metamath_lib
)