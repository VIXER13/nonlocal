cmake_minimum_required(VERSION 3.16)

project(finite_elements_tests)

add_library(finite_elements_test_lib OBJECT
    init_lagrangian_elements_1d.cpp
    init_lagrangian_elements_2d.cpp
    init_triangle_elements_2d.cpp
    init_serendipity_elements_2d.cpp
    quadratures_test.cpp
    element_1d_test.cpp
    element_2d_test.cpp
)

target_include_directories(finite_elements_test_lib PUBLIC ".")
target_link_libraries(finite_elements_test_lib
    Boost::ut
    metamath_lib
)