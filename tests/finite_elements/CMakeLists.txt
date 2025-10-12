cmake_minimum_required(VERSION 3.16)

project(finite_elements_tests)

set(INCLUDES $ENV{CPATH})

add_library(
    finite_elements_test_lib
    OBJECT
    quadratures_test.cpp
    element_1d_test.cpp
    element_2d_test.cpp
)

target_include_directories(
    finite_elements_test_lib
    PUBLIC #
        ${LIBRARIES_DIR}
        ${INCLUDES}
)
target_link_libraries(
    finite_elements_test_lib
    PRIVATE #
        metamath_lib
        Boost::ut
)
