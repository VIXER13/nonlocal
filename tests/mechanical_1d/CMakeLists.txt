cmake_minimum_required(VERSION 3.16)

project(mechanical_1d_test)


add_library(
    mechanical_1d_test_lib
    OBJECT
    stationary.cpp
)
target_include_directories(
    mechanical_1d_test_lib
    PUBLIC #
        ${LIBRARIES_DIR}
        
)

target_link_libraries(
    mechanical_1d_test_lib
    PRIVATE #
        finite_element_solver_1d_lib
        mechanical_equation_solver_1d_lib
        Boost::ut
)
