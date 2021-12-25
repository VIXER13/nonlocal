cmake_minimum_required(VERSION 3.17)
project(finite_element_solver_1d)

add_subdirectory(thermal)

set(FINITE_ELEMENT_SOLVER_1D_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_library(finite_element_solver_1d_lib INTERFACE)
target_sources(finite_element_solver_1d_lib INTERFACE
               finite_element_matrix_1d.hpp
               boundary_condition_1d.hpp
               right_part_1d.hpp
               parameters_1d.hpp
               influence_functions_1d.hpp)
target_include_directories(finite_element_solver_1d_lib INTERFACE ${FINITE_ELEMENT_SOLVER_1D_LIB_DIR})
target_link_libraries(finite_element_solver_1d_lib INTERFACE heat_equation_solver_1d_lib)