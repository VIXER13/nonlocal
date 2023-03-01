cmake_minimum_required(VERSION 3.16)
project(structural_solver_2d)

set(STRUCTURAL_SOLVER_2D_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_library(structural_solver_2d_lib INTERFACE)
target_sources(structural_solver_2d_lib INTERFACE
    equilibrium_equation_2d.hpp
    mechanical_boundary_conditions_2d.hpp
    mechanical_parameters_2d.hpp
    mechanical_solution_2d.hpp
    stiffness_matrix_2d.hpp
    temperature_condition_2d.hpp
)
target_include_directories(structural_solver_2d_lib INTERFACE ${STRUCTURAL_SOLVER_2D_LIB_DIR})