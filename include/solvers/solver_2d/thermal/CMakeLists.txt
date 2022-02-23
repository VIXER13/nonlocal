cmake_minimum_required(VERSION 3.17)
project(heat_equation_solver_2d)

set(HEAT_EQUATION_SOLVER_2D_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_library(heat_equation_solver_2d_lib INTERFACE)
target_sources(heat_equation_solver_2d_lib INTERFACE
               heat_equation_solver_2d.hpp
               heat_equation_solution_2d.hpp
               heat_equation_parameters_2d.hpp
               thermal_conductivity_matrix_2d.hpp
               heat_capacity_matrix_2d.hpp
               stationary_heat_equation_solver_2d.hpp
               nonstationary_heat_equation_solver_2d.hpp)
target_include_directories(heat_equation_solver_2d_lib INTERFACE ${HEAT_EQUATION_SOLVER_2D_LIB_DIR})