cmake_minimum_required(VERSION 3.17)
project(heat_equation_solver_1d)

set(HEAT_EQUATION_SOLVER_1D_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_library(heat_equation_solver_1d_lib INTERFACE)
target_sources(heat_equation_solver_1d_lib INTERFACE
               thermal_conductivity_matrix_1d.hpp
               heat_capacity_matrix_1d.hpp
               boundary_condition_third_kind_1d.hpp
               stationary_heat_equation_solver_1d.hpp
               nonstationary_heat_equation_solver.hpp)
target_include_directories(heat_equation_solver_1d_lib INTERFACE ${HEAT_EQUATION_SOLVER_2D_LIB_DIR})