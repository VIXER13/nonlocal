cmake_minimum_required(VERSION 3.16)
project(heat_equation_solver_1d)

set(HEAT_EQUATION_SOLVER_1D_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_library(heat_equation_solver_1d_lib INTERFACE)
target_sources(
    heat_equation_solver_1d_lib INTERFACE
    thermal_parameters_1d.hpp
    thermal_conductivity_matrix_assembler_1d.hpp
    heat_capacity_matrix_assembler_1d.hpp
    convection_condition_1d.hpp
    stationary_heat_equation_solver_1d.hpp
    nonstationary_heat_equation_solver_1d.hpp
    heat_equation_solution_1d.hpp
)
target_include_directories(heat_equation_solver_1d_lib INTERFACE 
    ${HEAT_EQUATION_SOLVER_2D_LIB_DIR}
)
target_link_libraries(heat_equation_solver_1d_lib INTERFACE
    mesh_1d_lib
    finite_element_solver_1d_base_lib
)