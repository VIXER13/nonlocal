cmake_minimum_required(VERSION 3.16)
project(heat_equation_solver_1d)

add_library(heat_equation_solver_1d_lib INTERFACE)
target_sources(
    heat_equation_solver_1d_lib INTERFACE
    thermal_parameters_1d.hpp
    thermal_conductivity_assembler.hpp
    heat_capacity_assembler.hpp
    init_problem_settings.hpp
    convection_condition_1d.hpp
    stationary_heat_equation_solver_1d.hpp
    nonstationary_heat_equation_solver_1d.hpp
    heat_equation_solution_1d.hpp
)
target_include_directories(heat_equation_solver_1d_lib INTERFACE ${LIBRARIES_DIR})
target_link_libraries(heat_equation_solver_1d_lib
    INTERFACE
    mesh_1d_lib
    finite_element_solver_1d_base_lib
)
