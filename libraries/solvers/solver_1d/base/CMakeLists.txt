cmake_minimum_required(VERSION 3.16)
project(finite_element_solver_1d_base)

add_library(finite_element_solver_1d_base_lib INTERFACE)
target_sources(finite_element_solver_1d_base_lib INTERFACE
    assebmler_base.hpp
    assemble_matrix_portrait.hpp
    boundary_condition_first_kind_1d.hpp
    boundary_condition_second_kind_1d.hpp
    boundary_conditions_1d.hpp
    finite_element_matrix.hpp
    problem_settings.hpp
    right_part_1d.hpp
    solution_1d.hpp
)
target_include_directories(finite_element_solver_1d_base_lib INTERFACE ${LIBRARIES_DIR})
target_link_libraries(finite_element_solver_1d_base_lib INTERFACE 
    assembler_base_lib
    mesh_1d_lib
)