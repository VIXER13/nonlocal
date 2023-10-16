cmake_minimum_required(VERSION 3.16)
project(finite_element_solver_2d_base)

set(FINITE_ELEMENT_SOLVER_2D_BASE_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_library(finite_element_solver_2d_base_lib INTERFACE)
target_sources(finite_element_solver_2d_base_lib INTERFACE
    boundary_condition_first_kind_2d.hpp
    boundary_condition_second_kind_2d.hpp
    boundary_conditions_2d.hpp
    finite_element_matrix_2d.hpp
    index_initializer.hpp
    indexator_base.hpp
    matrix_separator_base.hpp
    mesh_runner_types.hpp
    right_part_2d.hpp
    shift_initializer.hpp
    solution_2d.hpp
    solvers_utils.hpp
)
target_include_directories(finite_element_solver_2d_base_lib INTERFACE
    ${FINITE_ELEMENT_SOLVER_2D_BASE_LIB_DIR}
    ${CONAN_INCLUDE_DIRS_EIGEN}
)
target_link_libraries(finite_element_solver_2d_base_lib INTERFACE
    mesh_2d_lib
)