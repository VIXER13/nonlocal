cmake_minimum_required(VERSION 3.16)
project(finite_element_solver_2d_base)

add_library(finite_element_solver_2d_base_lib INTERFACE)
target_sources(
    finite_element_solver_2d_base_lib
    INTERFACE #
        boundary_condition_first_kind_2d.hpp
        boundary_condition_second_kind_2d.hpp
        boundary_conditions_2d.hpp
        finite_element_matrix.hpp
        indices_initializer.hpp
        integrator.hpp
        matrix_assembler.hpp
        problem_settings.hpp
        right_part_2d.hpp
        solution_2d.hpp
        solvers_utils.hpp
)
target_include_directories(
    finite_element_solver_2d_base_lib
    INTERFACE #
        ${LIBRARIES_DIR}
)
target_link_libraries(
    finite_element_solver_2d_base_lib
    INTERFACE #
        assembler_base_lib
        mesh_2d_lib
)

if(ENABLE_OMP)
    target_link_libraries(
        finite_element_solver_2d_base_lib
        INTERFACE #
            OpenMP::OpenMP_CXX
    )
endif()
