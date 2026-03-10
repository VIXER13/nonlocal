cmake_minimum_required(VERSION 3.16)
project(structural_solver_2d)

add_library(structural_solver_2d_lib INTERFACE)
target_sources(
    structural_solver_2d_lib
    INTERFACE #
        elastic_parameters.hpp
        equilibrium_equation_2d.hpp
        evaluate_mechanical_parameters.hpp
        init_problem_settings.hpp
        mechanical_boundary_conditions_2d.hpp
        mechanical_parameters_2d.hpp
        mechanical_solution_2d.hpp
        stiffness_matrix_2d.hpp
        temperature_condition_2d.hpp
        thermal_expansion_parameters.hpp
)
target_include_directories(
    structural_solver_2d_lib
    INTERFACE #
        ${LIBRARIES_DIR}
)
target_link_libraries(
    structural_solver_2d_lib
    INTERFACE #
        finite_element_solver_2d_base_lib
        mesh_2d_lib
)

if(ENABLE_OMP)
    target_link_libraries(
        structural_solver_2d_lib
        INTERFACE #
            OpenMP::OpenMP_CXX
    )
endif()
