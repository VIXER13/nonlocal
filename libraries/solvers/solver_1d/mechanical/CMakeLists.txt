cmake_minimum_required(VERSION 3.16)
project(mechanical_equation_solver_1d)

add_library(mechanical_equation_solver_1d_lib INTERFACE)
target_sources(
    mechanical_equation_solver_1d_lib
    INTERFACE #
		mechanical_parameters_1d.hpp
        stiffness_matrix_assembler.hpp
        mass_matrix_assembler.hpp
        init_problem_settings.hpp
        spring_condition_1d.hpp
        mechanical_boundary_conditions_1d.hpp
        stationary_mechanical_equation_solver_1d.hpp
        harmonic_mechanical_equation_solver_1d.hpp
        mechanical_equation_solution_1d.hpp
)
target_include_directories(
    mechanical_equation_solver_1d_lib
    INTERFACE #
        ${LIBRARIES_DIR}
)
target_link_libraries(
    mechanical_equation_solver_1d_lib
    INTERFACE #
        mesh_1d_lib
        finite_element_solver_1d_base_lib
)

if(ENABLE_OMP)
    target_link_libraries(
        mechanical_equation_solver_1d_lib
        INTERFACE #
            OpenMP::OpenMP_CXX
    )
endif()
