cmake_minimum_required(VERSION 3.16)
project(slae_solver)

add_library(slae_solver_lib INTERFACE)
target_sources(
    slae_solver_lib
    INTERFACE #
        preconditioner_base.hpp
        identity_preconditioner.hpp
        diagonal_preconditioner.hpp
        ildlt_preconditioner.hpp
        ilu0_preconditioner.hpp
        ldlt_solver.hpp
        lu_solver.hpp
        iterative_solver_base.hpp
        conjugate_gradient.hpp
        stable_biconjugate_gradient.hpp
        gmres.hpp
)
target_include_directories(
    slae_solver_lib
    INTERFACE #
        ${LIBRARIES_DIR}
)
target_link_libraries(
    slae_solver_lib
    INTERFACE #
        parallel_lib
        linear_lib
)

if(ENABLE_OMP)
    target_link_libraries(
        slae_solver_lib
        INTERFACE #
            OpenMP::OpenMP_CXX
    )
endif()
