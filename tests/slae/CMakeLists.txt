cmake_minimum_required(VERSION 3.16)

project(slae_test)

add_library(
    slae_test_lib
    OBJECT
        conjugate_gradient.cpp
        diagonal_preconditioner.cpp
        identity_preconditioner.cpp
        ildlt_preconditioner.cpp
        ilu0_preconditioner.cpp
        stable_biconjugate_gradient.cpp
        gmres.cpp
)

target_include_directories(
    slae_test_lib
    PUBLIC #
        ${LIBRARIES_DIR}
        
)
target_link_libraries(
    slae_test_lib
    PRIVATE #
        linear_lib
        slae_solver_lib
        Boost::ut
)
