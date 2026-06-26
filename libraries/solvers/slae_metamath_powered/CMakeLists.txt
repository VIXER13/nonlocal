cmake_minimum_required(VERSION 3.16)
project(slae_solver_metamath_powered)

add_library(slae_solver_metamath_powered_lib INTERFACE)
target_sources(
    slae_solver_metamath_powered_lib
    INTERFACE #
        conjugate_gradient.hpp
)
target_include_directories(
    slae_solver_metamath_powered_lib
    INTERFACE #
        ${LIBRARIES_DIR}
)
target_link_libraries(
    slae_solver_metamath_powered_lib
    INTERFACE #
        parallel_lib
        linear_lib
)

if(ENABLE_OMP)
    target_link_libraries(
        slae_solver_metamath_powered_lib
        INTERFACE #
            OpenMP::OpenMP_CXX
    )
endif()
