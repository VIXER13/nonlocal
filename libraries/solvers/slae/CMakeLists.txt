cmake_minimum_required(VERSION 3.16)
project(slae_solver_1d)

set(SLAE_SOLVER_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_library(slae_solver_lib INTERFACE)
target_sources(slae_solver_lib INTERFACE conjugate_gradient.hpp conjugate_gradient_old.hpp)
target_include_directories(slae_solver_lib INTERFACE ${SLAE_SOLVER_LIB_DIR})
target_link_libraries(slae_solver_lib INTERFACE 
    parallel_utils_lib
)