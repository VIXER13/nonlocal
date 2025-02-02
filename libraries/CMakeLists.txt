cmake_minimum_required(VERSION 3.16)

set(LIBRARIES_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(config)
add_subdirectory(constants)
add_subdirectory(logger)
add_subdirectory(mesh/mesh_1d)
add_subdirectory(mesh/mesh_2d)
add_subdirectory(metamath)
add_subdirectory(parallel)
add_subdirectory(solvers/slae)
add_subdirectory(solvers/solver_1d)
add_subdirectory(solvers/solver_2d)