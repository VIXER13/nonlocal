cmake_minimum_required(VERSION 3.16)

set(LIBRARIES_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(config)
add_subdirectory(constants)
add_subdirectory(logger)
add_subdirectory(math_expression)
add_subdirectory(mesh)
add_subdirectory(metamath)
add_subdirectory(parallel)
add_subdirectory(solvers)
