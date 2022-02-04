cmake_minimum_required(VERSION 3.10)

project(mechanical_2d)

add_executable(mechanical_quad mechanical_quad.cpp)
target_include_directories(mechanical_quad PUBLIC ".")
target_link_libraries(mechanical_quad metamath_lib mesh_2d_lib slae_solver_lib finite_element_solver_2d_lib MPI_utils_lib)