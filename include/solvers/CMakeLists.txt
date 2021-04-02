cmake_minimum_required(VERSION 3.17)
project(finite_element_solvers)

set(FINITE_ELEMENT_SOLVERS_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_library(finite_element_solvers_lib INTERFACE)
target_sources(finite_element_solvers_lib INTERFACE finite_element_solver_base.hpp
                                                    heat_equation_solver.hpp
                                                    structural_solver.hpp)
target_include_directories(finite_element_solvers_lib INTERFACE ${FINITE_ELEMENT_SOLVERS_LIB_DIR} mesh_lib)