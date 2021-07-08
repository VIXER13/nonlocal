cmake_minimum_required(VERSION 3.10)

project(thermal_1d)

add_executable(nonlocal_stationary stationary.cpp)
target_include_directories(nonlocal_stationary PUBLIC ".")
target_link_libraries(nonlocal_stationary metamath_lib mesh_lib finite_element_solvers_lib)

add_executable(nonlocal_nonstationary nonstationary.cpp)
target_include_directories(nonlocal_nonstationary PUBLIC ".")
target_link_libraries(nonlocal_nonstationary metamath_lib mesh_lib finite_element_solvers_lib)