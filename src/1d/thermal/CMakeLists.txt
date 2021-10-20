cmake_minimum_required(VERSION 3.10)

project(thermal_1d)

add_executable(test test.cpp)
target_include_directories(test PUBLIC ".")
target_link_libraries(test metamath_lib mesh_1d_lib finite_element_solver_1d_lib MPI_utils_lib)

add_executable(nonlocal_stationary thermal_stationary_1d.cpp)
target_include_directories(nonlocal_stationary PUBLIC ".")
target_link_libraries(nonlocal_stationary metamath_lib mesh_1d_lib finite_element_solver_1d_lib MPI_utils_lib)

add_executable(nonlocal_nonstationary thermal_nonstationary_1d.cpp)
target_include_directories(nonlocal_nonstationary PUBLIC ".")
target_link_libraries(nonlocal_nonstationary metamath_lib mesh_1d_lib finite_element_solver_1d_lib MPI_utils_lib)