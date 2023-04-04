cmake_minimum_required(VERSION 3.16)

project(thermal_2d)

add_executable(thermal_stationary_2d thermal_stationary_2d.cpp)
target_include_directories(thermal_stationary_2d PUBLIC ".")
target_link_libraries(thermal_stationary_2d 
    finite_element_solver_2d_lib
    parallel_utils_lib
    config_lib
)