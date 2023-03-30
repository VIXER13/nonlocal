cmake_minimum_required(VERSION 3.16)

project(thermal_1d)

add_executable(thermal_stationary_1d thermal_stationary_1d.cpp)
target_include_directories(thermal_stationary_1d PUBLIC ".")
target_link_libraries(thermal_stationary_1d 
    finite_element_solver_1d_lib
    config_lib
)

add_executable(thermal_nonstationary_1d thermal_nonstationary_1d.cpp)
target_include_directories(thermal_nonstationary_1d PUBLIC ".")
target_link_libraries(thermal_nonstationary_1d 
    finite_element_solver_1d_lib
    config_lib
)