cmake_minimum_required(VERSION 3.16)

project(mechanical_2d)

add_executable(thermomechanical thermomechanical.cpp)
target_include_directories(thermomechanical PUBLIC ".")
target_link_libraries(thermomechanical 
    finite_element_solver_2d_lib
    config_lib
)