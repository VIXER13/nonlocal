cmake_minimum_required(VERSION 3.16)

project(NonLocFEM)

add_executable(NonLocFEM 
    main.cpp
)

target_include_directories(NonLocFEM PUBLIC 
    "."
)

target_link_libraries(NonLocFEM
    config_lib
    finite_element_solver_1d_lib
    finite_element_solver_2d_lib
)