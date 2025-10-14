cmake_minimum_required(VERSION 3.16)

project(thermal_2d_test)


add_library(
    thermal_2d_test_lib
    OBJECT
    boundary_conditions_2d.cpp
    flux_stability.cpp
)
target_include_directories(
    thermal_2d_test_lib
    PUBLIC #
        ${LIBRARIES_DIR}
        
        ${PROJECT_ROOT}/build
)

target_link_libraries(
    thermal_2d_test_lib
    PRIVATE #
        finite_element_solver_2d_lib
        heat_equation_solver_2d_lib
        Boost::ut
)
