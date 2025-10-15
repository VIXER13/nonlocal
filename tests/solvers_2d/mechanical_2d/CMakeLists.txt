cmake_minimum_required(VERSION 3.16)

project(mechanical_2d_test)


add_library(mechanical_2d_test_lib OBJECT saint_venant.cpp)
target_include_directories(
    mechanical_2d_test_lib
    PUBLIC #
        ${LIBRARIES_DIR}
        
        ${PROJECT_ROOT}/build
)

target_link_libraries(
    mechanical_2d_test_lib
    PRIVATE #
        finite_element_solver_2d_lib
        structural_solver_2d_lib
        Boost::ut
)
