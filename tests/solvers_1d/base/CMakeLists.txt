cmake_minimum_required(VERSION 3.16)

project(solvers_1d_base_tests)

add_library(solvers_1d_base_test_lib OBJECT assemblers.cpp boundaries.cpp stationary.cpp linear_solver.cpp)

target_include_directories(solvers_1d_base_test_lib PRIVATE ${PROJECT_ROOT})

target_link_libraries(
    solvers_1d_base_test_lib
    PRIVATE #
        assembler_base_lib
        finite_element_solver_1d_lib
        Boost::ut
)
