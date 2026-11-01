cmake_minimum_required(VERSION 3.16)

project(solvers_base_tests)

add_library(solvers_base_test_lib OBJECT coefficient.cpp)

target_link_libraries(
    solvers_base_test_lib
    PRIVATE #
        assembler_base_lib
        Boost::ut
)
