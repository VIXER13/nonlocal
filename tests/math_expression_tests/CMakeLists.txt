cmake_minimum_required(VERSION 3.16)

project(math_expression_tests)


add_library(math_expression_test_lib OBJECT math_expression_test.cpp)
target_include_directories(
    math_expression_test_lib
    PUBLIC #
        ${LIBRARIES_DIR}
        
)
target_link_libraries(
    math_expression_test_lib
    PRIVATE #
        math_expression_lib
        Boost::ut
)
