cmake_minimum_required(VERSION 3.16)

project(math_expression)

add_library(math_expression_lib STATIC 
    math_expression.cpp
    utils.cpp
)
target_include_directories(math_expression_lib PUBLIC ${LIBRARIES_DIR})