cmake_minimum_required(VERSION 3.16)

project(math_expression)

set(INCLUDES $ENV{CPATH})
if (NOT DEFINED NO_CONAN)
    set(INCLUDES ${CONAN_INCLUDE_DIRS_EIGEN})
endif()

add_library(math_expression_lib STATIC 
    math_expression.cpp
    utils.cpp
)

target_include_directories(math_expression_lib PUBLIC 
    "." 
    ${INCLUDES}
)