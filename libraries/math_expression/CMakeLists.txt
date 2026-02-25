cmake_minimum_required(VERSION 3.16)

project(math_expression)

add_library(math_expression_lib STATIC utils.cpp tokenizer.cpp)
target_include_directories(
    math_expression_lib
    PUBLIC #
        ${LIBRARIES_DIR}
)
target_link_libraries(
    math_expression_lib
    PUBLIC #
)

if(ENABLE_OMP)
    target_link_libraries(
        math_expression_lib
        PUBLIC #
            OpenMP::OpenMP_CXX
    )
endif()
