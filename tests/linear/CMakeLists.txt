cmake_minimum_required(VERSION 3.16)

project(linear_test)

add_library(
    linear_test_lib
    OBJECT
        fixed_matrix.cpp
        norm.cpp
        operators.cpp
        scalar_product.cpp
        sparse_matrix.cpp
)

target_include_directories(
    linear_test_lib
    PUBLIC #
        ${LIBRARIES_DIR}
        
)
target_link_libraries(
    linear_test_lib
    PRIVATE #
        linear_lib
        Boost::ut
)
