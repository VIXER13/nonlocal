cmake_minimum_required(VERSION 3.16)
project(linear)

add_library(linear_lib INTERFACE)
target_sources(
    linear_lib
    INTERFACE #
        linear.hpp
)
target_include_directories(
    linear_lib
    INTERFACE #
        ${LIBRARIES_DIR}
)
