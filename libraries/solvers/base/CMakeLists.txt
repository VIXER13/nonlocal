cmake_minimum_required(VERSION 3.16)
project(assembler_base)

add_library(assembler_base_lib INTERFACE)
target_sources(assembler_base_lib INTERFACE
    equation_parameters.hpp
    utils.hpp
)
target_include_directories(assembler_base_lib INTERFACE ${LIBRARIES_DIR})
target_link_libraries(assembler_base_lib INTERFACE)