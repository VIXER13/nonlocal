cmake_minimum_required(VERSION 3.17)
project(symbolic_functions)

set(SYMBOLIC_FUNCTIONS_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_library(symbolic_functions_lib INTERFACE)
target_sources(symbolic_functions_lib INTERFACE symbolic_functions.hpp)
target_include_directories(symbolic_functions_lib INTERFACE ${SYMBOLIC_FUNCTIONS_LIB_DIR})