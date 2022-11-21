cmake_minimum_required(VERSION 3.16)
project(metamath_functions)

set(METAMATH_FUNCTIONS_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_library(metamath_functions_lib INTERFACE)
target_sources(metamath_functions_lib INTERFACE metamath_functions.hpp)
target_include_directories(metamath_functions_lib INTERFACE ${FUNCTIONS_LIB_DIR})