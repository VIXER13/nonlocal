cmake_minimum_required(VERSION 3.16)
project(constants)

set(CONSTANTS_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_library(constants_lib INTERFACE)
target_sources(constants_lib INTERFACE constants.hpp)
target_include_directories(constants_lib INTERFACE ${CONSTANTS_LIB_DIR})