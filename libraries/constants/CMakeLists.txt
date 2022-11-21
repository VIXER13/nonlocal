cmake_minimum_required(VERSION 3.17)
project(nonlocal_constants)

set(NONLOCAL_CONSTANTS_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_library(nonlocal_constants_lib INTERFACE)
target_sources(nonlocal_constants_lib INTERFACE nonlocal_constants.hpp)
target_include_directories(nonlocal_constants_lib INTERFACE ${NONLOCAL_CONSTANTS_LIB_DIR})