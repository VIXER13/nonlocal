cmake_minimum_required(VERSION 3.16)
project(symbolic_base)

set(SYMBOLIC_BASE_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_library(symbolic_base_lib INTERFACE)
target_sources(symbolic_base_lib INTERFACE symbolic_base.hpp)
target_include_directories(symbolic_base_lib INTERFACE ${SYMBOLIC_BASE_LIB_DIR})