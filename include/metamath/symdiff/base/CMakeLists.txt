cmake_minimum_required(VERSION 3.17)
project(symdiff_base)

set(SYMDIFF_BASE_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_library(symdiff_base_lib INTERFACE)
target_sources(symdiff_base_lib INTERFACE symdiff_base.hpp)
target_include_directories(symdiff_base_lib INTERFACE ${SYMDIFF_BASE_LIB_DIR})