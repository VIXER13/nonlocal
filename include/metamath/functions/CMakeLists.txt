cmake_minimum_required(VERSION 3.16)
project(functions)

set(FUNCTIONS_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_library(functions_lib INTERFACE)
target_sources(functions_lib INTERFACE functions.hpp)
target_include_directories(functions_lib INTERFACE ${FUNCTIONS_LIB_DIR})
target_link_libraries(functions_lib INTERFACE constants_lib)