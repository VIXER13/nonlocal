cmake_minimum_required(VERSION 3.16)
project(symbolic)

add_subdirectory(base)
add_subdirectory(functions)
add_subdirectory(utils)

set(SYMBOLIC_LIB_DIR symbolic)

add_library(symbolic_lib INTERFACE)
target_sources(symbolic_lib INTERFACE symbolic.hpp)
target_include_directories(symbolic_lib INTERFACE ${SYMBOLIC_LIB_DIR})
target_link_libraries(symbolic_lib INTERFACE symbolic_base_lib symbolic_functions_lib symbolic_utils_lib)
