cmake_minimum_required(VERSION 3.16)
project(symdiff)

add_subdirectory(base)
add_subdirectory(functions)

set(SYMDIFF_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_library(symdiff_lib INTERFACE)
target_sources(symdiff_lib INTERFACE symdiff.hpp)
target_include_directories(symdiff_lib INTERFACE ${SYMDIFF_LIB_DIR})
target_link_libraries(symdiff_lib INTERFACE symdiff_base_lib symdiff_functions_lib)