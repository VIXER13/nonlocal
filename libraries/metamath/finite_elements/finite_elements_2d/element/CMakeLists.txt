cmake_minimum_required(VERSION 3.16)
project(elements_2d)

set(ELEMENTS_2D_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_library(elements_2d_lib INTERFACE)
target_sources(elements_2d_lib INTERFACE
    element_2d.hpp
    element_2d_serendipity.hpp
    basis/basis_2d.hpp
)
target_include_directories(elements_2d_lib INTERFACE ${ELEMENTS_2D_LIB_DIR})