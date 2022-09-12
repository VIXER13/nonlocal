cmake_minimum_required(VERSION 3.16)
project(finite_elements_geometry_2d)

set(FINITE_ELEMENTS_GEOMETRY_2D_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_library(finite_elements_geometry_2d_lib INTERFACE)
target_sources(finite_elements_geometry_2d_lib INTERFACE geometry_2d.hpp
                                                         geometric_primitives/geometric_primitives_2d.hpp)
target_include_directories(finite_elements_geometry_2d_lib INTERFACE ${FINITE_ELEMENTS_GEOMETRY_2D_LIB_DIR})