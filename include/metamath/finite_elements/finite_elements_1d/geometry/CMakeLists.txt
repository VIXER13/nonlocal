cmake_minimum_required(VERSION 3.16)
project(finite_elements_geometry_1d)

set(FINITE_ELEMENTS_GEOMETRY_1D_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_library(finite_elements_geometry_1d_lib INTERFACE)
target_sources(finite_elements_geometry_1d_lib INTERFACE geometry_1d.hpp
                                                         geometry_primitives.hpp)
target_include_directories(finite_elements_geometry_1d_lib INTERFACE ${FINITE_ELEMENTS_GEOMETRY_1D_LIB_DIR})