cmake_minimum_required(VERSION 3.16)
project(finite_element_geometry)

set(FINITE_ELEMENT_GEOMETRY_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_library(finite_element_geometry_lib INTERFACE)
target_sources(finite_element_geometry_lib INTERFACE geometry_1d.hpp
                                                     geometry_2d.hpp)
target_include_directories(finite_element_geometry_lib INTERFACE ${FINITE_ELEMENT_GEOMETRY_LIB_DIR})
target_link_libraries(finite_element_geometry_lib INTERFACE symdiff_lib)