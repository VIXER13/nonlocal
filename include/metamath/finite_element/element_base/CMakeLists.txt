cmake_minimum_required(VERSION 3.16)
project(finite_element_base)

set(FINITE_ELEMENT_BASE_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_library(finite_element_base_lib INTERFACE)
target_sources(finite_element_base_lib INTERFACE element_base.hpp
                                                 element_integrate_base.hpp)
target_include_directories(finite_element_base_lib INTERFACE ${FINITE_ELEMENT_BASE_LIB_DIR})