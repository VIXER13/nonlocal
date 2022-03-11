cmake_minimum_required(VERSION 3.16)
project(finite_element)

add_subdirectory(geometry)
add_subdirectory(quadrature)
add_subdirectory(element_base)
add_subdirectory(element_1d)
add_subdirectory(element_2d)

set(FINITE_ELEMENT_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_library(finite_element_lib INTERFACE)
target_sources(finite_element_lib INTERFACE finite_element.hpp)
target_include_directories(finite_element_lib INTERFACE ${FINITE_ELEMENT_LIB_DIR})
target_link_libraries(finite_element_lib INTERFACE finite_element_1d_lib finite_element_2d_lib)