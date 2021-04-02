cmake_minimum_required(VERSION 3.17)
project(finite_element_2d)

set(FINITE_ELEMENT_2D_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_library(finite_element_2d_lib INTERFACE)
target_sources(finite_element_2d_lib INTERFACE element_2d_serendipity.hpp
                                               element_2d_integrate.hpp
                                               basis/basis.hpp)
target_include_directories(finite_element_2d_lib INTERFACE ${FINITE_ELEMENT_2D_LIB_DIR})
target_link_libraries(finite_element_2d_lib INTERFACE finite_element_base_lib
                                                      finite_element_quadrature_lib)