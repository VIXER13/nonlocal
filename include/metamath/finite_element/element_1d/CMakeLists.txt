cmake_minimum_required(VERSION 3.16)
project(finite_element_1d)

set(FINITE_ELEMENT_1D_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_library(finite_element_1d_lib INTERFACE)
target_sources(finite_element_1d_lib INTERFACE element_1d_integrate.hpp
                                               basis/basis.hpp)
target_include_directories(finite_element_1d_lib INTERFACE ${FINITE_ELEMENT_1D_LIB_DIR})
target_link_libraries(finite_element_1d_lib INTERFACE finite_element_base_lib
                                                      finite_element_quadrature_lib)