cmake_minimum_required(VERSION 3.16)
project(elements_1d)

set(ELEMENTS_1D_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_library(elements_1d_lib INTERFACE derivative_element_basis_1d.hpp)
target_sources(elements_1d_lib INTERFACE
    element_1d.hpp
    element_1d_integrate.hpp
    basis/lagrangian_elements_1d.hpp
)
target_include_directories(elements_1d_lib INTERFACE ${ELEMENTS_1D_LIB_DIR})