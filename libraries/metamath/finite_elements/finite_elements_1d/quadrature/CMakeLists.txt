cmake_minimum_required(VERSION 3.16)
project(quadrature_1d)

set(QUADRATURE_1D_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_library(quadrature_1d_lib INTERFACE)
target_sources(quadrature_1d_lib INTERFACE
    quadrature_1d.hpp
    gaussian_quadrature.hpp
)
target_include_directories(quadrature_1d_lib INTERFACE ${QUADRATURE_1D_LIB_DIR})