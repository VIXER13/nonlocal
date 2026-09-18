cmake_minimum_required(VERSION 3.16)
project(quadrature_2d)

add_library(quadrature_2d_lib INTERFACE)
target_sources(
    quadrature_2d_lib
    INTERFACE #
        make_quadrature_2d.hpp
)
target_include_directories(
    quadrature_2d_lib
    INTERFACE #
        ${LIBRARIES_DIR}
)
