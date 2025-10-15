cmake_minimum_required(VERSION 3.16)
project(finite_elements_geometry_1d)

add_library(finite_elements_geometry_1d_lib INTERFACE)
target_sources(
    finite_elements_geometry_1d_lib
    INTERFACE #
        geometry_1d.hpp
        geometry_primitives.hpp
)
target_include_directories(
    finite_elements_geometry_1d_lib
    INTERFACE #
        ${LIBRARIES_DIR}
)
