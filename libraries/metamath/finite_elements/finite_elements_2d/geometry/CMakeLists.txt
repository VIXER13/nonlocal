cmake_minimum_required(VERSION 3.16)
project(finite_elements_geometry_2d)

add_library(finite_elements_geometry_2d_lib INTERFACE)
target_sources(
    finite_elements_geometry_2d_lib
    INTERFACE #
        geometry_2d.hpp
        geometric_primitives/geometric_primitives_2d.hpp
)
target_include_directories(
    finite_elements_geometry_2d_lib
    INTERFACE #
        ${LIBRARIES_DIR}
)
