cmake_minimum_required(VERSION 3.16)
project(elements_1d)

add_library(elements_1d_lib INTERFACE)
target_sources(
    elements_1d_lib
    INTERFACE #
        element_1d.hpp
        element_1d_integrate.hpp
        basis/lagrangian_elements_1d.hpp
)
target_include_directories(
    elements_1d_lib
    INTERFACE #
        ${LIBRARIES_DIR}
)
