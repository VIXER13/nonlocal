cmake_minimum_required(VERSION 3.17)
project(finite_elements_1d)

add_subdirectory(geometry)
add_subdirectory(quadrature)
add_subdirectory(element)

set(FINITE_ELEMENTS_LIB_1D_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_library(finite_elements_1d_lib INTERFACE)
target_sources(finite_elements_1d_lib INTERFACE finite_elements_1d.hpp)
target_include_directories(finite_elements_1d_lib INTERFACE ${FINITE_ELEMENTS_LIB_1D_DIR})
target_link_libraries(finite_elements_1d_lib INTERFACE
    finite_elements_geometry_1d_lib
    quadrature_1d_lib
    elements_1d_lib
)