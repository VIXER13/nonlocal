cmake_minimum_required(VERSION 3.17)
project(finite_elements_2d)

add_subdirectory(geometry)
#add_subdirectory(quadrature)
add_subdirectory(element)

set(FINITE_ELEMENTS_LIB_2D_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_library(finite_elements_2d_lib INTERFACE)
target_sources(finite_elements_2d_lib INTERFACE finite_elements_2d.hpp)
target_include_directories(finite_elements_2d_lib INTERFACE ${FINITE_ELEMENTS_LIB_2D_DIR})
target_link_libraries(finite_elements_2d_lib INTERFACE
    finite_elements_geometry_2d_lib
    elements_2d_lib
)