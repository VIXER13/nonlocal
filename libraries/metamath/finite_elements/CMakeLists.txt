cmake_minimum_required(VERSION 3.16)
project(finite_elements)

add_subdirectory(base)
add_subdirectory(finite_elements_1d)
add_subdirectory(finite_elements_2d)

set(FINITE_ELEMENTS_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_library(finite_elements_lib INTERFACE)
target_sources(finite_elements_lib INTERFACE finite_elements.hpp)
target_include_directories(finite_elements_lib INTERFACE ${FINITE_ELEMENT_LIB_DIR})
target_link_libraries(finite_elements_lib INTERFACE
    finite_elements_base_lib
    finite_elements_1d_lib
    finite_elements_2d_lib
)