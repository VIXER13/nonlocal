cmake_minimum_required(VERSION 3.16)
project(finite_elements)

add_subdirectory(base)
add_subdirectory(finite_elements_1d)
add_subdirectory(finite_elements_2d)

add_library(finite_elements_lib INTERFACE)
target_sources(finite_elements_lib INTERFACE finite_elements.hpp)
target_include_directories(finite_elements_lib INTERFACE ${LIBRARIES_DIR})
target_link_libraries(finite_elements_lib INTERFACE
    finite_elements_base_lib
    finite_elements_1d_lib
    finite_elements_2d_lib
)