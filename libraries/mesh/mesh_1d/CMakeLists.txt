cmake_minimum_required(VERSION 3.16)
project(mesh_1d)

add_library(mesh_1d_lib INTERFACE)
target_sources(mesh_1d_lib INTERFACE 
    mesh_1d.hpp
    mesh_1d_utils.hpp
)
target_include_directories(mesh_1d_lib INTERFACE ${LIBRARIES_DIR})
target_link_libraries(mesh_1d_lib
    INTERFACE
    metamath_lib
    nonlocal_constants_lib
)