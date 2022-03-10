cmake_minimum_required(VERSION 3.16)
project(mesh_1d)

set(MESH_1D_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_library(mesh_1d_lib INTERFACE)
target_sources(mesh_1d_lib INTERFACE mesh_1d.hpp)
target_include_directories(mesh_1d_lib INTERFACE ${MESH_1D_LIB_DIR})