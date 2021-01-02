cmake_minimum_required(VERSION 3.17)
project(mesh)

set(MESH_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_library(mesh_lib INTERFACE)
target_sources(mesh_lib INTERFACE mesh.hpp)
target_include_directories(mesh_lib INTERFACE ${MESH_LIB_DIR} metamath_lib)