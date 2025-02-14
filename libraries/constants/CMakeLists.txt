cmake_minimum_required(VERSION 3.16)

project(nonlocal_constants)

add_library(nonlocal_constants_lib INTERFACE)
target_sources(nonlocal_constants_lib INTERFACE nonlocal_constants.hpp)
target_include_directories(nonlocal_constants_lib INTERFACE ${LIBRARIES_DIR})