cmake_minimum_required(VERSION 3.16)
project(metamath_types)

add_library(metamath_types_lib INTERFACE)
target_sources(metamath_types_lib INTERFACE metamath_types.hpp)
target_include_directories(metamath_types_lib INTERFACE ${LIBRARIES_DIR})