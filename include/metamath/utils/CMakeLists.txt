cmake_minimum_required(VERSION 3.16)
project(metamath_utils)

set(METAMATH_UTILS_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_library(metamath_utils_lib INTERFACE)
target_sources(metamath_utils_lib INTERFACE metamath_utils.hpp)
target_include_directories(metamath_utils_lib INTERFACE ${METAMATH_UTILS_LIB_DIR})