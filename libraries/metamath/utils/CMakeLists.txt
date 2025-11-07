cmake_minimum_required(VERSION 3.16)
project(metamath_utils)

add_library(metamath_utils_lib INTERFACE)
target_sources(
    metamath_utils_lib
    INTERFACE #
        metamath_utils.hpp
)
target_include_directories(
    metamath_utils_lib
    INTERFACE #
        ${LIBRARIES_DIR}
)
