cmake_minimum_required(VERSION 3.16)
project(symbolic_utils)

add_library(symbolic_utils_lib INTERFACE)
target_sources(
    symbolic_utils_lib
    INTERFACE #
        symbolic_utils.hpp
)
target_include_directories(
    symbolic_utils_lib
    INTERFACE #
        ${LIBRARIES_DIR}
)
