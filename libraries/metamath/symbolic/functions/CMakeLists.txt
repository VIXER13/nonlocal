cmake_minimum_required(VERSION 3.16)
project(symbolic_functions)

add_library(symbolic_functions_lib INTERFACE)
target_sources(
    symbolic_functions_lib
    INTERFACE #
        symbolic_functions.hpp
)
target_include_directories(
    symbolic_functions_lib
    INTERFACE #
        ${LIBRARIES_DIR}
)
