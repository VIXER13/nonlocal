cmake_minimum_required(VERSION 3.16)
project(symbolic)

set(SYMBOLIC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(base)
add_subdirectory(functions)
add_subdirectory(utils)

add_library(symbolic_lib INTERFACE)
target_sources(
    symbolic_lib
    INTERFACE #
        symbolic.hpp
)
target_include_directories(
    symbolic_lib
    INTERFACE #
        ${LIBRARIES_DIR}
)
target_link_libraries(
    symbolic_lib
    INTERFACE #
        symbolic_base_lib
        symbolic_functions_lib
        symbolic_utils_lib
)
