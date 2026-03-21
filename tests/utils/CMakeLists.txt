cmake_minimum_required(VERSION 3.16)
project(unit_tests_utils)

add_library(unit_tests_utils_lib INTERFACE)
target_sources(
    unit_tests_utils_lib
    INTERFACE #
        error.hpp
)
target_include_directories(
    unit_tests_utils_lib
    INTERFACE #
        ${LIBRARIES_DIR}
)
target_link_libraries(
    unit_tests_utils_lib
    INTERFACE #
        mesh_2d_lib
)
