cmake_minimum_required(VERSION 3.16)
project(metamath)

add_subdirectory(functions)
add_subdirectory(symbolic)
add_subdirectory(finite_elements)
add_subdirectory(utils)
add_subdirectory(types)

add_library(metamath_lib INTERFACE)
target_sources(
    metamath_lib
    INTERFACE #
        metamath.hpp
)
target_include_directories(
    metamath_lib
    INTERFACE #
        ${LIBRARIES_DIR}
)
target_link_libraries(
    metamath_lib
    INTERFACE
        metamath_types_lib
        metamath_utils_lib
        metamath_functions_lib
        symbolic_lib
        finite_elements_lib
)
