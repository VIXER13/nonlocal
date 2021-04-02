cmake_minimum_required(VERSION 3.17)
project(metamath)

add_subdirectory(constants)
add_subdirectory(functions)
add_subdirectory(symdiff)
add_subdirectory(finite_element)

set(METAMATH_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_library(metamath_lib INTERFACE)
target_sources(metamath_lib INTERFACE metamath.hpp)
target_include_directories(metamath_lib INTERFACE ${METAMATH_LIB_DIR})
target_link_libraries(metamath_lib INTERFACE constants_lib functions_lib symdiff_lib finite_element_lib)