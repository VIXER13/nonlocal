cmake_minimum_required(VERSION 3.16)

project(metamath_instantiations)

set(METAMATH_PCH "${CMAKE_CURRENT_SOURCE_DIR}/pch.hpp" CACHE FILEPATH "Path to the shared metamath PCH header")

add_library(metamath_instantiations_lib STATIC metamath_instantiations.cpp)
target_include_directories(metamath_instantiations_lib PUBLIC ${LIBRARIES_DIR})
target_link_libraries(metamath_instantiations_lib PUBLIC metamath_lib)
