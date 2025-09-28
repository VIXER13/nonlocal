cmake_minimum_required(VERSION 3.16)

project(thermal_1d_test)

set(INCLUDES $ENV{CPATH})
if (NOT DEFINED NO_CONAN)
    set(INCLUDES ${CONAN_INCLUDE_DIRS_BOOST-EXT-UT})
endif()

add_library(thermal_1d_test_lib OBJECT
    flux_relaxation.cpp
    nonstationary.cpp
    stationary.cpp
)
target_include_directories(thermal_1d_test_lib PUBLIC
    ${LIBRARIES_DIR}
    ${INCLUDES}
)

target_link_libraries(thermal_1d_test_lib
    finite_element_solver_1d_lib
    heat_equation_solver_1d_lib
)
