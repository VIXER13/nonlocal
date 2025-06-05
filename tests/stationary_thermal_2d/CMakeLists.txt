cmake_minimum_required(VERSION 3.16)

project(stationary_thermal_2d)

set(INCLUDES $ENV{CPATH})
if (NOT DEFINED NO_CONAN)
    set(INCLUDES ${CONAN_INCLUDE_DIRS_BOOST-EXT-UT})
endif()

add_library(stationary_thermal_2d_lib OBJECT 
    boundary_conditions_2d.cpp
    init_2d_mesh.cpp
)

target_include_directories(stationary_thermal_2d_lib PUBLIC
    ${LIBRARIES_DIR}
    ${INCLUDES}
)

target_link_libraries(stationary_thermal_2d_lib
    heat_equation_solver_2d_lib
)
