cmake_minimum_required(VERSION 3.16)

project(solvers_2d)

set(INCLUDES $ENV{CPATH})
if (NOT DEFINED NO_CONAN)
    set(INCLUDES ${CONAN_INCLUDE_DIRS_BOOST-EXT-UT})
endif()

FileEmbedAdd(./meshes/plate_10x1_h0.125.su2)

add_subdirectory(thermal_2d)

add_library(solvers_2d_test_lib INTERFACE)
target_include_directories(solvers_2d_test_lib INTERFACE
    ${LIBRARIES_DIR}
    ${INCLUDES}
    ${PROJECT_ROOT}/build
)
target_link_libraries(solvers_2d_test_lib
    INTERFACE
    thermal_2d_test_lib
)
