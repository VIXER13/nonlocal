cmake_minimum_required(VERSION 3.16)

project(solvers_2d)


FileEmbedAdd(./meshes/plate_10x1_h0.125.su2)
FileEmbedAdd(./meshes/square_1x1_16el_uniform_mesh.su2)

add_subdirectory(thermal_2d)
add_subdirectory(mechanical_2d)

add_library(solvers_2d_test_lib INTERFACE)
target_include_directories(
    solvers_2d_test_lib
    INTERFACE #
        ${LIBRARIES_DIR}
        
        ${PROJECT_ROOT}/build
)
target_link_libraries(
    solvers_2d_test_lib
    INTERFACE #
        thermal_2d_test_lib
        mechanical_2d_test_lib
)
