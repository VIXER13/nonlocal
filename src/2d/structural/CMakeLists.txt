cmake_minimum_required(VERSION 3.10)

project(structural_2d)

add_executable(nonlocal_saint_venant structural_saint_venant.cpp)
target_include_directories(nonlocal_saint_venant PUBLIC ".")
target_link_libraries(nonlocal_saint_venant metamath_lib mesh_lib finite_element_solvers_lib)
if (WITH_MPI)
    target_link_libraries(nonlocal_saint_venant MPI_utils_lib)
endif()

add_executable(nonlocal_ellipse structural_ellipse.cpp)
target_include_directories(nonlocal_ellipse PUBLIC ".")
target_link_libraries(nonlocal_ellipse metamath_lib mesh_lib finite_element_solvers_lib)