cmake_minimum_required(VERSION 3.16)
project(mesh_2d)

add_library(mesh_2d_lib INTERFACE)
target_sources(mesh_2d_lib INTERFACE 
    elements_set.hpp
    mesh_2d.hpp
    mesh_2d_utils.hpp
    mesh_container_2d_utils.hpp
    mesh_container_2d.hpp
    mesh_parser.hpp
    su2_parser.hpp
    vtk_elements_set.hpp
    indexator_base.hpp
)
target_include_directories(mesh_2d_lib INTERFACE ${LIBRARIES_DIR})
target_link_libraries(mesh_2d_lib
    INTERFACE
    metamath_lib
    parallel_lib
    nonlocal_constants_lib
)

if(ENABLE_OMP)
    target_link_libraries(mesh_2d_lib INTERFACE OpenMP::OpenMP_CXX)
endif()