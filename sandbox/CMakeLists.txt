cmake_minimum_required(VERSION 3.16)

project(nonlocal_sandbox)

set(PROJECT_ROOT ${CMAKE_SOURCE_DIR})
message(STATUS "PROJECT_ROOT: ${PROJECT_ROOT}")

add_executable(sandbox nonlocal_sandbox.cpp)
target_include_directories(
    sandbox
    PRIVATE #
        ${LIBRARIES_DIR}
)
target_link_libraries(
    sandbox
    PRIVATE #
        OpenMP::OpenMP_CXX
        Eigen3::Eigen
        slae_solver_lib
)
