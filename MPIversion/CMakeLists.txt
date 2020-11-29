cmake_minimum_required(VERSION 3.10)
project(nonlocalMPI)

find_package(MPI REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(PETSc REQUIRED IMPORTED_TARGET petsc)

add_subdirectory(include/metamath)
add_subdirectory(include/mesh)
add_subdirectory(include/solvers)

set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

set(CMAKE_CXX_STANDARD 17)

add_executable(nonlocal src/main.cpp)
target_include_directories(nonlocal PUBLIC ${MPI_C_INCLUDE_DIRS})
target_link_libraries(nonlocal ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES} ${PETSc_LIBRARIES} metamath_lib mesh_lib finite_element_solvers_lib)