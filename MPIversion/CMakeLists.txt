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

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_compile_options(-O2)

message(STATUS $ENV{PETSc_LIBRARIES})

add_executable(nonlocal src/structural_long_rect_with_T.cpp)
target_include_directories(nonlocal PUBLIC ${MPI_C_INCLUDE_DIRS})
target_link_libraries(nonlocal ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES} ${PETSc_LIBRARIES} metamath_lib mesh_lib finite_element_solvers_lib)
