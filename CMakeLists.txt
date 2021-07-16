cmake_minimum_required(VERSION 3.10)
project(nonlocalMPI)

set(INCLUDES)
set(LINKS)

add_subdirectory(include/metamath)
add_subdirectory(include/mesh)
add_subdirectory(include/solvers)

set(CMAKE_CXX_STANDARD 17)

find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_compile_options(-O2)

add_subdirectory(src/1d/thermal)
add_subdirectory(src/2d/thermal)
