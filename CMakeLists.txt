cmake_minimum_required(VERSION 3.16)

project(nonlocalMPI)

if (WITH_MPI)
    find_package(MPI)
    if(MPI_FOUND)
        add_definitions(-DMPI_USE=1)
#        set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
#        set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
    endif()
endif()

find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(CMAKE_CXX_STANDARD 20)
add_compile_options(-O2)

#libraries
add_subdirectory(include/parallel_utils)
add_subdirectory(include/metamath)
add_subdirectory(include/constants)
add_subdirectory(include/mesh/mesh_1d)
add_subdirectory(include/mesh/mesh_2d)
add_subdirectory(include/solvers/slae)
add_subdirectory(include/solvers/solver_1d)
add_subdirectory(include/solvers/solver_2d)

#tests
add_subdirectory(src/1d/thermal)
add_subdirectory(src/2d/thermal)
add_subdirectory(src/2d/mechanical)
add_subdirectory(src/2d/thermomechanical)