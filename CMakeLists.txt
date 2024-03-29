cmake_minimum_required(VERSION 3.16)

project(nonlocal)

if (WITH_MPI)
    find_package(MPI)
    if(MPI_FOUND)
        add_definitions(-DMPI_USED=1)
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

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

# libraries
add_subdirectory(libraries/logger)
add_subdirectory(libraries/parallel_utils)
add_subdirectory(libraries/metamath)
add_subdirectory(libraries/constants)
add_subdirectory(libraries/mesh/mesh_1d)
add_subdirectory(libraries/mesh/mesh_2d)
add_subdirectory(libraries/solvers/slae)
add_subdirectory(libraries/solvers/solver_1d)
add_subdirectory(libraries/solvers/solver_2d)
add_subdirectory(libraries/config)

# NonLocFEM
add_subdirectory(src)

# unit tests
include(CMakeTools/EmbeddedFiles.cmake)
add_subdirectory(tests)