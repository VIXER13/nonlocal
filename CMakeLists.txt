cmake_minimum_required(VERSION 3.16)

project(nonlocal)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

if (MPI_BUILD)
    find_package(MPI)
    if(MPI_FOUND)
        add_definitions(-DMPI_BUILD=1)
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
enable_testing()

if (NOT DEFINED NO_CONAN)
    if(NOT EXISTS "${CMAKE_BINARY_DIR}/conanbuildinfo.cmake")
        message(STATUS "Downloading conanbuildinfo.cmake from https://github.com/conan-io/cmake-conan")
        file(DOWNLOAD "https://raw.githubusercontent.com/conan-io/cmake-conan/master/conan.cmake" "${CMAKE_BINARY_DIR}/conanbuildinfo.cmake")
        include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
        conan_cmake_run(CONANFILE conanfile.txt BASIC_SETUP CMAKE_TARGETS KEEP_RPATHS)
    else()
        include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
    endif()
    conan_basic_setup()
endif()

# libraries
add_subdirectory(libraries/logger)
add_subdirectory(libraries/parallel)
add_subdirectory(libraries/metamath)
add_subdirectory(libraries/constants)
add_subdirectory(libraries/mesh/mesh_1d)
add_subdirectory(libraries/mesh/mesh_2d)
add_subdirectory(libraries/solvers/slae)
add_subdirectory(libraries/solvers/solver_1d)
add_subdirectory(libraries/solvers/solver_2d)
add_subdirectory(libraries/config)

# NonLocFEM
add_subdirectory(NonLocFEM)

# unit tests
include(CMakeTools/EmbeddedFiles.cmake)
add_subdirectory(tests)