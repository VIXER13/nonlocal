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
enable_testing()

### Conan 2 integration (CMakeDeps + CMakeToolchain) ###
# The toolchain file is passed from the "conan install" step (see README).
# After that we just rely on find_package with CONFIG mode to load dependencies.
find_package(Eigen3 REQUIRED CONFIG)
find_package(nlohmann_json REQUIRED CONFIG)
find_package(spectra REQUIRED CONFIG)
find_package(ut REQUIRED CONFIG)

# libraries
add_subdirectory(libraries)

# NonLocFEM
add_subdirectory(NonLocFEM)

# unit tests
include(CMakeTools/EmbeddedFiles.cmake)
add_subdirectory(tests)