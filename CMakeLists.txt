cmake_minimum_required(VERSION 3.16)

project(nonlocal)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

if(MPI_BUILD)
    find_package(MPI)
    if(MPI_FOUND)
        add_definitions(-DMPI_BUILD=1)
    endif()
endif()

option(ENABLE_OMP "Enable OpenMP parallel regions" ON)
if(ENABLE_OMP)
    # Required when enabled: fail fast if not available.
    find_package(OpenMP REQUIRED)
    message(STATUS "OpenMP enabled")
endif()

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
enable_testing()

# Conan 2 integration (CMakeDeps + CMakeToolchain)
# The toolchain file is passed from the "conan install" step (see README).
# After that we just rely on find_package with CONFIG mode to load dependencies.
find_package(Eigen3 REQUIRED CONFIG)
find_package(nlohmann_json REQUIRED CONFIG)
find_package(spectra REQUIRED CONFIG)
find_package(ut REQUIRED CONFIG)

add_compile_options(
    $<$<COMPILE_LANG_AND_ID:CXX,Clang,AppleClang>:-Wno-missing-template-arg-list-after-template-kw>
)

if(ENABLE_ANALYSIS)
    include(CheckCXXCompilerFlag)
    set(CMAKE_REQUIRED_LINK_OPTIONS "-fsanitize=address")
    check_cxx_compiler_flag(-fsanitize=address HAS_ASAN)
    unset(CMAKE_REQUIRED_LINK_OPTIONS)
    check_cxx_compiler_flag(-fsanitize=undefined HAS_UBSAN)

    if(NOT HAS_ASAN OR NOT HAS_UBSAN)
        message(WARNING "Compiler does not support requested sanitizers")
        return()
    endif()

    add_compile_options(
        -fsanitize=address,undefined
        -fno-omit-frame-pointer
        -g
        -fno-sanitize-recover=all
    )
    add_link_options(
        -fsanitize=address,undefined
    )

    message(STATUS "Sanitizers enabled: ASan + UBSan")
endif()

# libraries
add_subdirectory(libraries)

# NonLocFEM
add_subdirectory(NonLocFEM)

# unit tests
include(CMakeTools/EmbeddedFiles.cmake)
add_subdirectory(tests)

if(ENABLE_ANALYSIS)
    include(${CMAKE_SOURCE_DIR}/CMakeTools/PVS-Studio.cmake)
    pvs_studio_add_target(TARGET code.analyze ALL
                        OUTPUT FORMAT json
                        ANALYZE NonLocFEM unit_tests
                        MODE GA OP 64 MISRA
                        LOG ${CMAKE_SOURCE_DIR}/NonLocFEM.json
                        CONFIG ${CMAKE_SOURCE_DIR}/nonlocal.pvsconfig
                        ARGS -e *.conan*)
endif()