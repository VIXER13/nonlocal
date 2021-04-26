cmake_minimum_required(VERSION 3.10)
project(nonlocalMPI)

set(INCLUDES)
set(LINKS)

find_package(MPI)
if(MPI_FOUND)
    add_definitions(-DMPI_USE=1)
    set(INCLUDES "${MPI_C_INCLUDE_DIRS};${INCLUDES}")
#    set(LINKS "${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES};${LINKS}")

    find_package(PkgConfig)
    if(PkgConfig_FOUND)
        pkg_check_modules(PETSc IMPORTED_TARGET petsc)
        if(PETSc_FOUND)
            set(INCLUDES "${PETSc_INCLUDE_DIRS};${INCLUDES}")
#            set(LINKS "${PETSc_LIBRARIES};${LINKS}")
        endif()

        pkg_check_modules(MKL IMPORTED_TARGET mkl)
        if(MKL_FOUND)
            set(INCLUDES "${MKL_INCLUDE_DIRS};${INCLUDES}")
#            set(LINKS "${MKL_LINK_LIBRARIES} ${MKL_LDFLAGS};${LINKS}")
        endif()
    endif()
endif()

add_subdirectory(include/metamath)
add_subdirectory(include/mesh)
add_subdirectory(include/solvers)

if(MPI_FOUND)
    set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
    set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
endif()

set(CMAKE_CXX_STANDARD 17)

find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_compile_options(-O2)

add_executable(nonlocal src/thermal_examples/thermal_nonstationary_orthotropic_quad.cpp)
target_include_directories(nonlocal PUBLIC ${INCLUDES})
target_link_libraries(
    nonlocal
    ${MKL_LINK_LIBRARIES} ${MKL_LDFLAGS}
    ${PETSc_LIBRARIES}
    ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES}
    metamath_lib mesh_lib finite_element_solvers_lib)
