#cmake_minimum_required(VERSION 3.17)
#project(solver)
#
#add_subdirectory(solver_1d)
#add_subdirectory(solver_2d)
#
#set(SOLVER_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR})
#
#add_library(solver_lib INTERFACE)
#target_sources(solver_lib INTERFACE)
#target_include_directories(solver_lib INTERFACE ${SOLVER_LIB_DIR} metamath_lib mesh_lib)
#target_link_libraries(solver_lib INTERFACE finite_element_solver_1d_lib finite_element_solver_2d_lib)



cmake_minimum_required(VERSION 3.17)
project(finite_element_solvers)

add_subdirectory(solver_1d)
add_subdirectory(solver_2d)

set(FINITE_ELEMENT_SOLVERS_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_library(finite_element_solvers_lib INTERFACE)
target_sources(finite_element_solvers_lib INTERFACE)
#               solver_2d/finite_element_solver_base.hpp
#               solver_2d/thermal/heat_equation_solver_2d.hpp
#               solver_2d/structural/structural_solver.hpp)
target_include_directories(finite_element_solvers_lib INTERFACE ${FINITE_ELEMENT_SOLVERS_LIB_DIR} mesh_lib metamath_lib)
target_link_libraries(finite_element_solvers_lib INTERFACE
                      finite_element_solver_1d_lib
                      finite_element_solver_2d_lib
                     )




#find_package(MPI)
#if(MPI_FOUND)
#    add_definitions(-DMPI_USE=1)
#    set(INCLUDES "${MPI_C_INCLUDE_DIRS};${INCLUDES}")
##    set(LINKS "${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES};${LINKS}")
#
#    find_package(PkgConfig)
#    if(PkgConfig_FOUND)
#        pkg_check_modules(PETSc IMPORTED_TARGET petsc)
#        if(PETSc_FOUND)
#            set(INCLUDES "${PETSc_INCLUDE_DIRS};${INCLUDES}")
##            set(LINKS "${PETSc_LIBRARIES};${LINKS}")
#        endif()
#
#        pkg_check_modules(MKL IMPORTED_TARGET mkl)
#        if(MKL_FOUND)
#            set(INCLUDES "${MKL_INCLUDE_DIRS};${INCLUDES}")
##            set(LINKS "${MKL_LINK_LIBRARIES} ${MKL_LDFLAGS};${LINKS}")
#        endif()
#    endif()
#endif()
#
#if(MPI_FOUND)
#    set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
#    set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
#endif()



#add_executable(nonlocal src/1d/thermal/nonstationary.cpp)
#target_include_directories(nonlocal PUBLIC ${INCLUDES})
#target_link_libraries(
#    nonlocal
#    ${MKL_LINK_LIBRARIES} ${MKL_LDFLAGS}
#    ${PETSc_LIBRARIES} ${PETSc_LDFLAGS}
#    ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES}
#    metamath_lib mesh_lib finite_element_solvers_lib)