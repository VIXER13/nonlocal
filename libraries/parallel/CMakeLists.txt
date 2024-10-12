cmake_minimum_required(VERSION 3.16)

project(parallel)

set(INCLUDES $ENV{CPATH})
if (NOT DEFINED NO_CONAN)
    set(INCLUDES ${CONAN_INCLUDE_DIRS_EIGEN})
endif()

add_library(parallel_lib STATIC 
    uniform_ranges.cpp
    MPI_utils.cpp 
    OMP_utils.cpp
)
target_include_directories(parallel_lib PUBLIC "." 
    ${MPI_C_INCLUDE_DIRS} 
    ${MPI_CXX_INCLUDE_DIRS}
    ${INCLUDES}
)
target_link_libraries(parallel_lib 
    ${MPI_C_LIBRARIES} 
    ${MPI_CXX_LIBRARIES}
)