cmake_minimum_required(VERSION 3.16)

project(parallel_utils)

add_library(parallel_utils_lib STATIC 
    uniform_ranges.cpp
    MPI_utils.cpp 
    OMP_utils.cpp
)
target_include_directories(parallel_utils_lib PUBLIC "." 
    ${MPI_C_INCLUDE_DIRS} 
    ${MPI_CXX_INCLUDE_DIRS}
    ${CONAN_INCLUDE_DIRS_EIGEN}
)
target_link_libraries(parallel_utils_lib 
    ${MPI_C_LIBRARIES} 
    ${MPI_CXX_LIBRARIES}
)