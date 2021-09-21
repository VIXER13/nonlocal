cmake_minimum_required(VERSION 3.10)

project(MPI_utils)

add_library(MPI_utils_lib STATIC MPI_utils.cpp)
target_include_directories(MPI_utils_lib PUBLIC "." ${MPI_C_INCLUDE_DIRS} ${MPI_CXX_INCLUDE_DIRS})
target_link_libraries(MPI_utils_lib ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})