cmake_minimum_required(VERSION 3.10)

project(MPI_utils)

message(STATUS "CMAKE_C_COMPILER = ${CMAKE_C_COMPILER}")
message(STATUS "CMAKE_CXX_COMPILER = ${CMAKE_CXX_COMPILER}")
message(STATUS "MPI_C_INCLUDE_DIRS = ${MPI_C_INCLUDE_DIRS}")
message(STATUS "MPI_CXX_INCLUDE_DIRS = ${MPI_CXX_INCLUDE_DIRS}")
message(STATUS "MPI_C_LIBRARIES = ${MPI_C_LIBRARIES}")
message(STATUS "MPI_CXX_LIBRARIES = ${MPI_CXX_LIBRARIES}")

add_library(MPI_utils_lib STATIC MPI_utils.cpp)
target_include_directories(MPI_utils_lib PUBLIC "." ${MPI_C_INCLUDE_DIRS} ${MPI_CXX_INCLUDE_DIRS})
target_link_libraries(MPI_utils_lib ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})