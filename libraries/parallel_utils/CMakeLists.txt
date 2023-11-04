cmake_minimum_required(VERSION 3.16)

project(parallel_utils)

if (MPI_BUILD)
    find_package(MPI)
    if(MPI_FOUND)
        add_definitions(-DMPI_BUILD=1)
        set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
        set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
    endif()
endif()

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