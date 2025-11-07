cmake_minimum_required(VERSION 3.16)

project(parallel_tests)


add_library(parallel_test_lib OBJECT uniform_ranges_test.cpp)
target_include_directories(
    parallel_test_lib
    PUBLIC #
        ${LIBRARIES_DIR}
        
)
target_link_libraries(
    parallel_test_lib
    PRIVATE #
        parallel_lib
        Boost::ut
)
