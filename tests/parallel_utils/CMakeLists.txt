cmake_minimum_required(VERSION 3.16)

project(parallel_utils_tests)

add_library(parallel_utils_test_lib OBJECT 
    init_uniform_ranges_test.cpp
)
target_include_directories(parallel_utils_test_lib PUBLIC ".")
target_link_libraries(parallel_utils_test_lib
    Boost::ut
    parallel_utils_lib
)