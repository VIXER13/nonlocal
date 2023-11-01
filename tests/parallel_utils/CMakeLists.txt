cmake_minimum_required(VERSION 3.16)

project(parallel_utils_tests)

add_library(parallel_utils_test_lib OBJECT 
    uniform_ranges_test.cpp
)
target_include_directories(parallel_utils_test_lib PUBLIC
    "."
    ${CONAN_INCLUDE_DIRS_BOOST-EXT-UT}
)
target_link_libraries(parallel_utils_test_lib
    parallel_utils_lib
)