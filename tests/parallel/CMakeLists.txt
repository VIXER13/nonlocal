cmake_minimum_required(VERSION 3.16)

project(parallel_tests)

set(INCLUDES $ENV{CPATH})
if (NOT DEFINED NO_CONAN)
    set(INCLUDES ${CONAN_INCLUDE_DIRS_BOOST-EXT-UT})
endif()

add_library(parallel_test_lib OBJECT 
    uniform_ranges_test.cpp
)
target_include_directories(parallel_test_lib PUBLIC
    "."
    ${INCLUDES}
)
target_link_libraries(parallel_test_lib
    PRIVATE
    parallel_lib
)