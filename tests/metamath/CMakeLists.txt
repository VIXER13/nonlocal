cmake_minimum_required(VERSION 3.16)

project(metamath_tests)


add_library(metamath_test_lib OBJECT power.cpp)
target_include_directories(
    metamath_test_lib
    PUBLIC #
        ${LIBRARIES_DIR}
        
)
target_link_libraries(
    metamath_test_lib
    PRIVATE #
        metamath_lib
        Boost::ut
)
