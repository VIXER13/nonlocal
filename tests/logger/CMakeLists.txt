cmake_minimum_required(VERSION 3.16)

project(logger_tests)


add_library(logger_test_lib OBJECT logger_test.cpp)
target_include_directories(
    logger_test_lib
    PUBLIC #
        ${LIBRARIES_DIR}
        
)
target_link_libraries(
    logger_test_lib
    PRIVATE #
        logger_lib
        Boost::ut
)
