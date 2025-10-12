cmake_minimum_required(VERSION 3.16)

project(config_tests)

set(INCLUDES $ENV{CPATH})

FileEmbedAdd(test_exceptions.json)
FileEmbedAdd(test_utils.json)

add_library(config_test_lib OBJECT test_exceptions.cpp test_utils.cpp)

target_include_directories(
    config_test_lib
    PRIVATE #
        ${LIBRARIES_DIR}
        ${INCLUDES}
        ${PROJECT_ROOT}/build/embedded_files
)
target_link_libraries(
    config_test_lib
    PRIVATE #
        config_lib
        Boost::ut
)
