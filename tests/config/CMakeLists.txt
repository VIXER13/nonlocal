cmake_minimum_required(VERSION 3.16)

project(config_tests)

FileEmbedAdd(test_exceptions.json)
FileEmbedAdd(test_utils.json)

set(INCLUDES $ENV{CPATH})
if (NOT DEFINED NO_CONAN)
    set(INCLUDES ${CONAN_INCLUDE_DIRS_BOOST-EXT-UT})
endif()

add_library(config_test_lib OBJECT
    test_exceptions.cpp
    test_utils.cpp
)
target_include_directories(config_test_lib PRIVATE
    ${LIBRARIES_DIR}
    ${INCLUDES}
    embedded_files
)
target_link_libraries(config_test_lib 
    PRIVATE
    embedded_files
    config_lib
)
