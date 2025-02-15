cmake_minimum_required(VERSION 3.16)

project(parser_tests)

set(INCLUDES $ENV{CPATH})
if (NOT DEFINED NO_CONAN)
    set(INCLUDES ${CONAN_INCLUDE_DIRS_BOOST-EXT-UT})
endif()

add_library(parser_test_lib OBJECT 
    parser_t.cpp
)
target_include_directories(parser_test_lib PUBLIC
    ${LIBRARIES_DIR}
    ${INCLUDES}
)
target_link_libraries(parser_test_lib
    PRIVATE
    parser_lib
)


