cmake_minimum_required(VERSION 3.16)

project(math_expression_tests)

set(INCLUDES $ENV{CPATH})
if (NOT DEFINED NO_CONAN)
    set(INCLUDES ${CONAN_INCLUDE_DIRS_BOOST-EXT-UT})
endif()

add_library(math_expression_test_lib OBJECT 
    math_expression_t.cpp
)
target_include_directories(math_expression_test_lib PUBLIC
    ${LIBRARIES_DIR}
    ${INCLUDES}
)
target_link_libraries(math_expression_test_lib
    PRIVATE
    math_expression_lib
)


