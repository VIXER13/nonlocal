cmake_minimum_required(VERSION 3.16)

project(test_template)


set(INCLUDES $ENV{CPATH})
if (NOT DEFINED NO_CONAN)
    set(INCLUDES ${CONAN_INCLUDE_DIRS_BOOST-EXT-UT})
endif()

add_library(test_template_lib OBJECT 
    test_template.cpp
)
target_include_directories(test_template_lib PUBLIC
    "."
    ${INCLUDES}
)
