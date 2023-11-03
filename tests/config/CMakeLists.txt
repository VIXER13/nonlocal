cmake_minimum_required(VERSION 3.16)

project(config_tests)

FileEmbedAdd(required_fields.json)
FileEmbedAdd(double_conversion.json)

add_library(config_test_lib OBJECT 
    required_fields.cpp
    double_conversion.cpp
)
target_include_directories(config_test_lib PRIVATE
    "."
    ${CONAN_INCLUDE_DIRS_BOOST-EXT-UT}
    embedded_files
)
target_link_libraries(config_test_lib 
    PRIVATE
    embedded_files
    config_lib
)