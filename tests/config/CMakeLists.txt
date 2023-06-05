cmake_minimum_required(VERSION 3.16)

project(config_tests)

FileEmbedSetup()
FileEmbedAdd(required_fields.json)
FileEmbedAdd(reverse_conversion.json)

add_library(config_test_lib OBJECT 
    required_fields.cpp
    reverse_conversion.cpp
)
target_include_directories(config_test_lib PRIVATE ".")
target_link_libraries(config_test_lib 
    PRIVATE
    embedded_files

    PUBLIC
    Boost::ut
    config_lib
)