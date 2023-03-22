cmake_minimum_required(VERSION 3.16)

project(json_parser)

link_libraries(${JSONCPP_LIBRARIES})
add_library(json_parser_lib STATIC json_parser.cpp)
target_include_directories(json_parser_lib PUBLIC "." "/usr/include/json")
target_link_libraries(json_parser_lib ${JSONCPP_LIBRARIES})