cmake_minimum_required(VERSION 3.16)

project(json_parser)

link_libraries(${JSONCPP_LIBRARIES})
add_library(json_parser_lib STATIC 
    config_utils.cpp
    general_config_data.cpp
    thermal_config_data.cpp
)
target_include_directories(json_parser_lib PUBLIC "." ${JSONCPP_INCLUDE_DIRS})
target_link_libraries(json_parser_lib 
    ${JSONCPP_LIBRARIES}
    nonlocal_constants_lib
)