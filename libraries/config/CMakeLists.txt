cmake_minimum_required(VERSION 3.16)

project(config)

link_libraries(${JSONCPP_LIBRARIES})
add_library(config_lib STATIC 
    config_utils.cpp
    general_config_data.cpp
)
target_include_directories(config_lib PUBLIC "." ${JSONCPP_INCLUDE_DIRS})
target_link_libraries(config_lib 
    ${JSONCPP_LIBRARIES}
    nonlocal_constants_lib
)