cmake_minimum_required(VERSION 3.16)

project(config)

find_package(nlohmann_json REQUIRED)

link_libraries(${JSONCPP_LIBRARIES})
add_library(config_lib STATIC 
    config_utils.cpp
    # save_data.cpp
)
target_include_directories(config_lib PUBLIC ".")
target_link_libraries(config_lib 
    nlohmann_json::nlohmann_json
    nonlocal_constants_lib
)