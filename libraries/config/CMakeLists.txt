cmake_minimum_required(VERSION 3.16)

project(config)

link_libraries(${JSONCPP_LIBRARIES})
add_library(config_lib STATIC 
    config_utils.cpp
    save_data.cpp
)
target_include_directories(config_lib PUBLIC 
    "." 
    ${CONAN_INCLUDE_DIRS_NLOHMANN_JSON}
)
target_link_libraries(config_lib 
    nonlocal_constants_lib
)