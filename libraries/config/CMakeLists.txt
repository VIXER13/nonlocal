cmake_minimum_required(VERSION 3.16)

project(config)

set(INCLUDES $ENV{CPATH})
add_library(config_lib STATIC 
    config_utils.cpp
    read_mechanical_parameters.cpp
    read_mesh.cpp
    read_model.cpp
    save_data.cpp
    task_data.cpp
)
target_include_directories(config_lib PUBLIC
    ${LIBRARIES_DIR}
    ${INCLUDES}
)
target_link_libraries(config_lib
    PUBLIC
    nlohmann_json::nlohmann_json
    math_expression_lib
    logger_lib
    metamath_lib
    mesh_1d_lib
    mesh_2d_lib
)