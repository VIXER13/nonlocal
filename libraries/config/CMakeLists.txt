cmake_minimum_required(VERSION 3.16)

project(config)

set(INCLUDES $ENV{CPATH})
if (NOT DEFINED NO_CONAN)
    set(INCLUDES ${CONAN_INCLUDE_DIRS_NLOHMANN_JSON})
endif()

link_libraries(${JSONCPP_LIBRARIES})
add_library(config_lib STATIC 
    config_utils.cpp
    task_data.cpp
    save_data.cpp
    mesh_data.cpp
    materials_data.cpp
)
target_include_directories(config_lib PUBLIC 
    "." 
    ${INCLUDES}
)
target_link_libraries(config_lib 
    logger_lib
    nonlocal_constants_lib
)