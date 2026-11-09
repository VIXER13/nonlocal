cmake_minimum_required(VERSION 3.16)

project(config_tests)


FileEmbedAdd(test_exceptions.json)
FileEmbedAdd(test_utils.json)

add_library(config_test_lib
    OBJECT 
        test_exceptions.cpp 
        test_utils.cpp
        test_boundary_conditions_1d.cpp
        test_materials_1d.cpp
)

target_include_directories(
    config_test_lib
    PRIVATE #
        ${LIBRARIES_DIR}
        
        ${PROJECT_ROOT}/build/embedded_files
)
target_link_libraries(
    config_test_lib
    PRIVATE #
        config_lib
        Boost::ut
)
