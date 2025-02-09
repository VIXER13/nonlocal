cmake_minimum_required(VERSION 3.16)

project(finite_elements_tests)

set(INCLUDES $ENV{CPATH})
if (NOT DEFINED NO_CONAN)
    set(INCLUDES ${CONAN_INCLUDE_DIRS_BOOST-EXT-UT})
endif()

add_library(finite_elements_test_lib OBJECT
    init_lagrangian_elements_1d.cpp
    init_lagrangian_elements_2d.cpp
    init_triangle_elements_2d.cpp
    init_serendipity_elements_2d.cpp
    quadratures_test.cpp
    element_1d_test.cpp
    element_2d_test.cpp
)

target_include_directories(finite_elements_test_lib PUBLIC 
    ${LIBRARIES_DIR}
    ${INCLUDES}
)
target_link_libraries(finite_elements_test_lib
    PRIVATE
    metamath_lib
)