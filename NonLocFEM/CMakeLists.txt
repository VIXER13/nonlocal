cmake_minimum_required(VERSION 3.16)

project(NonLocFEM)

add_executable(NonLocFEM main.cpp determine_problem.cpp)
target_include_directories(
    NonLocFEM
    PRIVATE #
        ${LIBRARIES_DIR}
)
target_link_libraries(
    NonLocFEM
    PRIVATE #
        metamath_instantiations_lib
        config_lib
        finite_element_solver_1d_lib
        finite_element_solver_2d_lib
)
target_precompile_headers(NonLocFEM PRIVATE
    "$<BUILD_INTERFACE:${METAMATH_PCH}>")
