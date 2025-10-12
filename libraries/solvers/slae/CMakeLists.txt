cmake_minimum_required(VERSION 3.16)
project(slae_solver)

add_library(slae_solver_lib INTERFACE)
target_sources(slae_solver_lib INTERFACE conjugate_gradient.hpp)
target_include_directories(slae_solver_lib INTERFACE ${LIBRARIES_DIR})
target_link_libraries(slae_solver_lib
	INTERFACE
	parallel_lib
	Eigen3::Eigen
	Spectra::Spectra
)
