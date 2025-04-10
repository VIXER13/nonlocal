cmake_minimum_required(VERSION 3.16)
project(finite_elements_base)

add_library(finite_elements_base_lib INTERFACE)
target_sources(finite_elements_base_lib INTERFACE quadrature_base.hpp
                                                  finite_element_base.hpp
                                                  finite_element_integrate_base.hpp)
target_include_directories(finite_elements_base_lib INTERFACE ${LIBRARIES_DIR})