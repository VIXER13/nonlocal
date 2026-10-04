cmake_minimum_required(VERSION 3.16)

project(argparse)

add_library(
    argparse_lib
    STATIC
    cli.cpp
)
target_include_directories(
    argparse_lib
    PUBLIC #
        ${LIBRARIES_DIR}
)