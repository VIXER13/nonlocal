function(FileEmbedSetup)
    if (NOT EXISTS ${CMAKE_BINARY_DIR}/embedded_files)
        file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/embedded_files)
    endif ()

    if (NOT EXISTS ${CMAKE_BINARY_DIR}/embedded_files/embedded_files_empty.c)
        file(WRITE ${CMAKE_BINARY_DIR}/embedded_files/embedded_files_empty.c "")
    endif ()

    add_library(embedded_files ${CMAKE_BINARY_DIR}/embedded_files/embedded_files_empty.c)
    target_include_directories(embedded_files PUBLIC ${CMAKE_BINARY_DIR}/embedded_files)

endfunction()

function(FileEmbedAdd file)
    FileEmbedGenerate(${file} var)
    target_sources(embedded_files PUBLIC ${var})

    add_custom_command(
            OUTPUT ${var}
            COMMAND ${CMAKE_COMMAND}
            -DRUN_embedded_files_GENERATE=1
            -Dembedded_files_GENERATE_PATH=${file}
            -P ${CMAKE_SOURCE_DIR}/cmake/FileEmbed.cmake
            MAIN_DEPENDENCY ${file}
    )
endfunction()

function(FileEmbedGenerate file generated_c)

    get_filename_component(base_filename ${file} NAME)
    set(output_filename "${base_filename}.c")
    string(MAKE_C_IDENTIFIER ${base_filename} c_name)
    file(READ ${file} content HEX)

    # Separate into individual bytes.
    string(REGEX MATCHALL "([A-Fa-f0-9][A-Fa-f0-9])" SEPARATED_HEX ${content})

    set(output_c "")

    set(counter 0)
    foreach (hex IN LISTS SEPARATED_HEX)
        string(APPEND output_c "0x${hex},")
        MATH(EXPR counter "${counter}+1")
        if (counter GREATER 16)
            string(APPEND output_c "\n    ")
            set(counter 0)
        endif ()
    endforeach ()

    set(output_c 
"#include \"${c_name}.h\"
const char ${c_name}_data[] = {
    ${output_c}0
}\;
const size_t ${c_name}_size = sizeof(${c_name}_data)\;"
)

    set(output_h 
"#pragma once
#include <stddef.h>
extern const char ${c_name}_data[]\;
extern const size_t ${c_name}_size\;"
)

    if (NOT EXISTS ${CMAKE_BINARY_DIR}/embedded_files)
        file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}embedded_files)
    endif ()


    file(WRITE ${CMAKE_BINARY_DIR}/embedded_files/${c_name}.c
            ${output_c})

    file(WRITE ${CMAKE_BINARY_DIR}/embedded_files/${c_name}.h
            ${output_h})

    set(${generated_c} ${CMAKE_BINARY_DIR}/embedded_files/${c_name}.c PARENT_SCOPE)

endfunction()

if (RUN_embedded_files_GENERATE)
    FileEmbedGenerate(${embedded_files_GENERATE_PATH} var)
endif ()
