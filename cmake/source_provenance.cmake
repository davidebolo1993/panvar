# M15-K: the source revision actually built, captured at BUILD time (run on every build by the panvar_source_provenance
# target), so a commit after the last CMake configure is still recorded correctly. Usage:
#   cmake -DSOURCE_DIR=<checkout> -DOUTPUT=<header> -P source_provenance.cmake
# Writes PANVAR_SOURCE_COMMIT (`git rev-parse HEAD`, or "unavailable" outside a git checkout or without git) and
# PANVAR_SOURCE_TREE ("clean", "dirty" when tracked files differ from HEAD, "unknown" when that cannot be determined).
# Untracked files are not considered, as before. The header is replaced only when its content changes, so an unchanged
# revision recompiles nothing.
if(NOT SOURCE_DIR OR NOT OUTPUT)
    message(FATAL_ERROR "source_provenance.cmake needs -DSOURCE_DIR and -DOUTPUT")
endif()
set(commit "unavailable")
set(tree "unknown")
find_program(GIT_EXECUTABLE git)
if(GIT_EXECUTABLE)
    execute_process(COMMAND "${GIT_EXECUTABLE}" rev-parse --verify HEAD WORKING_DIRECTORY "${SOURCE_DIR}"
                    RESULT_VARIABLE rc OUTPUT_VARIABLE out OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    if(rc EQUAL 0 AND out MATCHES "^[0-9a-f]+$")
        set(commit "${out}")
        execute_process(COMMAND "${GIT_EXECUTABLE}" status --porcelain --untracked-files=no WORKING_DIRECTORY "${SOURCE_DIR}"
                        RESULT_VARIABLE src OUTPUT_VARIABLE status OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
        if(src EQUAL 0)
            if(status)
                set(tree "dirty")
            else()
                set(tree "clean")
            endif()
        endif()
    endif()
endif()
set(text "// Generated at build time by cmake/source_provenance.cmake; do not edit.\n#pragma once\n#define PANVAR_SOURCE_COMMIT \"${commit}\"\n#define PANVAR_SOURCE_TREE \"${tree}\"\n")
set(current "")
if(EXISTS "${OUTPUT}")
    file(READ "${OUTPUT}" current)
endif()
if(NOT current STREQUAL text)
    get_filename_component(dir "${OUTPUT}" DIRECTORY)
    file(MAKE_DIRECTORY "${dir}")
    file(WRITE "${OUTPUT}.tmp" "${text}")
    file(RENAME "${OUTPUT}.tmp" "${OUTPUT}")
endif()
