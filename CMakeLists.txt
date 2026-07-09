cmake_minimum_required(VERSION 3.16)

project(panvar VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(PANVAR_BUILD_TESTING "Build tests" ON)

find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)

set(MINIMAP2_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/minimap2)
set(MINIMAP2_LIB ${MINIMAP2_DIR}/libminimap2.a)
if(NOT EXISTS ${MINIMAP2_DIR}/minimap.h OR NOT EXISTS ${MINIMAP2_DIR}/Makefile)
    message(FATAL_ERROR
        "Missing minimap2 sources in external/minimap2. "
        "Run: git submodule update --init --recursive")
endif()
find_program(MINIMAP2_MAKE NAMES make gmake REQUIRED)
set(MINIMAP2_MAKE_ARGS)
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(arm64|aarch64|ARM64)$")
    list(APPEND MINIMAP2_MAKE_ARGS arm_neon=1 aarch64=1)
endif()
add_custom_target(minimap2_build
    COMMAND ${MINIMAP2_MAKE} -C ${MINIMAP2_DIR} libminimap2.a ${MINIMAP2_MAKE_ARGS}
    BYPRODUCTS ${MINIMAP2_LIB}
    COMMENT "Building minimap2 static library"
)
add_library(minimap2 STATIC IMPORTED GLOBAL)
set_target_properties(minimap2 PROPERTIES
    IMPORTED_LOCATION ${MINIMAP2_LIB}
)
add_dependencies(minimap2 minimap2_build)

# edlib (submodule, MIT) -- bit-parallel (Myers) banded edit distance behind fit_align.
# Built as its own target so the vendored source is exempt from our strict warning flags.
set(EDLIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/edlib/edlib)
if(NOT EXISTS ${EDLIB_DIR}/src/edlib.cpp OR NOT EXISTS ${EDLIB_DIR}/include/edlib.h)
    message(FATAL_ERROR
        "Missing edlib sources in external/edlib. "
        "Run: git submodule update --init --recursive")
endif()
add_library(edlib STATIC ${EDLIB_DIR}/src/edlib.cpp)
target_include_directories(edlib PUBLIC $<BUILD_INTERFACE:${EDLIB_DIR}/include>)
set_target_properties(edlib PROPERTIES POSITION_INDEPENDENT_CODE ON)

# Eigen (header-only) for the `associate` LMM eigendecomposition / linear algebra.
set(EIGEN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/eigen)
if(NOT EXISTS ${EIGEN_DIR}/Eigen/Dense)
    message(FATAL_ERROR
        "Missing Eigen headers in external/eigen. "
        "Run: git submodule update --init --recursive")
endif()

add_library(panvarlib
    src/gfa.cpp
    src/gfa_io.cpp
    src/align.cpp
    src/graph_utils.cpp
    src/graph_sort.cpp
    src/three_edge_cc.cpp
    src/union_find.cpp
    src/integrated_snarls.cpp
    src/bubble_path.cpp
    src/bubbles.cpp
    src/ref_path.cpp
    src/gtf.cpp
    src/gene_cn_kmer.cpp
    src/minimap2_align.cpp
    src/variant_call.cpp
    src/describe.cpp
    src/panphorte.cpp
    src/output.cpp
)

target_include_directories(panvarlib
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        ${MINIMAP2_DIR}
)

target_compile_options(panvarlib
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
)

# PRIVATE: edlib is an implementation detail of align.cpp, not part of panvarlib's
# public interface. For a static lib CMake still propagates the symbols to the final
# executable link, while keeping edlib out of the installed export set.
target_link_libraries(panvarlib PRIVATE edlib)

add_executable(panvar
    src/main.cpp
    src/cli_utils.cpp
    src/bubble_command.cpp
    src/inspect.cpp
    src/call_command.cpp
    src/benchmark_command.cpp
    src/describe_command.cpp
    src/panphorte_command.cpp
    src/associate_command.cpp
)
target_include_directories(panvar PRIVATE ${EIGEN_DIR})
target_link_libraries(panvar PRIVATE panvarlib minimap2 Threads::Threads ZLIB::ZLIB)
if(NOT WIN32)
    target_link_libraries(panvar PRIVATE m)
endif()
# GCC 8 and older keep std::filesystem symbols in a separate library.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9)
        target_link_libraries(panvar PRIVATE stdc++fs)
    endif()
endif()

target_compile_options(panvar
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
)

include(GNUInstallDirs)

install(TARGETS panvar panvarlib edlib
    EXPORT panvarTargets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(EXPORT panvarTargets
    FILE panvarTargets.cmake
    NAMESPACE panvar::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/panvar
)

if(PANVAR_BUILD_TESTING)
    include(CTest)
    if(BUILD_TESTING)
        # Fast, dependency-free smoke on a tiny hand-built graph that exercises every
        # call event type + confounders with exact-record assertions (no vg needed).
        set(PANVAR_SYN_SMOKE ${CMAKE_CURRENT_SOURCE_DIR}/tests/synthetic_smoke.sh)
        if(EXISTS ${PANVAR_SYN_SMOKE})
            add_test(
                NAME panvar_smoke_synthetic
                COMMAND /usr/bin/env bash ${PANVAR_SYN_SMOKE} $<TARGET_FILE:panvar> ${CMAKE_CURRENT_SOURCE_DIR}/tests/synthetic_data ${CMAKE_CURRENT_BINARY_DIR}/smoke_synthetic
            )
        endif()

        # Real-data integration smoke: the whole bubble -> inspect -> describe -> panphorte
        # -> call pipeline on a committed locus, asserting each stage's outputs. Registered on
        # C4 (the smallest graph, ~3s, no external tools) so the pipeline is exercised end-to-end
        # by ctest. Labelled "real" so it can be excluded with `ctest -LE real` if desired; the
        # larger loci (gstm1/lpa) stay manual: bash tests/real_smoke.sh ./build/panvar <gfa> <out>.
        set(PANVAR_REAL_SMOKE ${CMAKE_CURRENT_SOURCE_DIR}/tests/real_smoke.sh)
        set(PANVAR_REAL_SMOKE_GFA ${CMAKE_CURRENT_SOURCE_DIR}/tests/real_data/c4.gfa.gz)
        if(EXISTS ${PANVAR_REAL_SMOKE} AND EXISTS ${PANVAR_REAL_SMOKE_GFA})
            add_test(
                NAME panvar_smoke_real_c4
                COMMAND /usr/bin/env bash ${PANVAR_REAL_SMOKE} $<TARGET_FILE:panvar> ${PANVAR_REAL_SMOKE_GFA} ${CMAKE_CURRENT_BINARY_DIR}/smoke_real_c4
            )
            set_tests_properties(panvar_smoke_real_c4 PROPERTIES LABELS "real" TIMEOUT 300)
        endif()
    endif()
endif()
