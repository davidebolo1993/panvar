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

# abPOA (submodule + nested SIMDe, MIT) -- partial-order alignment behind `refine`'s residual MSA.
# Built as a static library via its own Makefile (like minimap2); its Makefile auto-detects arm64
# (uname -m) and selects the apple-m1/SIMDe SIMD path.
set(ABPOA_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/abpoa)
set(ABPOA_LIB ${ABPOA_DIR}/lib/libabpoa.a)
if(NOT EXISTS ${ABPOA_DIR}/include/abpoa.h OR NOT EXISTS ${ABPOA_DIR}/include/simde/simde/x86/sse2.h)
    message(FATAL_ERROR
        "Missing abPOA sources in external/abpoa (or its nested SIMDe). "
        "Run: git submodule update --init --recursive")
endif()
find_program(ABPOA_MAKE NAMES make gmake REQUIRED)
add_custom_target(abpoa_build
    COMMAND ${ABPOA_MAKE} -C ${ABPOA_DIR} lib/libabpoa.a
    BYPRODUCTS ${ABPOA_LIB}
    COMMENT "Building abPOA static library"
)
add_library(abpoa STATIC IMPORTED GLOBAL)
set_target_properties(abpoa PROPERTIES IMPORTED_LOCATION ${ABPOA_LIB})
add_dependencies(abpoa abpoa_build)

# minigraph (submodule, MIT) -- progressive graph generation (mg_ggen) + graph mapping (mg_map) behind
# `rebuild`. Its Makefile builds only the executable, so we build the objects and archive them ourselves.
#
# minimap2 vendored the same lh3 chaining/util code, so the two libraries export ~24 identical symbols
# (mg_lchain_dp, gfa_aux_*, seq_nt4_table, ksort/krmq instantiations) and cannot simply be linked side by
# side. Fix: partial-link (ld -r) all of minigraph's objects EXCEPT main.o into a single relocatable
# object, keeping only the handful of entry points `rebuild` calls global and localizing everything else.
# minigraph's internal calls then bind to its own copies, and nothing collides with minimap2.
set(MINIGRAPH_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/minigraph)
set(MINIGRAPH_LIB ${MINIGRAPH_DIR}/libminigraph.a)
set(MINIGRAPH_MERGED ${MINIGRAPH_DIR}/minigraph_merged.o)
if(NOT EXISTS ${MINIGRAPH_DIR}/minigraph.h OR NOT EXISTS ${MINIGRAPH_DIR}/Makefile)
    message(FATAL_ERROR
        "Missing minigraph sources in external/minigraph. "
        "Run: git submodule update --init --recursive")
endif()
find_program(MINIGRAPH_MAKE NAMES make gmake REQUIRED)
if(APPLE)
    set(MINIGRAPH_EXPORTS ${CMAKE_CURRENT_SOURCE_DIR}/cmake/minigraph_exports_macho.txt)
else()
    set(MINIGRAPH_EXPORTS ${CMAKE_CURRENT_SOURCE_DIR}/cmake/minigraph_exports_elf.txt)
endif()
add_custom_target(minigraph_build
    COMMAND ${MINIGRAPH_MAKE} -C ${MINIGRAPH_DIR} minigraph
    COMMAND ${CMAKE_COMMAND} -E rm -f ${MINIGRAPH_LIB} ${MINIGRAPH_MERGED}
    COMMAND /usr/bin/env bash ${CMAKE_CURRENT_SOURCE_DIR}/cmake/build_minigraph_lib.sh
            ${MINIGRAPH_DIR} ${MINIGRAPH_EXPORTS}
    BYPRODUCTS ${MINIGRAPH_LIB} ${MINIGRAPH_MERGED}
    COMMENT "Building minigraph static library (symbols localized against minimap2)"
)
add_library(minigraph STATIC IMPORTED GLOBAL)
set_target_properties(minigraph PROPERTIES IMPORTED_LOCATION ${MINIGRAPH_LIB})
add_dependencies(minigraph minigraph_build)

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
    src/bubble_alleles.cpp
    src/syncmer.cpp
    src/genotype.cpp
    src/genotype_blocks.cpp
    src/genotype_index.cpp
    src/genotype_markers.cpp
    src/genotype_reads.cpp
    src/pangenie_model.cpp
    src/bubble_path.cpp
    src/bubbles.cpp
    src/ref_path.cpp
    src/gtf.cpp
    src/gene_cn_kmer.cpp
    src/minimap2_align.cpp
    src/node_coverage.cpp
    src/variant_call.cpp
    src/describe.cpp
    src/panphorte.cpp
    src/poa.cpp
    src/refine.cpp
    src/rebuild.cpp
    src/output.cpp
)

target_include_directories(panvarlib
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        ${ABPOA_DIR}/include
)

# SYSTEM: minimap2's and minigraph's headers use C99 flexible array members, which our -Wpedantic flags
# in C++. The vendored sources are not ours to fix, so exempt them as the other vendored targets are.
target_include_directories(panvarlib SYSTEM PRIVATE ${MINIMAP2_DIR} ${MINIGRAPH_DIR})

target_compile_options(panvarlib
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
)

# PRIVATE: edlib is an implementation detail of align.cpp, not part of panvarlib's
# public interface. For a static lib CMake still propagates the symbols to the final
# executable link, while keeping edlib out of the installed export set.
target_link_libraries(panvarlib PRIVATE edlib abpoa minigraph ZLIB::ZLIB)

add_executable(panvar
    src/main.cpp
    src/cli_utils.cpp
    src/bubble_command.cpp
    src/inspect.cpp
    src/call_command.cpp
    src/genotype_command.cpp
    src/benchmark_command.cpp
    src/describe_command.cpp
    src/panphorte_command.cpp
    src/refine_command.cpp
    src/rebuild_command.cpp
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

        # Deterministic statistical assertions for `associate`: linear p against R's lm() (the normal
        # tail would be 83 orders of magnitude out at 9 df), the logistic separation contract, case/
        # control carrier counts, and GRM rejection. No RNG, no R at run time, ~1s.
        set(PANVAR_ASSOC_STATS ${CMAKE_CURRENT_SOURCE_DIR}/tests/associate_stats.sh)
        if(EXISTS ${PANVAR_ASSOC_STATS})
            add_test(
                NAME panvar_associate_stats
                COMMAND /usr/bin/env bash ${PANVAR_ASSOC_STATS} $<TARGET_FILE:panvar> ${CMAKE_CURRENT_BINARY_DIR}/associate_stats
            )
            set_tests_properties(panvar_associate_stats PROPERTIES TIMEOUT 120)
        endif()

        # bubble contract on hand-checkable graphs. CTest exercised the default internal cactus finder
        # nowhere: synthetic_smoke.sh supplies checked-in external snarls, and real_smoke.sh makes its
        # bubble assertions conditional on finding a first bubble -- and a header-only CSV is non-empty,
        # so a total discovery regression would escape.
        set(PANVAR_BUBBLE_STATS ${CMAKE_CURRENT_SOURCE_DIR}/tests/bubble_stats.sh)
        if(EXISTS ${PANVAR_BUBBLE_STATS})
            add_test(
                NAME panvar_bubble_stats
                COMMAND /usr/bin/env bash ${PANVAR_BUBBLE_STATS} $<TARGET_FILE:panvar> ${CMAKE_CURRENT_BINARY_DIR}/bubble_stats
            )
            set_tests_properties(panvar_bubble_stats PROPERTIES TIMEOUT 180)
        endif()

        # inspect contract. real_smoke.sh only checks that inspect's files exist and that the gzip is
        # readable, which passes whether or not the contents are right -- it would not have caught a
        # whole allele missing from the FASTA, a duplicate path name, or a bubbles CSV from a different
        # graph running to completion.
        set(PANVAR_INSPECT_STATS ${CMAKE_CURRENT_SOURCE_DIR}/tests/inspect_stats.sh)
        if(EXISTS ${PANVAR_INSPECT_STATS})
            add_test(
                NAME panvar_inspect_stats
                COMMAND /usr/bin/env bash ${PANVAR_INSPECT_STATS} $<TARGET_FILE:panvar> ${CMAKE_CURRENT_BINARY_DIR}/inspect_stats
            )
            set_tests_properties(panvar_inspect_stats PROPERTIES TIMEOUT 180)
        endif()

        # panphorte contract. The registered C4 smoke test passes with five normalized=no rows, so it
        # never exercises folding at all.
        set(PANVAR_PANPHORTE_STATS ${CMAKE_CURRENT_SOURCE_DIR}/tests/panphorte_stats.sh)
        if(EXISTS ${PANVAR_PANPHORTE_STATS})
            add_test(
                NAME panvar_panphorte_stats
                COMMAND /usr/bin/env bash ${PANVAR_PANPHORTE_STATS} $<TARGET_FILE:panvar> ${CMAKE_CURRENT_BINARY_DIR}/panphorte_stats
            )
            set_tests_properties(panvar_panphorte_stats PROPERTIES TIMEOUT 300)
        endif()

        # refine contract. real_smoke.sh only checks that refine's files exist and that `call` runs
        # afterwards, which passes whether or not the rewrite was correct.
        set(PANVAR_REFINE_STATS ${CMAKE_CURRENT_SOURCE_DIR}/tests/refine_stats.sh)
        if(EXISTS ${PANVAR_REFINE_STATS})
            add_test(
                NAME panvar_refine_stats
                COMMAND /usr/bin/env bash ${PANVAR_REFINE_STATS} $<TARGET_FILE:panvar> ${CMAKE_CURRENT_BINARY_DIR}/refine_stats
            )
            set_tests_properties(panvar_refine_stats PROPERTIES TIMEOUT 300)
        endif()

        # rebuild contract: the input is never destroyed, a healthy graph survives, a malformed one is
        # refused rather than silently repaired, and nothing is left behind. CTest previously exercised
        # no rebuild at all.
        set(PANVAR_REBUILD_STATS ${CMAKE_CURRENT_SOURCE_DIR}/tests/rebuild_stats.sh)
        if(EXISTS ${PANVAR_REBUILD_STATS})
            add_test(
                NAME panvar_rebuild_stats
                COMMAND /usr/bin/env bash ${PANVAR_REBUILD_STATS} $<TARGET_FILE:panvar> ${CMAKE_CURRENT_BINARY_DIR}/rebuild_stats ${CMAKE_CURRENT_SOURCE_DIR}/tests/real_data/c4.gfa.gz
            )
            # Real-data integration: several forced rebuilds of C4 (reproducibility across threads,
            # k either side of 31, the reference-ambiguity and verdict cases). Minutes, not seconds.
            set_tests_properties(panvar_rebuild_stats PROPERTIES TIMEOUT 900 LABELS "real")
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
