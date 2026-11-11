cmake_minimum_required(VERSION 3.16)

project(panvar VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(PANVAR_BUILD_TESTING "Build tests" ON)

# The slow real-data suites (rebuild drives minigraph and takes minutes) are registered but disabled
# unless this is on.
option(PANVAR_SLOW_TESTS "Enable the slow test label (the ~6 minute rebuild suite)" OFF)

find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)

set(MINIMAP2_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/minimap2)
set(MINIMAP2_LIB ${MINIMAP2_DIR}/libminimap2.a)
if(NOT EXISTS ${MINIMAP2_DIR}/minimap.h OR NOT EXISTS ${MINIMAP2_DIR}/Makefile)
    message(FATAL_ERROR
        "Missing minimap2 sources in external/minimap2. "
        "Run: git submodule update --init --recursive")
endif()
# find_program(... REQUIRED) needs CMake 3.18, so a missing make is reported by hand.
find_program(VENDORED_MAKE NAMES make gmake)
if(NOT VENDORED_MAKE)
    message(FATAL_ERROR "make (or gmake) is needed to build the vendored minimap2, abPOA and minigraph")
endif()
set(MINIMAP2_MAKE ${VENDORED_MAKE})
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

# SSHash (submodule + nested PTHash, MIT) -- the compressed k-mer dictionary `background --compact` writes.
# Compiled here as a plain static library from its sources, with panvar's flags rather than its own (which add
# -Werror, and -mavx2 -mbmi2 on x86-64): its x86 fast paths are guarded by __AVX2__, so the default build runs
# on any x86-64; -DPANVAR_SSHASH_AVX2=ON turns them on for machines that have AVX2.
set(SSHASH_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/sshash)
if(NOT EXISTS ${SSHASH_DIR}/include/dictionary.hpp OR NOT EXISTS ${SSHASH_DIR}/external/pthash/include/pthash.hpp)
    message(FATAL_ERROR
        "Missing SSHash sources in external/sshash (or its nested PTHash). "
        "Run: git submodule update --init --recursive")
endif()
option(PANVAR_SSHASH_AVX2 "Build SSHash with AVX2 and BMI2 (x86-64 machines that have them)" OFF)
add_library(sshash STATIC
    ${SSHASH_DIR}/external/gz/zip_stream.cpp
    ${SSHASH_DIR}/external/cityhash/cityhash.cpp
)   # its templated sources (src/*.cpp) are included by src/background_sshash.cpp, as SSHash's own tools do
target_include_directories(sshash SYSTEM PUBLIC
    ${SSHASH_DIR}/include
    ${SSHASH_DIR}
    ${SSHASH_DIR}/external/pthash/include
    ${SSHASH_DIR}/external/pthash/external/bits/include
    ${SSHASH_DIR}/external/pthash/external/fastmod
    ${SSHASH_DIR}/external/pthash/external/bits/external/essentials/include
    ${SSHASH_DIR}/external/pthash/external/xxHash
    ${SSHASH_DIR}/external/pthash/external/mm_file/include
)
target_compile_options(sshash PRIVATE -O3 -w)
if(PANVAR_SSHASH_AVX2 AND CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|amd64)$")
    target_compile_options(sshash PUBLIC -mavx2 -mbmi2)
endif()
target_link_libraries(sshash PUBLIC ZLIB::ZLIB Threads::Threads)

# htslib (submodule + nested htscodecs, MIT/BSD) -- SAM, BAM and CRAM input for `recruit`. Compiled here from
# its sources as a static library, with a generated config.h instead of its configure step: no network
# access (libcurl off), bzip2 and lzma CRAM codecs and libdeflate when the system has them. Its codec SIMD
# variants are compiled where the compiler enables them by default (NEON on 64-bit ARM); elsewhere CRAM uses
# the portable codecs.
set(HTSLIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/htslib)
if(NOT EXISTS ${HTSLIB_DIR}/hts.c OR NOT EXISTS ${HTSLIB_DIR}/htscodecs/htscodecs/htscodecs.c)
    message(FATAL_ERROR
        "Missing htslib sources in external/htslib (or its nested htscodecs). "
        "Run: git submodule update --init --recursive")
endif()
find_package(BZip2 QUIET)
find_package(LibLZMA QUIET)
find_library(LIBDEFLATE_LIBRARY NAMES deflate)
find_path(LIBDEFLATE_INCLUDE_DIR libdeflate.h)
set(HTSLIB_GEN ${CMAKE_CURRENT_BINARY_DIR}/htslib_generated)
set(HTSLIB_CONFIG "/* config.h for the vendored htslib, written by panvar's CMake */\n#ifndef _XOPEN_SOURCE\n#define _XOPEN_SOURCE 600\n#endif\n#define HAVE_DRAND48 1\n")
if(BZIP2_FOUND)
    string(APPEND HTSLIB_CONFIG "#define HAVE_LIBBZ2 1\n")
endif()
if(LIBLZMA_FOUND)
    string(APPEND HTSLIB_CONFIG "#define HAVE_LIBLZMA 1\n#define HAVE_LZMA_H 1\n")
endif()
if(LIBDEFLATE_LIBRARY AND LIBDEFLATE_INCLUDE_DIR)
    string(APPEND HTSLIB_CONFIG "#define HAVE_LIBDEFLATE 1\n")
endif()
file(WRITE ${HTSLIB_GEN}/config.h.new "${HTSLIB_CONFIG}")
file(WRITE ${HTSLIB_GEN}/version.h.new "#define HTS_VERSION_TEXT \"1.23\"\n#define HTSCODECS_VERSION_TEXT \"1.6.5\"\n")
file(WRITE ${HTSLIB_GEN}/config_vars.h.new
    "#define HTS_CC \"${CMAKE_C_COMPILER}\"\n#define HTS_CPPFLAGS \"\"\n#define HTS_CFLAGS \"\"\n#define HTS_LDFLAGS \"\"\n#define HTS_LIBS \"\"\n")
foreach(h config.h version.h config_vars.h)   # rewritten only when changed, so htslib is not rebuilt each time
    configure_file(${HTSLIB_GEN}/${h}.new ${HTSLIB_GEN}/${h} COPYONLY)
endforeach()
set(HTSLIB_SOURCES)
foreach(f kfunc kstring bcf_sr_sort bgzf errmod faidx header hfile hts hts_expr hts_os md5 multipart probaln realn
          regidx region sam sam_mods simd synced_bcf_reader vcf_sweep tbx textutils thread_pool vcf vcfutils
          cram/cram_codecs cram/cram_decode cram/cram_encode cram/cram_external cram/cram_index cram/cram_io
          cram/cram_stats cram/mFILE cram/open_trace_file cram/pooled_alloc cram/string_alloc
          htscodecs/htscodecs/arith_dynamic htscodecs/htscodecs/fqzcomp_qual htscodecs/htscodecs/htscodecs
          htscodecs/htscodecs/pack htscodecs/htscodecs/rANS_static4x16pr htscodecs/htscodecs/rANS_static32x16pr
          htscodecs/htscodecs/rANS_static32x16pr_neon htscodecs/htscodecs/rANS_static htscodecs/htscodecs/rle
          htscodecs/htscodecs/tokenise_name3 htscodecs/htscodecs/utils)
    list(APPEND HTSLIB_SOURCES ${HTSLIB_DIR}/${f}.c)
endforeach()
add_library(hts STATIC ${HTSLIB_SOURCES})
target_include_directories(hts PRIVATE ${HTSLIB_GEN} ${HTSLIB_DIR}/htscodecs)
target_include_directories(hts SYSTEM PUBLIC ${HTSLIB_DIR})
target_compile_options(hts PRIVATE -O2 -w)
set_target_properties(hts PROPERTIES POSITION_INDEPENDENT_CODE ON C_STANDARD 99)
target_link_libraries(hts PUBLIC ZLIB::ZLIB Threads::Threads)
if(NOT APPLE)
    target_link_libraries(hts PUBLIC m)
endif()
if(BZIP2_FOUND)
    target_link_libraries(hts PUBLIC BZip2::BZip2)
endif()
if(LIBLZMA_FOUND)
    target_link_libraries(hts PUBLIC LibLZMA::LibLZMA)
endif()
if(LIBDEFLATE_LIBRARY AND LIBDEFLATE_INCLUDE_DIR)
    target_include_directories(hts PRIVATE ${LIBDEFLATE_INCLUDE_DIR})
    target_link_libraries(hts PUBLIC ${LIBDEFLATE_LIBRARY})
endif()

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
set(ABPOA_MAKE ${VENDORED_MAKE})
set(ABPOA_MAKE_ARGS)
# On x86-64 abPOA compiles its aligner once per SIMD level and names each copy from the compiler's CPU
# macros. A compiler whose default target already has AVX2 names the SSE2 and SSE4.1 copies as AVX2, and
# the dispatcher's SSE symbols go missing at link time: build from the plain x86-64 baseline instead.
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|amd64)$")
    list(APPEND ABPOA_MAKE_ARGS
        "CFLAGS=-O3 -Wall -Wno-unused-function -Wno-misleading-indentation -DUSE_SIMDE -DSIMDE_ENABLE_NATIVE_ALIASES -march=x86-64")
endif()
add_custom_target(abpoa_build
    COMMAND ${ABPOA_MAKE} -C ${ABPOA_DIR} lib/libabpoa.a ${ABPOA_MAKE_ARGS}
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
set(MINIGRAPH_MAKE ${VENDORED_MAKE})
# minigraph's Makefile adds -msse4, which GCC rejects off x86; minigraph uses no SSE intrinsics, so
# elsewhere its flags are passed without it.
set(MINIGRAPH_MAKE_ARGS)
if(NOT CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|amd64|i[3-6]86|x86)$")
    list(APPEND MINIGRAPH_MAKE_ARGS "CFLAGS=-g -Wall -Wc++-compat -std=c99 -O3")
endif()
if(APPLE)
    set(MINIGRAPH_EXPORTS ${CMAKE_CURRENT_SOURCE_DIR}/cmake/minigraph_exports_macho.txt)
else()
    set(MINIGRAPH_EXPORTS ${CMAKE_CURRENT_SOURCE_DIR}/cmake/minigraph_exports_elf.txt)
endif()
add_custom_target(minigraph_build
    COMMAND ${MINIGRAPH_MAKE} -C ${MINIGRAPH_DIR} minigraph ${MINIGRAPH_MAKE_ARGS}
    COMMAND ${CMAKE_COMMAND} -E rm -f ${MINIGRAPH_LIB} ${MINIGRAPH_MERGED}
    COMMAND /usr/bin/env bash ${CMAKE_CURRENT_SOURCE_DIR}/cmake/build_minigraph_lib.sh
            ${MINIGRAPH_DIR} ${MINIGRAPH_EXPORTS}
    BYPRODUCTS ${MINIGRAPH_LIB} ${MINIGRAPH_MERGED}
    COMMENT "Building minigraph static library (symbols localized against minimap2)"
    VERBATIM
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
    src/variant_nodes.cpp
    src/graph_sort.cpp
    src/three_edge_cc.cpp
    src/union_find.cpp
    src/integrated_snarls.cpp
    src/bubble_alleles.cpp
    src/syncmer.cpp
    src/bubble_path.cpp
    src/bubbles.cpp
    src/ref_path.cpp
    src/gtf.cpp
    src/gene_cn_kmer.cpp
    src/minimap2_align.cpp
    src/variant_call.cpp
    src/cosigt.cpp
    src/path_inversion.cpp
    src/module_fold.cpp
    src/describe.cpp
    src/panphorte.cpp
    src/panphorte_refine.cpp
    src/poa.cpp
    src/refine.cpp
    src/rebuild.cpp
    src/output.cpp
    src/background.cpp
    src/background_store.cpp
    src/background_graph.cpp
    src/background_sshash.cpp
    src/anchors.cpp
    src/fasta_index.cpp
    src/locus_index.cpp
    src/recruit_reads.cpp
    src/read_input.cpp
    src/recruit_scoring.cpp
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
target_link_libraries(panvarlib PRIVATE edlib abpoa minigraph sshash hts ZLIB::ZLIB)

add_executable(panvar
    src/main.cpp
    src/cli_utils.cpp
    src/bubble_command.cpp
    src/inspect.cpp
    src/call_command.cpp
    src/benchmark_command.cpp
    src/describe_command.cpp
    src/panphorte_command.cpp
    src/refine_command.cpp
    src/rebuild_command.cpp
    src/associate_command.cpp
    src/background_command.cpp
    src/recruit_command.cpp
)
target_include_directories(panvar PRIVATE ${EIGEN_DIR})
target_compile_definitions(panvar PRIVATE PANVAR_VERSION="${PROJECT_VERSION}")
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

# panvar ships as a command-line tool, so only the executable is installed. panvarlib is not exported:
# it links IMPORTED vendored targets (abpoa, minigraph) that have no package config, so a consumer's
# find_package(panvar) could not resolve them.
install(TARGETS panvar RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

if(PANVAR_BUILD_TESTING)
    include(CTest)
    if(BUILD_TESTING)
        # Unit tests of the shared k-mer primitives: brute force on random sequences, and the existing
        # collector on the real loci's spelled paths.
        add_executable(panvar_kmer_test tests/unit/kmer_test.cpp)
        target_link_libraries(panvar_kmer_test PRIVATE panvarlib minimap2 Threads::Threads ZLIB::ZLIB)
        target_compile_options(panvar_kmer_test PRIVATE -Wall -Wextra -Wpedantic)
        add_test(
            NAME panvar_unit_kmer
            COMMAND panvar_kmer_test
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/real_data/paralogs/c4.gfa.gz
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/real_data/paralogs/gstm.gfa.gz
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/real_data/paralogs/cyp2d.gfa.gz
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/real_data/paralogs/acot.gfa.gz
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/real_data/tandem_repeat/lpa.gfa.gz
        )
        set_tests_properties(panvar_unit_kmer PROPERTIES TIMEOUT 300)
        add_executable(panvar_background_test tests/unit/background_test.cpp)
        target_link_libraries(panvar_background_test PRIVATE panvarlib minimap2 Threads::Threads ZLIB::ZLIB)
        target_compile_options(panvar_background_test PRIVATE -Wall -Wextra -Wpedantic)
        add_test(NAME panvar_unit_background
                 COMMAND panvar_background_test ${CMAKE_CURRENT_BINARY_DIR}/unit_background)
        set_tests_properties(panvar_unit_background PROPERTIES TIMEOUT 300)
        # background end to end: kills at every commit step, concurrent adds, segments attached, range
        # merges, compaction, anchoring after the fact, a full disk; dumps against an independent count.
        set(PANVAR_BG_CLI ${CMAKE_CURRENT_SOURCE_DIR}/tests/suites/background_cli.sh)
        if(EXISTS ${PANVAR_BG_CLI})
            add_test(NAME panvar_background_cli
                     COMMAND /usr/bin/env bash ${PANVAR_BG_CLI} $<TARGET_FILE:panvar> ${CMAKE_CURRENT_SOURCE_DIR}/tests/tools ${CMAKE_CURRENT_BINARY_DIR}/background_cli)
            set_tests_properties(panvar_background_cli PROPERTIES TIMEOUT 600)
        endif()
        # recruit --index-only on a synthetic locus with a paralog: outside counts against brute force,
        # classes, decoys, screen, reuse and rebuilds.
        set(PANVAR_RECRUIT_INDEX ${CMAKE_CURRENT_SOURCE_DIR}/tests/suites/recruit_index.sh)
        if(EXISTS ${PANVAR_RECRUIT_INDEX})
            add_test(NAME panvar_recruit_index
                     COMMAND /usr/bin/env bash ${PANVAR_RECRUIT_INDEX} $<TARGET_FILE:panvar> ${CMAKE_CURRENT_SOURCE_DIR}/tests/tools ${CMAKE_CURRENT_BINARY_DIR}/recruit_index)
            set_tests_properties(panvar_recruit_index PROPERTIES TIMEOUT 300)
            add_test(NAME panvar_recruit_input
                     COMMAND /usr/bin/env bash ${CMAKE_CURRENT_SOURCE_DIR}/tests/suites/recruit_input.sh $<TARGET_FILE:panvar> ${CMAKE_CURRENT_SOURCE_DIR}/tests/tools ${CMAKE_CURRENT_BINARY_DIR}/recruit_input)
            set_tests_properties(panvar_recruit_input PROPERTIES TIMEOUT 300)
        endif()
        # Fast, dependency-free smoke on a tiny hand-built graph that exercises every
        # call event type + confounders with exact-record assertions (no vg needed).
        set(PANVAR_SYN_SMOKE ${CMAKE_CURRENT_SOURCE_DIR}/tests/suites/synthetic_smoke.sh)
        if(EXISTS ${PANVAR_SYN_SMOKE})
            add_test(
                NAME panvar_smoke_synthetic
                COMMAND /usr/bin/env bash ${PANVAR_SYN_SMOKE} $<TARGET_FILE:panvar> ${CMAKE_CURRENT_SOURCE_DIR}/tests/synthetic_data ${CMAKE_CURRENT_BINARY_DIR}/smoke_synthetic
            )
        endif()

        # bubble -> panphorte -> call on a generated locus with four planted SVs: sites, positions,
        # carriers and copy numbers against the planted truth. Seconds.
        set(PANVAR_SYN_BENCH ${CMAKE_CURRENT_SOURCE_DIR}/tests/suites/synthetic_bench.sh)
        if(EXISTS ${PANVAR_SYN_BENCH})
            add_test(
                NAME panvar_synthetic_bench
                COMMAND /usr/bin/env bash ${PANVAR_SYN_BENCH} $<TARGET_FILE:panvar> ${CMAKE_CURRENT_BINARY_DIR}/synthetic_bench
            )
            set_tests_properties(panvar_synthetic_bench PROPERTIES TIMEOUT 120)
        endif()

        # Deterministic statistical assertions for `associate`: linear p against R's lm() (a t tail,
        # not a normal one, at small df), the logistic separation contract, case/control carrier
        # counts, and GRM rejection. No RNG, no R at run time, about a second.
        set(PANVAR_ASSOC_STATS ${CMAKE_CURRENT_SOURCE_DIR}/tests/suites/associate_stats.sh)
        if(EXISTS ${PANVAR_ASSOC_STATS})
            add_test(
                NAME panvar_associate_stats
                COMMAND /usr/bin/env bash ${PANVAR_ASSOC_STATS} $<TARGET_FILE:panvar> ${CMAKE_CURRENT_BINARY_DIR}/associate_stats
            )
            set_tests_properties(panvar_associate_stats PROPERTIES TIMEOUT 120)
        endif()

        # bubble contract on hand-checkable graphs, run through the default internal snarl finder and
        # through --snarls.
        set(PANVAR_BUBBLE_STATS ${CMAKE_CURRENT_SOURCE_DIR}/tests/suites/bubble_stats.sh)
        if(EXISTS ${PANVAR_BUBBLE_STATS})
            add_test(
                NAME panvar_bubble_stats
                COMMAND /usr/bin/env bash ${PANVAR_BUBBLE_STATS} $<TARGET_FILE:panvar> ${CMAKE_CURRENT_BINARY_DIR}/bubble_stats
            )
            set_tests_properties(panvar_bubble_stats PROPERTIES TIMEOUT 180)
        endif()

        # inspect contract: the contents of its outputs (FASTA alleles, path names, input checks), which
        # the smoke tests only check exist.
        set(PANVAR_INSPECT_STATS ${CMAKE_CURRENT_SOURCE_DIR}/tests/suites/inspect_stats.sh)
        if(EXISTS ${PANVAR_INSPECT_STATS})
            add_test(
                NAME panvar_inspect_stats
                COMMAND /usr/bin/env bash ${PANVAR_INSPECT_STATS} $<TARGET_FILE:panvar> ${CMAKE_CURRENT_BINARY_DIR}/inspect_stats
            )
            set_tests_properties(panvar_inspect_stats PROPERTIES TIMEOUT 180)
        endif()

        # Opt-in regressions against the committed real loci, for contracts only real data reaches
        # (overlapping sites, folding at hundreds of haplotypes, the POA guards, re-reading a rebuilt
        # GFA). A few minutes, so it carries the slow label.
        set(PANVAR_REAL_REGRESSIONS ${CMAKE_CURRENT_SOURCE_DIR}/tests/suites/real_regressions.sh)
        if(EXISTS ${PANVAR_REAL_REGRESSIONS})
            add_test(
                NAME panvar_real_regressions
                COMMAND /usr/bin/env bash ${PANVAR_REAL_REGRESSIONS} $<TARGET_FILE:panvar> ${CMAKE_CURRENT_BINARY_DIR}/real_regressions
            )
            set_tests_properties(panvar_real_regressions PROPERTIES TIMEOUT 900 LABELS "real;slow")
            if(NOT PANVAR_SLOW_TESTS)
                set_tests_properties(panvar_real_regressions PROPERTIES DISABLED TRUE)
            endif()
        endif()

        set(PANVAR_CALL_STATS ${CMAKE_CURRENT_SOURCE_DIR}/tests/suites/call_stats.sh)
        if(EXISTS ${PANVAR_CALL_STATS})
            add_test(
                NAME panvar_call_stats
                COMMAND /usr/bin/env bash ${PANVAR_CALL_STATS} $<TARGET_FILE:panvar> ${CMAKE_CURRENT_BINARY_DIR}/call_stats
            )
            set_tests_properties(panvar_call_stats PROPERTIES TIMEOUT 180)
        endif()

        # call invariants that the real loci never exercise: a bubble row against the reference
        # direction, P-line order, node renumbering, single-linkage chaining, EVENTID on a merged
        # substitution, a short self-loop REP. Synthetic and fast. Cases whose fix is a later plan item
        # print xfail/XPASS and do not fail the test.
        set(PANVAR_CALL_INVARIANTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/suites/call_invariants.sh)
        if(EXISTS ${PANVAR_CALL_INVARIANTS})
            add_test(
                NAME panvar_call_invariants
                COMMAND /usr/bin/env bash ${PANVAR_CALL_INVARIANTS} $<TARGET_FILE:panvar> ${CMAKE_CURRENT_BINARY_DIR}/call_invariants
            )
            set_tests_properties(panvar_call_invariants PROPERTIES TIMEOUT 180)
        endif()

        # ORIGINAL vs CURRENT call on the six real loci, field by field plus BIMBAM, with every
        # difference listed in tests/suites/call_regress.expected.tsv. Needs a pinned original build, e.g.
        # `cmake -DPANVAR_ORIG_BIN=/path/to/original/panvar -DPANVAR_SLOW_TESTS=ON`. Three input
        # sets (base, refined, altref) built from scratch take about two minutes, so it is gated like
        # panvar_real_regressions.
        set(PANVAR_ORIG_BIN "" CACHE FILEPATH
            "Pinned ORIGINAL panvar binary for panvar_call_regress (the test is disabled while unset)")
        set(PANVAR_CALL_REGRESS ${CMAKE_CURRENT_SOURCE_DIR}/tests/suites/call_regress.sh)
        if(EXISTS ${PANVAR_CALL_REGRESS})
            add_test(
                NAME panvar_call_regress
                COMMAND /usr/bin/env bash ${PANVAR_CALL_REGRESS} "${PANVAR_ORIG_BIN}" $<TARGET_FILE:panvar> ${CMAKE_CURRENT_BINARY_DIR}/call_regress
            )
            set_tests_properties(panvar_call_regress PROPERTIES TIMEOUT 900 LABELS "real;slow")
            if(NOT PANVAR_SLOW_TESTS OR NOT PANVAR_ORIG_BIN)
                set_tests_properties(panvar_call_regress PROPERTIES DISABLED TRUE)
            endif()
        endif()

        # benchmark contract: what the reported numbers mean on hand-checkable cases, which the smoke
        # test does not assert.
        set(PANVAR_BENCHMARK_STATS ${CMAKE_CURRENT_SOURCE_DIR}/tests/suites/benchmark_stats.sh)
        if(EXISTS ${PANVAR_BENCHMARK_STATS})
            add_test(
                NAME panvar_benchmark_stats
                COMMAND /usr/bin/env bash ${PANVAR_BENCHMARK_STATS} $<TARGET_FILE:panvar> ${CMAKE_CURRENT_BINARY_DIR}/benchmark_stats
            )
            set_tests_properties(panvar_benchmark_stats PROPERTIES TIMEOUT 180)
        endif()

        # describe contract: the feature dosages, which the smoke test does not assert.
        set(PANVAR_DESCRIBE_STATS ${CMAKE_CURRENT_SOURCE_DIR}/tests/suites/describe_stats.sh)
        if(EXISTS ${PANVAR_DESCRIBE_STATS})
            add_test(
                NAME panvar_describe_stats
                COMMAND /usr/bin/env bash ${PANVAR_DESCRIBE_STATS} $<TARGET_FILE:panvar> ${CMAKE_CURRENT_BINARY_DIR}/describe_stats
            )
            set_tests_properties(panvar_describe_stats PROPERTIES TIMEOUT 300)
        endif()

        # panphorte contract: folding itself, which the smoke test's graph does not exercise.
        set(PANVAR_PANPHORTE_STATS ${CMAKE_CURRENT_SOURCE_DIR}/tests/suites/panphorte_stats.sh)
        if(EXISTS ${PANVAR_PANPHORTE_STATS})
            add_test(
                NAME panvar_panphorte_stats
                COMMAND /usr/bin/env bash ${PANVAR_PANPHORTE_STATS} $<TARGET_FILE:panvar> ${CMAKE_CURRENT_BINARY_DIR}/panphorte_stats
            )
            set_tests_properties(panvar_panphorte_stats PROPERTIES TIMEOUT 300)
        endif()

        # refine contract: whether the rewrite is correct, which the smoke test does not assert.
        set(PANVAR_REFINE_STATS ${CMAKE_CURRENT_SOURCE_DIR}/tests/suites/refine_stats.sh)
        if(EXISTS ${PANVAR_REFINE_STATS})
            add_test(
                NAME panvar_refine_stats
                COMMAND /usr/bin/env bash ${PANVAR_REFINE_STATS} $<TARGET_FILE:panvar> ${CMAKE_CURRENT_BINARY_DIR}/refine_stats
            )
            set_tests_properties(panvar_refine_stats PROPERTIES TIMEOUT 300)
        endif()

        # rebuild's module fold on a synthetic copy-number module: every path exact, one pass per copy.
        # Seconds, so on by default (the real-data rebuild suite below is the slow one).
        set(PANVAR_MODULE_FOLD_STATS ${CMAKE_CURRENT_SOURCE_DIR}/tests/suites/module_fold_stats.sh)
        if(EXISTS ${PANVAR_MODULE_FOLD_STATS})
            add_test(
                NAME panvar_module_fold_stats
                COMMAND /usr/bin/env bash ${PANVAR_MODULE_FOLD_STATS} $<TARGET_FILE:panvar> ${CMAKE_CURRENT_BINARY_DIR}/module_fold_stats
            )
            set_tests_properties(panvar_module_fold_stats PROPERTIES TIMEOUT 300)
        endif()

        # rebuild contract: the input is never destroyed, a healthy graph survives, a malformed one is
        # refused rather than silently repaired, and nothing is left behind.
        set(PANVAR_REBUILD_STATS ${CMAKE_CURRENT_SOURCE_DIR}/tests/suites/rebuild_stats.sh)
        if(EXISTS ${PANVAR_REBUILD_STATS})
            add_test(
                NAME panvar_rebuild_stats
                COMMAND /usr/bin/env bash ${PANVAR_REBUILD_STATS} $<TARGET_FILE:panvar> ${CMAKE_CURRENT_BINARY_DIR}/rebuild_stats ${CMAKE_CURRENT_SOURCE_DIR}/tests/real_data/paralogs/c4.gfa.gz
            )
            # Several forced rebuilds of a real locus (reproducibility across threads, k either side
            # of 31, the reference-ambiguity and verdict cases). It drives minigraph and takes minutes,
            # so it is off by default: enable with `cmake -DPANVAR_SLOW_TESTS=ON`, then `ctest -L slow`.
            # Run it before anything that touches rebuild.
            set_tests_properties(panvar_rebuild_stats PROPERTIES TIMEOUT 900 LABELS "real;slow")
            if(NOT PANVAR_SLOW_TESTS)
                set_tests_properties(panvar_rebuild_stats PROPERTIES DISABLED TRUE)
            endif()
        endif()

        # Real-data integration smoke: the whole bubble -> inspect -> describe -> panphorte
        # -> call pipeline on the smallest committed locus (seconds, no external tools), asserting each
        # stage's outputs. Labelled "real" so `ctest -LE real` can skip it; other graphs run by hand:
        # bash tests/suites/real_smoke.sh ./build/panvar <gfa> <out>.
        set(PANVAR_REAL_SMOKE ${CMAKE_CURRENT_SOURCE_DIR}/tests/suites/real_smoke.sh)
        set(PANVAR_REAL_SMOKE_GFA ${CMAKE_CURRENT_SOURCE_DIR}/tests/real_data/paralogs/c4.gfa.gz)
        if(EXISTS ${PANVAR_REAL_SMOKE} AND EXISTS ${PANVAR_REAL_SMOKE_GFA})
            add_test(
                NAME panvar_smoke_real_c4
                COMMAND /usr/bin/env bash ${PANVAR_REAL_SMOKE} $<TARGET_FILE:panvar> ${PANVAR_REAL_SMOKE_GFA} ${CMAKE_CURRENT_BINARY_DIR}/smoke_real_c4
            )
            set_tests_properties(panvar_smoke_real_c4 PROPERTIES LABELS "real" TIMEOUT 300)
        endif()
    endif()
endif()
