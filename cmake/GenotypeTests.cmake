# Do not silently pass a release suite with missing scripts, oracles, or fixtures.
foreach(required IN ITEMS
    "tests/associate_null.sh"
    "tests/associate_stats.sh"
    "tests/benchmark_stats.sh"
    "tests/bubble_stats.sh"
    "tests/call_stats.sh"
    "tests/describe_stats.sh"
    "tests/fixtures/profiles/genotype_c4.tsv"
    "tests/fixtures/profiles/genotype_c4_v6_campaign_local.tsv"
    "tests/fixtures/profiles/genotype_experimental_sr_v1_acot.tsv"
    "tests/fixtures/profiles/genotype_experimental_sr_v1_c4.tsv"
    "tests/fixtures/profiles/genotype_experimental_sr_v1_cyp2d6.tsv"
    "tests/gwas/run_lpa_real.sh"
    "tests/gwas/validate_gemma.sh"
    "tests/genotype_smoke.sh"
    "tests/genotype_cli.sh"
    "tests/genotype_index.sh"
    "tests/genotype_assessment.sh"
    "tests/genotype_dosage.sh"
    "tests/genotype_resources.sh"
    "tests/genotype_exactness.sh"
    "tests/real_data/c4.gfa.gz"
    "tests/helpers/assessment/assessment_spec.py"
    "tests/helpers/assessment/compare.py"
    "tests/helpers/documentation_test.py"
    "tests/helpers/compound_annotation_test.cpp"
    "tests/helpers/genotype_assessment_main.cpp"
    "tests/helpers/genotype_batch_resource.cpp"
    "tests/helpers/genotype_block_delta.cpp"
    "tests/helpers/genotype_block_quality.cpp"
    "tests/helpers/genotype_block_support.cpp"
    "tests/helpers/genotype_blocked_fine.cpp"
    "tests/helpers/genotype_blocked_union.cpp"
    "tests/helpers/genotype_branch_lifetime.cpp"
    "tests/helpers/genotype_cli_contract.cpp"
    "tests/helpers/genotype_command_e2e.cpp"
    "tests/helpers/genotype_copying_transitions.cpp"
    "tests/helpers/genotype_copying_transitions_e2e.cpp"
    "tests/helpers/genotype_curation_contract.cpp"
    "tests/helpers/genotype_depth_anchors_na19240.inc"
    "tests/helpers/genotype_depth_anchors_na20355.inc"
    "tests/helpers/genotype_dosage.cpp"
    "tests/helpers/genotype_fragment_oracle.cpp"
    "tests/helpers/genotype_fragment_scorer.cpp"
    "tests/helpers/genotype_grouped_beam.cpp"
    "tests/helpers/genotype_headline.cpp"
    "tests/helpers/genotype_joint_search.cpp"
    "tests/helpers/genotype_joint_search_e2e.cpp"
    "tests/helpers/genotype_locus_adapter.cpp"
    "tests/helpers/genotype_mosaic_core.cpp"
    "tests/helpers/genotype_mosaic_proposal.cpp"
    "tests/helpers/genotype_prepared_locus_file.cpp"
    "tests/helpers/genotype_progressive_kernel.cpp"
    "tests/helpers/genotype_review_regressions.cpp"
    "tests/helpers/genotype_review_regressions.expected"
    "tests/helpers/genotype_scoring_context.cpp"
    "tests/helpers/genotype_seeding_exactness.cpp"
    "tests/helpers/genotype_sequence_report.cpp"
    "tests/helpers/genotype_union_placement_cache.cpp"
    "tests/helpers/junction_proposals_test.cpp"
    "tests/helpers/library_estimate.cpp"
    "tests/helpers/library_sensitivity.cpp"
    "tests/helpers/report_test.py"
    "tests/helpers/source_provenance.sh"
    "tests/helpers/source_tree_clean.sh"
    "tests/inspect_stats.sh"
    "tests/panphorte_stats.sh"
    "tests/real_regressions.sh"
    "tests/real_smoke.sh"
    "tests/rebuild_stats.sh"
    "tests/refine_stats.sh"
    "tests/synthetic_bench.sh"
    "tests/synthetic_smoke.sh"
)
    if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${required}")
        message(FATAL_ERROR "Required release test file is missing: ${required}")
    endif()
endforeach()

# Genotyper release regression tests. No campaign data or general research probe is required.
add_executable(genotype_dosage_test tests/helpers/genotype_dosage.cpp src/cli_utils.cpp)
target_link_libraries(genotype_dosage_test PRIVATE panvarlib)
target_compile_options(genotype_dosage_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME panvar_genotype_dosage COMMAND genotype_dosage_test)

add_executable(library_estimate_test tests/helpers/library_estimate.cpp src/cli_utils.cpp)
target_link_libraries(library_estimate_test PRIVATE panvarlib)
target_compile_options(library_estimate_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME panvar_library_estimate COMMAND library_estimate_test)

add_executable(genotype_mosaic_core_test tests/helpers/genotype_mosaic_core.cpp)
target_link_libraries(genotype_mosaic_core_test PRIVATE panvarlib)
target_compile_options(genotype_mosaic_core_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME panvar_genotype_mosaic_core COMMAND genotype_mosaic_core_test)
add_executable(genotype_joint_search_test tests/helpers/genotype_joint_search.cpp)
target_link_libraries(genotype_joint_search_test PRIVATE panvarlib)
target_compile_options(genotype_joint_search_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME panvar_genotype_joint_search COMMAND genotype_joint_search_test)
add_executable(genotype_block_support_test tests/helpers/genotype_block_support.cpp)
target_link_libraries(genotype_block_support_test PRIVATE panvarlib)
target_compile_options(genotype_block_support_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME panvar_genotype_block_support COMMAND genotype_block_support_test)
add_executable(genotype_joint_search_e2e_test tests/helpers/genotype_joint_search_e2e.cpp src/cli_utils.cpp)
target_link_libraries(genotype_joint_search_e2e_test PRIVATE panvarlib)
target_compile_options(genotype_joint_search_e2e_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME panvar_genotype_joint_search_e2e COMMAND genotype_joint_search_e2e_test)

add_executable(genotype_copying_transitions_test tests/helpers/genotype_copying_transitions.cpp src/cli_utils.cpp)
target_link_libraries(genotype_copying_transitions_test PRIVATE panvarlib)
target_compile_options(genotype_copying_transitions_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME panvar_genotype_copying_transitions COMMAND genotype_copying_transitions_test)
add_executable(genotype_copying_transitions_e2e_test tests/helpers/genotype_copying_transitions_e2e.cpp src/cli_utils.cpp)
target_link_libraries(genotype_copying_transitions_e2e_test PRIVATE panvarlib)
target_compile_options(genotype_copying_transitions_e2e_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME panvar_genotype_copying_transitions_e2e COMMAND genotype_copying_transitions_e2e_test)

add_executable(genotype_locus_adapter_test tests/helpers/genotype_locus_adapter.cpp src/cli_utils.cpp)
target_link_libraries(genotype_locus_adapter_test PRIVATE panvarlib)
target_compile_options(genotype_locus_adapter_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME panvar_genotype_locus_adapter
         COMMAND genotype_locus_adapter_test
                 ${CMAKE_CURRENT_SOURCE_DIR}/tests/real_data/c4.gfa.gz
                 ${CMAKE_CURRENT_BINARY_DIR}/genotype_locus_c4)

add_executable(genotype_command_e2e_test tests/helpers/genotype_command_e2e.cpp src/cli_utils.cpp)
target_link_libraries(genotype_command_e2e_test PRIVATE panvarlib)
target_compile_options(genotype_command_e2e_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME panvar_genotype_command_e2e
         COMMAND genotype_command_e2e_test ${CMAKE_CURRENT_BINARY_DIR}/genotype_e2e ${CMAKE_CURRENT_SOURCE_DIR})
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/genotype_e2e)
add_test(NAME panvar_genotype_blocked_union_route
         COMMAND genotype_command_e2e_test ${CMAKE_CURRENT_BINARY_DIR}/genotype_e2e_blocked ${CMAKE_CURRENT_SOURCE_DIR})
set_tests_properties(panvar_genotype_blocked_union_route PROPERTIES ENVIRONMENT "PANVAR_E2E_ONLY=blocked_union")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/genotype_e2e_blocked)

add_executable(genotype_progressive_kernel_test tests/helpers/genotype_progressive_kernel.cpp)
target_link_libraries(genotype_progressive_kernel_test PRIVATE panvarlib)
target_compile_options(genotype_progressive_kernel_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME panvar_genotype_progressive_kernel COMMAND genotype_progressive_kernel_test)
add_executable(genotype_blocked_union_test tests/helpers/genotype_blocked_union.cpp)
target_link_libraries(genotype_blocked_union_test PRIVATE panvarlib)
target_compile_options(genotype_blocked_union_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME panvar_genotype_blocked_union COMMAND genotype_blocked_union_test)
add_executable(genotype_blocked_fine_test tests/helpers/genotype_blocked_fine.cpp)
target_link_libraries(genotype_blocked_fine_test PRIVATE panvarlib)
target_compile_options(genotype_blocked_fine_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME panvar_genotype_blocked_fine COMMAND genotype_blocked_fine_test)

   # M10-C diagnostic

   # M10-B diagnostic

add_executable(genotype_review_regressions tests/helpers/genotype_review_regressions.cpp src/cli_utils.cpp)
target_link_libraries(genotype_review_regressions PRIVATE panvarlib)
target_compile_options(genotype_review_regressions PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME panvar_genotype_review_regressions
         COMMAND genotype_review_regressions ${CMAKE_CURRENT_SOURCE_DIR}/tests/helpers/genotype_review_regressions.expected)

   # M13-A diagnostic gates (not a ctest)

   # M12-A diagnostic fixtures (not a ctest)

add_executable(genotype_fragment_oracle_test tests/helpers/genotype_fragment_oracle.cpp)
target_link_libraries(genotype_fragment_oracle_test PRIVATE panvarlib)
target_compile_options(genotype_fragment_oracle_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME panvar_genotype_fragment_oracle COMMAND genotype_fragment_oracle_test)

add_executable(genotype_batch_resource_test tests/helpers/genotype_batch_resource.cpp src/cli_utils.cpp)
target_link_libraries(genotype_batch_resource_test PRIVATE panvarlib)
target_compile_options(genotype_batch_resource_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME panvar_genotype_batch_resource COMMAND genotype_batch_resource_test)
add_executable(compound_annotation_test tests/helpers/compound_annotation_test.cpp src/cli_utils.cpp)
target_link_libraries(compound_annotation_test PRIVATE panvarlib)
target_compile_options(compound_annotation_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME panvar_compound_annotation COMMAND compound_annotation_test)
add_executable(junction_proposals_test tests/helpers/junction_proposals_test.cpp)
target_link_libraries(junction_proposals_test PRIVATE panvarlib)
target_compile_options(junction_proposals_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME panvar_junction_proposals COMMAND junction_proposals_test)
add_executable(genotype_sequence_report_test tests/helpers/genotype_sequence_report.cpp src/cli_utils.cpp)
target_link_libraries(genotype_sequence_report_test PRIVATE panvarlib)
target_compile_options(genotype_sequence_report_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME panvar_genotype_sequence_report
         COMMAND genotype_sequence_report_test ${CMAKE_CURRENT_BINARY_DIR}/genotype_sequence_report_out)
add_executable(genotype_prepared_locus_file_test tests/helpers/genotype_prepared_locus_file.cpp src/cli_utils.cpp)
target_link_libraries(genotype_prepared_locus_file_test PRIVATE panvarlib)
target_compile_options(genotype_prepared_locus_file_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME panvar_genotype_prepared_locus_file
         COMMAND genotype_prepared_locus_file_test ${CMAKE_CURRENT_BINARY_DIR}/prepared_locus_file)

add_executable(genotype_fragment_scorer_test tests/helpers/genotype_fragment_scorer.cpp)
target_link_libraries(genotype_fragment_scorer_test PRIVATE panvarlib)
target_compile_options(genotype_fragment_scorer_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME panvar_genotype_fragment_scorer COMMAND genotype_fragment_scorer_test)

add_executable(genotype_scoring_context_test tests/helpers/genotype_scoring_context.cpp)
target_link_libraries(genotype_scoring_context_test PRIVATE panvarlib)
target_compile_options(genotype_scoring_context_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME panvar_genotype_scoring_context COMMAND genotype_scoring_context_test)

add_executable(genotype_block_delta_test tests/helpers/genotype_block_delta.cpp)
target_link_libraries(genotype_block_delta_test PRIVATE panvarlib)
target_compile_options(genotype_block_delta_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME panvar_genotype_block_delta COMMAND genotype_block_delta_test)
add_test(NAME panvar_genotype_partition_accounts COMMAND genotype_block_delta_test)
set_tests_properties(panvar_genotype_partition_accounts PROPERTIES ENVIRONMENT "PANVAR_BLOCK_DELTA_ONLY=partition_accounts")
add_test(NAME panvar_genotype_round1_admission COMMAND genotype_block_delta_test)
set_tests_properties(panvar_genotype_round1_admission PROPERTIES ENVIRONMENT "PANVAR_BLOCK_DELTA_ONLY=admission")
add_executable(genotype_union_placement_cache_test tests/helpers/genotype_union_placement_cache.cpp)
target_link_libraries(genotype_union_placement_cache_test PRIVATE panvarlib)
target_compile_options(genotype_union_placement_cache_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME panvar_genotype_union_placement_cache COMMAND genotype_union_placement_cache_test)
add_executable(genotype_seeding_exactness_test tests/helpers/genotype_seeding_exactness.cpp)
target_link_libraries(genotype_seeding_exactness_test PRIVATE panvarlib)
target_compile_options(genotype_seeding_exactness_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME panvar_genotype_seeding_exactness COMMAND genotype_seeding_exactness_test)
add_executable(genotype_headline_test tests/helpers/genotype_headline.cpp src/cli_utils.cpp)
target_link_libraries(genotype_headline_test PRIVATE panvarlib)
target_compile_options(genotype_headline_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME panvar_genotype_headline COMMAND genotype_headline_test)
add_executable(genotype_assessment_test tests/helpers/genotype_assessment_main.cpp src/cli_utils.cpp)
target_link_libraries(genotype_assessment_test PRIVATE panvarlib)
target_compile_options(genotype_assessment_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME panvar_genotype_assessment
         COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/tests/helpers/assessment/compare.py
                 $<TARGET_FILE:genotype_assessment_test> ${CMAKE_CURRENT_BINARY_DIR}/genotype_assessment)
add_executable(genotype_cli_contract_test tests/helpers/genotype_cli_contract.cpp src/cli_utils.cpp)
target_link_libraries(genotype_cli_contract_test PRIVATE panvarlib)
target_compile_options(genotype_cli_contract_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME panvar_source_provenance
         COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/tests/helpers/source_provenance.sh ${CMAKE_COMMAND}
                 ${CMAKE_CURRENT_SOURCE_DIR}/cmake/source_provenance.cmake ${CMAKE_CURRENT_BINARY_DIR}/source_provenance)
set_tests_properties(panvar_source_provenance PROPERTIES TIMEOUT 60)
add_test(NAME panvar_genotype_cli_contract
         COMMAND genotype_cli_contract_test ${CMAKE_CURRENT_BINARY_DIR}/genotype_cli_contract ${CMAKE_CURRENT_SOURCE_DIR})
add_executable(genotype_curation_contract_test tests/helpers/genotype_curation_contract.cpp src/cli_utils.cpp)
target_link_libraries(genotype_curation_contract_test PRIVATE panvarlib)
target_compile_options(genotype_curation_contract_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME panvar_genotype_curation_contract
         COMMAND genotype_curation_contract_test ${CMAKE_CURRENT_BINARY_DIR}/genotype_curation_contract)

add_executable(library_sensitivity_test tests/helpers/library_sensitivity.cpp src/cli_utils.cpp)
target_link_libraries(library_sensitivity_test PRIVATE panvarlib)
target_compile_options(library_sensitivity_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME panvar_library_sensitivity COMMAND library_sensitivity_test)

add_executable(genotype_block_quality_test tests/helpers/genotype_block_quality.cpp src/cli_utils.cpp)
target_link_libraries(genotype_block_quality_test PRIVATE panvarlib)
target_compile_options(genotype_block_quality_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME panvar_genotype_block_quality COMMAND genotype_block_quality_test)
add_executable(genotype_branch_lifetime_test tests/helpers/genotype_branch_lifetime.cpp)
target_link_libraries(genotype_branch_lifetime_test PRIVATE panvarlib)
target_compile_options(genotype_branch_lifetime_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME panvar_genotype_branch_lifetime COMMAND genotype_branch_lifetime_test)

add_executable(genotype_mosaic_proposal_test tests/helpers/genotype_mosaic_proposal.cpp)
target_link_libraries(genotype_mosaic_proposal_test PRIVATE panvarlib)
target_compile_options(genotype_mosaic_proposal_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME panvar_genotype_mosaic_proposal COMMAND genotype_mosaic_proposal_test)

add_executable(genotype_grouped_beam_test tests/helpers/genotype_grouped_beam.cpp)
target_link_libraries(genotype_grouped_beam_test PRIVATE panvarlib)
target_compile_options(genotype_grouped_beam_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME panvar_genotype_grouped_beam COMMAND genotype_grouped_beam_test)

# These are self-contained tests of the public reports, not historical campaign runners.
add_test(NAME panvar_genotype_report
         COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/tests/helpers/report_test.py ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
add_test(NAME panvar_genotype_source_tree_clean
         COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/tests/helpers/source_tree_clean.sh)
add_test(NAME panvar_genotype_documentation
         COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/tests/helpers/documentation_test.py ${CMAKE_CURRENT_SOURCE_DIR})

set_tests_properties(
    panvar_genotype_mosaic_core
    panvar_genotype_locus_adapter
    panvar_genotype_command_e2e
    panvar_genotype_blocked_union_route
    PROPERTIES LABELS "genotype;smoke")

set_tests_properties(
    panvar_genotype_cli_contract
    panvar_genotype_documentation
    panvar_genotype_curation_contract
    panvar_genotype_source_tree_clean
    PROPERTIES LABELS "genotype;cli")

set_tests_properties(
    panvar_genotype_prepared_locus_file
    PROPERTIES LABELS "genotype;index")

set_tests_properties(
    panvar_genotype_assessment
    panvar_genotype_headline
    panvar_genotype_block_quality
    panvar_genotype_block_support
    panvar_genotype_sequence_report
    panvar_genotype_report
    panvar_compound_annotation
    PROPERTIES LABELS "genotype;assessment")

set_tests_properties(
    panvar_genotype_dosage
    panvar_library_estimate
    panvar_library_sensitivity
    PROPERTIES LABELS "genotype;dosage")

set_tests_properties(
    panvar_genotype_batch_resource
    panvar_genotype_branch_lifetime
    panvar_genotype_partition_accounts
    panvar_genotype_round1_admission
    PROPERTIES LABELS "genotype;resources")

set_tests_properties(
    panvar_genotype_copying_transitions
    panvar_genotype_copying_transitions_e2e
    panvar_genotype_joint_search
    panvar_genotype_joint_search_e2e
    panvar_genotype_progressive_kernel
    panvar_genotype_blocked_union
    panvar_genotype_blocked_fine
    panvar_genotype_review_regressions
    panvar_genotype_fragment_oracle
    panvar_genotype_fragment_scorer
    panvar_genotype_scoring_context
    panvar_genotype_block_delta
    panvar_genotype_union_placement_cache
    panvar_genotype_seeding_exactness
    panvar_genotype_mosaic_proposal
    panvar_genotype_grouped_beam
    panvar_junction_proposals
    PROPERTIES LABELS "genotype;exactness")

set_tests_properties(panvar_source_provenance PROPERTIES LABELS "genotype;index")
