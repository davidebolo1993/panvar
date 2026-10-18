#!/usr/bin/env python3
"""S2-G: complete-run comparison between the ordinary caller run (bulk indexed union) and the integrated run
(union_fragment_blocks fine_core_4096_experimental), with EXPLICIT allowances only:
  * every output file byte-identical after path normalisation, except:
    - call.mosaic_status.tsv: keys may differ only if listed below (route provenance written by one route, timing,
      process resources, executed/served work); every other key -- scores, ranking, genotype, phase, filters, sequence,
      completion -- must be identical and present in both;
    - call.settings.tsv: identical except the integrated run's own union_fragment_blocks row;
    - call.repair_rounds.tsv: the scientific columns (0-27) identical (the rest are timing);
    - call.dosage_status.tsv: identical except seconds_* rows;
    - call.repair_memory.tsv: a memory/timing log of the refinement, not compared (declared).
  s2g_compare.py <control dir> <integrated dir> <out.tsv>
"""
import os
import re
import sys

INTEGRATED_ONLY_NAMES = {
    "union_route", "union_fragment_block_pairs", "union_origin_aggregate_policy", "union_evidence_layout", "union_store_layout",
    "union_origin_cache_bytes", "union_origin_whole_budget_bytes", "union_origin_fine_budget_bytes", "union_cap_physical_origins",
    "union_cap_state_classes", "union_cap_piece_hits", "union_cap_latent_alignment_pairs", "union_cap_mass_cells",
    "union_cap_cached_placements_consumed", "union_cap_cached_origins_consumed", "union_cap_cached_states_consumed",
    "union_cap_cached_prefilters_consumed", "union_caps_scope", "union_fragment_blocks_completed",
    "union_fragment_blocks_run_completed", "union_fragment_block_footprint_max_bytes", "seconds_union_block_combination",
    "footprint_checkpoints"}
INTEGRATED_ONLY_RE = [
    re.compile(r"^union_block[0-9]+_(fragments|executed_physical_origins|direct_physical_origins|served_physical_origins|"
               r"whole_origins_constructed|fine_origins_constructed|whole_entities_admitted|whole_entities_declined|"
               r"whole_held_bytes|fine_windows_eligible|fine_windows_admitted|fine_windows_declined|fine_charged_bytes|"
               r"fine_served_physical_origins|executed_state_classes|executed_dps|phys_footprint_after_bytes)$"),
    re.compile(r"^union_completed_blocks_(physical_origins|piece_hits|cached_placements_consumed|cached_origins_consumed|"
               r"cached_states_consumed|cached_prefilters_consumed|state_classes)$"),
    re.compile(r"^footprint_checkpoint[0-9]+_(stage|phys_footprint_bytes|phys_footprint_peak_bytes|seconds)$")]
ORDINARY_ONLY_RE = re.compile(r"^(indexed_union_|seconds_indexed_union_|worker_seconds_indexed_union_|threads_indexed_union_)")
TIMING_RE = re.compile(r"^(seconds_|worker_seconds_)|^elapsed_seconds$")
SHARED_NAMES = {
    # process resources
    "resources_phys_footprint_peak_bytes", "resources_rss_peak_bytes", "resources_footprint_at_repair_start_bytes",
    "union_end_phys_footprint_bytes", "threads_peak_active_workers",
    # the arms' profile identity (the integrated arm's extra profile carries the route row)
    "profile_sha256",
    # timing and resources of the S2-A early pass and checkpoint
    "checkpoint_seconds_first_evaluated_prediction", "checkpoint_seconds_committed", "early_incumbent_seconds_evaluation",
    "early_incumbent_seconds_preparation", "early_incumbent_seconds_selection", "early_incumbent_footprint_after_bytes",
    # the last scorer call's memory accounts
    "scorer_mass_matrix_bytes", "scorer_haplotype_sequence_bytes",
    # executed / served union work, summed over blocks (held bytes: the largest block)
    "candidate_bases_scanned", "start_sort_in_place", "start_sort_out_of_place", "start_sort_fallbacks",
    "threads_start_sort_peak_worker_bytes", "union_origin_cache_eligible_segments", "union_origin_cache_admitted_segments",
    "union_origin_cache_declined_segments", "union_origin_cache_hits", "union_origin_cache_origins_constructed",
    "union_origin_cache_payload_bytes", "union_origin_cache_peak_bytes", "union_cached_origins_consumed",
    "union_cached_states_consumed", "physical_origins", "statistical_state_classes", "latent_alignment_pairs", "valid_fr_states",
    "piece_hits", "full_read_verifications", "dp_evaluations", "proposed_starts", "unique_starts", "accepted_mate_placements",
    "reused_physical_origins", "reused_state_aggregates", "merged_state_classes",
    "union_prefilter_cache_window_contents", "union_prefilter_cache_window_index_bytes", "union_prefilter_cache_positions_indexed",
    "union_prefilter_cache_decision_table_bytes", "union_prefilter_cache_decisions_held",
    "union_prefilter_cache_decisions_not_retained", "union_cached_prefilters_consumed", "union_prefilter_cache_misses",
    # the union pass's emission-cache work (per block, summed; the DP count also varies with scheduling) and served work
    "emission_cache_dp_hits", "emission_cache_dp_misses", "emission_cache_dp_requests", "emission_cache_entries",
    "emission_cache_executed_dps", "emission_cache_executed_prefilters", "emission_cache_payload_bytes",
    "reused_latent_alignment_pairs", "reused_valid_fr_states",
}


def kv(p, d):
    out = {}
    for line in open(p, errors="replace").read().replace(d, "<D>").splitlines():
        a = line.split("\t", 1)
        if len(a) == 2:
            out[a[0]] = a[1]
    return out


def integrated_only(k):
    return k in INTEGRATED_ONLY_NAMES or any(r.match(k) for r in INTEGRATED_ONLY_RE)


def status_differences(a, b, da, db):
    x, y = kv(a, da), kv(b, db)
    bad, allowed = [], []
    for k in sorted(set(x) | set(y)):
        if k in x and k in y and x[k] == y[k]:
            continue
        ok = (k not in x and integrated_only(k)) or (k not in y and ORDINARY_ONLY_RE.search(k)) or \
             (k in x and k in y and (TIMING_RE.search(k) or k in SHARED_NAMES))
        (allowed if ok else bad).append(k)
    return bad, allowed


def main(a, b, out):
    skip = re.compile(r"^(profile\..*|run\.log)$")
    fa = {f for f in os.listdir(a) if not skip.match(f)}
    fb = {f for f in os.listdir(b) if not skip.match(f)}
    rows, bad = [], []
    if fa != fb:
        bad.append(f"file sets differ: {sorted(fa ^ fb)}")
    for f in sorted(fa & fb):
        pa, pb = os.path.join(a, f), os.path.join(b, f)
        if f == "call.repair_memory.tsv":
            rows.append(f"{f}\tnot_compared\tdeclared memory/timing log")
            continue
        if f == "call.mosaic_status.tsv":
            diff, allowed = status_differences(pa, pb, a, b)
            rows.append(f"{f}\tallowed_differing_keys\t{len(allowed)}")
            rows.append(f"{f}\tunexpected_differing_keys\t{','.join(diff) or 'none'}")
            if diff:
                bad.append("status keys outside the allowances: " + ",".join(diff))
            continue
        if f == "call.settings.tsv":
            x, y = kv(pa, a), kv(pb, b)
            y.pop("union_fragment_blocks", None)
            if x != y:
                bad.append("settings differ beyond the route row")
            continue
        if f == "call.repair_rounds.tsv":
            cut = lambda p: [line.split("\t")[:28] for line in open(p).read().splitlines()]  # noqa: E731
            if cut(pa) != cut(pb):
                bad.append("trajectory differs")
            continue
        if f == "call.dosage_status.tsv":
            x = {k: v for k, v in kv(pa, a).items() if not k.startswith("seconds_")}
            y = {k: v for k, v in kv(pb, b).items() if not k.startswith("seconds_")}
            if x != y:
                bad.append("dosage status differs")
            continue
        if open(pa, errors="replace").read().replace(a, "<D>") != open(pb, errors="replace").read().replace(b, "<D>"):
            bad.append(f"{f} differs")
    with open(out, "w") as o:
        o.write("key\tvalue\n")
        o.write(f"outputs_identical\t{'yes' if not bad else 'NO: ' + ' | '.join(bad)}\n")
        for r in rows:
            o.write(r + "\n")
    print(open(out).read(), end="")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main(os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2]), sys.argv[3]))
