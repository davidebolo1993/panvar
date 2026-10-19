#!/usr/bin/env python3
"""S2-L comparator v3.1 (Amendment 1): v3 with ONE correction -- union_fragment_block_footprint_max_bytes, which only the
blocked route writes (S2-G), is a route-only row (category route, one run only) instead of a both-runs operational row.
v3 (s2l_compare.py) and its three failed 5x results are unchanged.

S2-L comparator v3: complete-run scientific equivalence with EXPLICIT, categorised operational allowances.
Earlier comparators (s2g_compare, s2j_compare v1, s2j_compare_v2, compare_6gib) and their failed results are unchanged.

  s2l_compare.py <run A> <run B> <out.tsv> --declared-settings k1,k2,... [--categories c1,c2,...]

Protected (any difference fails): every scientific output file byte for byte after path normalisation (genotypes, block
genotypes, candidates, scores, equivalence set, locus genotype, sequence report, dosage files, library sensitivity), the
rounds table's scientific columns (0-27), every settings VALUE except the declared settings (and the key sets must be
equal), and every status key not named by an enabled category -- including configured caps, budgets and every
prefilter configuration row (budget, consumption cap, read and window length, on/off).

Categories (each an explicit list or anchored pattern):
  operational  timing (seconds_*, worker_seconds_*, elapsed_seconds), process resources, the last scorer call's memory
               accounts, scheduling-variable emission/DP counters, and the AUDITED union prefilter-cache counters only:
               union_cached_prefilters_consumed, union_completed_blocks_cached_prefilters_consumed,
               union_prefilter_cache_misses, union_prefilter_cache_decisions_held,
               union_prefilter_cache_decisions_not_retained (S2-K: the hit/miss split and table admission order depend on
               thread scheduling; decisions are pure functions of content, so science cannot change)
  profile      profile, profile_name, profile_sha256 -- accepted only if the settings check passed (the arms' effective
               settings are equal apart from the declared settings)
  route        the blocked union route's own rows (present in one run only) and the ordinary indexed batch rows (present
               in the other only), plus union work/cache rows that legitimately differ between bulk and blocked execution
  admission    the round-1 admission policy's rows (one run only), the refinement's executed-work totals and margins,
               and repair_origin_cache_bytes_end
  checkpoints  footprint_checkpoints and footprint_checkpoint<N>_* in either or both runs (indices shift between revisions)
  envelope     footprint_limit_bytes (the declared envelope)
"""
import argparse
import os
import re
import sys

TIMING_RE = re.compile(r"^(seconds_|worker_seconds_)|^elapsed_seconds$|^checkpoint_seconds_|^early_incumbent_seconds_")
OPERATIONAL = {
    "resources_phys_footprint_peak_bytes", "resources_rss_peak_bytes", "resources_footprint_at_repair_start_bytes",
    "union_end_phys_footprint_bytes", "threads_peak_active_workers", "early_incumbent_footprint_after_bytes",
    "scorer_mass_matrix_bytes", "scorer_haplotype_sequence_bytes",
    "emission_cache_dp_hits", "emission_cache_dp_misses", "emission_cache_dp_requests", "emission_cache_entries",
    "emission_cache_executed_dps", "emission_cache_executed_prefilters", "emission_cache_payload_bytes",
    "dp_evaluations", "full_read_verifications",
    "union_cached_prefilters_consumed", "union_completed_blocks_cached_prefilters_consumed", "union_prefilter_cache_misses",
    "union_prefilter_cache_decisions_held", "union_prefilter_cache_decisions_not_retained"}
OPERATIONAL_RE = re.compile(r"^union_block[0-9]+_(executed_dps|phys_footprint_after_bytes)$")
PROFILE = {"profile", "profile_name", "profile_sha256"}
ROUTE_ONLY = {
    "union_route", "union_fragment_block_pairs", "union_origin_aggregate_policy", "union_evidence_layout", "union_store_layout",
    "union_origin_cache_bytes", "union_origin_whole_budget_bytes", "union_origin_fine_budget_bytes", "union_cap_physical_origins",
    "union_cap_state_classes", "union_cap_piece_hits", "union_cap_latent_alignment_pairs", "union_cap_mass_cells",
    "union_cap_cached_placements_consumed", "union_cap_cached_origins_consumed", "union_cap_cached_states_consumed",
    "union_cap_cached_prefilters_consumed", "union_caps_scope", "union_fragment_blocks_completed",
    "union_fragment_blocks_run_completed", "union_fragment_block_footprint_max_bytes"}
ROUTE_ONLY_RE = re.compile(
    r"^union_block[0-9]+_(fragments|executed_physical_origins|direct_physical_origins|served_physical_origins|"
    r"whole_origins_constructed|fine_origins_constructed|whole_entities_admitted|whole_entities_declined|whole_held_bytes|"
    r"fine_windows_eligible|fine_windows_admitted|fine_windows_declined|fine_charged_bytes|fine_served_physical_origins|"
    r"executed_state_classes|executed_dps|phys_footprint_after_bytes)$|"
    r"^union_completed_blocks_(physical_origins|piece_hits|cached_placements_consumed|cached_origins_consumed|"
    r"cached_states_consumed|cached_prefilters_consumed|state_classes)$|^seconds_union_block_combination$")
ORDINARY_ONLY_RE = re.compile(r"^(indexed_union_|seconds_indexed_union_|worker_seconds_indexed_union_|threads_indexed_union_)")
ROUTE_WORK = {  # union execution work that legitimately differs between bulk and blocked execution (S2-G)
    "candidate_bases_scanned", "start_sort_in_place", "start_sort_out_of_place", "start_sort_fallbacks",
    "threads_start_sort_peak_worker_bytes", "union_origin_cache_eligible_segments", "union_origin_cache_admitted_segments",
    "union_origin_cache_declined_segments", "union_origin_cache_hits", "union_origin_cache_origins_constructed",
    "union_origin_cache_payload_bytes", "union_origin_cache_peak_bytes", "union_cached_origins_consumed",
    "union_cached_states_consumed", "physical_origins", "statistical_state_classes", "latent_alignment_pairs", "valid_fr_states",
    "piece_hits", "proposed_starts", "unique_starts", "accepted_mate_placements", "reused_physical_origins",
    "reused_state_aggregates", "merged_state_classes", "reused_latent_alignment_pairs", "reused_valid_fr_states",
    "union_prefilter_cache_window_contents", "union_prefilter_cache_window_index_bytes", "union_prefilter_cache_positions_indexed",
    "union_prefilter_cache_decision_table_bytes"}
ADMISSION_ONLY = {
    "repair_round1_origin_admission", "repair_admission_restricted_plans", "repair_admission_planned_windows",
    "repair_admission_singleton_windows", "repair_admission_eligible_windows", "repair_admission_misses_while_restricted",
    "repair_admission_singleton_misses_not_retained", "repair_admission_retained_while_restricted",
    "repair_admission_hits_while_restricted", "repair_admission_group_reuses_while_restricted"}
ADMISSION_WORK = {
    "repair_cumulative_latent_pairs", "repair_cumulative_physical_origins", "repair_cumulative_piece_hits",
    "repair_margin_latent_pairs", "repair_margin_physical_origins", "repair_margin_piece_hits", "repair_origin_cache_bytes_end"}
CHECKPOINT_RE = re.compile(r"^footprint_checkpoint([0-9]+_(stage|phys_footprint_bytes|phys_footprint_peak_bytes|seconds|"
                           r"allocator_in_use_bytes)|s)$")
ENVELOPE = {"footprint_limit_bytes"}
SCIENTIFIC_SKIP = re.compile(r"^(profile\..*|run\.log|call\.refinement_.*|call\.repair_memory\.tsv)$")


def kv(path, d):
    out = {}
    for line in open(path, errors="replace").read().replace(d, "<D>").splitlines():
        a = line.split("\t", 1)
        if len(a) == 2:
            out[a[0]] = a[1]
    return out


def classify(k, in_a, in_b, cats, settings_ok):
    """The category that allows key k to differ (or be one-sided), or None."""
    both = in_a and in_b
    if "operational" in cats and both and (TIMING_RE.search(k) or k in OPERATIONAL or OPERATIONAL_RE.match(k)):
        return "operational"
    if "profile" in cats and both and k in PROFILE and settings_ok:
        return "profile"
    if "route" in cats:
        if (ROUTE_ONLY_RE.match(k) or k in ROUTE_ONLY) and in_a != in_b:
            return "route"
        if ORDINARY_ONLY_RE.search(k) and in_a != in_b:
            return "route"
        if both and k in ROUTE_WORK:
            return "route"
    if "admission" in cats:
        if k in ADMISSION_ONLY and in_a != in_b:
            return "admission"
        if both and k in ADMISSION_WORK:
            return "admission"
    if "checkpoints" in cats and CHECKPOINT_RE.match(k):
        return "checkpoints"
    if "envelope" in cats and both and k in ENVELOPE:
        return "envelope"
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("out")
    ap.add_argument("--declared-settings", default="")
    ap.add_argument("--categories", default="operational,profile")
    args = ap.parse_args()
    a, b = os.path.abspath(args.a), os.path.abspath(args.b)
    declared = {s for s in args.declared_settings.split(",") if s}
    cats = {c for c in args.categories.split(",") if c}
    fa = {f for f in os.listdir(a) if not SCIENTIFIC_SKIP.match(f)}
    fb = {f for f in os.listdir(b) if not SCIENTIFIC_SKIP.match(f)}
    bad, rows = [], []
    if fa != fb:
        bad.append(f"file sets differ: {sorted(fa ^ fb)}")
    # settings first: values (not sources) must agree on every key except the declared settings, with equal key sets
    settings_ok = False
    if "call.settings.tsv" in fa & fb:
        sa = {k: v.split("\t")[0] for k, v in kv(os.path.join(a, "call.settings.tsv"), a).items()}
        sb = {k: v.split("\t")[0] for k, v in kv(os.path.join(b, "call.settings.tsv"), b).items()}
        undeclared = sorted(k for k in set(sa) | set(sb) if sa.get(k) != sb.get(k) and k not in declared)
        settings_ok = not undeclared
        rows.append(f"settings\tundeclared_differences\t{','.join(undeclared) or 'none'}")
        rows.append(f"settings\tdeclared_differences\t{','.join(sorted(k for k in declared if sa.get(k) != sb.get(k))) or 'none'}")
        if undeclared:
            bad.append("effective settings differ outside the declared settings: " + ",".join(undeclared))
    for f in sorted(fa & fb):
        pa, pb = os.path.join(a, f), os.path.join(b, f)
        if f == "call.settings.tsv":
            continue
        if f == "call.mosaic_status.tsv":
            x, y = kv(pa, a), kv(pb, b)
            used, unexpected = {}, []
            for k in sorted(set(x) | set(y)):
                if k in x and k in y and x[k] == y[k]:
                    continue
                cat = classify(k, k in x, k in y, cats, settings_ok)
                if cat is None:
                    unexpected.append(k + ("" if k in x and k in y else (" (only in A)" if k in x else " (only in B)")))
                else:
                    used[cat] = used.get(cat, 0) + 1
            rows.append(f"status\tallowed_by_category\t{','.join(f'{c}:{n}' for c, n in sorted(used.items())) or 'none'}")
            rows.append(f"status\tunexpected\t{','.join(unexpected) or 'none'}")
            if unexpected:
                bad.append("status keys outside the declared categories: " + ",".join(unexpected))
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
    with open(args.out, "w") as o:
        o.write(f"key\tvalue\ncomparator\ts2l_compare_v3.1\ncategories\t{','.join(sorted(cats))}\n"
                f"declared_settings\t{','.join(sorted(declared)) or 'none'}\n")
        o.write(f"science_identical\t{'yes' if not bad else 'NO: ' + ' | '.join(bad)}\n")
        for r in rows:
            o.write(r + "\n")
    print(open(args.out).read(), end="")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
