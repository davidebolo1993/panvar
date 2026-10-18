#!/usr/bin/env python3
"""S2-J: complete-run comparison, existing round-1 origin admission vs reuse_only_experimental, explicit allowances only:
  * every output file byte-identical after path normalisation, except:
    - call.mosaic_status.tsv: keys may differ only if listed below -- the policy's own rows (reuse-only run only), the
      refinement's executed work and cap margins (windows not retained in round 1 are enumerated again), timing, process
      resources and scheduling-variable emission/DP counters; every other key must be identical and present in both;
    - call.settings.tsv: identical except the policy row;
    - call.repair_rounds.tsv: the scientific columns (0-27) identical; the rest (work, cache, timing) are reported;
    - call.dosage_status.tsv: identical except seconds_* rows;
    - call.repair_memory.tsv: a memory/timing log (not compared).
  s2j_compare.py <existing dir> <reuse dir> <out.tsv>
"""
import os
import re
import sys

POLICY_ROWS = {"repair_round1_origin_admission", "repair_admission_restricted_plans", "repair_admission_planned_windows",
               "repair_admission_singleton_windows", "repair_admission_eligible_windows",
               "repair_admission_misses_while_restricted", "repair_admission_singleton_misses_not_retained",
               "repair_admission_retained_while_restricted", "repair_admission_hits_while_restricted",
               "repair_admission_group_reuses_while_restricted"}
WORK = {"repair_cumulative_latent_pairs", "repair_cumulative_physical_origins", "repair_cumulative_piece_hits",
        "repair_margin_latent_pairs", "repair_margin_physical_origins", "repair_margin_piece_hits"}
RESOURCES = {"resources_phys_footprint_peak_bytes", "resources_rss_peak_bytes", "resources_footprint_at_repair_start_bytes",
             "union_end_phys_footprint_bytes", "threads_peak_active_workers", "profile_sha256",
             "checkpoint_seconds_first_evaluated_prediction", "checkpoint_seconds_committed",
             "early_incumbent_seconds_evaluation", "early_incumbent_seconds_preparation", "early_incumbent_seconds_selection",
             "early_incumbent_footprint_after_bytes", "union_fragment_block_footprint_max_bytes", "scorer_mass_matrix_bytes",
             "scorer_haplotype_sequence_bytes"}
SCHEDULING = {"emission_cache_dp_hits", "emission_cache_dp_misses", "emission_cache_dp_requests", "emission_cache_entries",
              "emission_cache_executed_dps", "emission_cache_executed_prefilters", "emission_cache_payload_bytes",
              "dp_evaluations", "full_read_verifications"}
TIMING_RE = re.compile(r"^(seconds_|worker_seconds_)|^elapsed_seconds$")
RESOURCE_RE = re.compile(r"^footprint_checkpoint[0-9]+_(phys_footprint_bytes|phys_footprint_peak_bytes|seconds)$|"
                         r"^union_block[0-9]+_(phys_footprint_after_bytes|executed_dps)$")


def kv(p, d):
    out = {}
    for line in open(p, errors="replace").read().replace(d, "<D>").splitlines():
        a = line.split("\t", 1)
        if len(a) == 2:
            out[a[0]] = a[1]
    return out


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
            continue
        if f == "call.mosaic_status.tsv":
            x, y = kv(pa, a), kv(pb, b)
            unexpected, allowed = [], []
            for k in sorted(set(x) | set(y)):
                if k in x and k in y and x[k] == y[k]:
                    continue
                ok = (k not in x and k in POLICY_ROWS) or (k in x and k in y and (
                    k in WORK or k in RESOURCES or k in SCHEDULING or TIMING_RE.search(k) or RESOURCE_RE.match(k)))
                (allowed if ok else unexpected).append(k)
            rows.append(f"{f}\tallowed_differing_keys\t{','.join(allowed) or 'none'}")
            rows.append(f"{f}\tunexpected_differing_keys\t{','.join(unexpected) or 'none'}")
            for k in ("repair_cumulative_physical_origins", "repair_cumulative_latent_pairs", "seconds_repair", "seconds_total",
                      "resources_phys_footprint_peak_bytes"):
                rows.append(f"{f}\t{k}\t{x.get(k, 'NA')} -> {y.get(k, 'NA')}")
            if unexpected:
                bad.append("status keys outside the allowances: " + ",".join(unexpected))
            continue
        if f == "call.settings.tsv":
            x, y = kv(pa, a), kv(pb, b)
            y.pop("repair_round1_origin_admission", None)
            if x != y:
                bad.append("settings differ beyond the policy row")
            continue
        if f == "call.repair_rounds.tsv":
            cut = lambda p: [line.split("\t")[:28] for line in open(p).read().splitlines()]  # noqa: E731
            if cut(pa) != cut(pb):
                bad.append("trajectory differs")
            rows.append(f"{f}\trounds\t{len(cut(pa)) - 1} -> {len(cut(pb)) - 1}")
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
