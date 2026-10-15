#!/usr/bin/env python3
"""Bounded performance measurement of the existing exact caches (not a campaign; no truth is read).

For each (locus, coverage) cell, runs the production caller twice, one job at a time with 4 threads:
  default  the built-in cache budgets (emission and union caches 0; repair origin cache 384 MiB)
  c4v6     the exact-cache budgets of profiles/genotype_c4_v6_campaign_local.tsv (emission, union placement /
           origin / prefilter caches, start-sort scratch, repair origin cache 768 MiB, consumption caps)
Every other setting is identical between the two (the locus base profile). Each run is bounded by the profile's
time_limit_seconds and footprint_limit_bytes and by a process-group backstop; its complete/incomplete status and the
caller's own work counters and stage timings are recorded, never extrapolated. Then the two runs' scientific outputs
are compared (perf_compare): a performance-only setting must leave calls, ties, scores, GQ, support, filters and
FASTA byte-identical.

Usage: perf_bounded.py <out dir> <panvar> <cells tsv> [<cell id>...]
cells tsv: id, locus, coverage, graph, reference, curation (or none), reads1, reads2, library, base_profile
"""
import csv
import os
import re
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "cyp2d6_pilot"))
import pilot  # noqa: E402  (run_timed: process-group-safe timeouts)

C4V6_CACHE = {"emission_cache_bytes": 536870912, "repair_origin_cache_bytes": 805306368,
              "union_placement_cache_bytes": 536870912, "max_union_cached_placements_consumed": 2000000000,
              "union_origin_cache_bytes": 1073741824, "max_union_cached_origins_consumed": 2000000000,
              "max_union_cached_states_consumed": 2000000000, "union_prefilter_cache_bytes": 134217728,
              "max_union_cached_prefilters_consumed": 2000000000, "start_sort_scratch_bytes": 335544320}
STAGE_KEYS = ["seconds_model", "seconds_markers", "seconds_proposal", "seconds_spelling",
              "seconds_rescoring_fragment_preparation", "seconds_rescoring_placement_and_mass",
              "seconds_rescoring_diploid_combination", "seconds_union_placement_cache_fill",
              "seconds_union_origin_cache_fill", "seconds_repair", "seconds_repair_scoring", "seconds_tie_branching",
              "seconds_total"]
WORK_KEYS = ["status", "call_class", "no_call_reason", "piece_hits", "physical_origins", "statistical_state_classes",
             "union_placement_cache_hits", "union_placement_cache_misses", "union_cached_placements_consumed",
             "union_origin_cache_hits", "union_origin_cache_misses", "union_cached_origins_consumed",
             "union_cached_states_consumed", "union_cached_prefilters_consumed", "union_prefilter_cache_misses",
             "union_placement_cache_peak_bytes", "union_origin_cache_peak_bytes", "union_prefilter_cache_decision_table_bytes",
             "repair_rounds_run", "repair_cumulative_physical_origins", "repair_cumulative_state_classes",
             "pool_repair_scored", "resources_rss_peak_bytes", "resources_phys_footprint_peak_bytes"]


def kv(path):
    out = {}
    if os.path.exists(path):
        for line in open(path).read().splitlines():
            p = line.split("\t")
            if len(p) >= 2:
                out.setdefault(p[0], p[1])
    return out


def write_profile(path, base, extra):
    lines = [l.rstrip("\n") for l in open(base)]
    keys = {l.split("\t")[0] for l in lines if "\t" in l}
    clash = keys & set(extra)
    if clash:
        raise SystemExit(f"{base} already sets {sorted(clash)}")
    with open(path, "w") as f:
        f.write("\n".join(l if not l.startswith("profile\t") else l + ("-c4v6cache" if extra else "") for l in lines)
                + "\n" + "".join(f"{k}\t{v}\n" for k, v in extra.items()))


def run_cell(out, panvar, cell, mode):
    d = os.path.join(out, "runs", f"{cell['id']}.{mode}")
    os.makedirs(d, exist_ok=True)
    prefix = os.path.join(d, "call")
    if os.path.exists(prefix + ".run.tsv"):
        raise SystemExit(f"{prefix}.run.tsv exists; a finished run is never repeated into the same place")
    profile = os.path.join(d, "profile.tsv")
    write_profile(profile, cell["base_profile"], C4V6_CACHE if mode == "c4v6" else {})
    limit = pilot.profile_time_limit(profile)
    cmd = ["/usr/bin/time", "-l", panvar, "genotype", "--preset", "short-read", "--profile", profile, "--threads", "4",
           "--graph", cell["graph"], "--reference-path", cell["reference"], "--reads1", cell["reads1"],
           "--reads2", cell["reads2"], "--library", cell["library"], "--dosage-mode", "repair", "--out", prefix, "--quiet"]
    if cell["curation"] != "none":
        cmd[cmd.index("--reads1"):cmd.index("--reads1")] = ["--panel-curation", cell["curation"]]
    t0 = time.monotonic()
    with open(prefix + ".call.log", "w") as log:
        log.write("$ " + " ".join(cmd) + "\n")
        log.flush()
        code, _ = pilot.run_timed(cmd, log, limit + pilot.EXIT_ALLOWANCE_SECONDS)
    wall = time.monotonic() - t0
    text = open(prefix + ".call.log").read()
    m = re.search(r"(\d+)\s+maximum resident set size", text)
    st = kv(prefix + ".mosaic_status.tsv")
    rec = {"cell": cell["id"], "locus": cell["locus"], "coverage": cell["coverage"], "mode": mode, "exit": code,
           "complete": "yes" if code in (0, 1) else "no", "wall_seconds": f"{wall:.1f}",
           "peak_rss_bytes": m.group(1) if m else "NA", "time_limit_seconds": limit}
    rec.update({k: st.get(k, "NA") for k in STAGE_KEYS + WORK_KEYS})
    with open(prefix + ".run.tsv", "w") as f:
        f.write("\t".join(map(str, rec)) + "\n" + "\t".join(map(str, rec.values())) + "\n")
    print("\t".join(f"{k}={rec[k]}" for k in ("cell", "mode", "exit", "wall_seconds", "seconds_total", "status")), flush=True)
    return rec


def main(out, panvar, cells_tsv, chosen):
    os.makedirs(out, exist_ok=True)
    with open(cells_tsv) as f:
        cells = [c for c in csv.DictReader(f, delimiter="\t") if not chosen or c["id"] in chosen]
    for cell in cells:
        for mode in ("default", "c4v6"):
            rec = run_cell(out, panvar, cell, mode)
            path = os.path.join(out, "measurements.tsv")
            new = not os.path.exists(path)
            with open(path, "a") as f:
                if new:
                    f.write("\t".join(rec) + "\n")
                f.write("\t".join(map(str, rec.values())) + "\n")


if __name__ == "__main__":
    main(os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2]), sys.argv[3], sys.argv[4:])
