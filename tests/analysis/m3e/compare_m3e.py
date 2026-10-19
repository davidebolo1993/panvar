#!/usr/bin/env python3
"""M3-E integration comparator (frozen before execution).
1. Annotation-only: each run (block_support_rule experimental; joint off/on) against the matching M3-D run: identical
   output-file sets and normalised contents, with these allowances only -- the M3-D volatile set (timing, threads,
   resources, memory/footprint/rss, binary identity, and the 7 emission/prefilter cache work counters of M3-D Amendment 1);
   .block_genotypes.tsv compared on its first 29 (pre-existing) columns; the status rows block_support_* and
   phase_competitors_* and the settings row block_support_rule added. Rounds compared on their science columns.
2. Cross-check: each block's genotype_support_v1x against the independent saved-evidence audit (same label; margin within
   1e-3), except blocks the audit flags as near ties.
  compare_m3e.py <runs> <m3d runs> <audit blocks.tsv> <out prefix> cell...
"""
import csv
import os
import re
import sys

runs, m3d, audit_path, prefix, cells = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5:]
VOLATILE = re.compile(r"^(seconds_|threads_|resources_|time_|wall_|peak_|binary|panvar_|build_|command|hostname|pid|date|started|"
                      r"finished|write_|running_|.*footprint|.*rss|.*_seconds|.*memory|.*_bytes_end|allocator|checkpoint|"
                      r"emission_cache_dp_hits|emission_cache_dp_misses|emission_cache_executed_prefilters|emission_cache_executed_dps|"
                      r"scorer_emission_cache_arena_bytes|union_cached_prefilters_consumed|union_prefilter_cache_misses)", re.I)
ADDED = re.compile(r"^(block_support_|phase_competitors_)")


def norm(path, run_prefix, suffix):
    lines = open(path, errors="replace").read().replace(run_prefix, "<prefix>").split("\n")
    if suffix == ".repair_rounds.tsv" and lines:
        head = lines[0].split("\t")
        cut = head.index("seconds") if "seconds" in head else len(head)
        return [("\t".join(l.split("\t")[:cut])) for l in lines]
    if suffix == ".block_genotypes.tsv":
        return ["\t".join(l.split("\t")[:29]) for l in lines]
    out = []
    for l in lines:
        key = l.split("\t", 1)[0]
        if VOLATILE.match(key) or ADDED.match(key) or (suffix == ".settings.tsv" and key == "block_support_rule"):
            continue
        out.append(l)
    return out


checks = []
for cell in cells:
    for arm in ("off", "on"):
        new, old = f"{runs}/run.{cell}.{arm}", f"{m3d}/run.{cell}.{arm}"
        fn = sorted(f[4:] for f in os.listdir(new) if f.startswith("call."))
        fo = sorted(f[4:] for f in os.listdir(old) if f.startswith("call."))
        checks.append((cell, arm, "file_set_equal", fn == fo, f"{len(fn)} | {len(fo)}"))
        for suffix in sorted(set(fn) & set(fo)):
            if suffix.endswith((".repair_memory.tsv", ".time", ".log")):
                continue
            a, b = norm(f"{new}/call{suffix}", f"{new}/call", suffix), norm(f"{old}/call{suffix}", f"{old}/call", suffix)
            diff = "" if a == b else next((f"{x!r} | {y!r}" for x, y in zip(a, b) if x != y), "length differs")
            checks.append((cell, arm, f"unchanged {suffix}", a == b, diff[:240]))
audit = {}
for r in csv.DictReader(open(audit_path), delimiter="\t"):
    audit[(r["run"], r["block"])] = r
for cell in cells:
    for arm in ("off", "on"):
        for r in csv.DictReader(open(f"{runs}/run.{cell}.{arm}/call.block_genotypes.tsv"), delimiter="\t"):
            a = audit[(f"{cell}.{arm}", r["block_id"])]
            if a["near_tie"] == "True":
                checks.append((cell, arm, f"crosscheck {r['block_id']}", True, "near tie in the rounded audit: not compared"))
                continue
            same = r["genotype_support_v1x"] == a["genotype_support_v1x"]
            if same and a["margin_L"] not in ("-", "") and r["data_support_margin_v1x"] != "NA":
                same = abs(float(r["data_support_margin_v1x"]) - float(a["margin_L"])) <= 1e-3
            checks.append((cell, arm, f"crosscheck {r['block_id']}", same,
                           f"{r['genotype_support_v1x']} {r['data_support_margin_v1x']} | {a['genotype_support_v1x']} {a['margin_L']}"))
with open(prefix + ".checks.tsv", "w") as f:
    f.write("cell\tarm\tcheck\tpass\tdetail\n")
    for c in checks:
        f.write("\t".join(str(x) for x in c) + "\n")
failed = [c for c in checks if not c[3]]
print(f"{len(checks)} checks, {len(failed)} failed")
for c in failed:
    print("FAILED", *c)
sys.exit(1 if failed else 0)
