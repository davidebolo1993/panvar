#!/usr/bin/env python3
"""M3-D matched-run comparator (TRUTH-FREE; frozen before execution).

1. Option-off identity: for each cell, every output file of the new binary's off run must equal the baseline binary's,
   after normalising the run prefix and dropping volatile lines (timing, threads, resources, memory/footprint/rss, binary
   identity) and the volatile columns of .repair_rounds.tsv (columns from `seconds` on). The set of files must match.
   .repair_memory.tsv and .time are resource logs and are compared by presence only.
2. On vs off (reported): status, call class, called sequence digests, block genotypes, the final winner's score
   components, moves, ties (winner set size), joint-search rows, pool sizes, and the whole-caller wall time and peak RSS
   (ledger) -- the incremental overhead.
Writes <prefix>.identity.tsv, <prefix>.matched.tsv and <prefix>.selection.sha256 (the frozen selection).

  compare_m3d.py <runs dir> <ledger.tsv> <out prefix> cell...
"""
import csv
import hashlib
import os
import re
import sys

runs, ledger_path, prefix, cells = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
VOLATILE = re.compile(r"^(seconds_|threads_|resources_|time_|wall_|peak_|binary|panvar_|build_|command|hostname|pid|date|started|"
                      r"finished|write_|running_|.*footprint|.*rss|.*_seconds|.*memory|.*_bytes_end|allocator|checkpoint)", re.I)
# Amendment 1 (only when M3D_CACHE_COUNTERS_VOLATILE=1): the emission/prefilter cache work counters, shown to vary between two
# runs of the SAME baseline binary on the same inputs (scheduling-dependent cache races at 4 threads)
CACHE = re.compile(r"^(emission_cache_dp_hits|emission_cache_dp_misses|emission_cache_executed_prefilters|emission_cache_executed_dps|"
                   r"scorer_emission_cache_arena_bytes|union_cached_prefilters_consumed|union_prefilter_cache_misses)$")
CACHE_VOLATILE = os.environ.get("M3D_CACHE_COUNTERS_VOLATILE") == "1"
RESOURCE_FILES = (".repair_memory.tsv", ".time", ".log")
ledger = {r["unit"]: r for r in csv.DictReader(open(ledger_path), delimiter="\t")}


def files(d, name):
    return sorted(f[len(name):] for f in os.listdir(d) if f.startswith(name + ".") and not f.endswith((".partial.tsv",)))


def normalised(path, run_prefix, suffix):
    text = open(path, errors="replace").read().replace(run_prefix, "<prefix>")
    lines = text.split("\n")
    if suffix == ".repair_rounds.tsv" and lines:
        head = lines[0].split("\t")
        cut = head.index("seconds") if "seconds" in head else len(head)
        return "\n".join("\t".join(l.split("\t")[:cut]) for l in lines)
    keep = []
    for l in lines:
        key = l.split("\t", 1)[0]
        if VOLATILE.match(key) or (CACHE_VOLATILE and CACHE.match(key)):
            continue
        keep.append(l)
    return "\n".join(keep)


def status(d):
    out = {}
    p = f"{d}/call.mosaic_status.tsv"
    if os.path.exists(p):
        for l in open(p):
            k, _, v = l.rstrip("\n").partition("\t")
            out[k] = v
    return out


identity, matched = [], []
for cell in cells:
    base, off, on = (f"{runs}/run.{cell}.{arm}" for arm in ("baseline", "off", "on"))
    fb, fo = files(base, "call"), files(off, "call")
    identity.append((cell, "file_set_equal", fb == fo, f"{len(fb)} | {len(fo)}"))
    for suffix in sorted(set(fb) & set(fo)):
        if suffix.endswith(RESOURCE_FILES):
            continue
        a = normalised(f"{base}/call{suffix}", f"{base}/call", suffix)
        b = normalised(f"{off}/call{suffix}", f"{off}/call", suffix)
        diff = "" if a == b else next((f"{x!r} | {y!r}" for x, y in zip(a.split("\n"), b.split("\n")) if x != y), "length differs")
        identity.append((cell, f"identical {suffix}", a == b, diff[:300]))
    so, sn = status(off), status(on)
    def digest(d):
        p = f"{d}/call.mosaic_haplotypes.fa"
        if not os.path.exists(p):
            return "none"
        seqs = sorted("".join(l.strip() for l in chunk.split("\n")[1:]) for chunk in open(p).read().split(">")[1:])
        return hashlib.sha256("|".join(seqs).encode()).hexdigest()[:16]
    def blocks(d):
        p = f"{d}/call.block_genotypes.tsv"
        return hashlib.sha256(open(p, "rb").read()).hexdigest()[:16] if os.path.exists(p) else "none"
    for key in ["status", "call_class", "no_call_reason", "phase_support", "repair_rounds_run", "repair_moves", "repair_stop_reason",
                "repair_converged", "repair_final_incumbent_score", "repair_final_incumbent_hap1", "repair_final_incumbent_hap2",
                "winner_set_size", "winner_candidate_set_probability", "numeric_separation_gap", "fragment_term_gap_to_runner_up",
                "pool_identities", "pool_repair_scored", "repair_cumulative_physical_origins", "repair_cumulative_state_classes",
                "seconds_total", "joint_search_status", "joint_search_decision", "joint_search_move", "joint_search_single_alternatives",
                "joint_search_catalogue_raw", "joint_search_catalogue_unique", "joint_search_winner_set_size", "joint_search_certified_gap",
                "joint_search_origins_reused", "joint_search_origins_recomputed", "joint_search_full_placements",
                "joint_search_census_bytes_peak", "joint_search_transient_bytes_peak", "joint_search_transient_bytes_after",
                "joint_search_live_transient_after", "joint_search_physical_origins", "joint_search_state_classes", "joint_search_footprint_peak",
                "pool_joint_scored", "joint_continued", "joint_leg_rounds", "joint_leg_moves", "joint_leg_stop", "joint_leg_converged",
                "pool_joint_leg_scored", "seconds_joint_search", "seconds_joint_leg"]:
        matched.append((cell, key, so.get(key, "-"), sn.get(key, "-")))
    matched.append((cell, "called_sequence_digest", digest(off), digest(on)))
    matched.append((cell, "block_genotypes_digest", blocks(off), blocks(on)))
    for arm in ("baseline", "off", "on"):
        u = ledger.get(f"run.{cell}.{arm}", {})
        matched.append((cell, f"wall_seconds_{arm}", u.get("wall_seconds", "-"), u.get("peak_rss_gib", "-")))
with open(f"{prefix}.identity.tsv", "w") as f:
    f.write("cell\tcheck\tpass\tdetail\n")
    for r in identity:
        f.write("\t".join(str(x) for x in r) + "\n")
with open(f"{prefix}.matched.tsv", "w") as f:
    f.write("cell\tkey\toff\ton\n")
    for r in matched:
        f.write("\t".join(str(x) for x in r) + "\n")
digest = hashlib.sha256(open(f"{prefix}.matched.tsv", "rb").read()).hexdigest()
open(f"{prefix}.selection.sha256", "w").write(digest + "\n")
failed = [r for r in identity if not r[2]]
print(f"identity: {len(identity)} checks, {len(failed)} failed")
for r in failed:
    print("FAILED", r[0], r[1], r[3])
print("selection digest", digest)
sys.exit(1 if failed else 0)
