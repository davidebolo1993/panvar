#!/usr/bin/env python3
"""M4-B constant-mode identity: every output file of the new binary's constant run equals the baseline binary's (a728eb71),
after removing the run prefix, the volatile status/settings rows (timing, threads, resources, memory/footprint/rss, binary
identity) and M3-D Amendment 1's scheduling-dependent cache counters, and the timing columns of .repair_rounds.tsv. Resource
logs (.repair_memory.tsv, .time) are compared by presence. Exit 1 on any difference.

  identity.py <baseline prefix> <new prefix> <out tsv>
"""
import os
import re
import sys

a, b, out = sys.argv[1:4]
VOL = re.compile(r"^(seconds_|threads_|resources_|time_|wall_|peak_|binary|panvar_|build_|command|hostname|pid|date|started|finished|write_|"
                 r"running_|.*footprint|.*rss|.*_seconds|.*memory|.*_bytes_end|allocator|checkpoint|emission_cache_dp_hits|emission_cache_dp_misses|"
                 r"emission_cache_executed_prefilters|emission_cache_executed_dps|scorer_emission_cache_arena_bytes|union_cached_prefilters_consumed|"
                 r"union_prefilter_cache_misses|threads\t)", re.I)


def files(prefix):
    d, n = os.path.dirname(prefix), os.path.basename(prefix)
    return sorted(f[len(n):] for f in os.listdir(d) if f.startswith(n + ".") and "partial" not in f and not f.endswith(".sorted.gfa"))


def norm(path, prefix, suffix):
    lines = open(path, errors="replace").read().replace(prefix, "<prefix>").split("\n")
    if suffix == ".repair_rounds.tsv" and lines:
        head = lines[0].split("\t")
        cut = head.index("seconds") if "seconds" in head else len(head)
        return "\n".join("\t".join(l.split("\t")[:cut]) for l in lines)
    return "\n".join(l for l in lines if not VOL.match(l))


fa, fb = files(a), files(b)
rows = [("file_set_equal", fa == fb, f"{len(fa)} | {len(fb)}")]
for s in sorted(set(fa) & set(fb)):
    if s.endswith((".repair_memory.tsv", ".time", ".log", ".tie_branch_memory.tsv")):
        continue
    x, y = norm(a + s, a, s), norm(b + s, b, s)
    first = "" if x == y else next((f"{p!r} | {q!r}" for p, q in zip(x.split("\n"), y.split("\n")) if p != q), "length differs")
    rows.append((f"identical {s}", x == y, first[:300]))
open(out, "w").write("check\tpass\tdetail\n" + "".join(f"{c}\t{p}\t{d}\n" for c, p, d in rows))
bad = [r for r in rows if not r[1]]
print(f"identity: {len(rows)} checks, {len(bad)} failed")
for r in bad:
    print("FAILED", r[0], r[2])
sys.exit(1 if bad else 0)
