#!/usr/bin/env python3
"""M6-E Amendment 2 (panel only): measurements for the state-count estimate and the structural layer.

  1. length-class frequencies by physical copy, haplotype and donor (periods, tails)
  2. haplotype arrangements: copy count and length-class pattern per template, by haplotype and donor
  3. local state counts: sliding windows (400 bp and 800 bp, step 200 bp) over each type's reference coordinate; per
     window the SIMPLE SNV sites and the number of distinct local allele patterns among unique components
  4. non-simple intervals (from analyse.a1): the number of distinct local event signatures among unique components
     (an upper bound on distinct local sequences: co-optimal alignments can split one sequence into several signatures)

  extra.py <prepare dir> <events tsv> <analyse prefix> <out prefix>
"""
import builtins
import csv
import os
import sys
from collections import Counter, defaultdict

import numpy as np

prep, events_path, ana, prefix = sys.argv[1:5]
ALLOW = {os.path.abspath(p) for p in [events_path, f"{prep}/components.tsv", f"{prep}/copies.tsv", f"{prep}/prepare.tsv",
                                       f"{ana}.sites.tsv", f"{ana}.intervals.tsv"]}
_open = builtins.open
OUTDIR = os.path.dirname(os.path.abspath(prefix))


def guarded(path, mode="r", *a, **k):
    p = os.path.abspath(path)
    if "w" in mode or "a" in mode:
        assert p.startswith(OUTDIR), f"write outside the unit directory: {p}"
    else:
        assert p in ALLOW or p.startswith(OUTDIR), f"not an allowlisted panel-only input: {p}"
    return _open(path, mode, *a, **k)


builtins.open = guarded
T = lambda p: list(csv.DictReader(open(p), delimiter="\t"))
comps = {r["unique_record"]: r for r in T(f"{prep}/components.tsv")}
copies = T(f"{prep}/copies.tsv")
refs = {l.split("\t")[1]: l.split("\t")[2].split()[0] for l in open(f"{prep}/prepare.tsv") if l.startswith("reference\t")}
ref_bp = {k: int(l.split("(")[1].split()[0]) for l in open(f"{prep}/prepare.tsv") if l.startswith("reference\t") for k in [l.split("\t")[1]]}
ev = defaultdict(list)
for r in T(events_path):
    if r["event"] != "SUMMARY":
        ev[r["label"]].append(r)
sites = T(f"{ana}.sites.tsv")
intervals = T(f"{ana}.intervals.tsv")
S = []
say = lambda k, v: S.append((k, v))

# 1. length classes by copy / haplotype / donor
for kind in ("period", "tail"):
    cp = [c for c in copies if c["type"] == kind]
    for lc in ("long", "short"):
        sel = [c for c in cp if comps[c["unique_record"]]["length_class"] == lc]
        say(f"{kind}.length_class_{lc}", f"copies {len(sel)}/{len(cp)}; haplotypes {len({c['template'] for c in sel})}; donors {len({c['donor'] for c in sel})}")

# 2. haplotype arrangements (copy order: period1, period2, tail)
order = {"period1": 0, "period2": 1, "tail": 2}
hap = defaultdict(list)
for c in copies:
    if c["type"] in ("period", "tail"):
        hap[c["template"]].append((order[c["position"]], comps[c["unique_record"]]["length_class"][0].upper()))
arr = {t: "".join(x[1] for x in sorted(v)) for t, v in hap.items()}
by_h = Counter(arr.values())
by_d = Counter()
for pat in by_h:
    by_d[pat] = len({t.split("#")[0] for t, a in arr.items() if a == pat})
say("arrangement.copies_per_haplotype", ", ".join(f"{n}:{c}" for n, c in sorted(Counter(len(a) for a in arr.values()).items())))
say("arrangement.patterns_by_haplotype(donors)", ", ".join(f"{p}:{n}({by_d[p]})" for p, n in by_h.most_common()))
dip = Counter()
for d in {t.split("#")[0] for t in arr}:
    ts = [t for t in arr if t.split("#")[0] == d]
    if len(ts) == 2:
        dip[sum(len(arr[t]) for t in ts)] += 1
say("arrangement.diploid_total_copies_by_donor", ", ".join(f"{n}:{c}" for n, c in sorted(dip.items())))


# per-component alleles at SIMPLE SNV positions
def alleles(kind, label_prefix):
    simple = sorted(int(s["ref_pos"]) for s in sites if s["type"] == label_prefix and s["cls"] == "SIMPLE")
    names = sorted((k for k in comps if comps[k]["type"] == kind), key=lambda k: int(k.split("_u")[1]))
    snv = {k: {int(e["ref_pos"]): e["query_bases"] for e in ev[f"{label_prefix}:{k}"] if e["event"] == "SNV"} for k in names}
    return simple, names, snv


# 3. windows
for kind in ("period", "tail"):
    simple, names, snv = alleles(kind, kind)
    L = ref_bp[kind]
    for wbp in (400, 800):
        nsite, npat = [], []
        for w0 in range(0, max(1, L - wbp + 1), 200):
            pos = [p for p in simple if w0 <= p < w0 + wbp]
            nsite.append(len(pos))
            npat.append(len({tuple(snv[k].get(p, "") for p in pos) for k in names}))
        nsite, npat = np.array(nsite), np.array(npat)
        say(f"{kind}.windows_{wbp}bp", f"{len(npat)} windows; SIMPLE sites/window median {np.median(nsite):.0f} max {nsite.max()}; "
            f"distinct local patterns among {len(names)} unique components: median {np.median(npat):.0f}, 90th pct {np.percentile(npat, 90):.0f}, max {npat.max()}")

# 4. non-simple intervals: distinct local event signatures
for kind in ("period", "tail"):
    names = [k for k in comps if comps[k]["type"] == kind]
    sig_counts = []
    for iv in intervals:
        if iv["type"] != kind or iv["class"] == "SIMPLE":
            continue
        b, e = int(iv["begin"]), int(iv["end"])
        sigs = set()
        for k in names:
            sig = []
            for x in ev[f"{kind}:{k}"]:
                p, n = int(x["ref_pos"]), max(1, int(x["ref_len"]))
                if p < e and p + n > b:
                    sig.append((x["event"], p, x["ref_len"], x["query_len"], x["query_bases"]))
            sigs.add(tuple(sorted(sig)))
        sig_counts.append((e - b, len(sigs), iv["class"], b, e))
    big = sorted(sig_counts, key=lambda s: -s[0])[:5]
    say(f"{kind}.nonsimple_intervals_distinct_signatures", f"{len(sig_counts)} intervals; signatures median {np.median([s[1] for s in sig_counts]):.0f}, max {max(s[1] for s in sig_counts)}; "
        "largest " + ", ".join(f"{s[3]}-{s[4]}({s[0]} bp,{s[2]}):{s[1]}" for s in big))

with open(f"{prefix}.summary.tsv", "w") as f:
    for k, v in S:
        f.write(f"{k}\t{v}\n")
print(open(f"{prefix}.summary.tsv").read())
