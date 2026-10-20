#!/usr/bin/env python3
"""M8-A identifiability signal for the masked-emission alternative (no fitting): for every held-out haplotype (family
excluded from training), using only its sequence-supported blocks, the minimum number of template switches a copying
path over training templates needs (a min-switch dynamic programme; unsupported blocks are skipped, never projected).
0 switches everywhere would push any fitted redraw parameter to its lower boundary (unidentifiable from below).

  masked.py <support.full dir> <out dir>
"""
import csv
import os
import sys
from collections import Counter, defaultdict

src, d = sys.argv[1], sys.argv[2]
T = lambda p: list(csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t"))
fams = T(f"{src}/families.tsv")
sup = T(f"{src}/support.tsv")
rows, summ = [], []
for locus in ("C4", "CYP2D6", "ACOT"):
    tab = [l.rstrip("\n").split("\t") for l in open(f"{src}/{locus}.block_sequences.tsv")][1:]
    blocks = list(dict.fromkeys(r[0] for r in tab))
    seq = {}
    for b in blocks:
        name = None
        for line in open(f"{src}/{locus}.{b}.alleles.fa"):
            if line.startswith(">"):
                name = line[1:].split()[0]
                seq[(b, name)] = ""
            else:
                seq[(b, name)] += line.strip().upper()
    path = defaultdict(dict)
    for r in tab:
        if r[1] == "panel":
            for t in [x for x in r[5].split(",") if x and x != "-"]:
                path[t][b] = None
                path[t][r[0]] = r[3]
    fam_of = {}
    for f in fams:
        if f["locus"] == locus:
            for dn in f["donors"].split(","):
                fam_of[dn] = (f["family"], f["status"])
    for s in sup:
        if s["locus"] != locus or s["assessable"] != "yes":
            continue
        t = s["template"]
        fam_donors = {dn for dn, (fa, _) in fam_of.items() if fa == s["family"]}
        train = [x for x in path if x.split("#")[0] not in fam_donors]
        unsupported = set(filter(None, s["unsupported_blocks"].split(",")))
        obs = [b for b in blocks if b not in unsupported]
        compat = [[seq[(b, path[x][b])] == seq[(b, path[t][b])] for x in train] for b in obs]
        INF = 10 ** 9
        cost = [0 if c else INF for c in compat[0]] if obs else []
        for k in range(1, len(obs)):
            best = min(cost)
            cost = [min(cost[j], best + 1) if compat[k][j] else INF for j in range(len(train))]
        ms = min(cost) if obs else 0
        rows.append(dict(locus=locus, family=s["family"], status=s["status"], template=t, supported_blocks=len(obs), unsupported_blocks=len(unsupported),
                         min_switches_on_supported_blocks=ms if ms < INF else "infeasible"))
    for st in ("UNINSPECTED", "INSPECTED_DEVELOPMENT"):
        rs = [r for r in rows if r["locus"] == locus and r["status"] == st]
        if rs:
            c = Counter(r["min_switches_on_supported_blocks"] for r in rs)
            summ.append(dict(locus=locus, status=st, haplotypes=len(rs), supported_blocks_median=sorted(r["supported_blocks"] for r in rs)[len(rs) // 2],
                             min_switches_distribution=", ".join(f"{k}:{c[k]}" for k in sorted(c, key=str)),
                             needing_at_least_one_switch=sum(1 for r in rs if r["min_switches_on_supported_blocks"] not in (0, "infeasible"))))
os.makedirs(d, exist_ok=True)
for name, rr in (("masked", rows), ("masked_summary", summ)):
    keys = list(rr[0])
    with open(f"{d}/{name}.tsv", "w") as f:
        w = csv.DictWriter(f, fieldnames=keys, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rr)
print(open(f"{d}/masked_summary.tsv").read())
