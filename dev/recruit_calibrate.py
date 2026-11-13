#!/usr/bin/env python3
"""Calibrates the certain cut on a benchmark (dev/recruit_benchmark.sh): for each technology and each cut
(private sites, density floor, decoy share), recall on the locus sets and false certain on the negative
sets (near, decoy, genome, new, wgnew; a genome set's fragments count once, its false certain over all
loci), from the pooled evidence recruit writes; then the cut with the
best recall whose 95% upper bound on false certain per million fragments meets the target.

    recruit_calibrate.py <benchmark dir> [--target 1.0] [--arm main] [--out calibration.tsv]

The pooled evidence leaves out mate rescue, so recall here is a little below recruit's own; recruit's own
false certain on the negatives (rescue included) is printed beside the default cut.
"""
import argparse
import csv
import glob
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from recruit_score import poisson_upper  # noqa: E402

ERROR = {"short": 0.002, "illumina": 0.002, "ont": 0.03, "hifi": 0.001}
SITES = [1, 2, 3, 4]
DENSITY = [0.0, 0.02, 0.05, 0.1]
SHARE = [0.0, 0.05, 0.1]
NEGATIVE = ("near", "decoy", "genome", "new", "wgnew")


def evidence(table, locus):
    with open(table) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["locus"] != locus:
                continue
            yield tuple(int(r[c]) for c in ("sites_p", "sites_d", "runs_n", "seg_kmers", "kmers", "windows")), r["fate"]


def certain(v, e, m, dens, share):
    sp, sd, runs, seg, kmers, windows = v
    if windows > 0 and kmers * 2 < windows:
        return False
    lam = (e + 0.001) * seg
    if runs > lam + 3.1 * math.sqrt(lam) + 2:
        return False
    return sp >= m and sd <= share * sp and sp >= dens * seg / 31


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bench")
    ap.add_argument("--target", type=float, default=1.0, help="false certain per million fragments (95%% upper bound)")
    ap.add_argument("--arm", default="main")
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    summary = list(csv.DictReader(open(os.path.join(a.bench, "summary.tsv")), delimiter="\t"))
    cuts = [(m, d, s) for m in SITES for d in DENSITY for s in SHARE]
    rows = []
    own = {}       # recruit's own false certain (mate rescue included) on the negatives, per technology
    own_n = {}
    for tech in ERROR:
        e = ERROR[tech]
        pos_n = neg_n = 0
        pos = [0] * len(cuts)
        neg = [0] * len(cuts)
        counted = set()   # a genome set is one table for all loci: its fragments count once
        for r in summary:
            if r["tech"] != tech or r["arm"] != a.arm:
                continue
            src, locus = r["source"], r["locus"]
            if src == "locus":
                table = os.path.join(a.bench, f"{a.arm}.{locus}.{tech}.locus.tsv")
            elif src in NEGATIVE:
                table = os.path.join(a.bench, f"{a.arm}.{src}.{tech}.tsv" if src in ("genome", "new", "wgnew")
                                     else f"{a.arm}.{locus}.{tech}.{src}.tsv")
            else:
                continue
            if not os.path.exists(table):
                continue
            if src == "locus":
                counts = pos
                pos_n += int(r["fragments"])
            else:
                counts = neg
                if table not in counted:
                    counted.add(table)
                    neg_n += int(r["fragments"])
            for v, fate in evidence(table, locus):
                if counts is neg and fate == "certain":
                    own[tech] = own.get(tech, 0) + 1
                for i, (m, d, s) in enumerate(cuts):
                    if certain(v, e, m, d, s):
                        counts[i] += 1
        if pos_n == 0:
            continue
        own_n[tech] = neg_n
        for i, (m, d, s) in enumerate(cuts):
            up = 1e6 * poisson_upper(neg[i]) / neg_n if neg_n else float("inf")
            rows.append({"tech": tech, "min_sites": m, "min_density": d, "decoy_share": s, "recall": pos[i] / pos_n,
                         "false_certain": neg[i], "negatives": neg_n, "per_million": 1e6 * neg[i] / neg_n if neg_n else 0,
                         "upper95_per_million": up})
    if a.out:
        with open(a.out, "w") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(rows)
    print(f"target: 95% upper bound <= {a.target} false certain per million fragments ({a.arm} arm)")
    print("tech\tcut (sites, density, share)\trecall\tfalse certain / negatives\tupper95 per million\tmeets target")
    for tech in ERROR:
        t = [r for r in rows if r["tech"] == tech]
        if not t:
            continue
        ok = [r for r in t if r["upper95_per_million"] <= a.target]
        best = max(ok or t, key=lambda r: (r["recall"], -r["upper95_per_million"]))
        default = next(r for r in t if (r["min_sites"], r["min_density"], r["decoy_share"]) == (2, 0.05, 0.05))
        for label, r in (("best", best), ("default", default)):
            print(f"{tech}\t{label} ({r['min_sites']}, {r['min_density']}, {r['decoy_share']})\t{r['recall']:.4f}\t"
                  f"{r['false_certain']} / {r['negatives']}\t{r['upper95_per_million']:.2f}\t{'yes' if r['upper95_per_million'] <= a.target else 'no'}")
        k, n = own.get(tech, 0), own_n[tech]
        up = 1e6 * poisson_upper(k) / n
        print(f"{tech}\trecruit as run (rescue included)\t\t{k} / {n}\t{up:.2f}\t{'yes' if up <= a.target else 'no'}")


if __name__ == "__main__":
    main()
