#!/usr/bin/env python3
"""Scores recruitment tables against the truth of simulated read sets, for several decision rules applied
to the same evidence (each table row keeps it), so one pass over the reads serves every rule.

    recruit_score.py --out <summary.tsv> <arm> <locus> <tech> <source> <fragments> <table> [...]
Each set is six arguments: the arm (main: GRCh38 + CHM13 background; g38: GRCh38 only), the locus, the
technology (short, ont, hifi), the source (locus, locus1 = one mate only, near, decoy, boundary,
genome, new), the fragments simulated, and recruit's table.

Rules (a read is never certain when uncallable or divergent, as recruit decides):
    runs1, runs2, runs3  >= 1, 2, 3 independent private sites, no decoy site
    share5, share10      >= 2 private sites, decoy sites at most 5% / 10% of them
    dense5               share5, and at least 0.05 private sites per 31 callable bases (recruit's default)
    kmers                >= 3 private k-mers, no decoy k-mer
    fraction             >= half of the k-mers in the locus, one private
Divergence uses each technology's error rate (short 0.002, ont 0.03, hifi 0.001).

Columns: certain, ambiguous and unrecruited counts; recall (locus sources); recall on the hard subset
(reads whose k-mers are mostly shared or repeat); false certain per million fragments and per million
screened fragments, with a 95% upper bound on the former; for the 1x genome sets, the false-certain reads
expected in 30x of the simulated chromosomes.
"""
import argparse
import csv
import math
import sys

RULES = ["runs1", "runs2", "runs3", "share5", "share10", "dense5", "kmers", "fraction"]
ERROR = {"short": 0.002, "ont": 0.03, "hifi": 0.001}


def poisson_upper(k, conf=0.95):
    """Exact upper confidence bound on a Poisson mean after observing k events."""
    lo, hi = float(k), float(k) + 10 + 10 * math.sqrt(k + 1)
    target = 1 - conf
    def cdf(lam):
        term = s = math.exp(-lam)
        for i in range(1, k + 1):
            term *= lam / i
            s += term
        return s
    for _ in range(100):
        mid = (lo + hi) / 2
        if cdf(mid) > target:
            lo = mid
        else:
            hi = mid
    return hi


def flags(r, e):
    k, w, runs = int(r["kmers"]), int(r["windows"]), int(r["runs_n"])
    uncallable = w > 0 and k * 2 < w
    lam = (e + 0.001) * k
    divergent = runs > lam + 3.1 * math.sqrt(lam) + 2
    return uncallable or divergent


def certain(rule, r, e):
    if flags(r, e):
        return False
    sp, sd, n_p, n_d = int(r["sites_p"]), int(r["sites_d"]), int(r["n_p"]), int(r["n_d"])
    if rule.startswith("runs"):
        return sd == 0 and sp >= int(rule[4:])
    if rule.startswith("share"):
        return sp >= 2 and sd <= int(rule[5:]) / 100 * sp
    if rule == "dense5":   # recruit's default: share5, and private sites at least 0.05 per 31 callable bases
        return sp >= 2 and sd <= 0.05 * sp and sp >= 0.05 * int(r["kmers"]) / 31
    if rule == "kmers":
        return n_d == 0 and n_p >= 3
    k = int(r["kmers"])
    return k > 0 and n_p >= 1 and (n_p + int(r["n_s"]) + int(r["n_r"])) / k >= 0.5


def unrecruited(r):
    sp, sd, n_s = int(r["sites_p"]), int(r["sites_d"]), int(r["n_s"])
    return (sp == 0 and sd >= 2) or (sp == 0 and sd < 2 and n_s == 0)


def hard(r):
    k = int(r["kmers"])
    return k > 0 and (int(r["n_s"]) + int(r["n_r"])) / k > 0.5


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("sets", nargs="+")
    a = ap.parse_args()
    if len(a.sets) % 6:
        sys.exit("sets come in sixes: arm locus tech source fragments table")
    out = []
    for i in range(0, len(a.sets), 6):
        arm, locus, tech, source, total, table = a.sets[i:i + 6]
        total = int(total)
        e = ERROR[tech]
        rows = []
        with open(table) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                if r["locus"] == locus:
                    rows.append(r)
        screened = len({r["fragment"] for r in rows})
        n_hard = sum(hard(r) for r in rows)
        positive = source in ("locus", "locus1")
        for rule in RULES:
            c = sum(certain(rule, r, e) for r in rows)
            u = sum(1 for r in rows if not certain(rule, r, e) and unrecruited(r))
            amb = len(rows) - c - u
            ch = sum(1 for r in rows if hard(r) and certain(rule, r, e))
            row = {"arm": arm, "locus": locus, "tech": tech, "source": source, "rule": rule, "fragments": total,
                   "screened": screened, "certain": c, "ambiguous": amb, "unrecruited": u,
                   "recall": f"{c / total:.4f}" if positive else "",
                   "recruited": f"{(c + amb) / total:.4f}" if positive else "",
                   "hard_reads": n_hard if positive else "",
                   "hard_recall": f"{ch / n_hard:.4f}" if positive and n_hard else "",
                   "false_certain_per_million": "" if positive or source == "boundary" else f"{1e6 * c / total:.2f}",
                   "upper95_per_million": "" if positive or source == "boundary" else f"{1e6 * poisson_upper(c) / total:.2f}",
                   "per_million_screened": "" if positive or source == "boundary" or screened == 0 else f"{1e6 * c / screened:.1f}",
                   "per_30x": f"{30 * c}" if source in ("genome", "new") else ""}
            out.append(row)
    cols = list(out[0].keys())
    with open(a.out, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in out:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")


if __name__ == "__main__":
    main()
