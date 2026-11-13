#!/usr/bin/env python3
"""Scores recruitment tables against the truth of simulated read sets.

    recruit_score.py --sim <sim dir> --out <summary.tsv> <arm> <locus> <tech> <source> <reads> <table> [...]
Each set is six arguments: the arm (main: GRCh38 + CHM13 background; g38: GRCh38 only; wg: the whole
GRCh38), the locus, the technology (short, illumina, ont, hifi), the source (locus, locus1 = one mate only,
near, decoy, boundary, genome, new, wgnew), the reads (the first mates' FASTQ, which is counted and, for
boundary sets, holds each fragment's span) and recruit's table.

Columns: fragments, screened; certain (of which partial and rescued), ambiguous and unrecruited;
recall and recruited (certain or ambiguous) for the locus sources; false certain per million fragments
with its 95% upper bound for the negative sources, and for the genome sets the false-certain fragments
expected in 30x. Boundary sets are split by the true span: fragments wholly inside the locus, straddling
an edge, and wholly outside (where certain is false); for long reads straddling, the median difference
between the locus segment recruit reports and the part truly inside.
"""
import argparse
import gzip
import math
import re
import sys

DEPTH = {"genome": 1.0, "new": 1.0, "wgnew": None}
WGNEW_DEPTH = {"short": 0.5, "illumina": 0.5, "ont": 0.2, "hifi": 0.2}
NEGATIVE = {"near", "decoy", "genome", "new", "wgnew"}


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


def fragment_name(h):
    n = h[1:].split()[0]
    return n[:-2] if len(n) > 2 and n[-2] == "/" and n[-1] in "12" else n


def read_fastq(path):
    """(name, length) of every record."""
    with gzip.open(path, "rt") as f:
        while True:
            h = f.readline()
            if not h:
                return
            s = f.readline()
            f.readline()
            f.readline()
            yield fragment_name(h), len(s.rstrip("\n"))


SPAN = {   # the fragment's span in its boundary source sequence, 0-based end-exclusive
    "short": re.compile(r"boundary\|(\w+?)_(\d+)_(\d+)_(\d+)_(\d+)_\d+:"),       # wgsim: 1-based outer ends
    "illumina": re.compile(r"boundary\|(\w+?)_(\d+)_(\d+)_(\d+)_(\d+)_[+-]$"),
    "ont": re.compile(r"boundary\|(\w+?)_(\d+)_(\d+)_(\d+)_(\d+)_[+-]$"),
    "hifi": re.compile(r"boundary\|(\w+?)_(\d+)_(\d+)_(\d+)_(\d+)_[+-]$"),
}


def truth_span(tech, name, length):
    m = SPAN[tech].search(name)
    if not m:
        sys.exit(f"recruit_score: no span in {name}")
    s0 = int(m.group(2))
    a, b = int(m.group(4)), int(m.group(5))
    if tech == "short":
        a -= 1
    return s0 + a, s0 + b


def locus_interval(sim, locus):
    with open(f"{sim}/{locus}/{locus}.index.tsv") as f:
        for line in f:
            x = line.rstrip("\n").split("\t")
            if x[0] == "locus" and x[1] == "interval":
                chrom, rest = x[2].rsplit(":", 1)
                a, b = rest.split("-")
                return int(a), int(b)
    sys.exit(f"recruit_score: no interval for {locus}")


def segment_length(seg):
    if seg == ".":
        return 0
    a, b = seg.split("-")
    return int(b) - int(a)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sim", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("sets", nargs="+")
    a = ap.parse_args()
    if len(a.sets) % 6:
        sys.exit("sets come in sixes: arm locus tech source reads table")
    counted = {}
    out = []
    for i in range(0, len(a.sets), 6):
        arm, locus, tech, source, reads, table = a.sets[i:i + 6]
        if reads not in counted:
            counted[reads] = list(read_fastq(reads)) if source == "boundary" else sum(1 for _ in read_fastq(reads))
        total = len(counted[reads]) if source == "boundary" else counted[reads]
        rows = {}
        with open(table) as f:
            head = f.readline().rstrip("\n").split("\t")
            for line in f:
                r = dict(zip(head, line.rstrip("\n").split("\t")))
                if r["locus"] == locus:
                    rows[r["fragment"]] = r
        fate = [r["fate"] for r in rows.values()]
        c = fate.count("certain")
        amb = fate.count("ambiguous")
        row = {"arm": arm, "locus": locus, "tech": tech, "source": source, "fragments": total, "screened": len(rows),
               "certain": c, "partial": sum(1 for r in rows.values() if r["fate"] == "certain" and r["reason"] == "partial"),
               "rescued": sum(1 for r in rows.values() if r["fate"] == "certain" and r["reason"] == "rescued"),
               "ambiguous": amb, "unrecruited": fate.count("unrecruited"), "recall": "", "recruited": "",
               "false_certain_per_million": "", "upper95_per_million": "", "per_30x": "",
               "in": "", "in_certain": "", "straddle": "", "straddle_certain": "", "straddle_recruited": "",
               "out": "", "out_certain": "", "segment_error_median": ""}
        if source in ("locus", "locus1"):
            row["recall"] = f"{c / total:.4f}"
            row["recruited"] = f"{(c + amb) / total:.4f}"
        if source in NEGATIVE:
            row["false_certain_per_million"] = f"{1e6 * c / total:.2f}"
            row["upper95_per_million"] = f"{1e6 * poisson_upper(c) / total:.2f}"
            depth = DEPTH.get(source, 0) or (WGNEW_DEPTH[tech] if source == "wgnew" else 0)
            if depth:
                row["per_30x"] = f"{30 / depth * c:.0f}"
        if source == "boundary":
            ls, le = locus_interval(a.sim, locus)
            n = {"in": 0, "straddle": 0, "out": 0}
            cert = {"in": 0, "straddle": 0, "out": 0}
            recr = 0
            errors = []
            for name, length in counted[reads]:
                x, y = truth_span(tech, name, length)
                inside = max(0, min(y, le) - max(x, ls))
                kind = "in" if inside >= y - x else "out" if inside == 0 else "straddle"
                n[kind] += 1
                r = rows.get(name)
                if r is None:
                    continue
                if r["fate"] == "certain":
                    cert[kind] += 1
                    if kind == "straddle" and tech in ("ont", "hifi"):
                        errors.append(abs(segment_length(r["seg1"]) - inside))
                if kind == "straddle" and r["fate"] in ("certain", "ambiguous"):
                    recr += 1
            errors.sort()
            row.update({"in": n["in"], "in_certain": cert["in"], "straddle": n["straddle"], "straddle_certain": cert["straddle"],
                        "straddle_recruited": recr, "out": n["out"], "out_certain": cert["out"],
                        "segment_error_median": errors[len(errors) // 2] if errors else ""})
        out.append(row)
    cols = list(out[0].keys())
    with open(a.out, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in out:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")


if __name__ == "__main__":
    main()
