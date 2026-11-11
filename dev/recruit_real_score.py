#!/usr/bin/env python3
"""Compares recruitment of real pairs with where an aligner put them (dev/recruit_real.sh).

    recruit_real_score.py --out <summary.tsv> <locus> <interval> <regions.bed> <aligned.tsv> <reads_1.fq.gz> <table> [...]
The aligned table has one row per primary alignment: name, mate, contig, position, MAPQ and, when
known, the alignment's end (else position + 150). A mate is in the locus when its alignment overlaps it.
Each pair is classed by its two mates:
    locus_q    both in the locus, both MAPQ >= 30 (the aligner is sure: expect certain)
    locus_low  both in the locus, a MAPQ below 30 (paralogous for the aligner)
    edge       one mate in the locus, the other outside it
    flank_q    both in the 200 kb flanks, 2 kb or more from the locus, MAPQ >= 30 (expect never certain)
    decoy_q    both in decoy windows, MAPQ >= 30 (expect never certain)
    decoy_low  both in decoy windows, a MAPQ below 30
    other      anything else (a mate elsewhere; pairs with an unmapped mate or a mate outside the fetched
               regions are not fetched as pairs and are left out)
and counted by recruit's fate.
"""
import argparse
import bisect
import gzip
import sys

CLASSES = ["locus_q", "locus_low", "edge", "flank_q", "decoy_q", "decoy_low", "other"]


def names(path):
    with gzip.open(path, "rt") as f:
        for i, line in enumerate(f):
            if i % 4 == 0:
                n = line[1:].split()[0]
                yield n[:-2] if n[-2:] in ("/1", "/2") else n


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("sets", nargs="+")
    a = ap.parse_args()
    if len(a.sets) % 6:
        sys.exit("sets come in sixes: locus interval regions aligned reads table")
    out = []
    for i in range(0, len(a.sets), 6):
        locus, interval, bed, aligned, reads, table = a.sets[i:i + 6]
        chrom, span = interval.rsplit(":", 1)
        ls, le = (int(x) for x in span.split("-"))
        decoys = {}
        with open(bed) as f:
            next(f)   # the locus and its flanks
            for line in f:
                c, s, e = line.split()[:3]
                decoys.setdefault(c, []).append((int(s), int(e)))
        for c in decoys:
            decoys[c].sort()
        starts = {c: [s for s, _ in v] for c, v in decoys.items()}

        def place(c, pos, end):
            if c == chrom and pos < le and end > ls:
                return "locus"
            if c in decoys:
                j = bisect.bisect_right(starts[c], pos) - 1
                if j >= 0 and pos < decoys[c][j][1]:
                    return "decoy"
            if c == chrom and (ls - 200000 <= pos < ls - 2000 or le + 2000 <= pos < le + 200000):
                return "flank"
            return "elsewhere"

        mates = {}
        with open(aligned) as f:
            for line in f:
                x = line.rstrip("\n").split("\t")
                n, m, c, pos, q = x[:5]
                start = int(pos) - 1   # SAM positions are 1-based; the interval is 0-based
                end = int(x[5]) - 1 if len(x) > 5 else start + 150
                mates.setdefault(n, {})[m] = (place(c, start, end), int(q))
        fate = {}
        with open(table) as f:
            head = f.readline().rstrip("\n").split("\t")
            for line in f:
                r = dict(zip(head, line.rstrip("\n").split("\t")))
                fate[r["fragment"]] = r["fate"]
        count = {k: {"n": 0, "certain": 0, "ambiguous": 0, "unrecruited": 0} for k in CLASSES}
        for n in names(reads):
            m = mates.get(n, {})
            p1, q1 = m.get("1", ("elsewhere", 0))
            p2, q2 = m.get("2", ("elsewhere", 0))
            good = q1 >= 30 and q2 >= 30
            if p1 == p2 == "locus":
                k = "locus_q" if good else "locus_low"
            elif "locus" in (p1, p2):
                k = "edge"
            elif p1 == p2 == "flank" and good:
                k = "flank_q"
            elif p1 == p2 == "decoy":
                k = "decoy_q" if good else "decoy_low"
            else:
                k = "other"
            count[k]["n"] += 1
            count[k][fate.get(n, "unrecruited")] += 1
        for k in CLASSES:
            c = count[k]
            out.append({"locus": locus, "class": k, "pairs": c["n"], "certain": c["certain"], "ambiguous": c["ambiguous"],
                        "unrecruited": c["unrecruited"],
                        "certain_share": f"{c['certain'] / c['n']:.4f}" if c["n"] else "",
                        "recruited_share": f"{(c['certain'] + c['ambiguous']) / c['n']:.4f}" if c["n"] else ""})
    cols = list(out[0].keys())
    with open(a.out, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in out:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")


if __name__ == "__main__":
    main()
