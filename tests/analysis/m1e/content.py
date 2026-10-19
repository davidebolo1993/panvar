#!/usr/bin/env python3
"""M1-E: pooled block content (canonical k-mers, k = 31 and 151) of each candidate against the truth blocks, and, per patch
set, the content difference between the matched 2+2 and 1+3 arrangements (what their different joins create).

  content.py BLOCKS_FA OUT.tsv
"""
import sys
from collections import Counter

COMP = str.maketrans("ACGT", "TGCA")


def fasta(path):
    out, name = {}, None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            name = line[1:].split()[0]
            out[name] = []
        elif name:
            out[name].append(line)
    return {k: "".join(v) for k, v in out.items()}


def kmers(seq, k):
    c, rc, n = Counter(), seq.translate(COMP)[::-1], len(seq)
    for i in range(n - k + 1):
        f, r = seq[i:i + k], rc[n - k - i:n - i]
        c[hash(f if f < r else r)] += 1
    return c


def diff(a, b):
    return sum(abs(a[x] - b[x]) for x in set(a) | set(b)) / 2


def main(blocks_fa, out):
    b = fasta(blocks_fa)
    labels = sorted({n.rsplit(".", 1)[0] for n in b if not n.startswith("truth.")})
    lines = ["label\tk\tpooled_vs_truth\tvs_matched_arrangement\n"]
    for k in (31, 151):
        cache = {n: kmers(s, k) for n, s in b.items()}
        tp = cache["truth.1"] + cache["truth.2"]
        pooled = {l: cache[f"{l}.1"] + cache[f"{l}.2"] for l in labels}
        for l in labels:
            other = l.replace("a22_", "e13_") if l.startswith("a22_") else l.replace("e13_", "a22_") if l.startswith("e13_") else None
            lines.append(f"{l}\t{k}\t{diff(pooled[l], tp):g}\t{diff(pooled[l], pooled[other]) if other else 'NA':}\n")
    open(out, "w").write("".join(lines))
    print("".join(lines), end="")


if __name__ == "__main__":
    main(*sys.argv[1:3])
