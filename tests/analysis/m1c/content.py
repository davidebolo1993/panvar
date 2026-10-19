#!/usr/bin/env python3
"""M1-C: read-visible content of each declared alternative (diagnostic evaluator; truth allowed).

Over WHOLE haplotypes (so new junctions created by a reallocation are counted), canonical k-mer multisets, k = 31 and 151:
  vs_reported_pooled  = sum_k |c_alt1 + c_alt2 - c_rep1 - c_rep2| / 2   (what a reallocation changed, junctions included)
  vs_truth_pooled     = sum_k |c_alt1 + c_alt2 - c_truth1 - c_truth2| / 2   (allocation-blind distance to truth)
  vs_truth_per_hap    = min over assignments of the per-haplotype sums / 2     (allocation-aware distance to truth)
A reallocation preserves the unit inventory but not necessarily the pooled fragment content; this measures the difference.

  content.py --explicit EXPLICIT.tsv --out OUT_PREFIX
"""
import argparse
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
    c = Counter()
    rc = seq.translate(COMP)[::-1]
    n = len(seq)
    for i in range(n - k + 1):
        f = seq[i:i + k]
        r = rc[n - k - i:n - i]
        c[hash(f if f < r else r)] += 1
    return c


def diff(a, b):
    return sum(abs(a[x] - b[x]) for x in set(a) | set(b)) / 2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--explicit", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    alts, seqs = [], {}
    for line in open(a.explicit):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if f[1] not in seqs:
            seqs[f[1]] = fasta(f[1])
        alts.append((f[0], seqs[f[1]][f[2]], seqs[f[1]][f[3]]))
    by = {x[0]: x for x in alts}
    lines = ["label\tk\tvs_reported_pooled\tvs_truth_pooled\tvs_truth_per_hap\n"]
    for k in (31, 151):
        cache = {}

        def cnt(s):
            key = hash(s)
            if key not in cache:
                cache[key] = kmers(s, k)
            return cache[key]

        rep = cnt(by["reported"][1]) + cnt(by["reported"][2])
        t1, t2 = cnt(by["truth"][1]), cnt(by["truth"][2])
        tp = t1 + t2
        for label, s1, s2 in alts:
            c1, c2 = cnt(s1), cnt(s2)
            pooled = c1 + c2
            per_hap = min(diff(c1, t1) + diff(c2, t2), diff(c1, t2) + diff(c2, t1))
            lines.append(f"{label}\t{k}\t{diff(pooled, rep):g}\t{diff(pooled, tp):g}\t{per_hap:g}\n")
    open(a.out + ".content.tsv", "w").write("".join(lines))
    print("".join(lines), end="")


if __name__ == "__main__":
    main()
