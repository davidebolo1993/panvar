#!/usr/bin/env python3
"""M1-B Amendment 2 (diagnostic evaluator; truth allowed): is the target-block competition decided by content that short
fragments can see, or by arrangement they cannot?

For each declared alternative (Amendment 1 candidates file), take its two target-block alleles and compare their k-mer
content with the truth haplotypes' target-block sequences (canonical k-mers, k = 31 and k = 151):
  pooled   = sum_k |c_alt1(k) + c_alt2(k) - c_truth1(k) - c_truth2(k)| / 2   -- allocation-blind: homologues pooled
  per_hap  = min over the two assignments of (sum_k |c_alt_i - c_truth_a| + sum_k |c_alt_j - c_truth_b|) / 2
             -- allocation-aware
  missing / extra = pooled truth k-mer copies absent from the alternative / alternative copies absent from truth.
A diplotype can be allocation-wrong yet pooled-content-closest; a short-read likelihood that cannot see across a repeat
unit would then prefer it. Nothing here is fed back into any caller or score.

  content.py --sequences PREFIX --block ID --candidates FILE --out OUT_PREFIX
"""
import argparse
from collections import Counter

COMP = str.maketrans("ACGTN", "TGCAN")


def fasta(path):
    out, name = {}, None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            name = line[1:].split()[0]
            out[name] = []
        elif name:
            out[name].append(line)
    return {k: "".join(v).upper() for k, v in out.items()}


def kmers(seq, k):
    c = Counter()
    rc = seq.translate(COMP)[::-1]
    n = len(seq)
    for i in range(n - k + 1):
        f = seq[i:i + k]
        r = rc[n - k - i:n - i]
        c[f if f < r else r] += 1
    return c


def diff(a, b):
    keys = set(a) | set(b)
    missing = sum(max(b[k] - a[k], 0) for k in keys)   # in b (truth) but not in a
    extra = sum(max(a[k] - b[k], 0) for k in keys)
    return (missing + extra) / 2, missing, extra


def main():
    ap = argparse.ArgumentParser()
    for x in ("sequences", "block", "candidates", "out"):
        ap.add_argument("--" + x, required=True)
    a = ap.parse_args()
    alleles = fasta(f"{a.sequences}.{a.block}.alleles.fa")
    truth = fasta(f"{a.sequences}.{a.block}.truth.fa")
    t = int(next(l for l in open(a.candidates) if l.startswith("#")).split("block index ")[1].split(")")[0])
    rows = []
    for line in open(a.candidates):
        if line.startswith("#") or not line.strip():
            continue
        label, h1, h2 = line.rstrip("\n").split("\t")[:3]
        rows.append((label, h1.split(",")[t], h2.split(",")[t]))
    out = ["label\tallele1\tallele2\tk\tpooled\tpooled_missing\tpooled_extra\tper_hap\n"]
    cache = {}
    for k in (31, 151):
        t1, t2 = kmers(truth["truth1"], k), kmers(truth["truth2"], k)
        tp = t1 + t2
        for label, x, y in rows:
            for z in (x, y):
                if (z, k) not in cache:
                    cache[(z, k)] = kmers(alleles[f"allele_{z}"], k)
            cx, cy = cache[(x, k)], cache[(y, k)]
            pooled, miss, extra = diff(cx + cy, tp)
            per_hap = min(diff(cx, t1)[0] + diff(cy, t2)[0], diff(cx, t2)[0] + diff(cy, t1)[0])
            out.append(f"{label}\t{x}\t{y}\t{k}\t{pooled:g}\t{miss}\t{extra}\t{per_hap:g}\n")
    open(a.out + ".content.tsv", "w").write("".join(out))
    print("".join(out), end="")


if __name__ == "__main__":
    main()
