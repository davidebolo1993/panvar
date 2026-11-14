#!/usr/bin/env python3
"""Copies per haplotype of each gene in a pangene BED, against `panvar transcripts`' copies.

usage: transcripts_copies.py <prefix>.hst.tsv <pangene.bed> [BED_NAME=GTF_NAME ...]
A gene named differently in the BED and the GTF is given as an alias (STK19=WHR1). Per gene: haplotypes
compared, how many match, both totals, and the mismatches as truth>ours: count.
"""
import collections
import sys

tsv, bed = sys.argv[1], sys.argv[2]
alias = dict(a.split("=") for a in sys.argv[3:])
ours = {}
with open(tsv) as fh:
    ix = {k: i for i, k in enumerate(fh.readline().rstrip("\n").split("\t"))}
    for line in fh:
        r = line.rstrip("\n").split("\t")
        ours[(r[ix["gene_name"]], r[ix["haplotype"]])] = int(r[ix["copies"]])
truth = collections.Counter()
haps = set()
with open(bed) as fh:
    fh.readline()
    for line in fh:
        m, g = line.split("\t")[:2]
        truth[(alias.get(g, g), m)] += 1
        haps.add(m)
ours_haps = {h for _, h in ours}
common = sorted(haps & ours_haps)
print(f"haplotypes: ours {len(ours_haps)}, truth {len(haps)}, common {len(common)}")
print("gene\thaplotypes\tmatch\ttruth_total\tours_total\tmismatches")
for g in sorted({g for g, _ in truth}):
    mm = collections.Counter(f"{truth[(g, h)]}>{ours.get((g, h), 0)}" for h in common if truth[(g, h)] != ours.get((g, h), 0))
    print(f"{g}\t{len(common)}\t{len(common) - sum(mm.values())}\t{sum(truth[(g, h)] for h in common)}\t"
          f"{sum(ours.get((g, h), 0) for h in common)}\t{','.join(f'{k}:{n}' for k, n in sorted(mm.items())) or '.'}")
