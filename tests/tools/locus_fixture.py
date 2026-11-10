#!/usr/bin/env python3
"""Genomes and a locus graph for the locus index suite, the same on every platform: a reference whose
locus has a diverged paralog elsewhere (so some locus k-mers are shared and the paralog yields decoys),
two more genomes carrying the locus with small changes, and a GFA whose paths are the three locus
copies (one node each, so the graph is trivially valid).

    locus_fixture.py <out_dir>
"""
import os
import random
import sys

out = sys.argv[1]
os.makedirs(out, exist_ok=True)
rng = random.Random(11)


def seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def mutate(s, rate):
    b = list(s)
    for i in range(len(b)):
        if rng.random() < rate:
            b[i] = rng.choice([c for c in "ACGT" if c != b[i]])
    return "".join(b)


def write(name, contigs):
    with open(os.path.join(out, name + ".fa"), "w") as f:
        for cname, s in contigs:
            f.write(f">{cname}\n")
            for i in range(0, len(s), 60):
                f.write(s[i:i + 60] + "\n")


ls, le = 120000, 160000
ref = list(seq(300000))
paralog = mutate("".join(ref[125000:155000]), 0.02)          # a diverged copy of most of the locus
ref[230000:260000] = list(paralog)
ref = "".join(ref)
hap_a = mutate(ref[ls:le], 0.002)
hap_b = ref[ls:140000] + seq(500) + ref[140000:le]            # an insertion
write("ref", [("chr1", ref)])
write("ga", [("chr1", ref[:ls] + hap_a + ref[le:])])
write("gb", [("chr1", ref[:ls] + hap_b + ref[le:])])
with open(os.path.join(out, "locus.gfa"), "w") as f:
    f.write("H\tVN:Z:1.0\n")
    for i, s in enumerate([ref[ls:le], hap_a, hap_b], 1):
        f.write(f"S\t{i}\t{s}\n")
    f.write(f"P\tref#1#chr1:{ls}-{le}\t1+\t*\n")
    f.write(f"P\tga#1#chr1:{ls}-{le}\t2+\t*\n")
    f.write(f"P\tgb#1#chr1:{ls}-{le + 500}\t3+\t*\n")
