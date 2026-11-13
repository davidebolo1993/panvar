#!/usr/bin/env python3
"""Genomes for the background command suite: a reference and variants of it (a moved locus, a reversed
copy, a duplicate), each 300-400 kb, the same on every platform.

    background_fixture.py <out_dir>
"""
import os
import random
import sys

out = sys.argv[1]
os.makedirs(out, exist_ok=True)
rng = random.Random(7)


def seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def write(name, contigs):
    with open(os.path.join(out, name + ".fa"), "w") as f:
        for cname, s in contigs:
            f.write(f">{cname}\n")
            for i in range(0, len(s), 60):
                f.write(s[i:i + 60] + "\n")


ref = seq(300000)
ls, le = 120000, 160000                       # the locus
write("ref", [("chr1", ref), ("chr2", seq(30000))])
write("g1", [("chr1", ref[:50000] + seq(20000) + ref[50000:])])            # the locus moved by 20 kb
write("g2", [("a", seq(50000) + rc(ref[ls - 15000:le + 15000]) + seq(50000)), ("b", ref)])   # a reversed extra copy
write("g3", [("chr1", ref[:le] + ref[ls:le] + ref[le:])])                  # a tandem duplicate
write("g4", [("chr1", ref[:200000] + seq(5000) + ref[200000:])])
with open(os.path.join(out, "reads.fq"), "w") as f:
    f.write("@r1\n" + seq(100) + "\n+\n" + "I" * 100 + "\n")
