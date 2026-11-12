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

# split/: a locus holding a tandem array (near-identical units), and a genome whose longer array is broken
# across three contigs, the middle one array alone (no anchors reach it): the index must find that piece
# and subtract it; a paralogous array (1% diverged, no anchors either) in another genome must stay outside
split = os.path.join(out, "split")
os.makedirs(split, exist_ok=True)
unit = seq(5500)


def array(n):
    return "".join(mutate(unit, 0.005) for _ in range(n))


left, right = seq(20000), seq(20000)
ref2 = seq(120000) + left + array(4) + right + seq(120000)
ls2 = 120000
le2 = ls2 + len(left) + 4 * 5500 + len(right)
long_array = array(8)
cut1, cut2 = 3 * 5500 + 1700, 6 * 5500 + 2900
hap_s = left + long_array + right
write("split/ref", [("chr1", ref2)])
write("split/gs", [("chr1", ref2[:ls2] + left + long_array[:cut1]),
                   ("x", long_array[cut1:cut2]),
                   ("y", long_array[cut2:] + right + ref2[le2:])])
write("split/gp", [("chr1", seq(60000)), ("p", seq(3000) + "".join(mutate(unit, 0.01) for _ in range(3)) + seq(3000))])
with open(os.path.join(split, "locus.gfa"), "w") as f:
    f.write("H\tVN:Z:1.0\n")
    f.write(f"S\t1\t{ref2[ls2:le2]}\nS\t2\t{hap_s}\n")
    f.write(f"P\tref#1#chr1:{ls2}-{le2}\t1+\t*\n")
    f.write(f"P\tgs#1#chr1:{ls2}-{ls2 + len(hap_s)}\t2+\t*\n")
