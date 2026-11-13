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

# split_t/: a locus ending in its array, and a genome with a much longer array on one contig: no inside
# anchor follows the array and the flank sits too far past it, so the copy's end is estimated and stops
# inside the array; the rest of the array, on the same contig, must be found as a piece
left_t = seq(20000)
ref_t = seq(120000) + left_t + array(4) + seq(120000)
le_t = 120000 + len(left_t) + 4 * 5500
hap_t = left_t + array(16)
os.makedirs(os.path.join(out, "split_t"), exist_ok=True)
write("split_t/ref", [("chr1", ref_t)])
write("split_t/gl", [("chr1", ref_t[:120000] + hap_t + ref_t[le_t:])])
with open(os.path.join(out, "split_t", "locus.gfa"), "w") as f:
    f.write("H\tVN:Z:1.0\n")
    f.write(f"S\t1\t{ref_t[120000:le_t]}\nS\t2\t{hap_t}\n")
    f.write(f"P\tref#1#chr1:120000-{le_t}\t1+\t*\n")
    f.write(f"P\tgl#1#chr1:120000-{120000 + len(hap_t)}\t2+\t*\n")

# split_o/: the locus of split/, and a genome whose copy is whole but which holds more of its array (units its
# haplotype has) on a contig of its own (no copy is broken, no end estimated): the locus's tandem k-mers left with a few copies outside
# send the index looking in every genome
os.makedirs(os.path.join(out, "split_o"), exist_ok=True)
arr_o = array(6)
hap_o = left + arr_o + right
write("split_o/ref", [("chr1", ref2)])
write("split_o/go", [("chr1", ref2[:ls2] + hap_o + ref2[le2:]), ("orphan", arr_o[:4 * 5500])])
with open(os.path.join(out, "split_o", "locus.gfa"), "w") as f:
    f.write("H\tVN:Z:1.0\n")
    f.write(f"S\t1\t{ref2[ls2:le2]}\nS\t2\t{hap_o}\n")
    f.write(f"P\tref#1#chr1:{ls2}-{le2}\t1+\t*\n")
    f.write(f"P\tgo#1#chr1:{ls2}-{ls2 + len(hap_o)}\t2+\t*\n")

# split_e/: the locus of split/, and arrays that run into a contig end through a stub with errors (no window
# fits there, and its k-mers fall short of the density a window needs; the graph holds the array without the
# errors): in ge, a contig holding more of the array of a whole copy; in gt, a copy reaching its contig's end
os.makedirs(os.path.join(out, "split_e"), exist_ok=True)
arr_e, arr_g = array(6), array(7)
stub_e = mutate(arr_e[4 * 5500:4 * 5500 + 1700], 0.006)
stub_g = mutate(arr_g[5 * 5500:5 * 5500 + 2000], 0.01)
hap_e = left + arr_e + right
hap_g = left + arr_g + right
write("split_e/ref", [("chr1", ref2)])
write("split_e/ge", [("chr1", ref2[:ls2] + hap_e + ref2[le2:]), ("o", arr_e[:4 * 5500] + stub_e)])
write("split_e/gt", [("chr1", ref2[:ls2] + left + arr_g[:5 * 5500] + stub_g)])
with open(os.path.join(out, "split_e", "locus.gfa"), "w") as f:
    f.write("H\tVN:Z:1.0\n")
    f.write(f"S\t1\t{ref2[ls2:le2]}\nS\t2\t{hap_e}\nS\t3\t{hap_g}\n")
    f.write(f"P\tref#1#chr1:{ls2}-{le2}\t1+\t*\n")
    f.write(f"P\tge#1#chr1:{ls2}-{ls2 + len(hap_e)}\t2+\t*\n")
    f.write(f"P\tgt#1#chr1:{ls2}-{ls2 + len(hap_g)}\t3+\t*\n")

# split_c/: the locus of split/, and a genome with no copy of it: a paralog of the array (a haplotype's units
# the coordinate reference lacks, with 50 adjacent differences per 5 kb, dense enough for a window) between
# other sequence on a contig of its own, which must stay outside; and a contig of the array alone, which is
# a piece


def clustered(s, block=5000, run=50):
    b = list(s)
    for x in range(0, len(b) - block + 1, block):
        at = x + rng.randrange(block - run)
        for i in range(at, at + run):
            b[i] = rng.choice([c for c in "ACGT" if c != b[i]])
    return "".join(b)


os.makedirs(os.path.join(out, "split_c"), exist_ok=True)
write("split_c/gc", [("chr1", seq(60000)), ("p", seq(3000) + clustered(long_array[5500:4 * 5500]) + seq(3000)),
                     ("o", long_array[4 * 5500:7 * 5500])])
