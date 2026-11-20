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

# split_h/: the locus of split/, and a genome with no copy of it but a contig shorter than a piece window
# holding nothing but its array (a piece), and a short contig of other sequence
os.makedirs(os.path.join(out, "split_h"), exist_ok=True)
write("split_h/gh", [("chr1", seq(30000)), ("s", long_array[2000:6500]), ("n", seq(3000))])

# ends/: a locus whose two 10 kb ends are 30% diverged from the reference in another genome: no reference
# anchor holds there, so that copy is anchored in its middle only and its ends are uncertain. Its own k-mers
# there occur nowhere else, and must stay private: the copy is walked out through them. A copy of one end
# far away on the same contig stays outside.
os.makedirs(os.path.join(out, "ends"), exist_ok=True)
ref3 = seq(300000)
ls3, le3 = 120000, 160000
hap_d = mutate(ref3[ls3:ls3 + 10000], 0.3) + ref3[ls3 + 10000:le3 - 10000] + mutate(ref3[le3 - 10000:le3], 0.3)
write("ends/ref", [("chr1", ref3)])
write("ends/gd", [("chr1", ref3[:ls3] + hap_d + ref3[le3:] + seq(20000) + hap_d[2000:4000] + seq(20000))])
with open(os.path.join(out, "ends", "locus.gfa"), "w") as f:
    f.write("H\tVN:Z:1.0\n")
    f.write(f"S\t1\t{ref3[ls3:le3]}\nS\t2\t{hap_d}\n")
    f.write(f"P\tref#1#chr1:{ls3}-{le3}\t1+\t*\n")
    f.write(f"P\tgd#1#chr1:{ls3}-{le3}\t2+\t*\n")

# rep/: the locus of ends/ with a 3 kb repeat element (20 copies elsewhere) between its anchored middle and a 6 kb
# end that differs from the reference: the walk crosses the repeat to that end. The same element right past the
# locus's other end, in the flank, is not crossed: no locus sequence follows it.
os.makedirs(os.path.join(out, "rep"), exist_ok=True)
ref4 = list(seq(400000))
ls4, le4 = 120000, 160000
element = seq(3000)
ref4[ls4 + 6000:ls4 + 9000] = list(element)
for i in range(10):
    a = 200000 + i * 15000
    ref4[a:a + 3000] = list(element)
ref4 = "".join(ref4)
hap_r = mutate(ref4[ls4:ls4 + 6000], 0.3) + ref4[ls4 + 6000:le4 - 6000] + mutate(ref4[le4 - 6000:le4], 0.3)
far = "".join(seq(2000) + element for _ in range(10))
write("rep/ref", [("chr1", ref4)])
write("rep/gr", [("chr1", ref4[:ls4] + hap_r + element + ref4[le4:]), ("x", far)])
with open(os.path.join(out, "rep", "locus.gfa"), "w") as f:
    f.write("H\tVN:Z:1.0\n")
    f.write(f"S\t1\t{ref4[ls4:le4]}\nS\t2\t{hap_r}\n")
    f.write(f"P\tref#1#chr1:{ls4}-{le4}\t1+\t*\n")
    f.write(f"P\tgr#1#chr1:{ls4}-{le4}\t2+\t*\n")

# split_f/: the locus of split/, and a genome whose copy is whole and which holds more of its array mid-way
# along the same contig, far from the copy with other sequence on both sides (as a close paralog on a
# chromosome-long contig): the index looks in every genome, and that stretch is no piece
os.makedirs(os.path.join(out, "split_f"), exist_ok=True)
arr_f = array(6)
hap_f = left + arr_f + right
write("split_f/ref", [("chr1", ref2)])
write("split_f/gf", [("chr1", ref2[:ls2] + hap_f + ref2[le2:] + seq(150000) + arr_f[:4 * 5500] + seq(20000))])
with open(os.path.join(out, "split_f", "locus.gfa"), "w") as f:
    f.write("H\tVN:Z:1.0\n")
    f.write(f"S\t1\t{ref2[ls2:le2]}\nS\t2\t{hap_f}\n")
    f.write(f"P\tref#1#chr1:{ls2}-{le2}\t1+\t*\n")
    f.write(f"P\tgf#1#chr1:{ls2}-{ls2 + len(hap_f)}\t2+\t*\n")

# logo/: a locus with no paralog in the coordinate reference, and genomes holding the locus and one homolog each,
# held out in calibration: gu a 1% diverged paralog of 10 kb of it (its alleles there are the locus's, shared only
# with gu), gd an identical copy of 8 kb (shared k-mers only). Each homolog is cut in pieces of 400 bases (each with
# the next one's first 30, so every k-mer of the stretch is in one) laid in reverse order: no three anchors in order,
# so neither is located as a copy of the locus. A held-out haplotype gh is in the graph alone.
logo = os.path.join(out, "logo")
os.makedirs(logo, exist_ok=True)
ref5 = seq(300000)
ls5, le5 = 120000, 160000


def scrambled(a, b, rate):
    pieces = [mutate(ref5[x:min(x + 430, b)], rate) for x in range(a, b - 30, 400)]
    return "".join(p + seq(100) for p in reversed(pieces))


write("logo/ref", [("chr1", ref5)])
write("logo/gu", [("chr1", ref5), ("u", seq(20000) + scrambled(130000, 140000, 0.01) + seq(20000))])
write("logo/gd", [("chr1", ref5), ("d", seq(20000) + scrambled(145000, 153000, 0.0) + seq(20000))])
hap_h = mutate(ref5[ls5:le5], 0.002)
with open(os.path.join(logo, "locus.gfa"), "w") as f:
    f.write("H\tVN:Z:1.0\n")
    f.write(f"S\t1\t{ref5[ls5:le5]}\nS\t2\t{hap_h}\n")
    f.write(f"P\tref#1#chr1:{ls5}-{le5}\t1+\t*\n")
    f.write(f"P\tgh#1#chr1:{ls5}-{le5}\t2+\t*\n")

# gc: an identical copy of 8 kb of the locus, cut in pieces of 1 kb laid in reverse order with nothing between: no
# three anchors in order, so it is not located as a copy, yet nearly every k-mer of it is the locus's
write("logo/gc", [("chr1", ref5), ("c", seq(20000) + "".join(ref5[x:x + 1000] for x in reversed(range(121000, 129000, 1000))) + seq(20000))])
# ge: a 3% diverged paralog of the same 8 kb, in pieces like gu's: few of its k-mers are the locus's
write("logo/ge", [("chr1", ref5), ("e", seq(20000) + scrambled(121000, 129000, 0.03) + seq(20000))])
