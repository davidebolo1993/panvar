#!/usr/bin/env python3
"""Checks a locus index dump against brute force: each k-mer's outside count must equal its count over
the genomes minus its count inside the spans the index subtracted (its report's span rows); classes must
follow from it; decoys must be absent from the locus and present in the genomes.

    locus_outside.py <index.tsv report> <dump.tsv> <informative max> <genome.fa>...
Prints problems; exits 1 when any.
"""
import collections
import sys

report, dump, cmax = sys.argv[1], sys.argv[2], int(sys.argv[3])
genomes = sys.argv[4:]
cmax *= len(genomes)   # the cut is per background haplotype; each FASTA here is one
comp = str.maketrans("ACGT", "TGCA")


def canon(x):
    return min(x, x.translate(comp)[::-1])


def kmers(s):
    for piece in s.upper().replace("N", " ").split():
        for i in range(len(piece) - 30):
            yield canon(piece[i:i + 31])


def fasta(path):
    name, parts = None, []
    for line in open(path):
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(parts)
            name, parts = line[1:].split()[0], []
        else:
            parts.append(line.strip())
    if name is not None:
        yield name, "".join(parts)


spans = collections.defaultdict(list)   # genome -> (contig, start, end)
for line in open(report):
    f = line.rstrip("\n").split("\t")
    if f[0] == "span" and f[1] != "genome":
        spans[f[1]].append((f[2], int(f[3]), int(f[4])))
total, own = collections.Counter(), collections.Counter()
for g in genomes:
    name = g.split("/")[-1].rsplit(".fa", 1)[0]
    for contig, s in fasta(g):
        for k in kmers(s):
            total[k] += 1
        for c, a, b in spans.get(name, []):
            if c == contig:
                for k in kmers(s[a:b]):
                    own[k] += 1
bad = 0
locus = set()
decoys = []
n = {"P": 0, "S": 0, "R": 0}
for line in open(dump):
    f = line.rstrip("\n").split("\t")
    if f[0] == "kmer":
        continue
    if f[1] == "D":
        decoys.append(f[0])
        continue
    locus.add(f[0])
    n[f[1]] += 1
    expect = max(0, total[f[0]] - own[f[0]])
    out = int(f[2])
    if out != expect:
        bad += 1
        if bad <= 3:
            print(f"{f[0]}: outside {out}, brute force {expect}")
    want = "P" if out == 0 else "S" if out <= cmax else "R"
    if f[1] != want and not (want == "P" and f[1] == "S"):   # P may be demoted to S (low complexity, span)
        bad += 1
        print(f"{f[0]}: class {f[1]} with outside {out}")
    if f[5] == "1" and f[1] == "R":
        bad += 1
        print(f"{f[0]}: a repeat k-mer in the screen")
for d in decoys:
    if d in locus or total[d] == 0 or total[d] > cmax:
        bad += 1
        if bad <= 6:
            print(f"decoy {d}: in the locus or not 1..{cmax} times in the genomes ({total[d]})")
print(f"locus k-mers P {n['P']} S {n['S']} R {n['R']}, decoys {len(decoys)}, problems {bad}")
sys.exit(1 if bad else 0)
