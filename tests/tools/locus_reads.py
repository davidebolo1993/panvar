#!/usr/bin/env python3
"""Paired reads for the locus index suite, without a simulator: mates (150 bp) of fragments (400 bp) drawn
from a stretch of a FASTA, on either strand, error-free; the source is the read name's prefix.

    locus_reads.py <fasta> <contig> <start> <end> <pairs> <name> <seed> <out_prefix> [insert] [read length]
writes <out_prefix>_1.fq and <out_prefix>_2.fq.
"""
import random
import sys

fa, contig, start, end, pairs, name, seed, out = sys.argv[1:9]
start, end, pairs = int(start), int(end), int(pairs)
insert = int(sys.argv[9]) if len(sys.argv) > 9 else 400
length = int(sys.argv[10]) if len(sys.argv) > 10 else 150
rng = random.Random(int(seed))
seq, cur, parts = None, None, []
for line in open(fa):
    line = line.strip()
    if line.startswith(">"):
        if cur == contig:
            seq = "".join(parts)
        cur, parts = line[1:].split()[0], []
    else:
        parts.append(line)
if cur == contig:
    seq = "".join(parts)
comp = str.maketrans("ACGT", "TGCA")
with open(out + "_1.fq", "w") as f1, open(out + "_2.fq", "w") as f2:
    for i in range(pairs):
        p = rng.randrange(start, end - insert)
        frag = seq[p:p + insert]
        if rng.random() < 0.5:
            frag = frag.translate(comp)[::-1]
        m1, m2 = frag[:length], frag[-length:].translate(comp)[::-1]
        f1.write(f"@{name}{i}/1\n{m1}\n+\n{'I' * len(m1)}\n")
        f2.write(f"@{name}{i}/2\n{m2}\n+\n{'I' * len(m2)}\n")
