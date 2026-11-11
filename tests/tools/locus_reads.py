#!/usr/bin/env python3
"""Paired reads for the locus index suite, without a simulator: 150 bp mates of 400 bp fragments drawn
from a stretch of a FASTA, on either strand, error-free; the source is the read name's prefix.

    locus_reads.py <fasta> <contig> <start> <end> <pairs> <name> <seed> <out_prefix>
writes <out_prefix>_1.fq and <out_prefix>_2.fq.
"""
import random
import sys

fa, contig, start, end, pairs, name, seed, out = sys.argv[1:]
start, end, pairs = int(start), int(end), int(pairs)
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
        p = rng.randrange(start, end - 400)
        frag = seq[p:p + 400]
        if rng.random() < 0.5:
            frag = frag.translate(comp)[::-1]
        m1, m2 = frag[:150], frag[-150:].translate(comp)[::-1]
        f1.write(f"@{name}{i}/1\n{m1}\n+\n{'I' * 150}\n")
        f2.write(f"@{name}{i}/2\n{m2}\n+\n{'I' * 150}\n")
