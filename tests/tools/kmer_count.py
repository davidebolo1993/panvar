#!/usr/bin/env python3
"""Canonical 31-mer counts of FASTA files, summed, as sorted TSV (k-mer, count): an independent check of
`panvar background --dump` (sort its output to compare).

    kmer_count.py <fasta> [<fasta> ...] > counts.tsv
"""
import collections
import gzip
import sys

K = 31
COMP = str.maketrans("ACGT", "TGCA")


def contigs(path):
    opener = gzip.open if path.endswith(".gz") else open
    seq = []
    with opener(path, "rt") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq:
                    yield "".join(seq)
                seq = []
            elif line:
                seq.append(line.upper())
    if seq:
        yield "".join(seq)


counts = collections.Counter()
for path in sys.argv[1:]:
    for s in contigs(path):
        for piece in s.replace("N", " ").split():
            for i in range(len(piece) - K + 1):
                x = piece[i:i + K]
                if set(x) - set("ACGT"):
                    continue
                r = x.translate(COMP)[::-1]
                counts[min(x, r)] += 1
for kmer in sorted(counts):
    sys.stdout.write(f"{kmer}\t{counts[kmer]}\n")
