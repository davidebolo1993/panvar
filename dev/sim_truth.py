#!/usr/bin/env python3
"""Puts the true source span into simulated read names, so the read files alone carry their truth.

    sim_truth.py maf <reads.fq.gz> <reads.maf.gz> <sequence name> >> out.fq     (pbsim3, one sequence)
    sim_truth.py art <prefix> <out_1.fq.gz> <out_2.fq.gz>                       (ART paired, run with -sam)

Renamed reads end in _<sequence>_<start>_<end>_<strand>, 0-based and end-exclusive in the source sequence;
for a pair the span is the fragment's. Long reads that are mostly N are dropped.
"""
import gzip
import re
import sys


def fastq(f):
    while True:
        h = f.readline()
        if not h:
            return
        s, p, q = f.readline(), f.readline(), f.readline()
        yield h.rstrip("\n"), s.rstrip("\n"), p.rstrip("\n"), q.rstrip("\n")


def maf_spans(path):
    spans = {}
    with gzip.open(path, "rt") as f:
        lines = []
        for line in f:
            if line.startswith("s "):
                lines.append(line.split())
                if len(lines) == 2:
                    ref, read = lines
                    start, size = int(ref[2]), int(ref[3])
                    strand = "-" if "-" in (ref[4], read[4]) else "+"
                    spans[read[1]] = (start, start + size, strand)
                    lines = []
            elif line.startswith("a"):
                lines = []
    return spans


def do_maf(fq, maf, seqname):
    spans = maf_spans(maf)
    out = sys.stdout
    with gzip.open(fq, "rt") as f:
        for h, s, p, q in fastq(f):
            name = h[1:].split()[0]
            if name not in spans:
                sys.exit(f"sim_truth: {name} is not in {maf}")
            if 2 * sum(1 for c in s if c in "Nn") > len(s):
                continue
            a, b, strand = spans[name]
            out.write(f"@{name}_{seqname}_{a}_{b}_{strand}\n{s}\n+\n{q}\n")


def do_art(prefix, out1, out2):
    span = {}
    with open(prefix + ".sam") as f:
        for line in f:
            if line.startswith("@"):
                continue
            x = line.split("\t", 10)
            name, ref, pos, cigar = x[0], x[2], int(x[3]) - 1, x[5]
            length = sum(int(n) for n, op in re.findall(r"(\d+)([MDN=X])", cigar))
            a, b, _ = span.get(name, (pos, pos, ref))
            span[name] = (min(a, pos), max(b, pos + length), ref)
    with open(prefix + "1.fq") as f1, open(prefix + "2.fq") as f2, \
            gzip.open(out1, "wt", compresslevel=3) as o1, gzip.open(out2, "wt", compresslevel=3) as o2:
        for (h1, s1, _, q1), (h2, s2, _, q2) in zip(fastq(f1), fastq(f2)):
            name = h1[1:].split()[0]
            if name.endswith("/1"):
                name = name[:-2]
            a, b, ref = span[name]
            n = f"art{name.rsplit('-', 1)[-1]}_{ref}_{a}_{b}_+"
            o1.write(f"@{n}/1\n{s1}\n+\n{q1}\n")
            o2.write(f"@{n}/2\n{s2}\n+\n{q2}\n")


if __name__ == "__main__":
    if len(sys.argv) == 5 and sys.argv[1] == "maf":
        do_maf(*sys.argv[2:])
    elif len(sys.argv) == 5 and sys.argv[1] == "art":
        do_art(*sys.argv[2:])
    else:
        sys.exit(__doc__)
