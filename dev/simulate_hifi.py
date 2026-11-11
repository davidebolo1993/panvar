#!/usr/bin/env python3
"""HiFi-like reads: fragments of normal length (15 kb +- 3 kb by default) drawn uniformly over the
sequences (in proportion to their length), on either strand, with errors at a fixed rate (half
substitutions, a quarter insertions, a quarter deletions). pbsim3 makes HiFi only through PacBio's ccs,
which is not available everywhere, and caps a single pass at 99% accuracy.

    simulate_hifi.py <fasta> <depth> <out.fq.gz> [--seed N] [--error 0.001] [--length 15000] [--sd 3000] [--prefix name]
"""
import argparse
import gzip
import random


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fasta")
    ap.add_argument("depth", type=float)
    ap.add_argument("out")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--error", type=float, default=0.001)
    ap.add_argument("--length", type=int, default=15000)
    ap.add_argument("--sd", type=int, default=3000)
    ap.add_argument("--prefix", default="hifi")
    a = ap.parse_args()
    rng = random.Random(a.seed)
    seqs, name, parts = [], None, []
    for line in open(a.fasta):
        line = line.strip()
        if line.startswith(">"):
            if name is not None:
                seqs.append((name, "".join(parts).upper()))
            name, parts = line[1:].split()[0], []
        else:
            parts.append(line)
    if name is not None:
        seqs.append((name, "".join(parts).upper()))
    total = sum(len(s) for _, s in seqs)
    comp = str.maketrans("ACGTN", "TGCAN")
    n_reads = int(a.depth * total / a.length)
    weights = [len(s) for _, s in seqs]
    with gzip.open(a.out, "wt") as out:
        for i in range(n_reads):
            name, s = rng.choices(seqs, weights=weights)[0]
            length = min(len(s), max(1000, int(rng.gauss(a.length, a.sd))))
            start = rng.randrange(0, len(s) - length + 1)
            frag = s[start:start + length]
            if frag.count("N") > length // 2:
                continue
            strand = "+"
            if rng.random() < 0.5:
                frag = frag.translate(comp)[::-1]
                strand = "-"
            read = []
            for c in frag:
                x = rng.random()
                if x < a.error / 2:
                    read.append(rng.choice([b for b in "ACGT" if b != c]))
                elif x < a.error * 3 / 4:
                    read.append(c + rng.choice("ACGT"))
                elif x < a.error:
                    continue
                else:
                    read.append(c)
            r = "".join(read)
            out.write(f"@{a.prefix}{i}_{name}_{start}_{strand}\n{r}\n+\n{'I' * len(r)}\n")


if __name__ == "__main__":
    main()
