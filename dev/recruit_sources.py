#!/usr/bin/env python3
"""Sources of simulated reads for the recruitment benchmark, one locus at a time.

    recruit_sources.py held-out --gfa <locus.gfa[.gz]> --sample <name[,relative,...]> --out <held_out.gfa>
        the graph without the paths of the sample and its relatives (its index is built from this)
    recruit_sources.py locus --gfa <locus.gfa[.gz]> --sample <name> --out <locus.fa>
        the sample's haplotypes, spelled from the graph (truth: locus)
    recruit_sources.py near --fasta <genome.fa.gz> --region chr:start-end --flank <bp> --out <near.fa>
        the two flanks of the locus (truth: nearby, not locus)
    recruit_sources.py decoy --fasta <genome.fa.gz> --genome-name <name> --report <locus.index.tsv> [--context <bp>] --out <decoy.fa>
        the decoy windows of a locus index found in that genome, widened by --context on each side but
        never into the locus interval (truth: homolog)
    recruit_sources.py boundary --fasta <genome.fa.gz> --region chr:start-end --flank <bp> --out <boundary.fa>
        both edges of the locus, --flank on each side of each (truth: partly locus)
    recruit_sources.py genome --fasta <genome.fa.gz> --mask chr:start-end [...] --out <genome.fa>
        the genome with the loci masked (truth: elsewhere)
Intervals are 0-based, half-open; FASTA intervals are cut with samtools faidx.
"""
import argparse
import gzip
import subprocess
import sys


def lines(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            yield line.rstrip("\n")


def sample_of(fields):
    if fields[0] == "W":
        return fields[1]
    return fields[1].split("#")[0]


def write_fasta(out, name, seq):
    out.write(f">{name}\n")
    for i in range(0, len(seq), 80):
        out.write(seq[i:i + 80] + "\n")


def comp(s):
    return s.translate(str.maketrans("ACGTacgtNn", "TGCAtgcaNn"))[::-1]


def faidx(fasta, regions):
    """Sequences of 0-based half-open regions (contig, start, end) through samtools faidx."""
    if not regions:
        return []
    arg = [f"{c}:{s + 1}-{e}" for c, s, e in regions]
    out = subprocess.run(["samtools", "faidx", fasta] + arg, check=True, capture_output=True, text=True).stdout
    seqs, cur = [], None
    for line in out.splitlines():
        if line.startswith(">"):
            if cur is not None:
                seqs.append("".join(cur))
            cur = []
        else:
            cur.append(line)
    if cur is not None:
        seqs.append("".join(cur))
    return seqs


def region(text):
    contig, span = text.rsplit(":", 1)
    a, b = span.split("-")
    return contig, int(a), int(b)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("what", choices=["held-out", "locus", "near", "decoy", "boundary", "genome"])
    ap.add_argument("--genome-name", default="")
    ap.add_argument("--context", type=int, default=0)
    ap.add_argument("--gfa")
    ap.add_argument("--sample")
    ap.add_argument("--fasta")
    ap.add_argument("--region")
    ap.add_argument("--flank", type=int, default=2_000_000)
    ap.add_argument("--report")
    ap.add_argument("--mask", action="append", default=[])
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    if a.what == "held-out":
        n = 0
        drop = set(a.sample.split(","))
        with open(a.out, "w") as out:
            for line in lines(a.gfa):
                f = line.split("\t")
                if f[0] in ("P", "W") and sample_of(f) in drop:
                    n += 1
                    continue
                out.write(line + "\n")
        if n == 0:
            sys.exit(f"{a.gfa} has no path of {a.sample}")
        print(f"held out {n} paths of {a.sample}", file=sys.stderr)

    elif a.what == "locus":
        seq, walks = {}, []
        for line in lines(a.gfa):
            f = line.split("\t")
            if f[0] == "S":
                seq[f[1]] = f[2]
            elif f[0] == "P" and sample_of(f) == a.sample:
                walks.append((f[1], [(s[:-1], s[-1] == "-") for s in f[2].split(",")]))
            elif f[0] == "W" and f[1] == a.sample:
                steps, w, i = [], f[6], 0
                while i < len(w):
                    j = i + 1
                    while j < len(w) and w[j] not in "<>":
                        j += 1
                    steps.append((w[i + 1:j], w[i] == "<"))
                    i = j
                walks.append((f"{f[1]}#{f[2]}#{f[3]}", steps))
        if not walks:
            sys.exit(f"{a.gfa} has no path of {a.sample}")
        with open(a.out, "w") as out:
            for name, steps in walks:
                s = "".join(comp(seq[n]) if rev else seq[n] for n, rev in steps)
                write_fasta(out, "locus|" + name.replace(":", "_"), s)
        print(f"{len(walks)} haplotypes of {a.sample}", file=sys.stderr)

    elif a.what == "near":
        c, s, e = region(a.region)
        left = (c, max(0, s - a.flank), s)
        right = (c, e, e + a.flank)
        seqs = faidx(a.fasta, [left, right])
        with open(a.out, "w") as out:
            for (cc, ss, ee), sq in zip([left, right], seqs):
                write_fasta(out, f"near|{cc}_{ss}_{ee}", sq)

    elif a.what == "decoy":
        wins = []
        locus = None
        for line in open(a.report):
            f = line.rstrip("\n").split("\t")
            if f[0] == "locus" and f[1] == "interval":
                locus = region(f[2])
            if f[0] == "decoy_window" and f[1] != "genome":
                if a.genome_name and f[1] != a.genome_name:
                    continue
                c, s0, e0 = f[2], int(f[3]), int(f[4])
                s0, e0 = max(0, s0 - a.context), e0 + a.context
                if locus and c == locus[0] and s0 < locus[2] and e0 > locus[1]:   # never into the locus
                    if e0 - locus[2] > locus[1] - s0:
                        s0 = locus[2]
                    else:
                        e0 = locus[1]
                if e0 > s0:
                    wins.append((c, s0, e0))
        with open(a.out, "w") as out:
            for (c, s0, e0), sq in zip(wins, faidx(a.fasta, wins)):
                write_fasta(out, f"decoy|{c}_{s0}_{e0}", sq)
        print(f"{len(wins)} decoy windows", file=sys.stderr)

    elif a.what == "boundary":
        c, s0, e0 = region(a.region)
        edges = [(c, max(0, s0 - a.flank), s0 + a.flank), (c, max(0, e0 - a.flank), e0 + a.flank)]
        with open(a.out, "w") as out:
            for (cc, ss, ee), sq in zip(edges, faidx(a.fasta, edges)):
                write_fasta(out, f"boundary|{cc}_{ss}_{ee}", sq)

    elif a.what == "genome":
        masks = {}
        for m in a.mask:
            c, s, e = region(m)
            masks.setdefault(c, []).append((s, e))
        name, parts = None, []
        with open(a.out, "w") as out:
            def flush():
                if name is None:
                    return
                s = bytearray("".join(parts), "ascii")
                for ms, me in masks.get(name, []):
                    me = min(me, len(s))
                    s[ms:me] = b"N" * (me - ms)
                write_fasta(out, "genome|" + name, s.decode("ascii"))
            for line in lines(a.fasta):
                if line.startswith(">"):
                    flush()
                    name, parts = line[1:].split()[0], []
                else:
                    parts.append(line)
            flush()


if __name__ == "__main__":
    main()
