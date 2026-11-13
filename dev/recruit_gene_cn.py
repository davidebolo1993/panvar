#!/usr/bin/env python3
"""Per-gene copy number from a sample's certain reads (recruit --reads-out --certain-only), against the
truth tables of tests/real_data.

    recruit_gene_cn.py --gfa <locus graph> --sample <name> --reads <fq.gz>... [--genes <BED>] [--repeats <TSV>]

A unit is a gene of --genes (molecule, gene, start, end, strand: copies on each annotated haplotype, with
positions) or the repeat of --repeats (hapl, copies). Its k-mers: canonical 31-mers inside its copies
(any k-mer, for a repeat) whose count in at least 90% of the graph's annotated haplotypes equals that
haplotype's copies of the unit, so one k-mer reads one copy, and paralogs, which would add to the count,
are left out. Backbone k-mers: once in every annotated haplotype. A unit's copy number in the diploid
sample is 2 x (median depth of its k-mers in the reads) / (median depth of the backbone k-mers). The
truth: the sample's copies in the same table (its haplotypes need not be in the graph).
"""
import argparse
import collections
import gzip
import re
import statistics
import sys

K = 31
COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def kmers(seq):
    seq = seq.upper()
    rc = seq.translate(COMP)[::-1]
    n = len(seq) - K + 1
    out = []
    for i in range(max(0, n)):
        f = seq[i:i + K]
        if "N" in f:
            continue
        r = rc[n - 1 - i:n - 1 - i + K]
        out.append(f if f < r else r)
    return out


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_gfa(path):
    """Each path's sequence by name: P lines by their name, W lines as sample#haplotype#contig:start-end."""
    seqs, paths = {}, {}
    with opener(path) as f:
        for line in f:
            x = line.rstrip("\n").split("\t")
            if x[0] == "S":
                seqs[x[1]] = x[2]
            elif x[0] == "P":
                paths[x[1]] = [(st[:-1], st[-1] == "-") for st in x[2].split(",")]
            elif x[0] == "W":
                walk = re.findall(r"([<>])([^<>]+)", x[6])
                name = f"{x[1]}#{x[2]}#{x[3]}" + (f":{x[4]}-{x[5]}" if x[4] != "*" else "")
                paths[name] = [(node, d == "<") for d, node in walk]
    out = {}
    for name, steps in paths.items():
        out[name] = "".join(seqs[n].translate(COMP)[::-1] if rev else seqs[n] for n, rev in steps)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gfa", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--reads", nargs="+", required=True)
    ap.add_argument("--genes")
    ap.add_argument("--repeats")
    ap.add_argument("--agree", type=float, default=0.9)
    a = ap.parse_args()
    if not a.genes and not a.repeats:
        sys.exit("give --genes or --repeats")

    hap = read_gfa(a.gfa)
    # truth per haplotype: unit -> {haplotype: copies}; positions of each unit's copies
    copies = collections.defaultdict(lambda: collections.defaultdict(int))
    where = collections.defaultdict(list)   # unit -> [(haplotype, start, end)]
    annotated = set()
    if a.genes:
        with open(a.genes) as f:
            next(f)
            for line in f:
                mol, gene, s, e = line.rstrip("\n").split("\t")[:4]
                copies[gene][mol] += 1
                where[gene].append((mol, int(s), int(e)))
                annotated.add(mol)
    repeat_units = []
    if a.repeats:
        with open(a.repeats) as f:
            next(f)
            for line in f:
                mol, n = line.split("\t")[:2]
                copies["repeat"][mol] = int(n)
                annotated.add(mol)
        repeat_units.append("repeat")
    paths = sorted(p for p in hap if p in annotated)
    if not paths:
        sys.exit("no graph path is in the truth table")
    truth = {u: sum(n for m, n in c.items() if m.startswith(a.sample + "#")) for u, c in copies.items()}

    # candidates: k-mers inside a unit's copies (any k-mer of the haplotypes for a repeat)
    cand = collections.defaultdict(set)   # k-mer -> units
    counts = {}
    for p in paths:
        counts[p] = collections.Counter(kmers(hap[p]))
        for u, ws in where.items():
            for mol, s, e in ws:
                if mol == p:
                    for k in kmers(hap[p][s:e]):
                        cand[k].add(u)
    # agreement: haplotypes where the k-mer's count equals the unit's copies (absent counts as 0)
    zero = {u: sum(1 for p in paths if copies[u].get(p, 0) == 0) for u in copies}
    agree = collections.Counter()
    for p in paths:
        c = counts[p]
        for k, n in c.items():
            for u in list(cand.get(k, ())) + repeat_units:
                t = copies[u].get(p, 0)
                if n == t:
                    agree[(k, u)] += 1
                elif t == 0:
                    agree[(k, u)] -= 1   # present where the unit is absent: undo the absent-counts-as-agreeing below
    need = a.agree * len(paths)
    unit_kmers = collections.defaultdict(set)
    for (k, u), n in agree.items():
        if n + zero[u] >= need and any(copies[u].get(p, 0) > 0 for p in paths):
            unit_kmers[u].add(k)
    backbone = None
    for p in paths:
        once = {k for k, n in counts[p].items() if n == 1}
        backbone = once if backbone is None else backbone & once
    wanted = set(backbone)
    for s in unit_kmers.values():
        wanted |= s

    depth = collections.Counter()
    for path in a.reads:
        with opener(path) as f:
            for i, line in enumerate(f):
                if i % 4 == 1:
                    depth.update(k for k in kmers(line.strip()) if k in wanted)
    bb = statistics.median(depth[k] for k in backbone) if backbone else 0
    print(f"# {a.sample}: {len(paths)} annotated haplotypes in the graph; {len(backbone)} backbone k-mers, median depth {bb}", file=sys.stderr)
    print("\t".join(["sample", "unit", "kmers", "median_depth", "copy_number", "called", "truth", "match"]))
    for u in sorted(copies):
        ks = unit_kmers.get(u, set())
        if not ks or bb == 0:
            print("\t".join([a.sample, u, str(len(ks)), ".", ".", ".", str(truth[u]), "."]))
            continue
        med = statistics.median(depth[k] for k in ks)
        cn = 2 * med / bb
        called = int(cn + 0.5)
        print("\t".join([a.sample, u, str(len(ks)), f"{med:g}", f"{cn:.2f}", str(called), str(truth[u]), "yes" if called == truth[u] else "no"]))


if __name__ == "__main__":
    main()
