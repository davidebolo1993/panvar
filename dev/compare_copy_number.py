#!/usr/bin/env python3
"""compare_copy_number.py - one comparison of call's copy numbers with a truth set, haplotype by haplotype.

Usage:
  compare_copy_number.py --truth <file> --vcf <prefix.region.vcf> [--genes G1,G2] [--label L] [-o table.tsv]
  compare_copy_number.py --truth <file> --gene-table <prefix.dup_gene_cn.tsv> --genes G [--site ID] [--label L]
                         [-o table.tsv]

Truth, one of:
  a pangene BED (a header with molecule and gene columns): a haplotype's expected copy number is its rows of
  the --genes, 0 when it has none;
  a table with a header and one row per haplotype, its name then its copy number (e.g. lpa.repeats.tsv).

Called:
  --vcf         the copy numbers of the CNV record naming any of the --genes (the largest, by SVLEN or span),
                or of the largest CNV record when none names them: the module total
  --gene-table  the gene's own copy number from the per-gene table (one gene in --genes); a gene with rows at
                several sites needs --site, the bubble_id whose rows to take, as a copy number per site does
                not say how the sites combine

Prints the haplotypes compared, how many match exactly and the correlation. -o writes one row per haplotype
(label sample truth_cn called_cn is_reference level hybrid) for plot_cn_correlation.R.
"""
import argparse
import collections
import re
import sys

REFERENCE_NAMES = ("grch38", "chm13")


def match_truth(cn, exp):
    """Called haplotypes keyed by their truth name: the name itself when the truth has it, else without an _inv
    suffix (odgi flip adds one to a path it reverse-complemented). Two called paths for one truth name stop."""
    out, src = {}, {}
    for s, v in cn.items():
        t = re.sub(r"_inv$", "", s)
        if s in exp or t not in exp:
            t = s
        if t in out:
            sys.exit(f"called haplotypes {src[t]} and {s} both match truth haplotype {t}")
        out[t], src[t] = v, s
    return out


def load_truth(path, genes):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [line.rstrip("\n").split("\t") for line in fh if line.strip()]
    if "gene" in header:
        mi, gi = header.index("molecule") if "molecule" in header else 0, header.index("gene")
        if not genes:
            sys.exit("a pangene BED needs --genes")
        exp = {r[mi]: 0 for r in rows}
        for r in rows:
            if r[gi] in genes:
                exp[r[mi]] += 1
        return exp
    return {r[0]: int(r[1]) for r in rows if len(r) > 1 and r[1].lstrip("-").isdigit()}


def called_from_vcf(path, genes):
    samples, best, best_named = [], None, None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
            if info.get("SVTYPE") != "CNV":
                continue
            # a module record has no SVLEN: its carriers gain and lose copies; it is sized by its span
            size = abs(int(info["SVLEN"])) if "SVLEN" in info else int(info.get("END", f[1])) - int(f[1])
            if best is None or size > best[0]:
                best = (size, f)
            if set(info.get("GENES", "").split(",")) & genes and (best_named is None or size > best_named[0]):
                best_named = (size, f)
    rec = best_named or best
    if rec is None:
        sys.exit(f"no CNV record in {path}")
    cn = {}
    for s, g in zip(samples, rec[1][9:]):
        p = g.split(":")
        if len(p) > 1 and p[1].isdigit():
            cn[s] = int(p[1])
    return cn


def called_from_gene_table(path, gene, site):
    """The gene's copy number per haplotype, and which haplotypes carry a hybrid copy of it."""
    cn, hybrid, sites = {}, {}, set()
    with open(path) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        bi, si, gi = h.index("bubble_id"), h.index("sample"), h.index("genes")
        ci, ri = h.index("cn"), h.index("reliable")
        hi = h.index("hybrid_copies") if "hybrid_copies" in h else None
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[gi] != gene:
                continue
            sites.add(f[bi])
            if site and f[bi] != site:
                continue
            if f[ri] == "1" and f[ci].isdigit():
                cn[f[si]] = int(f[ci])
                hybrid[f[si]] = int(hi is not None and f[hi] not in (".", "0", ""))
    # One row per gene and haplotype at a site: rows at several sites are copy numbers of different modules.
    if not site and len(sites) > 1:
        sys.exit(f"{gene} has rows at sites {', '.join(sorted(sites, key=lambda b: (len(b), b)))} of {path}: give --site")
    if site and site not in sites:
        sys.exit(f"{gene} has no row at site {site} of {path}")
    return cn, hybrid


def pearson(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sx = sum((a - mx) ** 2 for a in x) ** 0.5
    sy = sum((b - my) ** 2 for b in y) ** 0.5
    return sxy / (sx * sy) if sx > 0 and sy > 0 else None


def main():
    ap = argparse.ArgumentParser(description="Call's copy numbers against a truth set, haplotype by haplotype.",
                                 formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=32, width=100))
    ap.add_argument("--truth", required=True, metavar="FILE", help="pangene BED, or haplotype and copy number")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--vcf", metavar="VCF", help="call's region VCF: the module total")
    src.add_argument("--gene-table", metavar="TSV", help="call's dup_gene_cn.tsv: one gene's copies")
    ap.add_argument("--genes", default="", metavar="G1,G2", help="genes counted in the truth")
    ap.add_argument("--site", default="", metavar="ID", help="--gene-table: the site, when the gene has several")
    ap.add_argument("--label", default="", help="name in the output")
    ap.add_argument("-o", "--out", default="", metavar="TSV", help="one row per haplotype, for plot_cn_correlation.R")
    args = ap.parse_args()

    genes = {g for g in args.genes.split(",") if g}
    exp = load_truth(args.truth, genes)
    if args.gene_table:
        if len(genes) != 1:
            sys.exit("--gene-table compares one gene: give it in --genes")
        (cn, hybrid), level = called_from_gene_table(args.gene_table, next(iter(genes)), args.site), "gene"
    else:
        if args.site:
            sys.exit("--site goes with --gene-table")
        cn, level = called_from_vcf(args.vcf, genes), "locus"
        hybrid = {}
    cn, hybrid = match_truth(cn, exp), match_truth(hybrid, exp)
    label = args.label or ",".join(sorted(genes)) or "copy number"

    shared = sorted(s for s in cn if s in exp)
    if not shared:
        sys.exit(f"{label}: no haplotype both called and in the truth (do the names match?)")
    x, y = [exp[s] for s in shared], [cn[s] for s in shared]
    exact = sum(a == b for a, b in zip(x, y))
    r = pearson(x, y)
    print(f"{label}: {len(shared)} haplotypes compared ({len(cn)} called, {len(exp)} in the truth)")
    print(f"  exact: {exact}/{len(shared)} ({100 * exact / len(shared):.1f}%)")
    print(f"  r: {'n/a (constant)' if r is None else f'{r:.3f}'}")
    off = collections.Counter((a, b) for a, b in zip(x, y) if a != b)
    if off:
        print("  differing (truth -> called: haplotypes): " +
              ", ".join(f"{a} -> {b}: {n}" for (a, b), n in sorted(off.items())))
    hy = [s for s in shared if hybrid.get(s)]
    if hy:
        print(f"  with a hybrid copy: {len(hy)}, exact {sum(exp[s] == cn[s] for s in hy)}")
    if args.out:
        with open(args.out, "w") as out:
            out.write("label\tsample\ttruth_cn\tcalled_cn\tis_reference\tlevel\thybrid\n")
            for s in shared:
                ref = int(any(k in s.lower() for k in REFERENCE_NAMES))
                out.write(f"{label}\t{s}\t{exp[s]}\t{cn[s]}\t{ref}\t{level}\t{hybrid.get(s, 0)}\n")


if __name__ == "__main__":
    main()
