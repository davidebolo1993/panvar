#!/usr/bin/env python3
"""Families from an associate example's cohort, for the comparison with GEMMA.

Pairs the example's samples into parents within each population and gives each couple children
that inherit one haplotype from each parent. A child's dosage in a row is the sum of its two
haplotypes' columns in describe's haplotype matrices (graph and variant). The phenotype is

    y = effect * z(causal dosage) + 0.1 * sex + family background + noise

where the background is shared as the kinship says (parents' values unrelated, a child's the mean
of its parents' plus half their variance), with --h2 of the variance outside the causal row.

Writes into <out>/:
  geno.bimbam.gz   children's dosages, each row rescaled to 0..2 (as describe --scale-dosage),
                   constant rows dropped; geno.raw.bimbam.gz the same rows unscaled
  samples.txt      the children, in matrix order
  phenotype.tsv    sample, phenotype, sex, PC1-PC5 (the parents' mean)
  kinship.tsv      n x n, no header: 1 on the diagonal, 0.5 between siblings
  gemma.pheno.txt, gemma.covar.txt   the same phenotype and covariates (with an intercept) for GEMMA

  make_families.py <example locus dir> <out> [--families 300] [--children 3] [--seed 1]
"""
import argparse
import collections
import gzip
import os
import random


def read_tsv(path):
    with open(path) as fh:
        hdr = next(fh).rstrip("\n").split("\t")
        return [dict(zip(hdr, line.rstrip("\n").split("\t"))) for line in fh if line.strip()]


def read_matrix(folder, substrate):
    """Haplotype columns and rows {id: [dosage or None per column]} of one describe matrix."""
    with gzip.open(os.path.join(folder, "samples.txt.gz"), "rt") as fh:
        cols = [line.strip() for line in fh if line.strip()]
    rows = {}
    with gzip.open(os.path.join(folder, f"bimbam_{substrate}.bimbam.gz"), "rt") as fh:
        for line in fh:
            f = [x.strip() for x in line.split(",")]
            rows[f[0]] = [None if x == "NA" else float(x) for x in f[3:]]
    return cols, rows


def main():
    ap = argparse.ArgumentParser(description="Families from an associate example cohort, for the GEMMA comparison.",
                                 formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=32, width=100))
    ap.add_argument("locus_dir", help="a locus of tests/synthetic_data/associate")
    ap.add_argument("out", help="output folder")
    ap.add_argument("--families", type=int, default=300, metavar="N", help="families (default 300)")
    ap.add_argument("--children", type=int, default=3, metavar="N", help="children per family (default 3)")
    ap.add_argument("--effect", type=float, default=0.4, metavar="F", help="effect per SD of the causal row (default 0.4)")
    ap.add_argument("--h2", type=float, default=0.4, metavar="F", help="family share of the other variance (default 0.4)")
    ap.add_argument("--causal", metavar="ID", help="causal row (default: the copy-number record)")
    ap.add_argument("--seed", type=int, default=1, metavar="N", help="random seed (default 1)")
    a = ap.parse_args()
    rng = random.Random(a.seed)
    d = a.locus_dir
    os.makedirs(a.out, exist_ok=True)

    haps = {r["sample"]: (r["hap1"], r["hap2"]) for r in read_tsv(os.path.join(d, "cohort.samples.tsv"))}
    covs = {r["sample"]: r for r in read_tsv(os.path.join(d, "cohort.phenotype.tsv"))}
    pop = {r["sample"]: r["population"] for r in read_tsv(os.path.join(d, "cohort.truth.tsv"))}

    # Parents: the example's samples paired within each population.
    by_pop = collections.defaultdict(list)
    for s in sorted(haps):
        by_pop[pop[s]].append(s)
    couples = []
    for p in sorted(by_pop):
        ss = by_pop[p][:]
        rng.shuffle(ss)
        couples += [(ss[i], ss[i + 1]) for i in range(0, len(ss) - 1, 2)]
    rng.shuffle(couples)
    if len(couples) < a.families:
        raise SystemExit(f"only {len(couples)} couples in the cohort, asked for {a.families} families")
    couples = couples[:a.families]

    rows, cols = {}, None
    for sub in ("variant", "graph"):
        c, r = read_matrix(os.path.join(d, "haplotype", sub), sub)
        if cols is not None and c != cols:
            raise SystemExit("the haplotype matrices list different columns")
        cols = c
        rows.update(r)
    col = {h: i for i, h in enumerate(cols)}
    causal = a.causal
    if causal is None:
        annot = os.path.join(d, "haplotype", "variant", "feature_annot.variant.tsv.gz")
        with gzip.open(annot, "rt") as fh:
            hdr = next(fh).rstrip("\n").split("\t")
            cnv = [l.split("\t")[0] for l in fh if dict(zip(hdr, l.rstrip("\n").split("\t"))).get("svtype") == "CNV"]
        if len(cnv) != 1:
            raise SystemExit(f"{len(cnv)} copy-number records; name the causal row with --causal")
        causal = cnv[0]

    # Children: one haplotype from each parent; background as the pedigree kinship says.
    sd_bg, sd_noise = a.h2 ** 0.5, (1 - a.h2) ** 0.5
    kids = []   # (name, family, hap columns, sex, PCs, background)
    for fi, (m, f) in enumerate(couples):
        bm, bf = rng.gauss(0, sd_bg), rng.gauss(0, sd_bg)
        pcs = [(float(covs[m][f"PC{k}"]) + float(covs[f][f"PC{k}"])) / 2 for k in range(1, 6)]
        for ci in range(a.children):
            h = (col[rng.choice(haps[m])], col[rng.choice(haps[f])])
            bg = (bm + bf) / 2 + rng.gauss(0, sd_bg / 2 ** 0.5)
            kids.append((f"fam{fi + 1:03d}_{ci + 1}", fi, h, rng.randint(1, 2), pcs, bg))

    def dosage(vals, h):
        x, y = vals[h[0]], vals[h[1]]
        return None if x is None or y is None else x + y

    raw, scaled = [], []
    for rid, vals in rows.items():
        g = [dosage(vals, k[2]) for k in kids]
        seen = [x for x in g if x is not None]
        lo, hi = min(seen), max(seen)
        if hi == lo:
            continue
        raw.append((rid, g))
        scaled.append((rid, [None if x is None else 2 * (x - lo) / (hi - lo) for x in g]))
    cg = dict(raw)[causal]
    if any(x is None for x in cg):
        raise SystemExit(f"the causal row {causal} has missing dosages")
    mu = sum(cg) / len(cg)
    sd = (sum((x - mu) ** 2 for x in cg) / len(cg)) ** 0.5
    y = [a.effect * (x - mu) / sd + 0.1 * k[3] + k[5] + rng.gauss(0, sd_noise) for x, k in zip(cg, kids)]

    def fmt(x):
        return "NA" if x is None else f"{x:.6g}"
    for name, rr in (("geno.raw.bimbam.gz", raw), ("geno.bimbam.gz", scaled)):
        with gzip.open(os.path.join(a.out, name), "wt") as fh:
            for rid, g in rr:
                fh.write(f"{rid}, A, B, " + ", ".join(fmt(x) for x in g) + "\n")
    with open(os.path.join(a.out, "samples.txt"), "w") as fh:
        fh.write("".join(k[0] + "\n" for k in kids))
    with open(os.path.join(a.out, "phenotype.tsv"), "w") as ph, \
            open(os.path.join(a.out, "gemma.pheno.txt"), "w") as gp, \
            open(os.path.join(a.out, "gemma.covar.txt"), "w") as gc:
        ph.write("sample\tphenotype\tsex\t" + "\t".join(f"PC{k}" for k in range(1, 6)) + "\n")
        for k, v in zip(kids, y):
            c = [str(k[3])] + [f"{p:.6g}" for p in k[4]]
            ph.write(f"{k[0]}\t{v:.6g}\t" + "\t".join(c) + "\n")
            gp.write(f"{v:.6g}\n")
            gc.write("1\t" + "\t".join(c) + "\n")
    with open(os.path.join(a.out, "kinship.tsv"), "w") as fh:
        for ki in kids:
            fh.write("\t".join("1" if ki is kj else ("0.5" if ki[1] == kj[1] else "0") for kj in kids) + "\n")
    print(f"{len(kids)} children in {a.families} families; {len(raw)} rows vary; causal {causal}")


if __name__ == "__main__":
    main()
