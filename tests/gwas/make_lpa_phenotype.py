#!/usr/bin/env python3
"""Build a diploid cohort + literature-based Lp(a) phenotype on the REAL LPA graph.

Reads panphorte's <prefix>.panphorte.copies.tsv (KIV-2 copies per real haplotype path), pairs
real haplotypes into synthetic diploid samples (cosigt-style), and assigns an Lp(a) phenotype
using the published INVERSE relationship between KIV-2 copy number and plasma Lp(a) (more KIV-2
repeats -> lower Lp(a)). The phenotype values are synthetic but the effect direction/size is
literature-plausible; the point is to check the GWAS recovers the KIV-2 locus on real topology.

  Lp(a) = BASE - SLOPE * (CN_A + CN_B) + N(0, SIGMA)        case = Lp(a) > median

Usage: make_lpa_phenotype.py <copies.tsv> <out_dir> [n_samples] [seed]
Writes <out_dir>/samples.tsv (sample, hap1, hap2) and <out_dir>/phenotypes.tsv (per sample).
"""
import os
import random
import statistics
import sys

copies_tsv = sys.argv[1]
out_dir = sys.argv[2]
N_SAMPLES = int(sys.argv[3]) if len(sys.argv) > 3 else 200
SEED = int(sys.argv[4]) if len(sys.argv) > 4 else 23
BASE, SLOPE, SIGMA = 150.0, 2.2, 14.0   # Lp(a)-like units; inverse KIV-2 effect

random.seed(SEED)
os.makedirs(out_dir, exist_ok=True)

# parse copies.tsv -> per-haplotype CN for the KIV-2 bubble (the one with the largest unit_bp)
rows = []
with open(copies_tsv) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    ci = {c: i for i, c in enumerate(hdr)}
    for line in fh:
        f = line.rstrip("\n").split("\t")
        rows.append((f[ci["path_name"]], int(f[ci["bubble_id"]]), int(f[ci["copies"]]),
                     int(f[ci["unit_bp"]])))
if not rows:
    sys.exit("no rows in copies.tsv (run panphorte --min-similarity 0.90 first)")
# KIV-2 bubble = the bubble id with the largest unit_bp
kiv2_bubble = max({(b, ub) for _, b, _, ub in rows}, key=lambda x: x[1])[0]
cn = {p: c for p, b, c, _ in rows if b == kiv2_bubble}
haps = sorted(cn)
print(f"KIV-2 bubble={kiv2_bubble}; {len(haps)} haplotypes; CN range {min(cn.values())}..{max(cn.values())}")

# diploid cohort: random pairs of real haplotypes. Covariates (Age, Sex, PC1-3) are synthetic; the
# phenotype carries the KIV-2 effect PLUS small covariate effects (so covariate adjustment matters),
# and ~5% of phenotypes are set to NA to exercise associate's missing-data filtering.
AGE_SLOPE, SEX_EFFECT, NA_FRAC = 0.25, 6.0, 0.05
samples, prows = [], []
for i in range(N_SAMPLES):
    a, b = random.choice(haps), random.choice(haps)
    dosage = cn[a] + cn[b]
    age = random.randint(35, 80)
    sex = random.randint(0, 1)
    pcs = [round(random.gauss(0, 1), 4) for _ in range(3)]
    lpa = (BASE - SLOPE * dosage + AGE_SLOPE * (age - 55) + SEX_EFFECT * sex
           + random.gauss(0, SIGMA))
    sid = f"ind{i:03d}"
    samples.append((sid, a, b))
    prows.append([sid, dosage, lpa, age, sex, pcs])

median_lpa = statistics.median(r[2] for r in prows)
for r in prows:
    r.append(1 if r[2] > median_lpa else 0)   # r[6] = case
na = {i for i in range(N_SAMPLES) if random.random() < NA_FRAC}  # samples with NA phenotype

with open(os.path.join(out_dir, "samples.tsv"), "w") as f:
    f.write("sample\thaplotype_1\thaplotype_2\n")
    for sid, a, b in samples:
        f.write(f"{sid}\t{a}\t{b}\n")
# legacy table (kiv2 truth dosage), kept for the copy-number sanity check
with open(os.path.join(out_dir, "phenotypes.tsv"), "w") as f:
    f.write("sample\tkiv2_dosage\tlpa_continuous\tcase_binary\n")
    for r in prows:
        f.write(f"{r[0]}\t{r[1]}\t{r[2]:.4f}\t{r[6]}\n")
# associate input: one phenotype column + covariates; some phenotypes NA (covariates kept).
def write_pheno(path, pheno_idx, fmt):
    with open(path, "w") as f:
        f.write("sample\tphenotype\tAge\tSex\tPC1\tPC2\tPC3\n")
        for i, r in enumerate(prows):
            val = "NA" if i in na else fmt(r[pheno_idx])
            f.write(f"{r[0]}\t{val}\t{r[3]}\t{r[4]}\t{r[5][0]}\t{r[5][1]}\t{r[5][2]}\n")
write_pheno(os.path.join(out_dir, "pheno.quant.tsv"), 2, lambda v: f"{v:.4f}")
write_pheno(os.path.join(out_dir, "pheno.binary.tsv"), 6, lambda v: str(v))
print(f"wrote {out_dir}/samples.tsv + phenotypes.tsv + pheno.{{quant,binary}}.tsv  "
      f"({len(samples)} diploid samples; dosage {min(r[1] for r in prows)}..{max(r[1] for r in prows)}; "
      f"cases {sum(r[6] for r in prows)}/{len(prows)}; {len(na)} NA phenotypes)")
