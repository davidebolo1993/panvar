#!/usr/bin/env python3
"""Build a structured diploid positive-control cohort on the real LPA graph.

Reads panphorte's <prefix>.panphorte.copies.tsv (KIV-2 copies per real haplotype path), partitions the
real haplotypes into a few simulated SUBPOPULATIONS, and samples diploid individuals by
drawing a subpopulation and then two haplotypes from it. This injects realistic POPULATION STRUCTURE:
subpopulations differ both in KIV-2 allele frequency and in a baseline Lp(a) offset that is NOT caused by
KIV-2 (an ancestry/environment confounder). A naive GWAS that ignores ancestry is therefore inflated
(genomic inflation lambda > 1); the supplied ancestry covariates are designed to remove it. The LMM
run is a separate negative control here because its locus-panel relationship matrix does not capture
the simulated structure. That contrast is the point of the example.

Phenotype (synthetic values, literature-plausible shape and effect direction):

  log10 Lp(a) = BASE - SLOPE*(CN_A+CN_B) + subpop_offset + age/sex effects + N(0, SIGMA)

so plasma Lp(a) is log-normal (right-skewed, ~0.3-300 mg/dL, median ~10-12; Coassin & Kronenberg 2022)
with the published INVERSE KIV-2 effect (more KIV-2 repeats -> smaller isoform -> lower Lp(a)). The binary
trait is the clinical high-risk cut Lp(a) > 50 mg/dL. The quantitative phenotype emitted is log10 Lp(a)
(so a linear model sees ~normal residuals). Biological caveat: in reality the SMALLER of the two isoforms
dominates Lp(a); we model the additive summed dosage for a clean, recoverable demo.

Outputs in <out_dir>:
  samples.tsv             sample <tab> haplotype_1 <tab> haplotype_2     (cosigt-style)
  pheno.quant.tsv         sample, phenotype(=log10 Lp(a)), Age, Sex, PC1..PC10   (with ~5% NA)
  pheno.binary.tsv        sample, phenotype(=high-risk 0/1), Age, Sex, PC1..PC10
  pheno.quant.nopc.tsv    same but WITHOUT the PCs (the naive/uncorrected analysis)
  pheno.binary.nopc.tsv
  phenotypes.tsv          legacy truth table (kiv2 dosage, raw Lp(a), subpop) for sanity checks
  kinship.tsv             n x n locus-panel relationship matrix from haplotype sharing
                          (only if --kinship-out and n is small enough; not a genome-wide GRM)

Usage:
  make_lpa_phenotype.py <copies.tsv> <out_dir> [--n N] [--seed S] [--subpops K] [--kinship-out PATH]
A bare positional 3rd arg is still accepted as N for backward compatibility.
"""
import argparse
import math
import os
import random
import statistics
import sys

# ---- literature-grounded constants (see docs/gwas/example.md "Data & literature resources") ----
BASE_LOG10 = 1.08      # 10**1.08 ~ 12 mg/dL median Lp(a) at the cohort-mean KIV-2 dosage
SLOPE_LOG10 = 0.022    # per summed-copy decrease in log10 Lp(a) (inverse KIV-2 effect)
SIGMA_LOG10 = 0.42     # residual sd on log10 scale (KIV-2 leaves substantial unexplained variance)
AGE_SLOPE_LOG10 = 0.0015   # mild age effect per year
SEX_EFFECT_LOG10 = 0.05    # mild sex effect
HIGH_RISK_MGDL = 50.0      # clinical high-risk Lp(a) threshold -> binary case
NA_FRAC = 0.05             # fraction of phenotypes set to NA (to exercise the missing-data filter)
N_PCS = 10                 # ancestry PCs written as covariate columns (top 2 carry the subpop signal)


def parse_args(argv):
    p = argparse.ArgumentParser(description="structured LPA cohort + Lp(a) phenotype")
    p.add_argument("copies_tsv")
    p.add_argument("out_dir")
    p.add_argument("n_pos", nargs="?", type=int, default=None, help="(legacy positional N)")
    p.add_argument("seed_pos", nargs="?", type=int, default=None, help="(legacy positional seed)")
    p.add_argument("--n", type=int, default=None, help="number of diploid individuals (default 200)")
    p.add_argument("--seed", type=int, default=None)
    p.add_argument("--subpops", type=int, default=3, help="number of subpopulations (default 3)")
    p.add_argument("--kinship-out", default=None,
                   help="write an n x n LPA-panel relationship matrix here (not a genome-wide GRM)")
    p.add_argument("--kinship-max-n", type=int, default=3000, help="skip the relationship matrix above this n")
    p.add_argument("--sim-markers", type=int, default=0,
                   help="also emit a synthetic many-null marker panel (causal KIV-2 + N "
                        "subpop-stratified null SNPs) for the structure-correction demo (needs numpy)")
    p.add_argument("--strat-sd", type=float, default=0.15,
                   help="per-subpop allele-frequency sd for the stratified null markers (default 0.15)")
    a = p.parse_args(argv)
    a.n = a.n if a.n is not None else (a.n_pos if a.n_pos is not None else 200)
    a.seed = a.seed if a.seed is not None else (a.seed_pos if a.seed_pos is not None else 23)
    return a


def main(argv):
    a = parse_args(argv)
    random.seed(a.seed)
    os.makedirs(a.out_dir, exist_ok=True)

    # parse copies.tsv -> per-haplotype CN for the KIV-2 bubble (the one with the largest unit_bp)
    # panphorte renamed this column to input_bubble_id, to say plainly that it is the id the site had
    # BEFORE folding: re-snarling the normalized graph reassigns ids, so it is not the bubble id the
    # downstream VCF uses. That is fine here -- this only needs a consistent key to group a haplotype's
    # copy rows by site -- but it must not be joined against a call record's BUBBLE_ID.
    rows = []
    with open(a.copies_tsv) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ci = {c: i for i, c in enumerate(hdr)}
        bid_col = next((c for c in ("input_bubble_id", "bubble_id") if c in ci), None)
        missing = [c for c in ("path_name", "copies", "unit_bp") if c not in ci]
        if bid_col is None or missing:
            sys.exit(f"{a.copies_tsv}: missing column(s) "
                     f"{missing + ([] if bid_col else ['input_bubble_id'])}; header is {hdr}")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rows.append((f[ci["path_name"]], int(f[ci[bid_col]]), int(f[ci["copies"]]),
                         int(f[ci["unit_bp"]])))
    if not rows:
        sys.exit("no rows in copies.tsv (run panphorte first)")
    kiv2_bubble = max({(b, ub) for _, b, _, ub in rows}, key=lambda x: x[1])[0]
    cn = {p: c for p, b, c, _ in rows if b == kiv2_bubble}
    haps = sorted(cn)
    print(f"KIV-2 bubble={kiv2_bubble}; {len(haps)} haplotypes; CN range "
          f"{min(cn.values())}..{max(cn.values())}")

    # subpopulations: each tilts the KIV-2 frequency and carries an Lp(a) offset independent of KIV-2 (the confounder)
    K = max(1, a.subpops)
    cn_lo, cn_hi = min(cn.values()), max(cn.values())
    cn_mid = 0.5 * (cn_lo + cn_hi)
    # per-subpop tilt of the CN sampling weights: subpop 0 favors LOW CN (-> high Lp(a)), last favors HIGH
    tilts = [(-1.0 + 2.0 * k / max(1, K - 1)) for k in range(K)] if K > 1 else [0.0]
    # per-subpop baseline log10 Lp(a) offset (the confounder), ~+-0.4 log units -> inflates a naive scan
    offsets = [(-0.4 + 0.8 * k / max(1, K - 1)) for k in range(K)] if K > 1 else [0.0]
    subpop_frac = [0.45, 0.35, 0.20] if K == 3 else [1.0 / K] * K

    def hap_weights(tilt):
        # exponential tilt over the CN range; tilt>0 favors high CN, tilt<0 favors low CN
        return [math.exp(tilt * (cn[h] - cn_mid) / max(1.0, cn_hi - cn_mid)) for h in haps]

    pop_weights = [hap_weights(t) for t in tilts]

    # pass 1: draw subpop, haplotypes, covariates (so we can center the KIV-2 effect on the cohort mean)
    samples, draws = [], []
    for i in range(a.n):
        sp = random.choices(range(K), weights=subpop_frac, k=1)[0]
        w = pop_weights[sp]
        h1 = random.choices(haps, weights=w, k=1)[0]
        h2 = random.choices(haps, weights=w, k=1)[0]
        age = max(35, min(85, int(round(random.gauss(55, 12)))))   # Moli-sani-like adult range
        sex = random.randint(0, 1)
        # PCs capture subpop ancestry: the top 2 are the subpop centroid on a circle (the real signal),
        # PC3..PC10 are individual noise (as in a real scree plot, only the leading PCs matter).
        ang = 2.0 * math.pi * sp / K
        pc = tuple([round(2.0 * math.cos(ang) + random.gauss(0, 0.5), 4),
                    round(2.0 * math.sin(ang) + random.gauss(0, 0.5), 4)]
                   + [round(random.gauss(0, 1), 4) for _ in range(N_PCS - 2)])
        sid = f"ind{i:05d}"
        samples.append((sid, h1, h2))
        draws.append(dict(sid=sid, sp=sp, dosage=cn[h1] + cn[h2], age=age, sex=sex, pc=pc))
    mean_dosage = statistics.mean(d["dosage"] for d in draws)  # center so median Lp(a) ~ 10**BASE

    # pass 2: phenotype = log-normal Lp(a), inverse KIV-2 effect (centered) + subpop confounder + noise
    prows = []
    for d in draws:
        log10lpa = (BASE_LOG10 - SLOPE_LOG10 * (d["dosage"] - mean_dosage) + offsets[d["sp"]]
                    + AGE_SLOPE_LOG10 * (d["age"] - 55) + SEX_EFFECT_LOG10 * d["sex"]
                    + random.gauss(0, SIGMA_LOG10))
        lpa = min(400.0, max(0.1, 10.0 ** log10lpa))            # mg/dL, log-normal, clipped
        prows.append(dict(sid=d["sid"], sp=d["sp"], dosage=d["dosage"], age=d["age"], sex=d["sex"],
                          pc=d["pc"], log10lpa=round(math.log10(lpa), 4),
                          lpa=round(lpa, 3), case=1 if lpa > HIGH_RISK_MGDL else 0))

    na = {i for i in range(a.n) if random.random() < NA_FRAC}
    ncase = sum(r["case"] for r in prows)
    print(f"{a.n} individuals across {K} subpops; KIV-2 dosage "
          f"{min(r['dosage'] for r in prows)}..{max(r['dosage'] for r in prows)}; "
          f"Lp(a) median {statistics.median(r['lpa'] for r in prows):.1f} mg/dL; "
          f"high-risk cases {ncase}/{a.n}; {len(na)} NA phenotypes")

    # ---- write tables ----
    with open(os.path.join(a.out_dir, "samples.tsv"), "w") as f:
        f.write("sample\thaplotype_1\thaplotype_2\n")
        for sid, h1, h2 in samples:
            f.write(f"{sid}\t{h1}\t{h2}\n")

    with open(os.path.join(a.out_dir, "phenotypes.tsv"), "w") as f:
        f.write("sample\tsubpop\tkiv2_dosage\tlpa_mgdl\tlog10_lpa\tcase_highrisk\n")
        for r in prows:
            f.write(f"{r['sid']}\t{r['sp']}\t{r['dosage']}\t{r['lpa']}\t{r['log10lpa']}\t{r['case']}\n")

    def write_pheno(path, key, fmt, with_pc):
        cols = ["sample", "phenotype", "Age", "Sex"] + ([f"PC{j}" for j in range(1, N_PCS + 1)] if with_pc else [])
        with open(path, "w") as f:
            f.write("\t".join(cols) + "\n")
            for i, r in enumerate(prows):
                val = "NA" if i in na else fmt(r[key])
                base = [r["sid"], val, str(r["age"]), str(r["sex"])]
                if with_pc:
                    base += [str(x) for x in r["pc"]]
                f.write("\t".join(base) + "\n")

    qfmt = lambda v: f"{v:.4f}"
    bfmt = lambda v: str(v)
    write_pheno(os.path.join(a.out_dir, "pheno.quant.tsv"), "log10lpa", qfmt, True)
    write_pheno(os.path.join(a.out_dir, "pheno.binary.tsv"), "case", bfmt, True)
    write_pheno(os.path.join(a.out_dir, "pheno.quant.nopc.tsv"), "log10lpa", qfmt, False)
    write_pheno(os.path.join(a.out_dir, "pheno.binary.nopc.tsv"), "case", bfmt, False)

    # optional locus-panel relationship matrix from LPA haplotype sharing; this is a test input for
    # the LMM path, not the external genome-wide GRM a real association study should supply
    if a.kinship_out:
        if a.n > a.kinship_max_n:
            print(f"  (skipping kinship: n={a.n} > --kinship-max-n={a.kinship_max_n}; "
                  f"use the PC covariate columns at this scale)")
        else:
            try:
                import numpy as np
            except ImportError:
                print("  (skipping kinship: numpy not available)")
            else:
                hidx = {h: j for j, h in enumerate(haps)}
                Zc = np.zeros((a.n, len(haps)))            # per-sample haplotype dosage (0/1/2)
                for i, (_, h1, h2) in enumerate(samples):
                    Zc[i, hidx[h1]] += 1.0
                    Zc[i, hidx[h2]] += 1.0
                mean = Zc.mean(axis=0)
                sd = Zc.std(axis=0)
                keep = sd > 0
                Zs = (Zc[:, keep] - mean[keep]) / sd[keep]
                Kmat = Zs @ Zs.T / Zs.shape[1]
                with open(a.kinship_out, "w") as f:
                    for i in range(a.n):
                        f.write("\t".join(f"{x:.6g}" for x in Kmat[i]) + "\n")
                print(f"  wrote kinship GRM {a.kinship_out} ({a.n}x{a.n})")

    # optional many-null synthetic panel for the structure-correction demo: real KIV-2 dosage +
    # subpop-stratified null SNPs. The nulls are spuriously associated under a naive scan (inflated
    # lambda) but corrected by the supplied PCs, while KIV-2 survives. The LMM run is retained as a
    # negative control for the locus-panel relationship matrix. See docs/gwas.md.
    if a.sim_markers > 0:
        try:
            import gzip
            import numpy as np
        except ImportError:
            print("  (skipping --sim-markers: numpy not available)")
        else:
            rng = np.random.default_rng(a.seed + 1)
            sp = np.array([r["sp"] for r in prows])
            M = a.sim_markers
            base = rng.uniform(0.1, 0.5, M)                       # baseline allele freq per marker
            dev = rng.normal(0.0, a.strat_sd, (M, K))             # per-subpop deviation (stratification)
            half = M // 2
            dev[:half] = 0.0                                      # half are unstratified (true nulls)
            freq = np.clip(base[:, None] + dev, 0.02, 0.98)       # M x K allele frequency
            geno = rng.binomial(2, freq[np.arange(M)[:, None], sp[None, :]])  # M x n null dosages (0/1/2)
            causal = np.array([r["dosage"] for r in prows], dtype=float)      # real KIV-2 summed CN

            gpath = os.path.join(a.out_dir, "geno.sim.bimbam.gz")
            apath = os.path.join(a.out_dir, "feature_annot.sim.tsv.gz")
            spath = os.path.join(a.out_dir, "sim.samples.txt")
            with open(spath, "w") as f:
                for r in prows:
                    f.write(r["sid"] + "\n")
            with gzip.open(gpath, "wt") as gf, gzip.open(apath, "wt") as af:
                af.write("feature_id\tlayer\tencoding\tbubbles\tnodes\n")
                # causal KIV-2 marker first (real node 4789, bubble 7)
                gf.write("KIV2, A, B, " + ", ".join(f"{x:g}" for x in causal) + "\n")
                af.write("KIV2\tsim\tcount\t7\t4789\n")
                for m in range(M):
                    fid = f"sim{m:05d}"
                    gf.write(fid + ", A, B, " + ", ".join(str(int(x)) for x in geno[m]) + "\n")
                    # fake node id = m+1 so the Manhattan spreads the nulls along x
                    af.write(f"{fid}\tsim\tpa\t.\t{m + 1}\n")
            print(f"  wrote synthetic structure-demo panel {gpath} ({M} null markers + 1 causal; "
                  f"{half} unstratified) + {os.path.basename(apath)} + {os.path.basename(spath)}")


if __name__ == "__main__":
    main(sys.argv[1:])
