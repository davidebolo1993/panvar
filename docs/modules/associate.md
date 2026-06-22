# Associate Module (Module 5 — GWAS)

CLI entrypoint:

- `panvar associate`

## What it does

`associate` runs a genome/region-wide association test of a phenotype against the genotypes produced by
`describe`. For each feature (a k-mer, or a node/edge dosage) it fits a generalized linear model

```
phenotype ~ genotype + covariates
```

with the genotype dosage as the tested term, reports the Wald test on the genotype coefficient, applies a
**minor-frequency (MAF) filter on the final genotypes**, and corrects for **multiple testing over the
features actually tested** (not a fixed genome-wide threshold, because panvar tests a *region*, not the
whole genome). It writes a per-feature results table plus a one-line summary with the significance
thresholds, ready for the Manhattan/QQ plotter.

The phenotype type is auto-detected: a **binary** (0/1) phenotype → **logistic** regression (effect
reported as `log_or`); anything else → **linear** regression (effect reported as `beta`). Override with
`--model`. For a quantitative trait you can also fit a **linear mixed model** (`--model lmm`) that absorbs
relatedness / population structure through a kinship matrix, or add the top kinship **PCs as covariates**
(`--pca N`) to the GLM — see [Population structure](#population-structure-pcs-and-the-lmm). Every run also
reports the **genomic-inflation factor λ** so you can see whether the test is well-calibrated.

New to GWAS? Start with the [GWAS primer](../gwas_primer.md), which explains genotypes/dosage, MAF,
multiple testing, kinship, λ, and LMM-vs-PCs in plain terms on this exact LPA example.

## Required inputs

- `--genotypes <bimbam.gz>` — a BIMBAM mean-genotype dosage matrix from `describe`
  (`bimbam_{kmers,graph}.bimbam.gz`, or the per-sample `bimbam_{kmers,graph}.samples.bimbam.gz` when the
  cohort is diploid). Format: `feature_id, A, B, dose_1, dose_2, …`; `NA` = a sample that does not
  traverse the feature's bubble (genuinely missing, dropped per feature), distinct from `0` (traverses but
  reference).
- `--samples <txt[.gz]>` — the sample (column) order, one per line (`describe`'s `bimbam.samples.txt.gz`
  or `bimbam.samples.samples.txt.gz`).
- `--phenotype <tsv>` — `sample <tab> phenotype [<tab> covariate1 …]`, header required. The first column is
  the sample id (must match `--samples`), the second is the phenotype, the rest are covariates (e.g. `Age`,
  `Sex`, `PC1…PCN`). Any cell may be `NA`/`.`/`-9`/empty; a sample with an **NA phenotype or any NA
  covariate is dropped** from the fit.
- `-o, --out-prefix <prefix>`

## Key options

```bash
panvar associate --genotypes <bimbam.gz> --samples <samples.txt> --phenotype <table.tsv> -o <prefix> [options]
```

- `--feature-annot <tsv.gz>` — `describe`'s `feature_annot.tsv.gz`; adds `layer` / `bubbles` / `nodes`
  (graph provenance) to the output so hits are traceable to the graph and, via `node_genes.tsv`, to genes.
- `--node-genes <tsv>` — `call`'s `node_genes.tsv` (from `call --gtf`); adds a **`gene` column** to the
  output by joining the node ids in each feature's `nodes` provenance to gene names. No GTF → omit it (the
  column is still written, valued `.`).
- `--min-maf <X>` — drop features whose **minor (non-modal) genotype frequency** is below `X` (default
  `0.01`). This is the "present in the graph but rare/invariant in the cohort" filter: it is computed from
  the **actual sample matrix**, so it removes features that carry no usable variation in *these* samples.
  It is well-defined for both presence/absence (0/1) and copy-number dosage.
- `--model <auto|linear|logistic|lmm>` — default `auto` (binary → logistic, else linear). `lmm` is a
  linear mixed model for a **quantitative** trait and needs a kinship source (`--kinship`/`--make-kinship`).
- `--kinship <path>` — a precomputed `n × n` genetic-relationship matrix (GRM), rows/cols in `--samples`
  order (whitespace/comma/tab-separated; subset to the used samples automatically). Used by `--model lmm`
  and `--pca`.
- `--make-kinship` — build the GRM from the genotype matrix itself (standardize each feature, `K = ZZᵀ/m`).
  **Caveat:** when the matrix is a single region this GRM is *proximally contaminated* by the very signal
  you test (it can over-correct); it is appropriate for a **genome-wide-like** marker panel. Prefer an
  external `--kinship` built from genome-wide markers when you have one.
- `--pca <N>` — add the top-N kinship **PCs as fixed covariates** to the GLM (cheap structure control; no
  mixed model). Needs `--kinship` or `--make-kinship`. Works for both linear and logistic.
- `-q, --quiet`

## Outputs

- `<prefix>.assoc.tsv` — one row per tested feature, sorted by `p`:
  `feature_id, layer, bubbles, nodes, n, minor_freq, beta|log_or, se, z, p, p_bonf, q_bh, gene`.
  - `n` — samples used for that feature (after NA-genotype/phenotype/covariate removal).
  - `p_bonf` — Bonferroni-adjusted p (`min(1, p · n_tests)`).
  - `q_bh` — Benjamini-Hochberg FDR q-value.
  - `gene` — gene name(s) behind the feature's nodes (`.` unless `--node-genes` is given).
- `<prefix>.summary.tsv` — `model`, `phenotype_type`, `covariates` (incl. any PCs), `pca_covariates`,
  `samples_used`, `features_tested`, `dropped_min_maf`, `dropped_fit`, the
  **`bonferroni_threshold` (`0.05 / n_tests`)**, the counts significant at Bonferroni and FDR<0.05, the
  genomic-inflation **`lambda_gc`**, and (LMM only) the fitted variance ratio `lmm_delta`. The plotter
  reads `features_tested` from here to draw the threshold line.

## Multiple testing — why not 5×10⁻⁸

The genome-wide `5e-8` threshold assumes ~10⁶ effectively independent tests across an entire genome.
`associate` tests only the features in the target region(s) — hundreds to a few thousand — so the honest
correction is over the **number of tests actually performed**:

- **Bonferroni** line `0.05 / n_tests` (a per-feature `p_bonf` is also reported), and
- **Benjamini-Hochberg** `q_bh` (FDR), the usual default for many correlated features.

Because k-mers from one node are highly correlated, `n_tests` over-counts the independent tests, so
Bonferroni here is *conservative*; an effective-number-of-tests / permutation threshold is a planned
refinement. The Manhattan plot draws the nominal (`0.05`) and Bonferroni lines, and highlights FDR<0.05
points, so "before vs after correction" is read directly off the figure.

## Population structure: PCs and the LMM

When a cohort contains **subpopulations** (ancestry groups) that differ both in allele frequencies and in
the trait mean, a naive test reports spurious associations at any marker whose frequency tracks ancestry —
the test is **inflated** (genomic-inflation `lambda_gc` > 1). `associate` offers two standard remedies:

- **PCs as covariates** (`--pca N`, or just include ancestry `PC1…PCN` columns in `--phenotype`): condition
  the regression on the top genetic principal components, so ancestry is regressed out. Cheap, works for
  linear and logistic.
- **Linear mixed model** (`--model lmm`): model relatedness as a random effect with covariance `σ²_g·K`
  (the kinship `K`) plus noise. `associate` uses the fast **EMMAX** approximation — eigendecompose `K`
  once, estimate the variance ratio `δ = σ²_e/σ²_g` once under the null, then test each feature by
  generalized least squares in the rotated space. Reported `lmm_delta` is that ratio.

Read calibration off `lambda_gc` (and the QQ plot): ~1 is well-calibrated, >1 flags residual structure.

> **Scope caveat (important).** `lambda_gc` is only meaningful when the matrix has **many independent null
> markers**. A single pangenome region (e.g. LPA alone) is essentially *one* correlated signal block with no
> nulls, so its λ is not interpretable and structure correction is neither needed nor diagnostic there — the
> operative controls for a single region are the **MAF filter** and **region-wide multiple testing** above.
> PCs/LMM earn their keep once you assemble a **genome-wide** panel of pangenome features (many regions). The
> [gwas_example.md](../gwas_example.md) shows both: the LPA region scan (KIV-2 recovery) *and* a
> genome-wide-like panel where a naive λ≈2–4 collapses to ≈1 under `--pca`/`--model lmm`.

## Genotype encoding

Genotype dosage is used **as-is** for the test (count/dosage), which is the right representation for
copy-number loci (a KIV-2 or CYP2D6 dosage is the signal). Presence/absence is simply the special case of a
0/1 dosage. The `feature_annot.tsv.gz` `encoding` column records `count` vs `pa` per feature for
downstream interpretation; binarization, if wanted, is applied at analysis time.

## Plotting

```bash
Rscript scripts/plot_associate.R \
  --assoc <prefix>.assoc.tsv --summary <prefix>.summary.tsv \
  --out <prefix> --title "my trait"
```

Writes `<prefix>.manhattan.png/pdf` (−log10 p along graph order, with the nominal and region-wide
Bonferroni threshold lines, FDR<0.05 / Bonferroni points highlighted) and `<prefix>.qq.png/pdf` (with the
genomic-inflation λ).

## Example

A full worked example on the real LPA graph — synthetic diploid cohort, covariates, and NA phenotypes —
is in [gwas_example.md](../gwas_example.md) and `tests/gwas/run_lpa_real.sh`. The KIV-2 copy-number locus
is recovered as the top hit with a **negative** effect (more KIV-2 repeats → lower Lp(a)), well past the
region-wide Bonferroni threshold, while the MAF filter drops the cohort-invariant features.
