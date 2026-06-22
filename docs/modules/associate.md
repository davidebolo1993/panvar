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
`--model`.

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
- `--min-maf <X>` — drop features whose **minor (non-modal) genotype frequency** is below `X` (default
  `0.01`). This is the "present in the graph but rare/invariant in the cohort" filter: it is computed from
  the **actual sample matrix**, so it removes features that carry no usable variation in *these* samples.
  It is well-defined for both presence/absence (0/1) and copy-number dosage.
- `--model <auto|linear|logistic>` — default `auto` (binary phenotype → logistic).
- `-q, --quiet`

## Outputs

- `<prefix>.assoc.tsv` — one row per tested feature, sorted by `p`:
  `feature_id, layer, bubbles, nodes, n, minor_freq, beta|log_or, se, z, p, p_bonf, q_bh`.
  - `n` — samples used for that feature (after NA-genotype/phenotype/covariate removal).
  - `p_bonf` — Bonferroni-adjusted p (`min(1, p · n_tests)`).
  - `q_bh` — Benjamini-Hochberg FDR q-value.
- `<prefix>.summary.tsv` — `model`, `phenotype_type`, `covariates`, `samples_used`, `features_tested`,
  `dropped_min_maf`, `dropped_fit`, the **`bonferroni_threshold` (`0.05 / n_tests`)**, and the counts of
  features significant at Bonferroni and at FDR<0.05. The plotter reads `features_tested` from here to draw
  the threshold line.

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
