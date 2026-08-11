# Module `associate`

CLI: `panvar associate`

## What it does

Tests a phenotype against the genotypes from `describe` across a pangenome region. For each unit it fits `phenotype ~ genotype + covariates`, reports a Wald test on the genotype term, applies a minor-allele-frequency (MAF) filter on the cohort genotypes, and corrects for the number of independent tests in the region. Phenotype type is auto-detected (binary or quantitative). The testable unit is chosen by `--unit`. Variant mode tests one genotype per structural-variant (SV) call (the `describe --variant-vcf` export) — the statistically honest unit, since the k-mers, nodes and edges within one variant are correlated rather than independent; correlated nearby variants are then collapsed by linkage-disequilibrium (LD) clumping, so an LD shadow is not counted as a separate hit. Feature mode keeps the fine-grained k-mer/node/edge tests but corrects with an effective number of independent tests (`Meff`, the distinct bubbles), because the raw feature count over-states how many independent tests were run. Both report Benjamini–Hochberg (BH) false-discovery-rate (FDR) control alongside the `Meff`-Bonferroni benchmark and the genomic-inflation `λ`. 
Beyond the threshold, `associate` also tests independence by conditioning: refitting each unit with the top signal(s) added as covariates, so a hit that merely tags a stronger nearby variant is exposed (its `p` collapses). The variant unit runs a forward-stepwise selection of jointly-independent signals (conditional-and-joint, COJO-style); the feature unit conditions on the single top feature with a within-bubble collinearity guard. 
For a quantitative trait it can also add the top kinship principal components (PCs) as covariates (`--pca N`), or fit a linear mixed model (LMM, `--model lmm`) against a kinship matrix, to control population structure.

Algorithm and worked trace: [algorithms/associate.md](../algorithms/associate.md).

## Required inputs

- `--genotypes <bimbam.gz>` — a `describe` BIMBAM matrix from one substrate folder: `<level>/<substrate>/bimbam_<substrate>.bimbam.gz`, where `<level>` is `haplotype` or `sample` (diploid cohort) and `<substrate>` is `kmers`/`graph` (feature unit) or `variant` (variant unit).
- `--samples <txt[.gz]>` — the sample (column) order (`describe`'s matching `*.samples[.samples].txt.gz`).
- `--phenotype <tsv>` — `sample <tab> phenotype [<tab> covariate…]`, header required; cells may be `NA` (a sample with NA phenotype or any NA covariate is dropped).
- `-o, --out-prefix <prefix>`.

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `--feature-annot <tsv.gz>` | the `feature_annot.<substrate>.tsv.gz` from the same folder as `--genotypes`; adds provenance and, for variants, `svtype`/`gene`/`AF`/`AN` | — |
| `--unit <auto\|variant\|feature>` | multiple-testing unit; `auto` picks `variant` when the `feature-annot` is the variant sidecar, else `feature` | `auto` |
| `--ld-r2 <X>` | variant unit: genotype r² above which a variant is an LD shadow of a lead (clumped, not counted in `Meff`) | `0.8` |
| `--min-ac <N>` | variant unit: flag low `AF` when the observed minority-genotype count < N (underpowered/unstable; such a variant also cannot anchor a clump) | `3` |
| `--cojo-p <X>` | variant unit: entry p for forward-stepwise conditional (COJO) signal selection | `0.05/Meff` |
| `--node-genes <tsv>` | `call`'s `node_genes.tsv` (from `--gtf`); adds a `gene` column | — |
| `--min-maf <X>` | drop features whose [minor non-modal frequency](../algorithms/associate.md#terms) < X, on the actual cohort | `0.01` |
| `--model <auto\|linear\|logistic\|lmm>` | `auto` = binary→logistic else linear; `lmm` = mixed model (quantitative; needs a kinship source) | `auto` |
| `--kinship <path>` | external (genome-wide) `n×n` genomic relationship matrix (GRM; rows/cols in `--samples` order) for `--model lmm` / `--pca`; panvar is local and does not build a GRM itself | — |
| `--pca <N>` | add the top-N kinship PCs as covariates to the generalized linear model (GLM; needs `--kinship`); usually you instead pass ancestry PCs as phenotype-table columns | off |
| `-q, --quiet` | less logging | off |

## Outputs

| file | contents |
|------|----------|
| `<prefix>.assoc.tsv` | one row per tested feature, sorted by `p` (columns below) |
| `<prefix>.summary.tsv` | run settings + diagnostics, one `key <tab> value` per line (below) |

`<prefix>.assoc.tsv` columns:

| column | meaning |
|--------|---------|
| `feature_id` | the unit tested — a variant id (variant unit) or a k-mer sequence / node id / edge key (feature unit) |
| `layer` | `variant`, `kmer`, or `graph` |
| `bubbles`, `nodes` | graph provenance (from `--feature-annot`): which bubble(s) / node(s) the unit comes from |
| `n` | number of samples used in this unit's fit |
| `minor_freq` | minor (non-modal) genotype frequency on the cohort (the MAF-filter quantity) |
| `beta` \| `log_or` | effect size on the genotype term — `beta` (linear) or `log_or` = log odds ratio (logistic) |
| `se` | standard error of the effect |
| `z` | Wald statistic, `effect / se` |
| `p` | Wald p-value |
| `p_bonf` | raw Bonferroni-adjusted p, `min(1, p · features_tested)` (over-conservative — kept for reference) |
| `p_bonf_meff` | effective-tests Bonferroni, `min(1, p · Meff)` — the honest correction |
| `q_bh` | Benjamini–Hochberg FDR q-value (the primary control) |
| `af`, `an` | (variant unit) allele frequency / traversing-haplotype count, carried from the VCF (Variant Call Format) |
| `low_af` | (variant unit) `1` when the minority-genotype count < `--min-ac` (underpowered/unstable), else `0`; `.` otherwise |
| `clump`, `is_lead` | (variant unit) LD-clump id and whether this is its lead variant (`1`); `.` in feature mode |
| `gene` | gene name (variant `GENES`, or via `--node-genes`), else `.` |
| `p_conditional` | conditional Wald p (variant: vs the COJO-selected set minus self; feature: vs the top lead); `NA` for a sole signal / the lead / a collinear feature |
| `cond_role` | variant: `signal` (independent) / `shadow`; feature: `lead` / `collinear` (same-event, r²>0.95, not scored) / `conditioned`; `.` otherwise |

`<prefix>.summary.tsv` keys:

| key | meaning |
|-----|---------|
| `model`, `phenotype_type` | model used (`linear`/`logistic`/`lmm`) and detected trait type |
| `covariates`, `pca_covariates` | covariate columns used / PCs added via `--pca` |
| `samples_used` | samples kept after dropping rows with NA phenotype or covariate |
| `features_tested` | units that passed the MAF filter and were tested |
| `unit` | `variant` or `feature` (the multiple-testing unit used) |
| `meff` | effective number of independent tests — LD-clump leads (variant) or distinct bubbles (feature) |
| `independent_variants` \| `distinct_bubbles` | the same `Meff` count, labelled for the active unit |
| `dropped_min_maf`, `dropped_fit` | units dropped by the MAF filter / by a failed model fit |
| `bonferroni_threshold` | raw region-wide threshold `0.05 / features_tested` (reference) |
| `bonferroni_threshold_meff` | the honest threshold `0.05 / Meff` |
| `significant_bonferroni`, `significant_bonferroni_meff`, `significant_fdr05` | counts passing raw Bonferroni / `Meff`-Bonferroni / BH FDR < 0.05 |
| `cojo_independent_signals` | (variant unit) number of jointly-independent signals from forward-stepwise conditioning |
| `lambda_gc` | genomic-inflation factor λ |
| `lmm_delta` | (LMM only) fitted variance ratio δ = σ²ₑ / σ²_g |

## Plotting

```bash
Rscript scripts/plot_associate.R \
  --assoc <prefix>.assoc.tsv \
  --summary <prefix>.summary.tsv \
  --out <prefix>
```

Writes `*.manhattan.{png,pdf}` — stacked panels: before correction (raw −log10 p with nominal and region-wide Bonferroni lines), after correction (Benjamini-Hochberg −log10 q with the q=0.05 line), and — when the `p_conditional`/`cond_role` columns are present — after conditioning (−log10 p_conditional, where shadows collapse below the line and only the conditioning signal(s) stay tall). x = node id (graph) or per-k-mer index ordered by node id (k-mers), with FDR/Bonferroni-significant genes flagged (ggrepel, from the `gene` column when `--node-genes` was passed) — and `*.qq.{png,pdf}` (with `λ`).

Script flags (need `Rscript` + `ggplot2`; `ggrepel` optional, for the gene labels):

- `--assoc <assoc.tsv>` — the association table (required).
- `--out <prefix>` — output prefix for the PNG/PDF files (required).
- `--summary <summary.tsv>` — read `features_tested` for the Bonferroni line (recommended).
- `--title <text>` — plot title.
- `--width` / `--height` / `--dpi` — Manhattan size (inches) and PNG resolution (defaults 10/7/150).


## Example

See the [GWAS example](../gwas.md) for a runnable association run on this locus, and the [LPA walkthrough](../walkthrough.md) for the full pipeline.

## Calibration

A p-value means nothing unless it is uniform when nothing is going on. `tests/associate_null.sh` permutes the phenotype table's sample labels — severing every genotype-phenotype link while leaving the genotype matrix, the missingness pattern and the phenotype/covariate joint distribution untouched — and reports type-I error, `lambda_GC`, and a per-feature uniformity test. Per feature, p-values across permutations are independent, which is what makes the uniformity test valid; pooled across features within one permutation they are not, so the pooled intervals it prints are optimistic and labelled as such.

Measured on the LPA cohort (6000 individuals; 492 cases / 5213 controls for the binary trait), 300 permutations:

| model | lambda_GC | features rejecting uniformity |
|-------|-----------|-------------------------------|
| linear, default `--min-maf` | 1.088 | 2/13 (95% bound is 2) |
| logistic, default `--min-maf` | 1.006 | 0/13 |
| logistic, `--min-maf 0` | 1.010 | 0/20 |

**Binary traits use a Rao score test.** The Wald test divides an estimate by its own standard error, and for a rare variant in an unbalanced case/control study the fit approaches separation: the coefficient grows, its standard error grows faster, and the statistic collapses. Measured here before the change, every feature below minor frequency 0.01 failed uniformity (`lambda_GC` 0.80, worst KS p 3e-12) while every feature above it passed — the default `--min-maf 0.01` was the only thing hiding it. The score test never fits the alternative, so it has no standard error to inflate.

Consequence for reading the output: `z` and `p` are the score test, while `log_or` and `se` remain the Wald maximum-likelihood effect size, so **p is not recoverable from `log_or`/`se`**. This is the same arrangement SAIGE and REGENIE use.

Residual limitation: at `p < 0.001` on very rare features the score test is still mildly anti-conservative (observed 0.0025 against a nominal 0.001). The saddlepoint approximation (SPA) is the standard remedy and is not implemented.
