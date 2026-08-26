# Module `associate`

CLI: `panvar associate`

## What it does

Tests a phenotype against `describe` dosages across a pangenome region. For each unit it fits `phenotype ~ genotype + covariates`, applies a cohort minor-frequency filter and reports a test on the genotype term: a Student-t Wald test for quantitative traits or a Rao score test for binary traits. Phenotype type is auto-detected.

### What is a unit

`--unit` chooses what a single test is about. Everything else in this page follows from that choice.

| unit | one row is | comes from | gets |
|------|-----------|------------|------|
| `variant` | a structural variant | `describe --variant-vcf` | LD clumping, forward-stepwise conditioning, `svtype` / `gene` / `AF` |
| `feature` | one k-mer, node or edge dosage | `describe`'s k-mer or graph substrate | conditioning on the single top feature, with a collinearity guard |

### Which number to read

Three multiple-testing summaries are reported, and they are not interchangeable:

| column | use it for | guarantee |
|--------|-----------|-----------|
| `q_bh` | the headline result | Benjamini–Hochberg FDR, under its usual dependence assumptions. Start here |
| `p_bonf` | a family-wise claim | genuine, and over-conservative in a region where units are correlated |
| `p_bonf_meff` | a regional guide only | none. `Meff` is a heuristic estimate of how many independent tests the region really holds |

Conditional analysis then separates independent signals from correlated shadows: variant mode runs forward-stepwise selection, feature mode conditions on its top feature. Quantitative analyses can additionally use kinship principal components (`--pca`) or an experimental linear mixed model (`--model lmm`), both of which need an external kinship matrix.

Algorithm and worked trace: [algorithms/associate.md](../algorithms/associate.md).

## Required inputs

- `--genotypes <bimbam.gz>` — a `describe` BIMBAM matrix from one substrate folder: `<level>/<substrate>/bimbam_<substrate>.bimbam.gz`, where `<level>` is `haplotype` or `sample` (diploid cohort) and `<substrate>` is `kmers`/`graph` (feature unit) or `variant` (variant unit).
- `--samples <txt[.gz]>` — the sample (column) order (`describe`'s matching `samples.txt.gz`).
- `--phenotype <tsv>` — `sample <tab> phenotype [<tab> covariate…]`, header required; cells may be `NA` (a sample with NA phenotype or any NA covariate is dropped).
- `-o, --out-prefix <prefix>`.

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `--feature-annot <tsv.gz>` | the `feature_annot.<substrate>.tsv.gz` from the same folder as `--genotypes`; adds provenance and, for variants, `svtype`/`gene`/`AF`/`AN` | — |
| `--unit <auto\|variant\|feature>` | multiple-testing unit; `auto` picks `variant` when the `feature-annot` is the variant sidecar, else `feature` | `auto` |
| `--ld-r2 <X>` | variant unit: genotype r² above which a variant is labelled an LD shadow of a p-seeded lead; clump count is diagnostic and is the `Meff` fallback when Li–Ji is unavailable | `0.8` |
| `--min-ac <N>` | variant unit: flag low `AF` when the observed minority-genotype count < N (underpowered/unstable; such a variant also cannot anchor a clump) | `3` |
| `--cojo-p <X>` | variant unit: entry p for forward-stepwise conditional (COJO) signal selection | `0.05/Meff` |
| `--node-genes <tsv>` | `call`'s `node_genes.tsv` (from `--gtf`); adds a `gene` column | — |
| `--min-maf <X>` | drop features whose [minor non-modal frequency](../algorithms/associate.md#1-filter-on-minor-frequency) < X, on the actual cohort | `0.01` |
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
| `n_conditional` | samples used in the conditional fit, which is the intersection of those observed for this unit and for everything conditioned on. It can be smaller than `n`, and where it is, the two p-values are not directly comparable |
| `minor_freq` | minor (non-modal) genotype frequency on the cohort (the MAF-filter quantity) |
| `beta` \| `log_or` | effect size on the genotype term — `beta` (linear) or `log_or` = log odds ratio (logistic) |
| `se` | standard error of the effect |
| `p_method` | which test produced `p`: `t`, `lmm`, `score`, `score_spa` (saddlepoint-corrected) or `score_exact` (evaluated exactly at the edge of the statistic's support) |
| `effect_status` | how the effect size was obtained, or `separation` where the maximum-likelihood fit diverged and the reported estimate is Firth's |
| `mac_case`, `mac_ctrl` | (binary trait) minor-allele carriers among cases and among controls. The total hides the split that governs reliability: one case against nineteen controls is far weaker evidence than ten against ten |
| `z` | test statistic. Linear/LMM: the Wald statistic, `effect / se`. Logistic: the Rao score statistic, which is not `log_or / se` (see [the algorithm](../algorithms/associate.md#2-fit-and-test-each-unit)) |
| `p` | two-sided p-value — Wald for linear/LMM, score for logistic |
| `p_bonf` | raw Bonferroni-adjusted p, `min(1, p · features_tested)` (over-conservative — kept for reference) |
| `p_bonf_meff` | regional effective-tests guide, `min(1, p · Meff)`; it has no formal family-wise guarantee |
| `q_bh` | Benjamini–Hochberg FDR q-value, the primary FDR summary |
| `af`, `an` | (variant unit) allele frequency / traversing-haplotype count, carried from the VCF (Variant Call Format) |
| `low_af` | (variant unit) `1` when the minority-genotype count < `--min-ac` (underpowered/unstable), else `0`; `.` otherwise |
| `clump`, `is_lead` | (variant unit) LD-clump id and whether this is its lead variant (`1`); `.` in feature mode |
| `gene` | gene name (variant `GENES`, or via `--node-genes`), else `.` |
| `p_conditional` | conditional p (same test as `p`) (variant: vs the COJO-selected set minus self; feature: vs the top lead); `NA` for a sole signal / the lead / a collinear feature |
| `cond_role` | variant: `signal` (independent) / `shadow`; feature: `lead` / `collinear` (same-event, r²>0.95, not scored) / `conditioned`; `.` otherwise |

`<prefix>.summary.tsv` keys:

| key | meaning |
|-----|---------|
| `model`, `phenotype_type` | model used (`linear`/`logistic`/`lmm`) and detected trait type |
| `covariates`, `pca_covariates` | covariate columns used / PCs added via `--pca` |
| `samples_used` | samples kept after dropping rows with NA phenotype or covariate |
| `features_tested` | units that passed the MAF filter and were tested |
| `unit` | `variant` or `feature` (the multiple-testing unit used) |
| `meff`, `meff_method` | regional effective-test estimate and the method that supplied it |
| `meff_eigen`, `meff_ld_clumping` | variant-unit Li–Ji and p-seeded clump estimates, reported separately |
| `independent_variants` \| `distinct_bubbles` | the selected `Meff`, labelled for the active unit |
| `dropped_min_maf`, `dropped_fit` | units dropped by the MAF filter / by a failed model fit |
| `bonferroni_threshold` | raw region-wide threshold `0.05 / features_tested` (reference) |
| `bonferroni_threshold_meff` | regional guide `0.05 / Meff` |
| `significant_bonferroni`, `significant_bonferroni_meff`, `significant_fdr05` | counts passing raw Bonferroni / `Meff`-Bonferroni / BH FDR < 0.05 |
| `cojo_independent_signals` | (variant unit) number of jointly-independent signals from forward-stepwise conditioning |
| `lambda_gc` | genomic-inflation factor λ |
| `lmm_delta` | (LMM only) fitted variance ratio δ = σ²ₑ / σ²_g |

`scripts/plot_associate.R` draws a Manhattan plot of the association table, marking the features that survive each correction and, where a conditional stage ran, the ones that remain independent. It documents its own options under `--help`.

## Limitations

- Variant `Meff` normally uses the phenotype-blind Li–Ji eigenvalue estimate. P-seeded LD clumps provide lead/shadow labels and a fallback estimate; because that count depends on the phenotype, it must not be interpreted as a formal effective-test count. Feature `Meff` is a biological bubble grouping rather than a statistical estimate. Read `p_bonf_meff` as a regional guide; raw Bonferroni is the family-wise reference, and BH is the FDR summary under its usual dependence assumptions.
- Both halves of the variant tier are region-scale. LD clumping compares every pair of features and the eigenvalue estimator forms a correlation matrix over all of them, so neither is genome-scale; above its cap the estimator falls back to clumping or to raw Bonferroni and says so.
- `lambda_gc` assumes most tests are null, which a single locus carrying one large effect violates. There it measures the effect rather than inflation, and the run says which situation it is in rather than reporting a number that reads as inflation either way.
- This is common single-variant association. Firth and the saddlepoint correction make a rare single-variant test better behaved, but there is no burden, collapsing or variance-component test, and rare binary p-values in the far tail remain around 1.7 times nominal.
- The linear mixed model is quantitative-trait only and experimental. Its only external check is a correlation against an established implementation, which cannot detect a systematic difference in effect size, standard error or p-value; the conditional stage is not implemented for LMM rows.

## Example

See the [GWAS example](../gwas.md) for a runnable association run, and the [walkthrough](../walkthrough.md) for the full pipeline.
