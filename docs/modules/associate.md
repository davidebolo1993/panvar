# Associate Module (Module 5 — GWAS)

CLI: `panvar associate`

## What it does

Tests a phenotype against the genotypes from `describe` across a pangenome region. For each unit it fits
`phenotype ~ genotype + covariates`, reports a Wald test on the genotype term, applies a MAF filter on the
cohort genotypes, and corrects for the number of independent tests in the region. Phenotype type is
auto-detected: binary → logistic (`log_or`), else linear (`beta`).

The testable unit is chosen by `--unit`. **Variant** mode tests one genotype per SV call (the
`describe --variant-vcf` export) — the statistically honest unit, since the k-mers, nodes and edges within one
variant are correlated rather than independent. Correlated nearby variants are then collapsed by LD-clumping,
so an LD shadow is not counted as a separate hit. **Feature** mode keeps the fine-grained k-mer/node/edge
tests but corrects with an effective number of independent tests (`Meff`, the distinct bubbles), because the
raw feature count over-states how many independent tests were run. Both report Benjamini–Hochberg FDR (the
primary control) alongside the `Meff`-Bonferroni benchmark and the genomic-inflation λ.

For a quantitative trait it can also add the top kinship PCs as covariates (`--pca N`), or fit a linear mixed
model (`--model lmm`) against a kinship matrix, to control population structure.

Algorithm and worked trace: [algorithms/associate.md](../algorithms/associate.md).

## Required inputs

- `--genotypes <bimbam.gz>` — a `describe` BIMBAM matrix: `bimbam_{kmers,graph}.bimbam.gz` (feature unit)
  or `bimbam_variant.bimbam.gz` (variant unit), or the per-sample `*.samples.bimbam.gz` for a diploid cohort.
- `--samples <txt[.gz]>` — the sample (column) order (`describe`'s matching `*.samples[.samples].txt.gz`).
- `--phenotype <tsv>` — `sample <tab> phenotype [<tab> covariate…]`, header required; cells may be `NA` (a
  sample with NA phenotype or any NA covariate is dropped).
- `-o, --out-prefix <prefix>`.

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `--feature-annot <tsv.gz>` | `describe`'s `feature_annot.tsv.gz` (feature unit) or `feature_annot.variant.tsv.gz` (variant unit); adds provenance and, for variants, `svtype`/`gene`/`AF`/`AN` | — |
| `--unit <auto\|variant\|feature>` | multiple-testing unit; `auto` picks `variant` when the feature_annot is the variant sidecar, else `feature` | `auto` |
| `--ld-r2 <X>` | variant unit: genotype r² above which a variant is an LD shadow of a lead (clumped, not counted in `Meff`) | `0.8` |
| `--min-ac <N>` | variant unit: flag `low_af` when the observed minority-genotype count < N (underpowered/unstable) | `3` |
| `--node-genes <tsv>` | `call`'s `node_genes.tsv` (from `--gtf`); adds a `gene` column | — |
| `--min-maf <X>` | drop features whose [minor non-modal frequency](../algorithms/associate.md#worked-trace--one-quantitative-feature) < X, on the actual cohort | `0.01` |
| `--model <auto\|linear\|logistic\|lmm>` | `auto` = binary→logistic else linear; `lmm` = mixed model (quantitative; needs a kinship source) | `auto` |
| `--kinship <path>` | external (genome-wide) `n×n` GRM (rows/cols in `--samples` order) for `--model lmm` / `--pca`; panvar is local and does not build a GRM itself | — |
| `--pca <N>` | add the top-N kinship PCs as covariates to the GLM (needs `--kinship`); usually you instead pass ancestry PCs as phenotype-table columns | off |
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
| `af`, `an` | (variant unit) allele frequency / traversing-haplotype count, carried from the VCF |
| `low_af` | (variant unit) `1` when the minority-genotype count < `--min-ac` (underpowered/unstable), else `0`; `.` otherwise |
| `clump`, `is_lead` | (variant unit) LD-clump id and whether this is its lead variant (`1`); `.` in feature mode |
| `gene` | gene name (variant `GENES`, or via `--node-genes`), else `.` |

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
| `lambda_gc` | genomic-inflation factor λ |
| `lmm_delta` | (LMM only) fitted variance ratio δ = σ²ₑ / σ²_g |

## Plotting

```bash
Rscript scripts/plot_associate.R --assoc <prefix>.assoc.tsv --summary <prefix>.summary.tsv \
  --out <prefix> --title "my trait"
```

Writes `*.manhattan.{png,pdf}` — two stacked panels, before correction (raw −log10 p with nominal +
region-wide Bonferroni lines) and after correction (Benjamini-Hochberg −log10 q with the q=0.05 line); x =
node id (graph) or per-k-mer index ordered by node id (k-mers), with FDR/Bonferroni-significant genes
flagged (ggrepel, from the `gene` column when `--node-genes` was passed) — and `*.qq.{png,pdf}` (with λ).

Script flags (needs `Rscript` + `ggplot2`; `ggrepel` optional, for the gene labels):

- `--assoc <assoc.tsv>` — the association table (required).
- `--out <prefix>` — output prefix for the PNG/PDF files (required).
- `--summary <summary.tsv>` — read `features_tested` for the Bonferroni line (recommended).
- `--title <text>` — plot title.
- `--width` / `--height` / `--dpi` — Manhattan size (inches) and PNG resolution (defaults 10 / 7 / 150).

## Example

Uses the committed example phenotype (`tests/gwas/lpa/`, simulated) over the BIMBAM matrices that
`describe` wrote under `results/real_data/lpa/gwas/desc/` (see the [GWAS example](../gwas/example.md) for the
full pipeline):

```bash
# region scan, feature unit (PC1..PC10 are covariate columns in the phenotype table)
./build/panvar associate \
  --genotypes results/real_data/lpa/gwas/desc/bimbam_graph.samples.bimbam.gz \
  --samples   results/real_data/lpa/gwas/desc/bimbam.samples.samples.txt.gz \
  --feature-annot results/real_data/lpa/gwas/desc/feature_annot.samples.tsv.gz \
  --node-genes results/real_data/lpa/call/call.node_genes.tsv \
  --phenotype tests/gwas/lpa/pheno.quant.tsv --min-maf 0.02 \
  -o results/real_data/lpa/gwas/assoc_graph_quant

# variant unit: test the SV calls directly (describe --variant-vcf export); --unit auto-detects it
./build/panvar associate \
  --genotypes results/real_data/lpa/gwas/desc/bimbam_variant.samples.bimbam.gz \
  --samples   results/real_data/lpa/gwas/desc/bimbam_variant.samples.samples.txt.gz \
  --feature-annot results/real_data/lpa/gwas/desc/feature_annot.variant.tsv.gz \
  --phenotype tests/gwas/lpa/pheno.quant.tsv -o results/real_data/lpa/gwas/assoc_variant

# optional: control structure with an EXTERNAL genome-wide GRM (LMM); panvar does not build one itself
./build/panvar associate --genotypes <panel.bimbam.gz> --samples <…> --phenotype <…> \
  --model lmm --kinship <genome_wide_grm.tsv> -o <…>_lmm
```

A worked end-to-end run with interpretation is in the [GWAS example](../gwas/example.md).
