# Associate Module (Module 5 — GWAS)

CLI: `panvar associate`

## What it does

Runs a region-wide association test of a phenotype against the genotypes from `describe`. Per feature (a
k-mer, or a node/edge dosage) it fits `phenotype ~ genotype + covariates` and reports a Wald test on the
genotype term, applies a MAF filter on the cohort genotypes, and corrects for multiple testing over
the features actually tested (region-wide Bonferroni + Benjamini–Hochberg FDR — not the genome-wide
`5e-8`). Phenotype type is auto-detected: binary → logistic (`log_or`), else linear (`beta`).

For a quantitative trait it can also fit a linear mixed model (`--model lmm`) using a kinship matrix, or
add the top kinship PCs as covariates (`--pca N`), to control population structure; every run reports
the genomic-inflation λ. Math + a worked trace: [algorithms/associate.md](../algorithms/associate.md).

## Required inputs

- `--genotypes <bimbam.gz>` — a `describe` BIMBAM matrix (`bimbam_{kmers,graph}.bimbam.gz`, or the
  per-sample `*.samples.bimbam.gz` for a diploid cohort).
- `--samples <txt[.gz]>` — the sample (column) order (`describe`'s `bimbam.samples[.samples].txt.gz`).
- `--phenotype <tsv>` — `sample <tab> phenotype [<tab> covariate…]`, header required; cells may be `NA` (a
  sample with NA phenotype or any NA covariate is dropped).
- `-o, --out-prefix <prefix>`.

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `--feature-annot <tsv.gz>` | `describe`'s `feature_annot.tsv.gz`; adds `layer`/`bubbles`/`nodes` provenance | — |
| `--node-genes <tsv>` | `call`'s `node_genes.tsv` (from `--gtf`); adds a `gene` column | — |
| `--min-maf <X>` | drop features whose [minor non-modal frequency](../algorithms/associate.md#worked-trace--one-quantitative-feature) < X, on the actual cohort | `0.01` |
| `--model <auto\|linear\|logistic\|lmm>` | `auto` = binary→logistic else linear; `lmm` = mixed model (quantitative; needs a kinship source) | `auto` |
| `--kinship <path>` | precomputed `n×n` GRM (rows/cols in `--samples` order) for `--model lmm` / `--pca` | — |
| `--make-kinship` | build the GRM from the genotype matrix (only valid for a genome-wide-like panel; region-only is proximally contaminated) | off |
| `--pca <N>` | add the top-N kinship [PCs as covariates](../gwas/primer.md#two-ways-to-correct-structure-pcs-and-the-lmm) to the GLM | off |
| `-q, --quiet` | less logging | off |

## Outputs

| file | contents |
|------|----------|
| `<prefix>.assoc.tsv` | one row per tested feature, sorted by `p` (columns below) |
| `<prefix>.summary.tsv` | run settings + diagnostics, one `key <tab> value` per line (below) |

`<prefix>.assoc.tsv` columns:

| column | meaning |
|--------|---------|
| `feature_id` | the feature tested — k-mer sequence (k-mer substrate) or node id / edge key (graph substrate) |
| `layer` | `kmer` or `graph` |
| `bubbles`, `nodes` | graph provenance (from `--feature-annot`): which bubble(s) / node(s) the feature comes from |
| `n` | number of samples used in this feature's fit |
| `minor_freq` | minor (non-modal) genotype frequency on the cohort (the MAF-filter quantity) |
| `beta` \| `log_or` | effect size on the genotype term — `beta` (linear) or `log_or` = log odds ratio (logistic) |
| `se` | standard error of the effect |
| `z` | Wald statistic, `effect / se` |
| `p` | Wald p-value |
| `p_bonf` | Bonferroni-adjusted p, `min(1, p · features_tested)` |
| `q_bh` | Benjamini–Hochberg FDR q-value |
| `gene` | gene name when `--node-genes` is given, else `.` |

`<prefix>.summary.tsv` keys:

| key | meaning |
|-----|---------|
| `model`, `phenotype_type` | model used (`linear`/`logistic`/`lmm`) and detected trait type |
| `covariates`, `pca_covariates` | covariate columns used / PCs added via `--pca` |
| `samples_used` | samples kept after dropping rows with NA phenotype or covariate |
| `features_tested` | features that passed the MAF filter and were tested (the multiple-testing denominator) |
| `dropped_min_maf`, `dropped_fit` | features dropped by the MAF filter / by a failed model fit |
| `bonferroni_threshold` | the region-wide threshold `0.05 / features_tested` |
| `significant_bonferroni`, `significant_fdr05` | counts passing Bonferroni / BH FDR < 0.05 |
| `lambda_gc` | genomic-inflation factor λ |
| `lmm_delta` | (LMM only) fitted variance ratio δ = σ²ₑ / σ²_g |

The plotter reads `features_tested` to draw the threshold line.

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

Matches `tests/gwas/run_lpa_real.sh` (region scan, PC-adjusted) and its structure-correction demo:

```bash
# region scan (PCs are covariate columns in the phenotype table)
./build/panvar associate \
  --genotypes results/real_data/lpa/gwas/desc/bimbam_graph.samples.bimbam.gz \
  --samples   results/real_data/lpa/gwas/desc/bimbam.samples.samples.txt.gz \
  --feature-annot results/real_data/lpa/gwas/desc/feature_annot.samples.tsv.gz \
  --node-genes results/real_data/lpa/call/call.node_genes.tsv \
  --phenotype results/real_data/lpa/gwas/real/pheno.quant.tsv --min-maf 0.02 \
  -o results/real_data/lpa/gwas/assoc_graph_quant

# structure correction on a genome-wide-like panel: --model lmm --kinship  (or --pca N)
./build/panvar associate --genotypes <panel.bimbam.gz> --samples <…> --phenotype <…> \
  --model lmm --kinship results/real_data/lpa/gwas/real/kinship.tsv -o <…>_lmm
```

Worked end-to-end run (region scan + naive→PC→LMM λ, GEMMA validation): [gwas/example.md](../gwas/example.md).
Math & worked trace: [algorithms/associate.md](../algorithms/associate.md). References:
[references.md](../references.md#associate).
