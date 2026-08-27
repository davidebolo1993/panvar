# Module `associate`

CLI: `panvar associate`

## What it does

`associate` tests whether a dosage produced by `describe` is associated with a phenotype after accounting for supplied covariates. The phenotype type is detected automatically, or the model can be selected explicitly.

### What is tested

`--unit` determines what one row and one statistical test represent:

| unit | one test is | input |
|------|-------------|-------|
| `variant` | one structural-variant dosage | `describe --variant-vcf` output |
| `feature` | one k-mer, node or edge dosage | `describe` k-mer or graph output |

Variant mode also reports LD clumps and forward-stepwise conditional signals. Feature mode conditions on the top feature and flags nearly collinear features.

### Which statistical test is used

| model | phenotype | p-value test | effect reported |
|-------|-----------|--------------|-----------------|
| `linear` | quantitative | Student-t Wald test | `beta`, `se` |
| `logistic` | binary | Rao score test, with saddlepoint or exact boundary handling when needed | `log_or`, `se`; Firth estimate under separation |
| `lmm` | quantitative with relatedness | EMMAX-style mixed-model Wald test | `beta`, `se` |

For logistic rows, `p` comes from the score test while `log_or` and `se` describe the fitted effect. Therefore `z` is not necessarily `log_or / se`. The `p_method` column records the exact method used.

### Which significance column to read

| column | meaning |
|--------|---------|
| `p` | marginal p-value for this unit |
| `p_conditional` | p-value after accounting for selected lead signal(s) |
| `q_bh` | Benjamini-Hochberg FDR result; the main multiple-testing summary |
| `p_bonf` | Bonferroni adjustment using every tested unit; conservative family-wise reference |
| `p_bonf_meff` | adjustment using an estimated effective test count; regional guide only |

Algorithm and worked trace: [algorithms/associate.md](../algorithms/associate.md).

## Required inputs

- `--genotypes <bimbam.gz>` — a `describe` BIMBAM matrix.
- `--samples <txt[.gz]>` — matching sample order from `describe`.
- `--phenotype <tsv>` — header plus `sample`, `phenotype` and optional covariate columns. A sample with any required `NA` is dropped.
- `-o, --out-prefix <prefix>`.

The genotype matrix and sample file must come from the same `describe` folder. Use the `variant` substrate for variant tests, or `kmers` / `graph` for feature tests.

## Key options

| flag | purpose | default |
|------|---------|---------|
| `--feature-annot <tsv.gz>` | add graph/variant provenance and help auto-detect the unit | — |
| `--unit <auto\|variant\|feature>` | choose the testing unit | `auto` |
| `--model <auto\|linear\|logistic\|lmm>` | choose the phenotype model | `auto` |
| `--min-maf <X>` | drop units with minor non-modal dosage frequency below X | `0.01` |
| `--ld-r2 <X>` | variant mode: r² threshold for lead/shadow clumps | `0.8` |
| `--min-ac <N>` | variant mode: flag very small minority-genotype counts and prevent them seeding a clump | `3` |
| `--cojo-p <X>` | variant mode: conditional entry threshold | `0.05 / Meff` |
| `--kinship <path>` | external genomic relationship matrix for `lmm` or `--pca` | — |
| `--pca <N>` | add the top N kinship PCs to the linear/logistic model | off |
| `--node-genes <tsv>` | add gene labels from `call`'s `node_genes.tsv` | — |

An external, preferably genome-wide kinship matrix must use the same sample order as `--samples`. Alternatively, ancestry PCs can be supplied directly as phenotype-table covariates.

## Outputs

| file | contents |
|------|----------|
| `<prefix>.assoc.tsv` | one row per tested unit, sorted by `p` |
| `<prefix>.summary.tsv` | run settings, filters and regional diagnostics |

### Main association fields

| group | columns | interpretation |
|-------|---------|----------------|
| identity | `feature_id`, `layer`, `bubbles`, `nodes`, `gene` | what was tested and where it came from |
| sample/filter | `n`, `minor_freq`, `af`, `an`, `low_af` | observations and frequency information |
| effect | `beta` or `log_or`, `se`, `effect_status` | fitted effect and how it was obtained |
| test | `z`, `p`, `p_method` | marginal test statistic and p-value |
| correction | `p_bonf`, `p_bonf_meff`, `q_bh` | multiple-testing summaries |
| independence | `p_conditional`, `n_conditional`, `cond_role` | result after conditioning |
| variant LD | `clump`, `is_lead` | variant-mode lead/shadow grouping |
| binary counts | `mac_case`, `mac_ctrl` | minority-genotype carriers in cases and controls |

`cond_role` is `signal` or `shadow` in variant mode, and `lead`, `conditioned` or `collinear` in feature mode. `n_conditional` may be smaller than `n` because all conditioned dosages must be observed.

The summary records the chosen model and unit, samples and features retained, dropped fits, `Meff` and its method, significance counts, the conditional signal count, and `lambda_gc`. LMM runs also report `lmm_delta`.

### Plotting

`scripts/plot_associate.R` plots marginal, FDR and conditional results. Run it with `--help` for its arguments.

## Limitations

- `p_bonf_meff` is descriptive. Variant `Meff` normally uses the phenotype-blind Li-Ji eigenvalue estimate; feature `Meff` groups annotated bubbles and is biological rather than a formal statistical estimate.
- This module performs common single-unit association. It does not implement burden, collapsing, SKAT or other rare-variant aggregate tests. SPA and Firth improve rare single-variant behaviour but do not turn it into rare-variant support.
- `lambda_gc` is hard to interpret within one strongly associated locus because the true signal itself can move it.
- The quantitative-trait LMM is experimental and has not yet been cross-checked against an external implementation with absolute tolerances on effect, standard error and p-value. Conditional analysis is not available for LMM rows.

## Example

See the runnable [GWAS example](../gwas.md) and the full [walkthrough](../walkthrough.md).
