# Module `associate`

CLI: `panvar associate`

## What it does

Asks, for every row of a `describe` matrix, whether the dosage goes with a phenotype: do samples carrying more copies of a k-mer, node or variant have a higher trait value, or are they more often cases? Each row gets an effect size and a p-value, with covariates such as sex or ancestry components taken into account.

| model | phenotype | effect |
|-------|-----------|--------|
| `linear` | a measurement | `beta`: change in phenotype per extra copy |
| `logistic` | case/control (0/1, or PLINK's 1/2) | `log_or`: change in log odds per extra copy |
| `lmm` | a measurement, in related samples (with a kinship matrix) | `beta` |

`--model auto` picks logistic for a 0/1 phenotype and linear otherwise. On top of the tests of each row:

- correction for the number of tests, and which significant rows are independent signals rather than echoes of a stronger one;
- with `--unit variant`, fine-mapping: each record's probability of being a cause, grouped into credible sets; `--unit feature` is for k-mers, nodes and edges;
- with `--group`, rare records tested together per gene or site, as one record alone has too few carriers;
- with `--haplotypes` and `--encoding`, the smaller or larger of a sample's two haplotype values besides their sum, for where one allele sets the phenotype;
- several matrices in one run, such as `call`'s records beside other callers' SNVs and SVs, to say which source holds the strongest signal and whether the others add to it.

## How it works

The tests are standard ones, implemented in `associate` with [Eigen](https://gitlab.com/libeigen/eigen): least squares; the logistic score test, with exact and [saddlepoint](https://doi.org/10.1016/j.ajhg.2017.05.014) tails and [Firth](https://doi.org/10.1093/biomet/80.1.27)'s fit where the ordinary one diverges; a REML mixed model as in [GEMMA](https://doi.org/10.1038/ng.2310), or with the variance fitted once as in [EMMAX](https://doi.org/10.1038/ng.548); [Benjamini-Hochberg](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x) and [Li-Ji](https://doi.org/10.1038/sj.hdy.6800717) corrections; [Freedman-Lane](https://doi.org/10.1080/07350015.1983.10509354) permutations; [COJO](https://doi.org/10.1038/ng.2213)-style selection; [SuSiE](https://doi.org/10.1111/rssb.12388) fine-mapping; and for rare records [SKAT](https://doi.org/10.1016/j.ajhg.2011.05.029) and [SKAT-O](https://doi.org/10.1016/j.ajhg.2012.06.007), with [Imhof](https://doi.org/10.1093/biomet/48.3-4.419)'s inversion, [SAIGE-GENE+](https://doi.org/10.1038/s41588-022-01178-w)'s variance adjustment and [ACAT](https://doi.org/10.1016/j.ajhg.2019.01.002). [GEMMA](https://github.com/genetics-statistics/GEMMA) gives the same linear and mixed-model results on the same BIMBAM files.

1. Each `--genotypes` matrix is read with its `--samples` file. With several, samples are matched by name and only those in every matrix are kept; each row keeps its `layer` (`--feature-annot`), so its source shows in every output. With `--haplotypes`, the columns are haplotypes, and `--encoding` forms each sample's value: `sum`, the dosage `describe --samples` writes; `min` and `max`, the smaller and the larger of its two haplotype values (rows `<id>:min`, `<id>:max`). For a present/absent record, `max` is the dominant coding and `min` the recessive one. Encoded rows are rows like any other: where the effect is additive they shadow the sum.
2. A sample enters only with a phenotype and every covariate. A phenotype of only 0 and 1 (or 1 and 2, as PLINK writes) is case/control. Inputs that cannot be tested are refused before any row: a constant phenotype, covariates that depend on one another, or covariates that reproduce a measured phenotype or separate cases from controls. With `--rint`, a measured phenotype is replaced by the normal quantile of its rank, so a few extreme values cannot drive the fits.
3. A row where nearly every sample has the same dosage cannot show an association: rows with `minor_freq` (the share of samples off the most common dosage, which works for copy numbers as well as for presence) under `--min-maf` are not tested.
4. For a measurement, the model is fitted by least squares: `beta`, its `se`, and a Student-t test, which allows for the noise being estimated from the same samples. A sample with no dosage in a row is left out of that row's fit. Each row is also tested with its dosage classes (0, 1, 2, … copies) each given its own mean: `p_genotypic` against the covariates alone, and `p_nonadditive` against the additive line, which flags a threshold, dominant or U-shaped response.
5. For case/control, the p-value comes from a score test, which needs only the model without the dosage: the usual Wald test fails for rare dosages, where the fit pushes the odds ratio towards infinity. With up to `--exact-carriers` carriers (with covariates, also at most 5% of the samples), the tail is counted exactly over the carriers' case/control outcomes (Fisher's exact test when there are no covariates); otherwise, past |z| = 2, it is corrected by a saddlepoint approximation. `log_or` and `se` come from the ordinary fit, or Firth's where it diverges (`effect_status` `separation`).
6. For related samples (`--model lmm`), the phenotype is modelled as dosage + covariates + a background shared in proportion to kinship + noise. The kinship's eigenvectors turn the data into coordinates where the samples are independent. `lmm_h2` is the background's share without any dosage; each row refits it by REML with the row in the model, as GEMMA does (`--lmm-fast` keeps the null share, as EMMAX). A missing dosage takes the row's mean. `--pca N` is a lighter alternative: the kinship's top N eigenvectors as covariates.
7. Testing many rows gives small p-values by chance. `p_bonf` multiplies `p` by the number of tests, which is strict as rows of one locus are correlated; `q_bh` controls the share of false hits among the significant rows; `p_bonf_meff` multiplies by the number of independent records (`Meff`, Li and Ji), reported only, as it undercounts on dense, correlated records.
8. `--permutations N` gives the family-wise threshold that holds however correlated the rows are: the run is repeated N times on phenotypes with no association but the same covariates, keeping each time the smallest p. The phenotype is shuffled, not the genotypes, so the rows keep their link to the covariates (residuals for a measurement, Freedman-Lane; phenotypes drawn from the covariates-only fit for case/control; rotated residuals for the LMM). A row's `p_perm` is (1 + k) / (N + 1), k the shuffles whose smallest p was at or below its own; `permutation_threshold` is the p below which `p_perm` is under 0.05, and needs 20 shuffles or more.
9. A strong variant makes every row correlated with it significant too. In the variant unit, records are selected forward, as in GCTA-COJO: the record with the smallest p given those already selected enters while it is below the entry threshold (`--cojo-p`, else the permutation threshold, else Bonferroni), and a selected record no longer below it given the others is dropped. Of two records that test the same given those selected, as a sample's `min` and `max` once the sum is in (max = sum − min), the one stronger on its own enters. Selected records are `signal`s; every other record is tested again with the signals as covariates, and is a `shadow`, or `collinear` when the signals reproduce it. A conditioning record's missing dosages take its mean, with an indicator of where it is missing: a SNV a deletion removes is missing in the deletion's carriers, and that carries the deletion's effect. Where a conditioning record separates cases from controls, the test runs on the other samples (`n_conditional`). Records are also grouped into clumps of r² above `--ld-r2`.
10. In the feature unit, the strongest feature is the `lead`; features with r² above 0.95 to it are `collinear`, the others `conditioned` on it. This tells whether anything beyond the lead is associated, not how many signals there are, which needs the variant unit. It is not run under the LMM.
11. Per source (`layer`), `<prefix>.layers.tsv` gives the strongest row, the share of the phenotype it explains (`best_r2`), its r² with the run's top row, and its role after conditioning. Each source's best row is tested against the top row in both directions (`p_best_given_top`, `p_top_given_best`). For the linear model in the variant unit, `delta_r2` compares two sources out of sample, over 1,000 bootstrap resamples: a source with more rows otherwise looks better just by picking its best on the same data. An interval spanning 0 means the two cannot be ranked.
12. In the variant unit, the records are fine-mapped with SuSiE, as susieR fits it: up to `--finemap-l` effects, each one record of unknown identity. A record's `pip` is the chance that some effect is it; an effect's credible set is the fewest correlated records holding `--cs-coverage` of its probability. The probabilities assume one of the tested records is the cause: near-copies of one event, such as the same deletion from two callers, share it.
13. With `--group`, the rare records (`minor_freq` at most `--rare-maf`) are tested together per gene, site, source or all at once, as the SKAT package does: each weighted by the Beta(1, 25) density of its `minor_freq`, so the rarest weigh most. The burden test suits records acting in one direction, SKAT records raising and lowering the phenotype, SKAT-O the best mixture of the two, corrected for the choice; `acat_p` combines them. For case/control, a group with few carriers is counted exactly as one record is, and otherwise its variances are adjusted as in SAIGE-GENE+, so few cases do not give far too small p-values. Rare records often ride on a common signal, so SKAT-O is also computed given the conditional selection's signals (`skato_p_conditional`). `--group` needs `sum` among the encodings.

## Required inputs

- `--genotypes <bimbam.gz>`: a BIMBAM matrix, from `describe` or another tool, as a regular file (it is read more than once). Repeat for several.
- `--samples <txt[.gz]>`: its column order, one name per line, such as the `samples.txt.gz` beside a `describe` matrix. One per `--genotypes`, in the same order.
- `--phenotype <tsv>`: a header, then one row per sample: sample, phenotype, then any covariates. A missing value is empty, `NA`, `NaN`, `.` or `null`, and a phenotype of -9 is missing too (PLINK's convention). A row must have as many columns as the header. A sample missing any value is left out.

A matrix and its samples file must come from the same folder (`describe`'s `haplotype/` or `sample/`); with `--haplotypes`, every matrix is from `haplotype/`. Row ids must be unique across matrices, and across `--feature-annot` files.

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `-o, --out-prefix <prefix>` | output prefix | `associate` |
| `--feature-annot <tsv.gz>` | `describe`'s sidecar, or a table like it: id, layer, encoding, sites and nodes, then, found by name, any of `svtype`, `gene`, `AF`, `AN`, `chrom`, `pos`, `end`. Repeat for several | — |
| `--haplotypes <tsv>` | the sample table given to `describe --samples` (sample, then one column per haplotype, contigs comma-separated): the matrices are per haplotype, and each sample's value is formed here | — |
| `--encoding <list>` | with `--haplotypes`: `sum`, `min` and/or `max` of a sample's haplotype values, comma-separated | `sum` |
| `--node-genes <tsv>` | `call`'s `node_genes.tsv` (`call --gtf`): genes from each row's nodes | — |
| `--model <auto\|linear\|logistic\|lmm>` | the model | `auto` |
| `--kinship <path>` | a genome-wide kinship matrix, n × n in the order of the first `--samples` (with `--haplotypes`, of its samples), no header; for `lmm` and `--pca` | — |
| `--lmm-fast` | LMM: the null variance share for every row (EMMAX) instead of each row's own | off |
| `--rint` | rank-based inverse-normal transform of a measured phenotype | off |
| `--pca <N>` | add the kinship's top N eigenvectors as covariates (refused when the N-th and next eigenvalues are equal, since the vectors are then not unique) | `0` |
| `--min-maf <X>` | skip rows where fewer than this share of samples differ from the most common dosage | `0.01` |
| `--unit <auto\|variant\|feature>` | the testing unit; `auto` is `variant` when most matrix rows are records; k-mers, nodes, edges and rows without an annotation count against | `auto` |
| `--permutations <N>` | N null replicates (shuffled residuals, or drawn phenotypes for case/control) for a family-wise threshold and `p_perm`; recommended in both units | `0` |
| `--seed <N>` | seed of the shuffles | `1` |
| `--ld-r2 <X>` | correlation (r²) above which a record joins a stronger one's clump (variant unit), and above which a conditioning missingness indicator that copies the tested row is left out (items 9 and 10) | `0.8` |
| `--min-ac <N>` | variant unit: flag `low_af` below N samples off the most common dosage | `3` |
| `--cojo-p <X>` | variant unit: p a record needs to enter the conditional selection | the permutation threshold with 20 or more shuffles, else 0.05 / tests |
| `--finemap-l <N>` | variant unit: fine-map with up to N effects (0 to 50); 0 turns fine-mapping off | `10` |
| `--cs-coverage <X>` | probability a credible set holds | `0.95` |
| `--group <gene\|bubble\|layer\|all>` | variant unit: test the rare records together per gene, per site, per source or all at once | — |
| `--rare-maf <X>` | with `--group`: a record is rare at `minor_freq` at most X | `0.01` |
| `--exact-carriers <N>` | logistic: compute the tail exactly when at most N samples differ from the most common dosage (0 to 20; with covariates, also at most 5% of the samples, else the saddlepoint at any z) | `16` |
| `-q, --quiet` | no progress logs | off |

## Outputs

| file | contents |
|------|----------|
| `<prefix>.assoc.tsv` | one row per tested row, sorted by `p` |
| `<prefix>.summary.tsv` | the run's settings, counts, thresholds and diagnostics, as key and value |
| `<prefix>.layers.tsv` | one row per layer: its strongest row and how many of its rows pass each threshold |
| `<prefix>.cs.tsv` | variant unit: one row per credible set (only the header otherwise) |
| `<prefix>.groups.tsv` | with `--group`: one row per group, sorted by `skato_p` |

Association columns:

| column | meaning |
|--------|---------|
| `feature_id`, `layer`, `bubbles`, `nodes` | the row, its source, and its sites and nodes from `--feature-annot` |
| `chrom`, `pos`, `end` | where the row lies on the reference, from `--feature-annot`; `.` when it gives none |
| `n`, `n_conditional` | samples in the test, and in the conditional test (for case/control, less those a separating conditioning record accounts for, item 9); under the LMM, `n` counts the samples with a dosage, while both tests use every sample |
| `minor_freq` | share of samples whose dosage is not the most common one |
| `beta` or `log_or`, `se` | the effect per extra copy and its standard error |
| `stat` | the test statistic: Student-t for `linear` and `lmm`, the score z for `logistic` |
| `p`, `p_method` | the p-value and how it was computed: `t`, `score`, `score_exact_carriers`, `score_spa`, `score_exact`, `score_normal_fallback`, `lmm` |
| `effect_status` | `ok`, or `separation` when the logistic fit diverged (the effect then comes from Firth's method; the p is unaffected) |
| `mac_case`, `mac_ctrl` | logistic: samples off the most common dosage among cases and among controls |
| `p_bonf` | `p` times the number of tests |
| `p_bonf_meff` | `p` times `Meff`; `.` in the feature unit |
| `q_bh` | Benjamini-Hochberg q-value |
| `p_perm` | (1 + shuffles whose smallest p was at or below this row's) / (shuffles + 1): its `p`, or under the LMM its p at the null share (item 6); `.` without `--permutations` |
| `af`, `an`, `low_af` | variant unit: the record's frequency and allele number from `--feature-annot` (over the panel's haplotypes for a haplotype-level sidecar, and over the samples' haplotypes for a sample-level one), and whether few samples are off the most common dosage |
| `clump`, `is_lead` | variant unit: the correlation group and whether the row leads it |
| `gene` | the variant's genes, else those of its nodes (`--node-genes`), `;`-separated |
| `p_conditional`, `cond_role` | the conditional p and the row's role: `signal`, `shadow` or `collinear` in the variant unit, `lead`, `conditioned` or `collinear` in the feature unit, `.` when nothing was conditioned on |
| `pip`, `cs` | variant unit: the chance the row is a cause, and its credible set (`.` in none) |
| `p_genotypic`, `p_nonadditive` | linear model: the dosage classes as a general effect, against the covariates alone and against the additive line; `NA` with fewer than 2 (3) classes |

Layer columns:

| column | meaning |
|--------|---------|
| `layer`, `rows_tested` | the source and its rows tested |
| `significant_bonferroni`, `significant_bonferroni_meff`, `significant_permutation`, `significant_fdr05` | its rows passing each threshold (`NA` where the threshold was not computed) |
| `best_feature`, `best_beta` or `best_log_or`, `best_p` | its strongest row (the smallest p, ties to the larger \|`stat`\|) |
| `best_r2` | the share of the phenotype that row explains beyond the covariates, t² / (t² + degrees of freedom), for the linear model; `NA` for the logistic model and the LMM, whose t is on a rotated scale. Rows with missing dosages are fitted on fewer samples, so compare layers on similar `n` |
| `r2_with_top` | that row's dosage correlation (r²) with the run's strongest row |
| `best_p_conditional`, `best_cond_role` | that row's conditional p and role |
| `cs_prob` | the layer's probability summed over the credible sets, each holding at least `--cs-coverage`: about 1 for each signal it holds alone |
| `p_best_given_top`, `p_top_given_best` | that row given the run's top row, and the top row given it, on the samples both measure; `NA` for the top layer |
| `delta_r2`, `delta_r2_low`, `delta_r2_high` | linear model, variant unit: this layer's explained share minus the top layer's, best rows re-picked and scored out of sample over bootstrap resamples: the median and a 95% interval |

Credible-set columns:

| column | meaning |
|--------|---------|
| `cs`, `rows` | the set, numbered from the strongest, and its records |
| `coverage`, `purity` | the probability it holds, and the smallest \|r\| between its records |
| `top_feature`, `top_prob` | its most probable record and that record's probability |
| `layers` | the set's probability per layer, `;`-separated |
| `features` | its records, by decreasing probability |

Group columns:

| column | meaning |
|--------|---------|
| `group`, `rows`, `carriers` | the gene, site, layer or `all`; its rare records; the samples carrying any of them |
| `burden_p`, `skat_p`, `skato_p` | the three tests' p-values |
| `skato_rho` | the SKAT-O weight of the burden part at the best p: 0 is SKAT, 1 the burden test |
| `p_bonf` | `skato_p` times the number of groups |
| `skato_p_conditional` | SKAT-O with the selected signals as covariates; `NA` when none was selected, or nothing is left to test |
| `acat_p` | burden, SKAT and SKAT-O under two weightings, combined (ACAT-O) |
| `features` | its records |

The summary holds:

- the model, phenotype type, covariates (`pca_covariates` among them), samples used and rows read;
- `features_tested`, `dropped_min_maf` and `dropped_fit`;
- the unit, `meff` with `meff_method`, `meff_eigen` and `meff_ld_clumping` (`NA` in the feature unit);
- the thresholds `bonferroni_threshold`, `bonferroni_threshold_meff` and `permutation_threshold` (a row is family-wise significant when its p is below it, under the LMM its p at the null share), with `permutations`;
- the counts `significant_bonferroni`, `significant_bonferroni_meff`, `significant_permutation` and `significant_fdr05`;
- `score_normal_fallback`, `lambda_gc`;
- in the variant unit, `cojo_entry_threshold` with `cojo_entry_source` (`--cojo-p`, `permutation` or `bonferroni`) and `cojo_independent_signals`; and `finemap_effects`, `credible_sets` and `finemap_converged`;
- with `--group`, `group_unit`, `rare_maf`, `rare_rows`, `groups_tested` and `groups_conditioned_on` (the signals);
- for the LMM, `lmm_delta`, `lmm_h2` (the share of phenotype variance the kinship explains) and `lmm_delta_at_boundary`.

### Visualization

`scripts/plot_associate.R` draws runs along the reference:

- One run is drawn over its region in Mb. Several (loci) share one row, as in a genome-wide plot: each locus a block named under it with its selection threshold, on its own y scale, so its hits are placed by its own values.
- Each row's −log10 p sits at its place (`chrom`, `pos`, `end`); a record is a bar over its span. Rows are coloured by source; credible-set members are ringed.
- The strongest row is named with its genes, the next ones with their ids (`plot_associate.R --top`: 5 for one run, 1 for several).
- The lines are Bonferroni, the permutation threshold when the run had `--permutations`, and, in dark red, the one the conditional selection used (step 9).
- In the variant unit, a second row shows −log10 p given the signals: what each row adds once they are accounted for.
- When a few rows stand far above the rest, the axis is broken. A linear or mixed model's p at the output's floor (1e-300) is placed from its statistic, so such rows keep their order (a case/control score test's stays at the floor); rows with no place on the reference are left out and counted.

```
Rscript scripts/plot_associate.R [<label>=]<associate prefix> [...]
```

## Limitations

- A dosage is tested as additive. Where one allele sets the phenotype, `--encoding min,max` tests the smaller and larger allele too (it needs haplotype-level matrices and adds rows, so the threshold is stricter); an effect of a particular pair of alleles is not tested.
- The exact and saddlepoint tails make a rare row's p right, not its power higher; `--group` pools rare records. With very few cases, a saddlepoint tail can still be several times too small (`--exact-carriers 20` counts more rows exactly), and with 1% to 2% cases SKAT and SKAT-O can exceed their level; the burden test does not.
- Sources are compared as their callers report them: one event called as several records divides its signal among them, so the source's best row underrates it.
- A SNV inside a deletion reads `NA` in the deletion's carriers, so it does not stand in for the deletion.
- In the feature unit Bonferroni is far too strict: use `--permutations`. Shuffles assume samples exchangeable given the covariates; relatives call for the LMM.
- The LMM's shuffles, conditional selection and fine-mapping use the null variance share. There is no LMM for case/control.
- Fine-mapping fills a missing dosage with the record's mean and fits case/control as a linear model, so its probabilities are approximate; `pip` summed over a source overstates it, `cs_prob` does not.
- The kinship and the PCs must come from outside the tested chromosome, the PCs from LD-pruned SNVs: built from the locus, they absorb the signal they adjust for.
- Calls whose accuracy differs between cases and controls, batches or coverage create associations no test removes: give batch and coverage as covariates.

## Example

`scripts/associate_examples.sh` tests `call`'s records beside SNVs and SVs called from the same haplotypes, on simulated cohorts of 1,000 samples at two loci, in `tests/examples/associate`. The inputs are in `tests/synthetic_data/associate`; how they are made is in [dev/README.md](../dev/README.md). Per locus, `associate` tests the three matrices together, with 1,000 permutations, four times:

- with sex and the PCs, each source's rare records also tested together (`--group layer`);
- with sex only, where the ancestry offset confounds the test;
- with the PCs, the simulated genetic part and the population offset, where any hit is a false one;
- as the first, on the haplotype-level matrices with `--encoding sum,min,max`.

`tests/examples/associate/sources.tsv` collects every run's layers table, and the first run of both loci is drawn with `scripts/plot_associate.R`, in `tests/examples/associate/associate.region.png`.

Where the phenotype follows a gene's copies, `call`'s copy-number record measures the simulated cause itself, so the run asks how closely the SNVs and SVs recover it; the repeat array, whose effect is not additive in its summed copies, is the harder case for every source. There the run on haplotypes keeps the smaller allele's copies beside the sum, and the rare records that are significant given the sum alone are not given both.
