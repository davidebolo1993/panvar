# Module `associate`

CLI: `panvar associate`

## What it does

Tests a phenotype against the genotypes from `describe` across a pangenome region. For each unit it fits `phenotype ~ genotype + covariates`, reports a test on the genotype term (Wald for quantitative traits, a Rao score test for binary ones — see [Calibration](#calibration)), applies a minor-allele-frequency (MAF) filter on the cohort genotypes, and corrects for the number of independent tests in the region. Phenotype type is auto-detected (binary or quantitative). The testable unit is chosen by `--unit`. Variant mode tests one genotype per structural-variant (SV) call (the `describe --variant-vcf` export) — the statistically honest unit, since the k-mers, nodes and edges within one variant are correlated rather than independent; correlated nearby variants are then collapsed by linkage-disequilibrium (LD) clumping, so an LD shadow is not counted as a separate hit. Feature mode keeps the fine-grained k-mer/node/edge tests but corrects with an effective number of independent tests (`Meff`, the distinct bubbles), because the raw feature count over-states how many independent tests were run. Both report Benjamini–Hochberg (BH) false-discovery-rate (FDR) control alongside the `Meff`-Bonferroni benchmark and the genomic-inflation `λ`. 
Beyond the threshold, `associate` also tests independence by conditioning: refitting each unit with the top signal(s) added as covariates, so a hit that merely tags a stronger nearby variant is exposed (its `p` collapses). The variant unit runs a forward-stepwise selection of jointly-independent signals (conditional-and-joint, COJO-style); the feature unit conditions on the single top feature with a within-bubble collinearity guard. 
For a quantitative trait it can also add the top kinship principal components (PCs) as covariates (`--pca N`), or fit a linear mixed model (LMM, `--model lmm`) against a kinship matrix, to control population structure.

Algorithm and worked trace: [algorithms/associate.md](../algorithms/associate.md).

## Required inputs

- `--genotypes <bimbam.gz>` — a `describe` BIMBAM matrix from one substrate folder: `<level>/<substrate>/bimbam_<substrate>.bimbam.gz`, where `<level>` is `haplotype` or `sample` (diploid cohort) and `<substrate>` is `kmers`/`graph` (feature unit) or `variant` (variant unit).
- `--samples <txt[.gz]>` — the sample (column) order (`describe`'s matching `.samples[.samples].txt.gz`).
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
| `z` | test statistic. Linear: the Wald statistic, `effect / se`. Logistic: the Rao score statistic, which is not `log_or / se` (see [Calibration](#calibration)) |
| `p` | two-sided p-value for `z` — Wald for linear, score for logistic |
| `p_bonf` | raw Bonferroni-adjusted p, `min(1, p · features_tested)` (over-conservative — kept for reference) |
| `p_bonf_meff` | effective-tests Bonferroni, `min(1, p · Meff)` — the honest correction |
| `q_bh` | Benjamini–Hochberg FDR q-value (the primary control) |
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

`scripts/plot_associate.R` draws a Manhattan plot of the association table, marking the features that survive each correction and, where a conditional stage ran, the ones that remain independent. It documents its own options under `--help`.

## Example

See the [GWAS example](../gwas.md) for a runnable association run, and the [walkthrough](../walkthrough.md) for the full pipeline.

## Effective tests

`Meff` is the denominator the Bonferroni threshold uses, so every tested feature has to be in it. In the variant tier a low-AF variant is barred from anchoring an LD clump — its `r^2` is unstable, so it must not claim shadows — but it is still tested, and it used to end up in no clump at all: its own threshold was then computed from a set it was not part of. Each such feature therefore gets a singleton clump: it counts toward the denominator but claims no shadows. The same rule applies in the feature tier, where a feature with no bubble annotation belongs to no block.

`Meff` is an LD-clumping heuristic, not an effective-test count. Clumps are seeded in p-value order, so the phenotype changes how many there are: in a correlation chain A—B—C with A and C uncorrelated, seeding at B gives one clump and seeding at A gives two. It is therefore not phenotype-blind and carries no formal family-wise guarantee. `p_bonf` (raw `0.05/n_tests`) and `q_bh` are the formal corrections and are what should be quoted; read `p_bonf_meff` as a regional guide. In the variant tier `Meff` is therefore now the phenotype-blind Li & Ji (2005) eigenvalue estimator on the genotype correlation matrix — `Meff = Σ [I(λᵢ ≥ 1) + frac(λᵢ)]` — which never looks at the phenotype, so the threshold it implies cannot be circular. Both are reported (`meff_eigen`, `meff_ld_clumping`) with `meff_method` naming which drove the threshold; `tests/associate_null.sh` also reports the empirical min-p 5% quantile, a maxT regional threshold that assumes nothing about independence (0.0031 there, against a Bonferroni 0.05/n of 0.0025). Feature-mode `Meff` (the number of distinct bubbles) is likewise a biological grouping, not a statistical one.

## Reading `lambda_gc`

`lambda_gc` is a genome-wide diagnostic: it reads the median chi-square on the assumption that most tests are null. `panvar` tests one locus, where a real signal and everything in linkage disequilibrium with it can be most of the tests — `lambda` then measures the signal, not inflation. At a locus carrying one large effect it measures that effect. The run summary now says which situation it is in, and only calls it an inflation estimate when there are at least 100 tests and fewer than a quarter of them are significant.

## Calibration

A p-value means nothing unless it is uniform when nothing is going on. `tests/associate_null.sh` permutes the phenotype table's sample labels — severing every genotype-phenotype link while leaving the genotype matrix, the missingness pattern and the phenotype/covariate joint distribution untouched — and reports type-I error, `lambda_gc`, and a per-feature uniformity test. Per feature, p-values across permutations are independent, which is what makes the uniformity test valid; pooled across features within one permutation they are not, so the pooled intervals it prints are optimistic and labelled as such.

On a cohort of a few thousand individuals over 300 permutations:

| model | lambda_gc | features rejecting uniformity |
|-------|-----------|-------------------------------|
| linear, default `--min-maf` | 1.088 | 2/13 (95% bound is 2) |
| logistic, default `--min-maf` | 1.006 | 0/13 |
| logistic, `--min-maf 0` | 1.010 | 0/20 |

Quantitative traits use a Student-t tail on `n - p` degrees of freedom, not the normal one: the residual variance is estimated from the same data, so the statistic is t-distributed. The difference is invisible at GWAS sample sizes but not at small ones — at 3 degrees of freedom, `t = 1.96` is `p = 0.145`, where the normal would say `0.050`. The implementation is validated against R's `pt()` to a relative error below 4e-12 over df 3–5000.

Binary traits use a Rao score test. The Wald test divides an estimate by its own standard error, and for a rare variant in an unbalanced case/control study the fit approaches separation: the coefficient grows, its standard error grows faster, and the statistic collapses. Measured here before the change, every feature below minor frequency 0.01 failed uniformity (`lambda_gc` 0.80, worst KS p 3e-12) while every feature above it passed — the default `--min-maf 0.01` was the only thing hiding it. The score test never fits the alternative, so it has no standard error to inflate.

Consequence for reading the output: `z` and `p` are the score test, while `log_or` and `se` remain a Wald-style effect size, so p is not recoverable from `log_or`/`se`. That effect size is the maximum-likelihood estimate normally, and Firth's penalised-likelihood estimate where the maximum-likelihood fit diverges under separation. This is the same arrangement SAIGE and REGENIE use.

Separation contract. Under near-complete separation the maximum-likelihood logistic fit diverges, so IRLS hitting the iteration cap returns failure rather than its last iterate. The score test does not fit the alternative and is unaffected, so the feature is still reported with a valid `p`, and `effect_status=separation` says what happened. The effect size comes from Firth's penalised likelihood (Jeffreys prior), which is finite and first-order unbiased where maximum likelihood diverges — so `log_or` and `se` are present unless Firth itself fails, in which case they are `NA`. A feature is only dropped (counted in the `(fit)` term of the run summary) when the score test also fails. There is never a silent fall back to a Wald p.

Scope: this is common single-variant association. `panvar associate` runs one test per feature. Firth and the saddlepoint correction make the single-variant rare test better behaved, but they are not rare-variant support: there is no burden test, no SKAT-style variance component and no collapsing of rare features into a gene or bubble unit. Rare features are tested individually and are simply underpowered.

Rare binary features are exploratory. The score test is well calibrated in the body of the distribution but still mildly anti-conservative in the far tail for very rare features, and a regional Bonferroni threshold can sit exactly there (`0.05/20 = 0.0025`, against a measured 0.0025 at a nominal 0.001). `associate` warns when features below 1% minor frequency, or with fewer than 10 minor-allele carriers among cases, are tested. `mac_case` and `mac_ctrl` are reported for binary traits because total MAC hides the split that governs reliability — 1 case / 19 controls is far weaker than 10 / 10. Keep `--min-maf 0.01` unless you accept that sub-threshold results are exploratory.

The kinship matrix is validated before use — row count and row width (a ragged matrix previously indexed past the end of a short row), finiteness, symmetry, and positive semi-definiteness. A matrix failing any of these is not a GRM, and the variance ratio the LMM derives from it would not mean anything.

The score test is saddlepoint-corrected past `|z| > 2`. The normal approximation matches only the first two cumulants of the score, which is why it drifts in the far tail exactly when the terms are skewed — a rare variant under case/control imbalance. The score is a sum of independent bounded terms, so its cumulant generating function is available exactly and the saddlepoint (Lugannani–Rice) expands about the point where the tilted distribution is centred on the observed value. An exact CGF does not make the result exact — Lugannani–Rice is still an approximation, and on the enumerated reference below it is about 1.4× the exact tail where the normal is 4.7× off in the opposite, dangerous direction. Below the cutoff the two agree to within printing precision and the normal is cheaper — the same gate SAIGE uses. The `p_method` column says which produced each p (`t`, `score`, `score_spa`, `lmm`).

At the edge of the score's support — a perfectly separating feature — the saddlepoint is at infinity and the expansion does not apply. Falling back to the normal tail there would use an approximation exactly where it is least trustworthy, so instead the probability is written down exactly: the boundary is a single Bernoulli configuration, `P(S = S_max) = Π_{g̃>0} μᵢ · Π_{g̃<0} (1−μᵢ)`, computed in logs. On the separation fixture that gives `0.00048828125 = 2·0.5¹²` where the normal tail said `0.000532`. Those rows are labelled `p_method=score_exact`. Without the boundary guard a root search finds a spurious solution far out and Lugannani–Rice returns 1.

Measured effect on the far tail, rare features included, over 1500 permutation replicates: type-I at `p < 0.001` falls from 0.0025 to 0.0017. Residual anti-conservatism remains; see below.

Reference case, from convolving the exact null distribution over 18 independent Bernoulli terms (4 carriers, 3 cases / 15 controls, `z = 3.55`):

| | two-sided p | ratio to exact |
|---|---|---|
| exact | 0.00182437858954 | 1.00 |
| saddlepoint | 0.00132 | 0.72 |
| normal | 0.000386 | 0.21 |

The two tails are evaluated independently, because an asymmetric score distribution puts them in different regimes: a threshold can be outside the support (probability exactly zero), on the boundary atom, or interior. Only when both tails are exact or provably zero is the result labelled `score_exact`; a mixture is reported as `score_spa`.

The saddlepoint is not uniformly better than the normal. It is an asymptotic expansion, and at very small `n` it need not win: on the asymmetric 16-sample fixture in the test suite the exact tail is `0.000184813`, the normal gives `0.000139` and the saddlepoint `0.000117` — further out, not closer. panvar's value there matches an independent Lugannani–Rice implementation to six significant figures, so that is the approximation's own error, not an implementation defect. The tests therefore pin the numerical result against that reference rather than asserting any superiority.

Residual limitations, measured. With 1500 replicates the tail is 1.7× nominal at `p < 0.001`, down from 2.5× before SPA but not at nominal. Under the parametric null — simulating the phenotype from the fitted covariate-only model, which unlike permutation preserves the genotype–covariate relationship — the rarest features (~18 carriers of 5705) show a genuinely non-uniform p-value distribution (median `p` 0.32 rather than 0.5), while their type-I error at 0.05 and 0.01 stays near nominal (0.033–0.060 and 0.003–0.017). So the distortion is in the middle of the distribution rather than at the decision threshold. A binary outcome always gives a discrete score distribution, so this is not proof of a continuous one; what the 300 distinct replicate p-values do rule out is that it is explained by obvious coarse ties.
