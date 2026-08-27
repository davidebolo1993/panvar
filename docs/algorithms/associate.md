# Module `associate` - algorithm

Mechanism for `associate`. For commands and outputs see [modules/associate.md](../modules/associate.md); references are listed in [references.md](../references.md#associate).

## How it works

### 1. Filter on minor frequency

Each row of the genotype matrix is one testing unit. Dosages are rounded, the most frequent dosage is identified, and the minor non-modal frequency is:

```text
minor_freq = 1 - frequency_of_modal_dosage
```

Units below `--min-maf` are not tested. This definition also works for multi-copy dosages, where a conventional biallelic MAF is not available.

### 2. Fit and test each unit

Every model tests the genotype coefficient in:

```text
phenotype ~ genotype + covariates
```

| model | fit and test | reason |
|-------|--------------|--------|
| linear | ordinary least squares; Student-t Wald p-value with residual degrees of freedom | residual variance is estimated from the cohort |
| logistic | covariate-adjusted Rao score p-value | remains usable when a rare genotype approaches separation |
| LMM | EMMAX-style generalized least squares after rotating by the kinship matrix | adjusts a quantitative trait for relatedness |

For binary traits, the score distribution receives extra tail handling:

- ordinary score p-value in the central range (`p_method=score`);
- saddlepoint approximation in a sufficiently extreme tail (`score_spa`);
- direct boundary probability when the observed score is at the edge of its support (`score_exact`).

The saddlepoint result is still an asymptotic approximation. At small sample sizes it is not guaranteed to be closer than the normal approximation.

Logistic effect estimation is separate from hypothesis testing. Maximum likelihood supplies `log_or` and `se` when it converges; Firth's penalised fit supplies a finite effect under separation. In both cases, `p` and `z` remain score-test quantities. A failed score test drops the unit rather than silently replacing it with a Wald p-value.

### 3. Correct for multiple testing

Three summaries are computed from the marginal p-values:

```text
p_bonf      = min(1, p * number_of_tests)
p_bonf_meff = min(1, p * Meff)
q_bh        = Benjamini-Hochberg FDR q-value
```

`p_bonf` is the conservative family-wise reference and `q_bh` is the main FDR summary. `p_bonf_meff` is a regional guide whose effective test count depends on the unit:

| unit | `Meff` |
|------|--------|
| variant | phenotype-blind Li-Ji estimate from genotype-correlation eigenvalues; LD-clump count is a fallback |
| feature | number of annotated bubble groups plus unannotated singleton features |

Variant LD clumps are still reported for interpretation. Variants are considered in p-value order; a non-low-frequency lead claims correlated variants above `--ld-r2` as shadows. Because lead selection uses phenotype p-values, the clump count is not used as the primary statistical estimate when Li-Ji is available.

### 4. Distinguish independent signals from shadows

Multiple-testing correction controls the number of claims; it does not show whether several hits represent the same underlying signal. `associate` therefore performs a separate conditional stage.

Variant mode uses forward-stepwise conditional analysis:

1. start with no selected variants;
2. test each remaining variant while including selected dosages as covariates;
3. add the smallest conditional p-value if it passes `--cojo-p`;
4. repeat until no variant enters.

Selected rows receive `cond_role=signal`; other rows are `shadow`. Each reported `p_conditional` conditions on the selected set excluding the target itself.

Feature mode uses the top marginal feature as a single lead. A feature with r² above 0.95 to that lead is labelled `collinear` and is not assigned a degenerate conditional p-value. Other features are tested with the lead dosage as an added covariate and labelled `conditioned`.

The conditional fit uses complete cases shared by the target and conditioning set, so `n_conditional` can differ from the marginal `n`.

### 5. Mixed-model adjustment

For `--model lmm`, the external kinship matrix `K` is validated, eigendecomposed once and used to rotate the phenotype, covariates and each genotype. A null-model variance ratio is estimated once; every feature is then tested by weighted least squares in the rotated space. This gives one cubic eigendecomposition followed by inexpensive per-feature fits.

`--pca N` is a lighter alternative: the top N eigenvectors of the same matrix are added as ordinary covariates, without fitting a random effect.

## Worked trace

The steps below follow the five above, one for one. Variant mode, a quantitative phenotype, five
variants over eight samples. Dosages are diploid alt counts:

```text
sample   V1  V2  V3  V4  V5   phenotype
s1        0   0   0   1   0        1.9
s2        0   0   1   0   0        2.1
s3        1   1   0   0   0        4.0
s4        1   1   1   1   0        4.2
s5        2   2   0   0   0        6.1
s6        2   2   1   0   0        6.0
s7        0   0   1   1   0        2.0
s8        1   1   0   0   1        4.1
```

`V2` tracks `V1` exactly. `V3` and `V4` vary independently of both. `V5` is carried by one sample.

1. Filter on minor frequency. `V5`'s minor genotype appears in 1 of 8 samples, a frequency of 0.0625
   in dosage terms but 1 carrier in count terms; under `--min-maf 0.01` it is retained, and under a
   stricter threshold it would be dropped and counted in `dropped_min_maf`. The other four pass.

2. Fit and test each unit. The phenotype is continuous, so `--model auto` selects linear and each
   variant gets `phenotype ~ dosage`, a Wald `t` test, and `beta` with its `se`. `V1` tracks the
   phenotype almost perfectly and is strongest; `V2` is identical to `V1` and scores the same; `V3`
   and `V4` are close to null.

3. Correct for multiple testing. Four variants were tested, so raw Bonferroni multiplies by 4 and its
   threshold is `0.05 / 4 = 0.0125`. Because `V1` and `V2` are perfectly correlated, Li-Ji estimates
   `meff = 3` rather than 4, and `p_bonf_meff` multiplies by 3 instead. `q_bh` is the
   Benjamini-Hochberg value and is the summary to read first.

4. Distinguish independent signals from shadows. Forward selection enters `V1`, then re-tests
   everything conditional on it. `V2` adds nothing once `V1` is in the model, so it is a `shadow` with
   a conditional p near 1. `V3` is unaffected by conditioning and stays a `signal`. Marginal
   significance and independence are visibly different columns.

5. Mixed-model adjustment. This trace used the linear model. With `--model lmm` and a kinship matrix,
   relatedness enters the fit in step 2 while every column above keeps its meaning; `p_conditional` is
   `NA`, because conditional analysis is not implemented for mixed-model rows.

The written rows:

```text
feature_id  layer    n  minor_freq   beta     se       z      p        p_bonf   p_bonf_meff  q_bh     p_conditional  cond_role
V1          variant  8       0.5000  2.100  0.061  34.43  3.4e-08  1.4e-07      1.0e-07  1.4e-07         2.9e-08  signal
V2          variant  8       0.5000  2.100  0.061  34.43  3.4e-08  1.4e-07      1.0e-07  1.4e-07            0.94  shadow
V3          variant  8       0.5000 -0.075  0.802  -0.09     0.93     1.00         1.00     0.93            0.91  signal
V4          variant  8       0.3750 -0.100  0.889  -0.11     0.91     1.00         1.00     0.93            0.88  signal
V5          variant  8       0.1250     .      .      .        .        .            .        .               .  .
```

`V2` shows why the conditional columns exist: it is as significant as `V1` on every marginal summary
and carries no independent evidence at all.
