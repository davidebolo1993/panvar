# Module `associate` - algorithm

Mechanism for the `associate` module. For usage/flags see [modules/associate.md](../modules/associate.md); References in [references.md](../references.md#associate).

`associate` tests each genotype feature against a phenotype and then corrects the result for the two problems a pangenome locus creates: many features read the same underlying variation, so the raw test count over-states the multiple-testing burden, and a strong signal casts marginal shadows onto features that merely tag it. The tested thing is the unit — a feature (a single k-mer or node/edge dosage) or a variant (one structural-variant, SV, call) — and the unit sets how both corrections work. Correction runs in two layers: a phenotype-blind burden layer that counts how many effectively independent tests the region holds, and a phenotype-aware independence layer that refits each hit conditioned on the lead signal(s).

## How it works

### 1. Filter on minor frequency

A unit is dropped when its minor (non-modal) frequency — the fraction of the cohort not carrying the most common rounded dosage — is below `--min-maf`. This is the minor-allele-frequency (MAF) quantity generalized to a multi-valued copy-number dosage, where a two-allele MAF is not defined: round each dosage, take the modal value's share `f`, and the minor frequency is `1 − f`.

### 2. Fit and test each unit

Per unit, `associate` fits `phenotype ~ genotype + covariates` and tests the genotype term. The reported effect size is the maximum-likelihood estimate, except where that fit diverges under separation, when Firth's penalised likelihood supplies a finite one instead. Quantitative traits get a Wald test with a Student-t tail on `n − p` degrees of freedom — `σ²` is estimated from the same residuals, so the statistic is t-distributed, and the normal tail is anti-conservative (at 9 df a `z` of 20.4 is `p = 7.6e-09` under t and `1.3e-92` under the normal). Binary traits get a covariate-adjusted Rao score test instead, because the Wald statistic collapses under near-separation; see [modules/associate.md](../modules/associate.md#calibration). For those, `z` and `p` come from the score test while `log_or` and `se` remain the maximum-likelihood effect size, so `p` is not recoverable from them.

- linear (quantitative trait) — ordinary least squares (OLS); `β = (XᵀX)⁻¹Xᵀy`, `Var(β) = σ²·(XᵀX)⁻¹` with `σ² = RSS/(n−p)` (RSS is the residual sum of squares).
- logistic (binary trait) — iteratively reweighted least squares (IRLS) Newton–Raphson; `β`, `Var(β) = (Xᵀ W X)⁻¹` at convergence.
- p — `z = β_genotype / se`, and for a quantitative trait `p = 2·P(T_{n−p} > |z|)` via the regularized incomplete beta (validated against R's `pt()` to a relative error below 3.3e-12). The normal form `erfc(|z|/√2)` — the upper tail computed directly to avoid the cancellation of `1 − Φ(|z|)` at large `|z|` — is what the score test uses for its chi-square-on-1-df tail. Floored at `1e-300`.

Which test is used, and why it is not always the obvious one:

- Quantitative traits use a Student-t tail on `n − p` degrees of freedom rather than the normal one. The residual variance is estimated from the same data, so the statistic is t-distributed; the difference is invisible at large sample sizes and substantial at small ones, where the normal tail is anticonservative by a wide margin.
- Binary traits use a Rao score test rather than a Wald test. A Wald statistic divides an estimate by its own standard error, and for a rare variant in an unbalanced case/control study the fit approaches separation: the coefficient grows and its standard error grows faster, so the statistic collapses toward zero exactly where the evidence is strongest. The score test never fits the alternative, so it has no standard error to inflate.
- Past `|z| > 2` the binary score is saddlepoint-corrected. The normal approximation matches only the first two cumulants, so it drifts in the far tail precisely when the terms are skewed, which is the rare-variant-under-imbalance case. The score is a sum of independent bounded terms, so its cumulant generating function is available exactly and the saddlepoint expansion is built from it. That does not make the result exact: the expansion is still an approximation, and at very small sample sizes it is not guaranteed to beat the normal.
- At the edge of the score's support, where a feature separates the outcome perfectly, the saddlepoint is at infinity and the expansion does not apply. Falling back to the normal tail there would use an approximation exactly where it is least trustworthy, so the probability is instead written down exactly: the boundary is a single configuration whose probability is a product over the samples, computed in logs. The two tails are evaluated independently, since an asymmetric score distribution can put them in different regimes, and the result is labelled exact only when both are.

Under separation the maximum-likelihood fit diverges, so an iteration limit is treated as failure rather than returning the last iterate. The score test is unaffected and still yields a valid p, and the effect size comes from Firth's penalised likelihood, which is finite where maximum likelihood is not. A feature is dropped only when the score test itself fails; there is never a silent fall back to a Wald p.

Consequence for reading a row: `z` and `p` come from the score test while `log_or` and `se` are a Wald-style effect size, so `p` cannot be recovered from them.

The unit is auto-detected from the `layer` column of `feature_annot` — `variant` when the majority of rows are variant-level, otherwise `feature` (override with `--unit`). Across all tested units the genomic-inflation factor `λ = median(z²)/0.4549` (the observed median z² over the null median of a χ²₁, `0.4549`) is reported in the summary: `λ ≈ 1` means the test is calibrated, `λ > 1` flags residual structure.

### 3. Correct the testing burden — `Meff` (phenotype-blind)

### 4. Adjust the p-values

`p_bonf` is the raw Bonferroni correction over every feature tested, `p_bonf_meff` scales by `Meff` instead, and `q_bh` is the Benjamini-Hochberg false-discovery rate. The first and last carry formal guarantees; the middle one is a regional guide, since `Meff` is a heuristic.


`Meff` is the effective number of independent tests in the region, always ≤ the raw feature count; it replaces the raw `features_tested` in the Bonferroni benchmark. How it is derived depends on the unit — the two are corrected by different mechanisms:

| unit | substrates | how `Meff` is computed | LD-clumping? | conditioning | per-row columns |
|------|-----------|------------------------|--------------|--------------|-----------------|
| feature | k-mers and graph (nodes/edges) | `Meff` = number of distinct bubbles the tested features map to (the `bubbles` column, split on `;`, de-duplicated) | no | single-lead + collinearity guard | `p_conditional`, `cond_role` |
| variant | SV calls (`describe --variant-vcf`) | `Meff` = Li–Ji eigenvalue estimate (phenotype-blind) | yes (genotype r²) | forward-stepwise (COJO) | `clump`, `is_lead`, `low_af`, `p_conditional`, `cond_role` |

Feature unit — membership-based `Meff`. k-mer counts and node/edge dosages are correlated by construction: every feature carries the `bubbles` it came from, and all the features inside one bubble are reads of the same local variation. So no genotype correlations are needed — each bubble collapses to one effective test and `Meff` = the number of distinct bubbles. There is no clumping and no r² in feature mode; it is a cheap set-count over the `bubbles` column, so it needs no genotype matrix. Both the k-mer run and the graph run use exactly this rule.

Variant unit — genotype-r² LD-clumping. Distinct SV calls can still be in linkage disequilibrium (LD) across different bubbles, which bubble-membership cannot see, so here correlation is measured directly. Greedy clump: sort variants by p ascending; the best unassigned variant becomes a lead (`is_lead=1`, new `clump` id); every still-unassigned variant whose genotype r² (squared Pearson on mean-imputed dosage) with that lead exceeds `--ld-r2` (default 0.8) is marked its shadow (`is_lead=0`, same `clump`); repeat over the remaining variants. `Meff` is the number of leads. Only the variant tier retains the full dosage vectors needed for r², so clumping is intrinsically variant-only. `low_af` flags a variant whose observed minor-allele count (`minor_freq · n`) is below `--min-ac` (default 3); a low-AF variant has an unstable r² and a fragile asymptotic p, so it is barred from anchoring a clump (and thus cannot inflate `Meff`) — it can still be claimed as a shadow of a genuine lead.

`Meff` only rescales the Bonferroni benchmark; Benjamini–Hochberg (BH) false-discovery-rate (FDR) control is unchanged and stays primary:

```text
p_bonf      = min(1, p · n_tests)   # raw, over-conservative (ignores correlation)
p_bonf_meff = min(1, p · Meff)      # the honest effective-tests Bonferroni
q_bh                                # Benjamini–Hochberg FDR — primary control, no Meff
```

The summary reports `unit`, `meff`, both Bonferroni thresholds, and a unit-named alias for the effective count: `independent_variants` (variant) or `distinct_bubbles` (feature).

### 5. Establish independence — conditioning (phenotype-aware)

Clumping and `Meff` fix the threshold; they do not establish that a hit is independent. A marginal test of a variant only weakly correlated with an extremely strong locus still comes out genome-wide significant — even when its r² is far below `--ld-r2`, so clumping never groups it. A conditional refit is the only honest way to separate a true signal from such a shadow, so `associate` runs one on every tier and reports `p_conditional` (a unit's p — same test as the marginal one — after the lead signal(s) enter the model as covariates, restricted to samples observed for the target and every conditioning feature, reported as `n_conditional`) plus `cond_role` (its status). This step is deliberately separate from the phenotype-blind `Meff`, so the threshold cannot become circular.

Variant tier — forward-stepwise conditional-and-joint (COJO-style) analysis. Select a set of jointly-independent signals:

- compute each variant's p conditioned on the currently-selected set (its dosage(s) added as covariates; the genotype of interest stays the target of the same test used marginally);
- add the variant with the smallest conditional p if it clears the entry threshold `--cojo-p` (default `0.05/Meff`);
- repeat until nothing new clears the bar.

Then every variant is reported conditioned on the selected set minus itself: a true signal stays significant and is tagged `cond_role=signal`; a shadow's `p_conditional` inflates (e.g. a neighbour that only tags the lead goes from genome-wide significant to ~null) and is tagged `cond_role=shadow`. The sole signal of a single-signal locus has an empty conditioning set, so its `p_conditional` is `NA`. The summary reports `cojo_independent_signals`. Linear/logistic only; linear-mixed-model (LMM) conditioning (which needs the rotation) is not done.

Feature tier (k-mer / graph) — single-lead, with a collinearity guard. Features inside one bubble are ~perfectly correlated (a node, its self-loop edge, and the repeat-unit k-mers all read the same copy number), so conditioning them on the lead feature is numerically degenerate. `associate` therefore conditions on the single top feature and:

- a feature with genotype r² > 0.95 against that lead is same-event redundant → `cond_role=collinear`, `p_conditional=NA` (flagged, not scored — its degenerate p would be meaningless);
- every other (cross-bubble) feature gets a real `p_conditional` and `cond_role=conditioned` — if it merely tagged the lead's variant, it collapses;
- the lead itself is `cond_role=lead`, `p_conditional=NA`.

Full forward-stepwise across thousands of correlated features would be unstable, so COJO stays variant-level; the feature tier uses single-lead conditioning purely to expose cross-bubble shadows. Both passes stream the genotype file twice so no feature × sample matrix is held in memory.

The two correction layers in order:

```text
1. MAF FILTER drop minor (non-modal) freq < --min-maf
2. TEST       per feature/variant: phenotype ~ genotype + covariates(+PCs)  -> Student-t (quantitative) / score (binary) -> marginal p
── Layer A: multiple-testing BURDEN (phenotype-blind) ───────────────────────────────────
3. Meff       variant: LD-clump by genotype r² > --ld-r2 (low-AF can't lead) -> Meff = #leads
              feature: Meff = #distinct bubbles  (no clumping, no r²)
4. p_adj      p_bonf = p·n_tests | p_bonf_meff = p·Meff | q_bh = BH-FDR (primary)
── Layer B: INDEPENDENCE (phenotype-aware) ──────────────────────────────────────────────
5. CONDITION  variant: forward-stepwise COJO (entry --cojo-p, default 0.05/Meff) -> signal/shadow
              feature: condition on top lead; r²>0.95 vs lead -> collinear, else conditioned
              -> p_conditional, cond_role  (cojo_independent_signals in summary)
```

`scripts/plot_associate.R` draws three Manhattan panels — raw `-log10(p)`, BH `-log10(q)`, and (when the columns are present) `-log10(p_conditional)`, where shadows collapse below the line and only the conditioning signal(s) stay tall.

## Worked trace

Five samples cross one feature; genotype dosage `g` (a copy number), phenotype `y` on a log scale. Covariates are omitted so the arithmetic stays legible — they only add columns to `X`:

```text
sample   g   y
s1       2   1.3
s2       2   1.6
s3       4   1.0
s4       6   0.8
s5       6   0.5
```

1. Filter on minor frequency. Round the dosages to `{2,2,4,6,6}`; the modal value occurs twice of five, so the minor (non-modal) frequency is `1 − 2/5 = 0.60`. That clears `--min-maf`, so the feature is tested.
2. Fit and test the unit. Fit `y ≈ a + b·g` by ordinary least squares. With `ḡ = 4` and `ȳ = 1.04`, the sums of squares are `Sgg = Σ(g−ḡ)² = 16` and `Sgy = Σ(g−ḡ)(y−ȳ) = −3.20`, giving slope `b = Sgy/Sgg = −0.20` and intercept `a = ȳ − b·ḡ = 1.84`.
3. Test the slope. The fitted values `1.44, 1.44, 1.04, 0.64, 0.64` leave residuals `−0.14, 0.16, −0.04, 0.16, −0.14`, so `RSS = 0.092`, `σ² = RSS/(n−2) = 0.0307`, `Var(b) = σ²/Sgg = 0.00192` and `se = 0.0438`. Then `z = b/se = −4.57`, and on `n − p = 3` degrees of freedom `p = 2·P(T₃ > 4.57) ≈ 0.020` — the normal tail would say `5e-6`, which at this sample size is simply wrong. Covariates only widen `X` and change `Var(b_g) = σ²·[(XᵀX)⁻¹]_gg`; a binary trait swaps this OLS for the IRLS-reweighted logistic fit and reports a score test rather than a Wald one.

Result for this feature:

```text
minor_freq=0.60  beta=−0.20  se=0.0438  z=−4.57  p≈0.020  p_method=t
```

4. Correct the testing burden. One feature is not a region. Across, say, four tested features `p = [5e-6, 2e-3, 0.03, 0.40]`:

```text
Bonferroni threshold = 0.05/4 = 0.0125    significant: 5e-6, 2e-3            (2 features)
p_bonf = min(1, p·4)                      = [2e-5, 8e-3, 0.12, 1.0]
BH q (sorted, q_i = p_i·m/rank, monotone) = [2e-5, 0.004, 0.04, 0.40]
  FDR < 0.05                              significant: 5e-6, 2e-3, 0.03      (3 features)
```

5. Adjust the p-values and establish independence. BH recovers one more than Bonferroni — the expected conservative-vs-FDR trade-off. But when the four features tag the same underlying variation, the raw `n_tests = 4` over-counts and the region-wide Bonferroni is too conservative; `Meff` (step 3) replaces that 4 with the effective test count — the number of distinct bubbles for a feature run, or the number of LD-clump leads for a variant run — and `p_bonf_meff` rescales accordingly.

## LMM (EMMAX) — the fast mixed model

The kinship matrix is validated before any of this runs: row count and row width, finiteness, symmetry, and positive semi-definiteness. A matrix failing any of them is not a covariance, and the variance ratio derived from it would not mean anything.

For `--model lmm`, relatedness is a random effect with covariance `σ²_g·K` (the kinship `K`). The fixed-effect rotation is done once (EMMAX), then each feature is a cheap generalized least squares (GLS):

- Eigendecompose `K = U diag(d) Uᵀ` (one symmetric eigendecomposition, via Eigen).
- Rotate the phenotype and covariates: `ỹ = Uᵀy`, `X̃ = UᵀX_cov`.
- Estimate the variance ratio `δ = σ²_e/σ²_g` once under the null by maximizing the restricted-maximum-likelihood (REML) profile likelihood over `δ` (1-D grid + golden-section), where row `i` has variance `(d_i + δ)`.
- Per feature, rotate the genotype `g̃ = Uᵀg`, append it to `X̃`, and solve a weighted least squares with weights `1/(d_i + δ)`; the Wald test on the genotype coefficient gives `β, se, z, p` exactly as above. `lmm_delta` in the summary is the fitted `δ`.

This costs one `O(n³)` eigendecomposition plus `O(n·p²)` per feature — orders of magnitude cheaper than re-fitting a full mixed model per feature. `--pca N` is the lighter alternative: it adds the top-N eigenvectors of `K` (principal components, PCs) as fixed covariates to the generalized linear model (GLM), with no variance-component step. Both need an external genome-wide `K` (`--kinship`); panvar is local and does not build one from its own region genotypes.
