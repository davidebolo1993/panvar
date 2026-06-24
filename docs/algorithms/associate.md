# associate — algorithm & worked example

Mechanism and a hand-traced example for **Module 5**. For usage/flags see [modules/associate.md](../modules/associate.md);
for the GWAS concepts (MAF, multiple testing, λ, kinship, PCs vs LMM) explained in context see
[gwas/example.md](../gwas/example.md); citations: [references.md](../references.md#associate).

## What it computes

Per feature, `associate` fits `phenotype ~ genotype + covariates` and reports a Wald test on the genotype
term:

- **linear** (quantitative trait) — ordinary least squares; `β = (XᵀX)⁻¹Xᵀy`,
  `Var(β) = σ²·(XᵀX)⁻¹` with `σ² = RSS/(n−p)`.
- **logistic** (binary trait) — IRLS Newton–Raphson; `β`, `Var(β) = (Xᵀ W X)⁻¹` at convergence.
- **Wald p** — `z = β_genotype / se`, `p = erfc(|z|/√2)` (the upper tail computed directly via `erfc` to
  avoid the catastrophic cancellation of `1 − Φ(|z|)` at large `|z|`), floored at `1e-300`.

Then: a **minor (non-modal) frequency** MAF filter on the cohort genotypes, **Bonferroni** (`0.05/n_tests`)
and **Benjamini–Hochberg** FDR over the tested features, and the **genomic-inflation factor**
`λ = median(z²)/0.4549`. Code: [`src/associate_command.cpp`](../../src/associate_command.cpp).

## Worked trace — one quantitative feature

Five samples; genotype dosage `g`, phenotype `y = log10 Lp(a)` (covariates omitted here for arithmetic — they
just add columns to `X`):

```text
sample   g   y
s1       2   1.3
s2       2   1.6
s3       4   1.0
s4       6   0.8
s5       6   0.5
```

**MAF (minor non-modal frequency).** Rounded dosages `{2,2,4,6,6}`: the most common value (mode) occurs 2×,
so `minor_freq = 1 − 2/5 = 0.60` ≥ `--min-maf` → tested.

**OLS fit** (`y ≈ a + b·g`): `ḡ=4`, `ȳ=1.04`; `Sgg=Σ(g−ḡ)²=16`, `Sgy=Σ(g−ḡ)(y−ȳ)=−3.20`.

```text
b = Sgy/Sgg = −3.20/16 = −0.20          a = ȳ − b·ḡ = 1.04 + 0.80 = 1.84
fitted ŷ = 1.84 − 0.20·g  → 1.44, 1.44, 1.04, 0.64, 0.64
residuals → −0.14, 0.16, −0.04, 0.16, −0.14   RSS = 0.092
σ² = RSS/(n−2) = 0.092/3 = 0.0307     Var(b) = σ²/Sgg = 0.0307/16 = 0.00192     se = 0.0438
z = b/se = −0.20/0.0438 = −4.57        p = erfc(4.57/√2) = erfc(3.23) ≈ 5.0e-6
```

So this feature: `β = −0.20` (more copies → lower log10 Lp(a), the inverse KIV-2 direction), `p ≈ 5e-6`.
With covariates the only change is `X = [1, g, Age, Sex, …]` and `Var(β_g) = σ²·[(XᵀX)⁻¹]_gg`; logistic
replaces OLS with the IRLS-reweighted version.

**Multiple testing** across, say, 4 tested features `p = [5e-6, 2e-3, 0.03, 0.40]`:

```text
Bonferroni threshold = 0.05/4 = 0.0125  → significant: 5e-6, 2e-3        (2 features)
p_bonf = min(1, p·4)                    = [2e-5, 8e-3, 0.12, 1.0]
BH q (sorted, q_i = p_i·m/rank, monotone) = [2e-5, 0.004, 0.04, 0.40]
  FDR < 0.05 → significant: 5e-6, 2e-3, 0.03                              (3 features)
```

BH recovers one more than Bonferroni — the expected conservative-vs-FDR trade-off. But the bigger problem in
a pangenome locus is that the tests are not independent: many of them tag the *same* underlying variant, so
the raw `n_tests` over-counts and the region-wide Bonferroni is too conservative. `associate` corrects this
with an effective number of independent tests, `Meff`, and how `Meff` is derived depends on the testing
**unit** — auto-detected from the `layer` column of `feature_annot` (`variant` when the majority of rows are
variant-level, otherwise `feature`; override with `--unit`).

### How the local correction is applied

The two units are corrected by **different** mechanisms — this is the key distinction:

| unit | substrates | how `Meff` is computed | LD-clumping? | conditioning | per-row columns |
|------|-----------|------------------------|--------------|--------------|-----------------|
| **feature** | **k-mers *and* graph** (nodes/edges) | `Meff` = number of **distinct bubbles** the tested features map to (the `bubbles` column, split on `;`, de-duplicated) | **no** | single-lead + collinearity guard | `p_conditional`, `cond_role` |
| **variant** | SV calls (`describe --variant-vcf`) | `Meff` = number of **LD-clump leads** | **yes** (genotype r²) | forward-stepwise (COJO) | `clump`, `is_lead`, `low_af`, `p_conditional`, `cond_role` |

**Feature unit — membership-based `Meff` (this covers both k-mers and graph features, identically).** k-mer
counts and node/edge dosages are correlated *by construction*: every feature carries the `bubbles` it came
from, and all the features inside one bubble are reads of the same local variation. So we don't need genotype
correlations at all — we collapse each bubble to one effective test and set `Meff` = the number of distinct
bubbles. There is **no clumping and no r² computation** in feature mode; it is a cheap set-count over the
`bubbles` column, so it needs no genotype matrix. Both the k-mer run and the graph run use exactly this rule.

**Variant unit — genotype-r² LD-clumping.** Distinct SV calls can still be in linkage disequilibrium across
*different* bubbles, which bubble-membership can't see, so here we measure correlation directly. Greedy clump:
sort variants by p ascending; the best unassigned variant becomes a **lead** (`is_lead=1`, new `clump` id);
every still-unassigned variant whose genotype r² (squared Pearson on mean-imputed dosage) with that lead
exceeds `--ld-r2` (default 0.8) is marked its **shadow** (`is_lead=0`, same `clump`); repeat over the
remaining variants. `Meff` is the number of leads. Only the variant tier retains the full dosage vectors
needed for r², so clumping is intrinsically variant-only. `low_af` flags a variant whose observed
minor-allele count (`minor_freq · n`) is below `--min-ac` (default 3); a low-AF variant has an unstable r²
and a fragile asymptotic p, so it is **barred from anchoring a clump** (and thus cannot inflate `Meff`) — it
can still be claimed as a shadow of a genuine lead.

**What the threshold then does (both units).** `Meff` only rescales the Bonferroni benchmark; BH-FDR is
unchanged and stays primary:

```text
p_bonf      = min(1, p · n_tests)   # raw, over-conservative (ignores correlation)
p_bonf_meff = min(1, p · Meff)      # the honest effective-tests Bonferroni
q_bh                                # Benjamini–Hochberg FDR — primary control, no Meff
```

The summary reports `unit`, `meff`, both Bonferroni thresholds, and a unit-named alias for the effective
count: `independent_variants` (variant) or `distinct_bubbles` (feature).

### Conditional / joint analysis (independence) — beyond the threshold

Clumping and `Meff` fix the *threshold*; they do **not** establish that a hit is independent. The reason is
decisive: a marginal test of a variant that is only weakly correlated with an *extremely* strong locus still
comes out genome-wide significant — even when its r² is far below `--ld-r2`, so clumping never groups it. A
conditional refit is the only honest way to separate a true signal from such a shadow, so `associate` runs one
on every tier and reports `p_conditional` (the conditional Wald p) plus `cond_role` (its status). This is a
**phenotype-aware** step, deliberately separate from the phenotype-blind `Meff` (so the threshold can't become
circular).

**Variant tier — forward-stepwise (COJO-style).** Select a set of jointly-independent signals:

1. compute each variant's p **conditioned on the currently-selected set** (its dosage(s) added as covariates;
   the genotype of interest stays the Wald target);
2. add the variant with the smallest conditional p **if** it clears the entry threshold `--cojo-p`
   (default `0.05/Meff`);
3. repeat until nothing new clears the bar.

Then every variant is reported conditioned on the selected set **minus itself**: a true signal stays
significant and is tagged `cond_role=signal`; a shadow's `p_conditional` inflates (e.g. a neighbour that only
tags the lead goes from genome-wide significant to ~null) and is tagged `cond_role=shadow`. The sole signal of
a single-signal locus has an empty conditioning set, so its `p_conditional` is `NA`. The summary reports
`cojo_independent_signals`. Linear/logistic only; LMM conditioning (which needs the rotation) is not done.

**Feature tier (k-mer / graph) — single-lead, with a collinearity guard.** Features inside one bubble are
~perfectly correlated (node `4789`, its self-loop edge, and the repeat-unit k-mers all read the same copy
number), so conditioning them on the lead feature is numerically degenerate. We therefore condition on the
single top feature and:

- a feature with genotype r² **> 0.95** against that lead is **same-event redundant** → `cond_role=collinear`,
  `p_conditional=NA` (flagged, not scored — its degenerate p would be meaningless);
- every other (cross-bubble) feature gets a real `p_conditional` and `cond_role=conditioned` — if it merely
  tagged the lead's variant, it collapses;
- the lead itself is `cond_role=lead`, `p_conditional=NA`.

Full forward-stepwise across thousands of correlated features would be unstable, so COJO stays variant-level;
the feature tier uses single-lead conditioning purely to expose cross-bubble shadows. Both passes stream the
genotype file twice so no feature × sample matrix is held in memory.

### Step-by-step (the two correction layers, in order)

```text
1. TEST       per feature/variant: phenotype ~ genotype + covariates(+PCs)  -> Wald -> marginal p
2. MAF FILTER drop minor (non-modal) freq < --min-maf
── Layer A: multiple-testing BURDEN (phenotype-blind) ───────────────────────────────────
3. Meff       variant: LD-clump by genotype r² > --ld-r2 (low-AF can't lead) -> Meff = #leads
              feature: Meff = #distinct bubbles  (no clumping, no r²)
4. p_adj      p_bonf = p·n_tests | p_bonf_meff = p·Meff | q_bh = BH-FDR (primary)
── Layer B: INDEPENDENCE (phenotype-aware) ──────────────────────────────────────────────
5. CONDITION  variant: forward-stepwise COJO (entry --cojo-p, default 0.05/Meff) -> signal/shadow
              feature: condition on top lead; r²>0.95 vs lead -> collinear, else conditioned
              -> p_conditional, cond_role  (cojo_independent_signals in summary)
```

`scripts/plot_associate.R` draws three Manhattan panels — raw `-log10(p)`, BH `-log10(q)`, and (when the
columns are present) `-log10(p_conditional)`, where shadows collapse below the line and only the conditioning
signal(s) stay tall.

## LMM (EMMAX) — the fast mixed model

For `--model lmm`, relatedness is a random effect with covariance `σ²_g·K` (the kinship `K`). The
fixed-effect rotation is done once (EMMAX), then each feature is a cheap GLS:

1. **Eigendecompose** `K = U diag(d) Uᵀ` (one symmetric eigendecomposition, via Eigen).
2. **Rotate** the phenotype and covariates: `ỹ = Uᵀy`, `X̃ = UᵀX_cov`.
3. **Estimate the variance ratio** `δ = σ²_e/σ²_g` **once under the null** by maximizing the REML profile
   likelihood over `δ` (1-D grid + golden-section), where row `i` has variance `(d_i + δ)`.
4. **Per feature**, rotate the genotype `g̃ = Uᵀg`, append it to `X̃`, and solve a weighted least squares
   with weights `1/(d_i + δ)`; the Wald test on the genotype coefficient gives `β, se, z, p` exactly as
   above. `lmm_delta` in the summary is the fitted `δ`.

This costs one `O(n³)` eigendecomposition plus `O(n·p²)` per feature — orders of magnitude cheaper than
re-fitting a full mixed model per feature. `--pca N` is the lighter alternative: it adds the top-N
eigenvectors of `K` as fixed covariates to the GLM, with no variance-component step. Both need an external
genome-wide `K` (`--kinship`); panvar is local and does not build one from its own region genotypes.
