# associate — algorithm & worked example

Mechanism and a hand-traced example for **Module 5**. For usage/flags see [modules/associate.md](../modules/associate.md);
for the GWAS concepts (MAF, multiple testing, λ, kinship, PCs vs LMM) in plain terms see
[gwas/primer.md](../gwas/primer.md); citations: [references.md](../references.md#associate).

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

BH recovers one more than Bonferroni — the expected conservative-vs-FDR trade-off. Because correlated
features (k-mers from one node) over-count the independent tests, the region-wide Bonferroni here is
conservative; this is *not* the genome-wide `5e-8` (see [gwas/primer.md](../gwas/primer.md#5-multiple-testing-why-not-p005-and-why-not-5e-8-either)).

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
eigenvectors of `K` as fixed covariates to the GLM, with no variance-component step. See
[gwas/primer.md](../gwas/primer.md#7-two-ways-to-correct-structure-pcs-and-the-lmm) for when to use which.
