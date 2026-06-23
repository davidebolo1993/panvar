# A GWAS primer (for panvar `associate`)

This page explains, from scratch and with small worked examples, the ideas behind a genome-wide association
study (GWAS) and how `panvar associate` implements them. It is meant to be read on its own: every concept is
introduced and then shown on a concrete case. The running thread is **LPA / KIV-2**, the copy-number trait;
[example.md](example.md) is the companion manual that runs the commands and reads the real numbers.

## The one question a GWAS asks

A GWAS repeats a single question, once per genetic marker:

> across a cohort of individuals, does the **genotype** at this marker track the **phenotype**?

The unit of analysis is the **individual (sample)** — one person. The **phenotype** is the trait, either
*quantitative* (a number, e.g. plasma Lp(a)) or *binary* (case/control). For each marker `associate` fits a
regression of phenotype on genotype, optionally adjusting for covariates:

```
phenotype  ~  genotype  +  Age + Sex + PC1 + PC2 + ...
```

and reports, for the **genotype** term, an effect size, its standard error, a z-statistic, and one p-value
(a Wald test). Quantitative traits use **linear** regression (the effect is `beta`); binary traits use
**logistic** regression (the effect is `log_or`, the log odds ratio). `associate --model auto` picks: a 0/1
phenotype triggers logistic, anything else linear.

## Genotype = dosage, and why the count matters

The genotype fed to the regression is a **dosage**: a number per individual per marker. For a SNP it is
0/1/2 — the number of copies of one allele. For a **copy-number** marker like KIV-2 it is simply the
**count**: how many repeat copies the individual carries (e.g. 28).

This distinction is the whole point for CNVs. A **presence/absence** test asks only "does this individual
carry the feature, yes/no?". That is fine for a SNP, but it **fails for KIV-2**: the repeat unit is present
in *everyone* (copy number ≥ 1), so presence/absence has frequency 1 — there is no contrast and nothing to
test, and the marker gets filtered out. All the signal is in the count. panvar therefore carries the true
per-individual dosage all the way through, and `associate` tests it directly; presence/absence is just the
special case where the dosage happens to be 0/1.

## Covariates: what you hold fixed

A covariate is something that affects the phenotype but isn't the marker you're testing — `Age`, `Sex`,
batch, ancestry. Putting it in the model estimates the genotype effect *holding the covariate fixed*. In
`associate` the covariates are simply the extra columns of the `--phenotype` table after the phenotype
column. Any sample with a missing (`NA`) phenotype **or** any missing covariate is dropped (`samples_used`
reflects this). The most important covariates are the ancestry principal components (see *population
structure* below).

## MAF filter: drop markers with no usable variation

If (almost) everyone has the same genotype, a marker carries (almost) no information and produces unstable
estimates, so GWAS routinely drops low-frequency markers. `associate` measures a **minor (non-modal)
frequency**: one minus the fraction of samples at the most common rounded dosage.

Example — rounded dosages across 10 samples `{2,2,2,2,2,2,3,3,4,1}`: the mode (2) occurs 6×, so
`minor_freq = 1 − 6/10 = 0.40`. With `--min-maf 0.02` this marker is kept; a marker where 99% share one
value would be dropped. Crucially this is computed on **your actual cohort**, not on the pangenome — a
variant can exist in the graph yet be invariant in *these* samples, and that is exactly what gets removed
(counted as `dropped_min_maf`).

## Multiple testing: not p<0.05, and not 5e-8 either

Run thousands of tests at raw `p < 0.05` and you expect thousands of false positives (5% of all true nulls).
You must correct for the number of tests. `associate` reports two standard corrections:

- **Bonferroni** — call significance at `p < 0.05 / n_tests`. Controls the chance of *any* false positive;
  conservative. Written as `bonferroni_threshold`, with a per-marker `p_bonf`.
- **Benjamini–Hochberg (BH) FDR** — the `q_bh` column. Controls the *expected fraction* of false positives
  among the markers you call significant. Less conservative; the usual choice when markers are correlated.

Example with 4 tested p-values `[5e-6, 2e-3, 0.03, 0.40]`: Bonferroni threshold `0.05/4 = 0.0125` keeps the
first two; BH q-values work out to `[2e-5, 0.004, 0.04, 0.40]`, so FDR<0.05 keeps the first **three** — BH
recovers one more, the expected conservative-vs-FDR trade-off.

The famous **5×10⁻⁸** is just Bonferroni for a *whole human genome* (~10⁶ independent common variants).
panvar tests a **region** (a handful of pangenome loci → hundreds to thousands of features), so `5e-8` is the
wrong scale; the right one is `0.05 / n_tests` over the features actually tested. Because many features
(k-mers from one node) are highly correlated, `n_tests` over-counts the independent tests, making this
Bonferroni conservative — a safe default.

## Population structure and genomic inflation

Suppose the cohort mixes two ancestry groups, and group B happens to have both a higher frequency of allele
X at many unrelated markers *and* a higher trait value, for reasons unrelated to X (diet, drift). A naive
GWAS then reports X — and every ancestry-tracking marker — as "associated", because genotype and phenotype
are correlated *through ancestry*. This is **confounding by population structure**, and it inflates results
everywhere.

You detect it with the **genomic-inflation factor λ**: compare the median test statistic to what the null
predicts. **λ ≈ 1** is well-calibrated; **λ > 1** signals inflation (structure, relatedness, miscalibration).
`associate` reports `lambda_gc`, and the QQ plot shows the same thing — points drifting above the diagonal.

> **λ needs nulls.** λ is only meaningful when most markers are *not* associated. A single pangenome region —
> LPA on its own — is essentially one correlated block of signal with no nulls, so its λ is **not**
> interpretable and structure correction is neither needed nor diagnostic there. For a single region the real
> safeguards are the MAF filter and region-wide multiple testing; λ, PCs and the LMM come into their own only
> on a **genome-wide** panel (many regions).

## Two ways to correct structure: PCs and the LMM

**Principal components (PCs).** Run a PCA on a genome-wide genotype matrix; the top PCs are axes of ancestry
(PC1 might separate the two groups). Add them as covariates and the regression conditions ancestry out. In
`associate`, either include `PC1…PCN` columns in the `--phenotype` table, or derive them with `--pca N` from
a kinship matrix. Cheap; works for linear and logistic.

**Linear mixed model (LMM).** Instead of a few PC axes, model relatedness *continuously* with a **kinship
matrix** `K`: the trait is `phenotype = fixed effects + u + noise`, where `u` is a random effect with
covariance `σ²_g · K`, so more-related individuals are expected to be more similar. `associate --model lmm`
fits this for a quantitative trait using the **EMMAX** speed trick — eigendecompose `K` once, estimate the
variance ratio `δ = σ²_e/σ²_g` once under the null, then test each marker by fast generalized least squares
in the rotated space (this is why panvar depends on Eigen). LMM handles gradual ancestry *and* cryptic
relatedness and is the genome-wide standard; PCs are the lighter alternative. You can use both at once.

Rule of thumb: **PCs** for discrete ancestry and speed; **LMM** for relatedness/continuous structure, or when
PCs leave λ above ~1.05.

## Kinship (the GRM)

A **kinship matrix** (genetic-relationship matrix, GRM) is an `n × n` table whose entry `(i, j)` measures how
genetically similar individuals `i` and `j` are, from genome-wide markers: standardize each marker into `Z`
(subtract mean, divide by sd), then `K = Z Zᵀ / m`. The diagonal is self-similarity; large off-diagonals are
relatives or same-ancestry pairs. The LMM uses `K` to know who should be correlated. `associate` gets it via
`--kinship <file>` (a precomputed GRM, ideally genome-wide — best practice) or `--make-kinship` (build it
from the genotype matrix being tested — convenient, but **only valid genome-wide**: on a single region `K` is
contaminated by the very signal you test and the LMM over-corrects; the tool warns).

## Reading the outputs

`<prefix>.assoc.tsv` — one row per tested feature, sorted by p; the `bubbles`/`nodes`/`gene` columns trace a
hit back to the graph and (with `--node-genes`) to a gene name. `<prefix>.summary.tsv` — the run's settings
and diagnostics; check **`lambda_gc`** first. The **Manhattan** plot shows −log10 p (before correction) and
−log10 q (after), so significance survives or collapses visibly; the **QQ** plot shows calibration with λ.
Column-by-column definitions are in [modules/associate.md](../modules/associate.md#outputs).

## The shortest possible recipe

```bash
# 1) genotypes from describe (per-sample dosage for a diploid cohort)
panvar describe -i graph.gfa --bubble-prefix-in pan --out-dir desc \
  --variant-nodes call.variant_nodes.tsv --samples cohort.samples.tsv --no-wide-matrix

# 2) association: quantitative trait, MAF filter, ancestry PCs already in the phenotype table
panvar associate \
  --genotypes desc/bimbam_graph.samples.bimbam.gz \
  --samples   desc/bimbam.samples.samples.txt.gz \
  --feature-annot desc/feature_annot.samples.tsv.gz \
  --node-genes call.node_genes.tsv \
  --phenotype pheno.tsv --min-maf 0.02 -o out

# 3) (genome-wide-like panel) control structure with an LMM, or with PCs
panvar associate ... --model lmm --kinship grm.tsv -o out_lmm
panvar associate ... --pca 5 --make-kinship       -o out_pca

# 4) plots
Rscript scripts/plot_associate.R --assoc out.assoc.tsv --summary out.summary.tsv --out out --title "my trait"
```

See [example.md](example.md) for the full LPA walk-through with real numbers, and
[modules/associate.md](../modules/associate.md) for the exhaustive flag and column reference.
