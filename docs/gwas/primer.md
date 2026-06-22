# A GWAS primer (for panvar `associate`)

This page explains, from scratch, the ideas behind a genome-wide association study (GWAS) and how
`panvar associate` implements them. It is written for someone who has the genotypes from `describe` and a
phenotype table and wants to *understand* what each knob does — not just run it. The running example is
**LPA / KIV-2** (the same one in [example.md](example.md)); read that for the worked numbers.

---

## 1. What a GWAS is

A GWAS asks one question, repeated once per genetic marker:

> across a cohort of individuals, does the **genotype** at this marker correlate with the **phenotype**?

- **Individual (sample):** one person. The unit of analysis.
- **Phenotype:** the trait. Either **quantitative** (a number, e.g. plasma Lp(a)) or **binary** (case/control).
- **Genotype / dosage:** a number per individual per marker. For a SNP it is 0/1/2 (copies of an allele); for
  a **copy-number** marker like KIV-2 it is the **count** (how many repeat copies the individual carries).
  panvar keeps the *raw count* — this is the whole point for CNVs (see [§3](#3-dosage-vs-presenceabsence)).

For each marker `associate` fits a regression of phenotype on genotype, with optional covariates:

```
phenotype  ~  genotype  +  Age + Sex + PC1 + PC2 + ...
```

and reports, for the **genotype** term: the effect size (`beta` for quantitative, `log_or` = log odds ratio
for binary), its standard error, a z-statistic, and a p-value (Wald test). One p-value per marker.

**Quantitative → linear regression. Binary → logistic regression.** `associate` picks automatically
(`--model auto`); a binary 0/1 phenotype triggers logistic, anything else linear.

---

## 2. Covariates: things you adjust for

A covariate is a variable that affects the phenotype but isn't the marker you're testing — `Age`, `Sex`,
ancestry components, batch. You put them in the model so the genotype effect is estimated *holding them
fixed*. In `associate` the covariates are simply the extra columns of the `--phenotype` table after the
phenotype column. A sample with a missing (`NA`) phenotype **or** any missing covariate is dropped from the
fit (you'll see `samples_used` reflect this).

The most important covariates in a GWAS are the **ancestry principal components** (PC1, PC2, …) — see
[§6](#6-population-structure-the-quiet-confounder).

---

## 3. Dosage vs presence/absence

A **presence/absence** test asks only "does this individual carry the feature, yes/no". That is fine for a
SNP or a rare indel. It **fails for a copy-number variant** like KIV-2: the repeat unit is present in
*everyone* (copy number ≥ 1), so presence/absence has frequency 1 — no contrast, nothing to test, and the
marker is filtered out. All the information is in the **count**. panvar therefore carries the true
per-individual dosage all the way through, and `associate` tests it directly. Presence/absence is just the
special case where the dosage happens to be 0/1.

---

## 4. MAF filter: drop markers with no usable variation

**Minor allele frequency (MAF)** is the frequency of the less-common genotype. A marker where (almost)
everyone has the same genotype carries (almost) no information and produces unstable, untrustworthy
estimates — so GWAS routinely drops low-MAF markers. `associate` uses a **minor (non-modal) frequency**:
`1 − (count of the most common rounded dosage) / n`. It is computed on **your actual cohort**, not on the
graph — a variant can exist in the pangenome but be invariant in *these* samples, and that is exactly what
this removes. `--min-maf 0.02` (say) drops anything where the minor genotype is below 2%. Reported as
`dropped_min_maf`.

---

## 5. Multiple testing: why not p<0.05, and why not 5e-8 either

You run thousands of tests, so a raw `p < 0.05` would give thousands of false positives (5% of all nulls).
You must correct for the number of tests. Two standard corrections, both reported:

- **Bonferroni:** declare significance at `p < 0.05 / n_tests`. Controls the chance of *any* false positive;
  conservative. `associate` writes this threshold as `bonferroni_threshold` and a per-marker `p_bonf`.
- **Benjamini–Hochberg (BH) FDR:** the `q_bh` column. Controls the *expected fraction* of false positives
  among the hits you call. Less conservative; the usual default when markers are correlated.

The famous **5×10⁻⁸** threshold is *Bonferroni for a whole human genome* (~10⁶ independent common variants).
panvar tests a **region** (one or a few pangenome loci → hundreds–thousands of features), so `5e-8` is the
wrong scale: use `0.05 / n_tests` over the features actually tested. Because many features (e.g. k-mers from
one node) are highly correlated, `n_tests` over-counts the *independent* tests, so this Bonferroni is
**conservative** — a safe default.

---

## 6. Population structure: the quiet confounder

Imagine your cohort mixes two ancestry groups. Group B happens to have, on average, both (a) a higher
frequency of allele X at many unrelated markers and (b) a higher trait value — for reasons that have nothing
to do with allele X (diet, environment, drift). A naive GWAS will report allele X (and every other
ancestry-tracking marker) as "associated", because genotype and phenotype are correlated *through ancestry*.
This is **confounding by population structure**, and it inflates results everywhere.

You detect it with the **genomic-inflation factor λ (lambda)**: take all the test statistics, and compare
their median to what you'd expect under the null. λ ≈ 1 means well-calibrated; **λ > 1 means inflation**
(structure, relatedness, or other miscalibration). `associate` reports `lambda_gc` in the summary, and the
QQ plot shows the same thing visually (points drifting above the diagonal).

> **λ needs nulls.** λ is only meaningful when most markers are *not* associated (true nulls). A single
> pangenome region — e.g. LPA on its own — is essentially **one correlated block of signal with no nulls**,
> so its λ is not interpretable, and structure correction is neither needed nor diagnostic there. For a
> single region, the real safeguards are the **MAF filter** and **region-wide multiple testing** above. λ,
> PCs, and the LMM come into their own once you test a **genome-wide** panel of features (many regions).

---

## 7. Two ways to correct structure: PCs and the LMM

**Principal components (PCs).** Run a PCA on a genome-wide genotype matrix; the top PCs are axes of ancestry
variation (PC1 might separate the two groups). Add them as **covariates** and the regression conditions
ancestry out. In `associate`: either put `PC1…PCN` columns in the `--phenotype` table (if you computed them
elsewhere), or let the tool derive them with `--pca N` from a kinship matrix. Cheap; works for linear and
logistic.

**Linear mixed model (LMM).** Instead of a few PC axes, model relatedness *continuously* with a **kinship
matrix** `K` (next section). The trait is `phenotype = fixed effects + u + noise`, where `u` is a random
effect with covariance `σ²_g · K`: individuals who are more related are expected to be more similar.
`associate --model lmm` fits this for a quantitative trait using the standard **EMMAX** speed trick:
eigendecompose `K` **once**, estimate the variance ratio `δ = σ²_e/σ²_g` **once** under the null, then test
every marker by fast generalized least squares in the rotated space. (This is why we depend on the Eigen
linear-algebra library.) LMM handles both gradual ancestry and cryptic relatedness, and is the genome-wide
standard; PCs are the lighter-weight alternative. You can use both at once (PCs as covariates *and* an LMM).

Rule of thumb: **PCs** for discrete ancestry and speed; **LMM** for relatedness/continuous structure or when
PCs leave λ above ~1.05.

---

## 8. Kinship (the GRM)

A **kinship matrix** (a.k.a. genetic-relationship matrix, GRM) is an `n × n` table where entry `(i, j)`
measures how genetically similar individuals `i` and `j` are, estimated from genome-wide markers: standardize
each marker (subtract mean, divide by sd) into `Z`, then `K = Z Zᵀ / m`. Diagonal ≈ self-similarity, large
off-diagonals ≈ relatives or same-ancestry pairs. The LMM uses `K` to know who should be correlated.

`associate` gets `K` in one of two ways:

- `--kinship <file>` — a precomputed GRM (ideally from **genome-wide** markers). Best practice.
- `--make-kinship` — build it from the genotype matrix you're testing. Convenient, but **only valid when that
  matrix is genome-wide-like**: if it's a single region, `K` is contaminated by the very signal you test (it
  "knows" the causal genotype) and the LMM will over-correct. The tool warns; prefer an external GRM.

---

## 9. Reading the outputs

`<prefix>.assoc.tsv` — one row per tested feature, sorted by p:
`feature_id, layer, bubbles, nodes, n, minor_freq, beta|log_or, se, z, p, p_bonf, q_bh, gene`. The
`bubbles`/`nodes`/`gene` columns trace a hit back to the graph and (with `--node-genes`) to a gene name.

`<prefix>.summary.tsv` — the run's settings and diagnostics: `model`, `samples_used`, `features_tested`,
`dropped_min_maf`, `bonferroni_threshold`, significant counts, and **`lambda_gc`** (and `lmm_delta` for an
LMM). Check `lambda_gc` first.

**Manhattan plot** (`scripts/plot_associate.R`): x = graph/genomic order, y = −log10 p; the nominal and
Bonferroni lines are drawn and FDR/Bonferroni-significant points highlighted — so "before vs after
correction" is one figure. **QQ plot**: observed vs expected p-values; on-diagonal = calibrated, an early
lift-off = inflation, annotated with λ.

---

## 10. The shortest possible recipe

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

See [example.md](example.md) for the full LPA walk-through (with real numbers), and
[modules/associate.md](../modules/associate.md) for the exhaustive flag reference.
