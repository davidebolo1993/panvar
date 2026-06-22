# Pangenome association with panvar (LPA copy-number example)

This page is a self-contained example of running a **genome-wide association study (GWAS)** on
panvar's output, and — importantly — of *why* the **multiplicity** (copy count) panvar emits matters.

The running example is **LPA / KIV-2**, the classic copy-number trait: the *LPA* gene contains a tandem
repeat (the KIV-2 VNTR) whose **copy number is inversely associated with plasma Lp(a)** (more KIV-2 copies
→ lower Lp(a) → lower cardiovascular risk). It is the perfect case to show multiplicity, because the repeat
unit is present in **every** individual (copy number ≥ 1) — so the signal is *purely* in the count.

---

## 1. What a GWAS actually does

A GWAS asks, across a cohort of **individuals (samples)**:

> *does the genotype at position X correlate with the phenotype?*

For each variant it runs **one statistical test** — phenotype vs genotype across all samples — and collects
the p-values. Two things are fixed:

- the **unit of analysis is the sample** (an individual);
- the **genotype is a dosage**: 0/1/2 copies of an allele for a diploid SNP, or a **count** for a copy-number variant (CNV).

The output is a list of p-values, one per variant, usually:

- corrected for multiple testing (here we report **Benjamini-Hochberg q-values** and Bonferroni), and
- drawn as a **Manhattan plot** (x = genomic position, y = −log₁₀ q; peaks = associated loci) plus a **QQ plot** (observed vs expected p-values; the genomic-inflation factor λ should be ≈ 1 if the test is well-calibrated and there is no hidden population structure).

## 2. How this maps onto a pangenome

**Why bother with a pangenome here?** A SNP-array / linear-reference GWAS can only test what the single
reference represents — it is blind to sequence that is **absent from the reference** (large insertions,
paralog content) and, crucially, to **VNTR/STR copy number**. That is exactly where Lp(a) lives: the KIV-2
copy number is the major determinant, and a SNP panel cannot genotype it. A pangenome represents **all**
haplotypes, so this structural variation becomes a testable genotype — copy number read directly as marker
**multiplicity**, with every hit tracing back to a graph node and the variant `call` made there.

`panvar` turns the graph into association features as follows:

1. **`cosigt` genotypes a sequenced sample** → it assigns the sample a **pair of haplotype paths** from the pangenome (the two alleles it best matches at the locus).
2. **`describe`** counts features per **haplotype path** (per bubble) on two substrates: **k-mers** (the primary, read-queryable layer) and **node/edge dosage** (a complementary graph-local layer). A feature localizes to the graph nodes it comes from, which trace back to a called variant.
3. The **sample genotype** at a feature is the **aggregate over its two assigned haplotypes**: `sample_count = count(hapA) + count(hapB)`. For KIV-2 this is exactly the **diploid copy number** (CN_A + CN_B). `describe --samples <cosigt.tsv>` does this aggregation and writes a sample-level file per substrate.
4. **GWAS** = phenotype ~ sample genotype, one test per feature, across samples.

Three consequences worth internalizing:

- **It is variant-level, not whole-haplotype.** Every feature points back to a bubble/variant; we never test "the whole haplotype" as one blob.
- **Testing is reference-free.** counts exist for every haplotype whether or not the reference traverses that sequence. A sample whose haplotypes carry a bubble **not spanned by the reference** is still genotyped and tested.
- **The reference is only for plotting.** A Manhattan needs an x-axis; `plot_associate.R` uses **graph/node order** (pangenome-native, places *every* variant including reference-disjoint ones), with position-less k-mers parked at the left.

## 3. Why multiplicity (counts), not just presence/absence

A **presence/absence** association asks only "does this sample contain the feature (yes/no)". That works for
SNP- and indel-like variants. It **fails for CNVs** like KIV-2: the repeat-unit feature is present in
*everyone* (CN ≥ 1), so presence/absence has **no contrast** (frequency = 1) and the variant is filtered out
before testing. The information is entirely in the **count**. panvar carries the true per-sample dosage
through to the genotype, so a **dosage** test recovers KIV-2 that a presence/absence test cannot.

> **Engine.** `panvar associate` ([associate](modules/associate.md)) tests the **dosage** directly
> (`phenotype ~ genotype + covariates`), so copy-number loci are first-class. It can also adjust for
> population structure with ancestry **PCs** (`--pca`) or a **linear mixed model** (`--model lmm` + a
> kinship matrix). The same genotypes are exported by `describe` as **BIMBAM** mean-genotype dosage, which
> is also **GEMMA-ready**, so an external dosage-based tool can be used instead if desired. New to any of
> these terms? The [GWAS primer](gwas_primer.md) explains them from scratch on this example.

## 4. The pipeline

The post-`panphorte` modules consume the **panphorte-normalized/sorted graph** and the **panphorte** prefix.

- **`describe --samples <cosigt.tsv>`** (one sample per line: `sample <tab> hap1 <tab> hap2 …`; haplotype
  names must match graph path names) writes, per substrate, a **per-sample BIMBAM dosage** matrix —
  `bimbam_kmers.samples.bimbam.gz` (k-mers) and `bimbam_graph.samples.bimbam.gz` (node/edge) — whose value
  is the summed dosage over the sample's assigned haplotypes (a haplotype listed twice → counted twice =
  homozygous), with `bimbam.samples.samples.txt.gz` (column order) and `feature_annot.samples.tsv.gz`
  (layer/bubbles/nodes). See [describe](modules/describe.md).
- **`panvar associate`** ([associate](modules/associate.md)) fits `phenotype ~ genotype + covariates` per
  feature (linear for a quantitative trait, logistic for a binary one), filters features by **minor
  (non-modal) genotype frequency** computed on the actual cohort, and corrects for multiple testing
  **over the features actually tested** (region-wide Bonferroni `0.05/n_tests` + Benjamini-Hochberg FDR).
- **`scripts/plot_associate.R`** draws the Manhattan (−log10 p along graph order, with the nominal and
  Bonferroni threshold lines, FDR points highlighted — i.e. before/after correction in one figure) and the
  QQ with the genomic-inflation λ.

## 5. Data & literature resources (how the synthetic cohort is grounded)

The phenotypes are **simulated**, but their shape and effect sizes follow the Lp(a)/KIV-2 literature so the
example behaves like a real Southern-European study (we have no consented individual-level cohort to ship):

- **KIV-2 copy number** ranges ~**1 to >40** per allele, >95% heterozygous, and explains **40–70%** of Lp(a)
  variance; carriers of *small* isoforms (≤22 KIV repeats) have up to **~5×** higher Lp(a) than large-isoform
  carriers — the **inverse** dose effect. (Schmidt et al. 2016; Coassin & Kronenberg 2022.) Our real LPA
  graph carries copies **1–32** across 466 haplotypes.
- **Plasma Lp(a)** is strongly **right-skewed (log-normal)**, spanning ~**0.3–300 mg/dL** with a population
  **median ~10–12 mg/dL**, and is somewhat higher in Southern than Northern Europe; **>50 mg/dL** is the
  common clinical high-risk cut. (Coassin & Kronenberg 2022; the Italian **Moli-sani** cohort, Frontiers
  Cardiovasc Med 2025.) We use these as the simulation's `BASE`/spread and the binary case threshold.
- **Covariates** (`Age` ~ adult 35–85, `Sex`, ancestry `PC1–PC3`) mirror an adult cohort design.

Full citations are in [references.md](references.md). `tests/gwas/make_lpa_phenotype.py` encodes these as
explicit constants (`BASE_LOG10`, `SLOPE_LOG10`, `SIGMA_LOG10`, `HIGH_RISK_MGDL`).

## 6. Example — real LPA graph, structured cohort, two results

`tests/gwas/make_lpa_phenotype.py` reads each real haplotype's KIV-2 copy number, partitions the haplotypes
into **3 Southern-European-like subpopulations** (differing KIV-2 frequencies **and** a baseline Lp(a)
offset = an ancestry confounder), and samples diploid individuals. Each gets a **log-normal Lp(a)** with the
inverse KIV-2 effect plus covariates (`Age, Sex, PC1–3`) and ~5% `NA` phenotypes. It also emits a
**synthetic genome-wide-like panel** (the real KIV-2 dosage + many subpop-stratified null SNPs) for the
structure-correction demo. One driver runs everything:

```bash
# numpy-capable python for the cohort sim; Rscript needs ggplot2 (conda activate base).
bash tests/gwas/run_lpa_real.sh build/panvar results/real_data/lpa/gwas python3 Rscript        # small/fast
bash tests/gwas/run_lpa_real.sh build/panvar results/real_data/lpa/gwas python3 Rscript --big   # ~10k cohort
```

It runs `bubble → panphorte → call → describe --samples → associate → plot_associate.R` and a self-check,
producing **two** results:

### 6a. Region scan — recover KIV-2 on the real pangenome

`associate` on the real KIV-2 features (graph & k-mer substrates, quantitative & binary), PC-adjusted, MAF
`0.02`. Every run lands on the called KIV-2 locus (**bubble 7**) with the correct **negative** effect (more
KIV-2 → lower Lp(a)), far past the region-wide Bonferroni line:

- **node/edge dosage** localizes the signal to the repeat node `4789` and its **self-loop edge**
  `4789+>4789+` — the most direct graph read of copy number;
- **k-mers** carry the same signal across the repeat-unit markers;
- with `--node-genes`, the hit is named **`LPA`** in the `gene` column.

The MAF filter drops the cohort-invariant features. Note that a *single region* is one correlated signal
block with no null markers, so its λ is **not** interpretable (see the primer) — the operative controls here
are the MAF filter and region-wide multiple testing, and the job is recovering KIV-2, which it does.

### 6b. Structure-correction demo — what PCs and the LMM buy you

On the synthetic genome-wide-like panel (causal KIV-2 + subpop-stratified nulls), the subpopulation
confounder inflates a naive scan and **buries** the causal signal; adding ancestry PCs, or running the LMM
with a panel-derived GRM, restores calibration and surfaces KIV-2 (representative run, n≈2 000):

| analysis | covariates | genomic inflation λ | KIV-2 rank | KIV-2 p |
|----------|-----------|---------------------|-----------|---------|
| naive            | Age, Sex            | **≈ 4** (inflated) | buried (~450) | recoverable but lost in nulls |
| + ancestry PCs   | Age, Sex, PC1–3     | **≈ 1.0–1.2**      | **1**         | ~1e-47 |
| LMM (`--make-kinship`) | Age, Sex (+ random effect) | **≈ 1.0** | **1**  | ~1e-49 |

This is the textbook before/after: a naive λ≫1 collapses to ≈1 once structure is modeled, while the true
KIV-2 signal rises to the top. (Exact numbers depend on `N`/seed; the driver asserts λ drops after
correction.)

![LPA structure demo — naive (inflated) vs PC-adjusted Manhattan](../results/real_data/lpa/gwas/sim_naive.manhattan.png)

## 7. Reading the outputs

- **`assoc_<sub>_<mode>.assoc.tsv`** — one row per tested feature, sorted by `p`:
  `feature_id, layer, bubbles, nodes, n, minor_freq, beta|log_or, se, z, p, p_bonf, q_bh, gene`. `feature_id`
  is the k-mer sequence (k-mer substrate) or the node id / edge key (graph substrate); `bubbles`/`nodes` give
  the graph provenance (from `feature_annot`), so a hit traces back through `call`'s `variant_nodes.tsv` to
  the KIV-2 `DUP`; `gene` names it (`LPA`) when `--node-genes` is supplied.
- **`assoc_<sub>_<mode>.summary.tsv`** — model, samples used, `features_tested`, `dropped_min_maf`, the
  `bonferroni_threshold` (`0.05/n_tests`), significant counts (Bonferroni and FDR<0.05), the
  genomic-inflation **`lambda_gc`**, and (LMM) `lmm_delta`.
- **Manhattan / QQ** — `*.manhattan.{png,pdf}` (raw −log10 p with the nominal `0.05` and Bonferroni lines;
  FDR<0.05 / Bonferroni points highlighted) and `*.qq.{png,pdf}` (genomic-inflation λ).
- **Structure demo** — `sim_{naive,pc,lmm}.{assoc.tsv,summary.tsv,manhattan.png,qq.png}`: compare the
  `lambda_gc` across the three.

## 8. Caveats (honest limits)

- The phenotypes are **simulated** (log-normal Lp(a) with the literature-grounded inverse KIV-2 effect +
  age/sex + a subpopulation offset + noise; high-risk binary at 50 mg/dL). This demonstrates *recovery of a
  planted signal under realistic structure*, **not** a real Lp(a) study; the subpopulations/PCs are
  synthetic ancestry, not real components.
- The **structure-correction panel is synthetic** (a genome-wide-like set of stratified null SNPs alongside
  the real KIV-2 dosage). It exists *because* a single pangenome region has no null markers — so λ and the
  PC/LMM correction are only meaningful at genome-wide-like scale. The **region scan** (§6a) is the result on
  the actual pangenome; its job is KIV-2 recovery via the MAF filter + region-wide multiple testing, not λ.
- The **LMM** here uses `--make-kinship` on the synthetic panel (a valid genome-wide-like GRM). On a real
  study, build the kinship from genome-wide markers and pass it via `--kinship`; a region-only `--make-kinship`
  is proximally contaminated and can over-correct.
- Markers are **bubble-local**: the per-sample dosage is read from the known graph traversal (`cosigt`
  assigns the haplotypes; provenance is recorded), so genome-wide uniqueness is irrelevant here — it would
  only matter if you screened raw reads with these features.
- Counts are the *graph-derived* per-haplotype multiplicity; with real reads you would carry read-depth
  uncertainty from `cosigt` into the dosage.

## See also
- [GWAS primer](gwas_primer.md) — the concepts (dosage, MAF, multiple testing, kinship, λ, LMM vs PCs) from scratch.
- [associate](modules/associate.md) — the GWAS engine, inputs, multiple-testing, and outputs.
- [describe](modules/describe.md) — the BIMBAM exports, `--samples` / `--variant-nodes` / `--variant-flank-bp`.
- [call](modules/call.md) — `node_track.tsv` / `variant_nodes.tsv` used for provenance/traceback.
