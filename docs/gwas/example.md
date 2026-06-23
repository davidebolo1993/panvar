# Pangenome association with panvar (LPA copy-number example)

A runnable, end-to-end example of a **GWAS on panvar's output**, and of why the **multiplicity** (copy
count) panvar emits matters. This page is a manual — commands, outputs, and how to read them. For the
concepts behind each step (dosage, MAF, multiple testing, kinship, λ, PCs vs LMM) see the
[GWAS primer](primer.md).

The trait is **LPA / KIV-2**: the *LPA* gene carries a tandem repeat (the KIV-2 VNTR) whose **copy number is
inversely associated with plasma Lp(a)** (more copies → lower Lp(a) → lower cardiovascular risk). It is the
ideal copy-number example because the repeat unit is present in **every** individual, so the entire signal is
in the count, not in presence/absence (see [primer §genotype = dosage](primer.md#genotype--dosage-and-why-the-count-matters)).

## 1. How panvar maps a GWAS onto a pangenome

A SNP-array or linear-reference GWAS can only test what the single reference represents; it is blind to
sequence **absent from the reference** and to **VNTR/STR copy number** — which is exactly where Lp(a) lives.
A pangenome represents all haplotypes, so this structural variation becomes a testable genotype. panvar turns
the graph into association features like so:

1. **`cosigt` genotypes a sequenced sample** → assigns it a **pair of haplotype paths** from the pangenome
   (the two alleles it best matches at the locus).
2. **`describe`** counts features per **haplotype path** on two substrates — **k-mers** (read-queryable) and
   **node/edge dosage** (graph-local) — each feature localizing to the graph nodes it comes from.
3. The **sample genotype** at a feature is the sum over its two haplotypes: `sample_count = count(hapA) +
   count(hapB)`. For KIV-2 that is the diploid copy number `CN_A + CN_B`. `describe --samples <cosigt.tsv>`
   writes this per-sample file per substrate.
4. **`associate`** runs `phenotype ~ sample genotype`, one test per feature, across samples.

Three properties follow: it is **variant-level** (every feature points back to a bubble/variant, never "the
whole haplotype" as one blob); it is **reference-free** (a haplotype carrying a bubble the reference never
traverses is still genotyped and tested); and the **reference is only for plotting** (the Manhattan x-axis
uses graph/node order for graph features and a per-k-mer index for k-mers).

## 2. The pipeline

The post-`panphorte` modules consume the panphorte-normalized/sorted graph and the panphorte prefix.
`describe --samples <cosigt.tsv>` (one sample per line: `sample <tab> hap1 <tab> hap2 …`, names matching graph
path names) writes per-sample BIMBAM dosage — `bimbam_{kmers,graph}.samples.bimbam.gz` — summed over each
sample's haplotypes, with `bimbam.samples.samples.txt.gz` (column order) and `feature_annot.samples.tsv.gz`
(provenance). `panvar associate` then fits the per-feature regression, applies the MAF filter, and corrects
for multiple testing; `scripts/plot_associate.R` draws the two-panel Manhattan + QQ. See
[describe](../modules/describe.md) and [associate](../modules/associate.md) for inputs and columns.

One driver runs everything on the real LPA graph and a structured synthetic cohort:

```bash
# numpy-capable python for the cohort sim; Rscript needs ggplot2 (conda activate base).
bash tests/gwas/run_lpa_real.sh build/panvar results/real_data/lpa/gwas python3 Rscript        # small/fast
bash tests/gwas/run_lpa_real.sh build/panvar results/real_data/lpa/gwas python3 Rscript --big   # ~10k cohort
```

`tests/gwas/make_lpa_phenotype.py` reads each real haplotype's KIV-2 copy number, partitions haplotypes into
**3 Southern-European-like subpopulations** (differing KIV-2 frequencies **and** a baseline Lp(a) offset = an
ancestry confounder), and samples diploid individuals with a **log-normal Lp(a)** carrying the inverse KIV-2
effect plus covariates (`Age, Sex, PC1–3`) and ~5% `NA` phenotypes. It also emits a **synthetic
genome-wide-like panel** (real KIV-2 dosage + many subpop-stratified null SNPs) for the structure-correction
demo. The script runs `bubble → panphorte → call → describe --samples → associate → plot_associate.R` and a
self-check, producing two results.

## 3. Result A — region scan: recover KIV-2 on the real pangenome

`associate` on the real KIV-2 features (graph & k-mer substrates, quantitative & binary), PC-adjusted, MAF
`0.02`. Every run lands on the called KIV-2 locus (**bubble 7**) with the correct **negative** effect, far
past the region-wide Bonferroni line:

- **node/edge dosage** localizes the signal to the repeat node `4789` and its self-loop edge `4789+>4789+` —
  the most direct graph read of copy number;
- **k-mers** carry the same signal across the repeat-unit markers;
- with `--node-genes`, the hit is named **`LPA`** in the `gene` column.

The MAF filter drops the cohort-invariant features. A single region is one correlated signal block with no
null markers, so its λ is **not** interpretable here (see [primer §λ needs nulls](primer.md#population-structure-and-genomic-inflation));
the operative controls are the MAF filter and region-wide multiple testing, and the job — recovering KIV-2 —
is done.

## 4. Result B — structure-correction demo: what PCs and the LMM buy you

On the synthetic genome-wide-like panel (causal KIV-2 + subpop-stratified nulls), the subpopulation
confounder inflates a naive scan and **buries** the causal signal; ancestry PCs, or the LMM with a
panel-derived GRM, restore calibration and surface KIV-2 (representative run, n≈2 000):

| analysis | covariates | genomic inflation λ | KIV-2 rank | KIV-2 p |
|----------|-----------|---------------------|-----------|---------|
| naive            | Age, Sex            | **≈ 4** (inflated) | buried (~450) | recoverable but lost in nulls |
| + ancestry PCs   | Age, Sex, PC1–3     | **≈ 1.0–1.2**      | **1**         | ~1e-47 |
| LMM (`--make-kinship`) | Age, Sex (+ random effect) | **≈ 1.0** | **1**  | ~1e-49 |

The textbook before/after: a naive λ≫1 collapses to ≈1 once structure is modeled, while the true KIV-2 signal
rises to the top. (Exact numbers depend on `N`/seed; the driver asserts λ drops after correction.)

![LPA structure demo — naive (inflated) vs PC-adjusted Manhattan](../../results/real_data/lpa/gwas/sim_naive.manhattan.png)

## 5. Reading the outputs

- **`assoc_<sub>_<mode>.assoc.tsv`** — one row per tested feature, sorted by `p`. `feature_id` is the k-mer
  sequence or node/edge key; `bubbles`/`nodes` give the graph provenance (from `feature_annot`), so a hit
  traces back through `call`'s `variant_nodes.tsv` to the KIV-2 `DUP`; `gene` names it (`LPA`) with
  `--node-genes`. Full column meanings: [associate outputs](../modules/associate.md#outputs).
- **`assoc_<sub>_<mode>.summary.tsv`** — model, `samples_used`, `features_tested`, `dropped_min_maf`, the
  `bonferroni_threshold`, significant counts, `lambda_gc`, and (LMM) `lmm_delta`.
- **Manhattan / QQ** — `*.manhattan.{png,pdf}`: two stacked panels (before correction = raw −log10 p with the
  nominal + Bonferroni lines; after = BH −log10 q with the q=0.05 line), points coloured by verdict, genes
  flagged. `*.qq.{png,pdf}` carries λ.
- **Structure demo** — `sim_{naive,pc,lmm}.{assoc.tsv,summary.tsv,manhattan.png,qq.png}`: compare `lambda_gc`
  across the three.

## 6. Validation against GEMMA

`panvar associate` is a from-scratch implementation, so we check it against **GEMMA** (the reference
mixed-model GWAS tool) on the *same* BIMBAM panel + phenotype/covariates — BIMBAM is GEMMA's native
mean-genotype format, so the genotypes load unchanged. `tests/gwas/validate_gemma.sh` runs both and reports
the Pearson correlation of effect size and of −log10 p over the shared features:

| comparison | r(β) | r(−log10 p) | top hit |
|------------|------|-------------|---------|
| linear (`associate --model linear` vs GEMMA `-lm`) | **1.0000** | **1.0000** | — |
| mixed (`associate --model lmm --kinship` vs GEMMA `-lmm`) | **0.9997** | **0.9997** | **match** |

The statistics are effectively identical — the engine is correct. **One revealing difference:** GEMMA
analyzes 2000 of the 2001 features and **silently drops the KIV-2 copy-number marker**, because its
allele-frequency model assumes a diploid 0–2 dosage and KIV-2's count (here 14–61) has "allele frequency" ≫
1, so its built-in MAF filter discards it. `panvar associate` has no such assumption and tests it directly
(KIV-2 p ≈ 1e-22). That is precisely why panvar carries the **raw count** through BIMBAM and tests
multiplicity itself — a generic 0–2 GWAS tool cannot genotype the very locus this example is about.

## 7. Caveats (honest limits)

- The phenotypes are **simulated** (log-normal Lp(a) with a literature-grounded inverse KIV-2 effect +
  age/sex + a subpopulation offset + noise; high-risk binary at 50 mg/dL). This demonstrates *recovery of a
  planted signal under realistic structure*, not a real Lp(a) study; the subpopulations/PCs are synthetic
  ancestry. The literature this is grounded in (Schmidt 2016; Coassin & Kronenberg 2022; Moli-sani 2025) is
  in [references.md](../references.md).
- The **structure-correction panel is synthetic** — it exists *because* a single pangenome region has no null
  markers, so λ and the PC/LMM correction are only meaningful at genome-wide-like scale. The **region scan**
  (§3) is the result on the actual pangenome.
- The **LMM** here uses `--make-kinship` on the synthetic panel (a valid genome-wide-like GRM). On a real
  study, build kinship from genome-wide markers and pass `--kinship`; a region-only `--make-kinship` is
  proximally contaminated and can over-correct.
- Counts are the *graph-derived* per-haplotype multiplicity; with real reads you would carry read-depth
  uncertainty from `cosigt` into the dosage.

## See also
- [GWAS primer](primer.md) — the concepts (dosage, MAF, multiple testing, kinship, λ, LMM vs PCs) from scratch.
- [associate](../modules/associate.md) — the GWAS engine, inputs, multiple-testing, and outputs.
- [describe](../modules/describe.md) — the BIMBAM exports, `--samples` / `--variant-nodes` / `--variant-flank-bp`.
- [call](../modules/call.md) — `variant_nodes.tsv` / `variant_paths.tsv` used for provenance/traceback.
