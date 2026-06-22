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
> (`phenotype ~ genotype + covariates`), so copy-number loci are first-class. The same genotypes are also
> exported as **BIMBAM** (GEMMA-ready) and **fsm-lite** (pyseer) by `describe`, so an external count- or
> presence/absence-based tool can be used instead if desired.

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

## 5. Example — real LPA graph, literature-based phenotype + covariates

`tests/gwas/make_lpa_phenotype.py` reads each **real** LPA haplotype's KIV-2 copy number from `panphorte`
(`copies.tsv`; KIV-2 unit ≈ 5.5 kb, copies **1–32** across 466 haplotypes), pairs real haplotypes into 200
synthetic diploid individuals, and assigns `Lp(a) = base − slope·(CN_A+CN_B) + age/sex effects + noise`
(synthetic values, literature-plausible **inverse** KIV-2 effect; case/control by median split). It also
emits **covariates** (`Age, Sex, PC1–PC3`) and sets ~5% of phenotypes to `NA` to exercise the filters. One
driver runs the whole thing across **both substrates** and **both trait codings**:

```bash
# Rscript needs ggplot2 (conda activate base). Also run by scripts/regen_results.sh.
bash tests/gwas/run_lpa_real.sh build/panvar results/real_data/lpa/gwas python3 Rscript
```

It runs `bubble → panphorte → call → describe --samples → panvar associate → plot_associate.R` and a
self-check. Observed top hits (193 of 200 samples used; 7 NA phenotypes dropped; 5 covariates):

| run | top feature | bubble | effect | p | p_bonf |
|-----|-------------|--------|--------|---|--------|
| graph / quant   | `4789+>4789+` (KIV-2 self-loop) | 7 | β = **−2.18** | 1.5e-52 | 2.6e-49 |
| kmers / quant   | a KIV-2 k-mer | 7 | β = **−2.18** | 5.4e-53 | 1.1e-49 |
| graph / binary  | node `4789` | 7 | log OR = **−0.27** | 1.9e-11 | 3.3e-08 |
| kmers / binary  | a KIV-2 k-mer | 7 | log OR = **−0.26** | 1.7e-11 | 3.3e-08 |

Every run lands on the called KIV-2 locus (**bubble 7**) with the correct **negative** effect (more KIV-2
→ lower Lp(a)), far past the region-wide Bonferroni line:

- **node/edge dosage** localizes the signal to the repeat node `4789` and its **self-loop edge**
  `4789+>4789+` — the most direct graph read of copy number;
- **k-mers** carry the same signal across the repeat-unit markers.

The MAF filter (`--min-maf 0.02`) drops the cohort-invariant features (e.g. 1413 of 3180 graph features),
and the Manhattan/QQ show a clean peak with λ near 1.

![LPA Lp(a) ~ graph dosage Manhattan](../results/real_data/lpa/gwas/assoc_graph_quant.manhattan.png)

## 6. Reading the outputs

- **`assoc_<sub>_<mode>.assoc.tsv`** — one row per tested feature, sorted by `p`:
  `feature_id, layer, bubbles, nodes, n, minor_freq, beta|log_or, se, z, p, p_bonf, q_bh`. `feature_id` is
  the k-mer sequence (k-mer substrate) or the node id / edge key (graph substrate); `bubbles`/`nodes` give
  the graph provenance (from `feature_annot`), so a hit traces back through `call`'s `variant_nodes.tsv` to
  the KIV-2 `DUP` and, with `call --gtf`, to the `LPA` gene.
- **`assoc_<sub>_<mode>.summary.tsv`** — model, samples used, `features_tested`, `dropped_min_maf`, the
  `bonferroni_threshold` (`0.05/n_tests`), and significant counts (Bonferroni and FDR<0.05).
- **Manhattan / QQ** — `*.manhattan.{png,pdf}` (raw −log10 p with the nominal `0.05` and Bonferroni lines;
  FDR<0.05 / Bonferroni points highlighted) and `*.qq.{png,pdf}` (genomic-inflation λ).

## 7. Caveats (honest limits)

- The phenotypes here are **simulated** from the genotype (inverse KIV-2→Lp(a) effect + small age/sex
  effects + Gaussian noise; binary trait by median split). This demonstrates *recovery of a planted
  signal*, **not** a real Lp(a) study; the PCs are random, not real ancestry components.
- `associate` uses a **fixed-effect GLM with covariates**. It does **not** yet model relatedness /
  population structure with a kinship matrix or mixed model (watch the QQ λ); an **LMM backend** is the
  planned next step. Multiple-testing is **region-wide** (Bonferroni over the tested features + BH FDR),
  which is the right scale here — and conservative, since correlated k-mers over-count the independent tests.
- Markers are **bubble-local**, and that is all the association needs: the per-sample dosage is read from
  the known graph traversal (`cosigt` assigns the haplotypes; provenance is recorded), so genome-wide
  uniqueness is irrelevant here — it would only matter if you screened raw reads with these features.
- Counts are the *graph-derived* per-haplotype multiplicity; with real reads you would carry read-depth
  uncertainty from `cosigt` into the dosage.

## See also
- [associate](modules/associate.md) — the GWAS engine, inputs, multiple-testing, and outputs.
- [describe](modules/describe.md) — the BIMBAM exports, `--samples` / `--variant-nodes` / `--variant-flank-bp`.
- [call](modules/call.md) — `node_track.tsv` / `variant_nodes.tsv` used for provenance/traceback.
