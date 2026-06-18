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
- **The reference is only for plotting.** A Manhattan needs an x-axis. panvar offers two (`plot_gwas.R --x ref|nodes`): reference bp (familiar; novel insertions plot at their bubble anchor) or **graph/node order** (pangenome-native; places *every* variant, including reference-disjoint ones).

## 3. Why multiplicity (counts), not just presence/absence

A **presence/absence** association asks only "does this sample contain the feature (yes/no)". That works for
SNP- and indel-like variants. It **fails for CNVs** like KIV-2: the repeat-unit feature is present in
*everyone* (CN ≥ 1), so presence/absence has **no contrast** (frequency = 1) and the variant is filtered out
before testing. The information is entirely in the **count**. panvar's sample-level files carry the true
per-sample count, so a **count-based** test recovers KIV-2 while a **presence/absence-based** test cannot.
The example below demonstrates exactly this contrast — on both feature substrates.

> **Harness.** `scripts/gwas_demo.py` is a small, transparent **scipy** harness that runs *both* a
> presence/absence test and a count test on the same file, so the comparison is apples-to-apples and
> reproducible. It stands in for any count-based or presence/absence-based association method; `panvar`
> emits the standard fsm-lite file (`<feature> | sample:count`) that such tools also read, and does not run
> the association itself.

## 4. The pipeline

The post-`panphorte` modules consume the **panphorte-normalized/sorted graph** and the **panphorte** prefix.

- **`describe --samples <cosigt.tsv>`** (one sample per line: `sample <tab> hap1 <tab> hap2 …`; haplotype
  names must match graph path names) writes a **sample-level file for each substrate** —
  `fsm_kmers.samples.txt.gz` (k-mers) and `fsm_graph.samples.txt.gz` (node/edge) — whose value is the summed
  dosage over the sample's assigned haplotypes (a haplotype listed twice → counted twice = homozygous). See
  [describe](modules/describe.md).
- **`scripts/gwas_demo.py`** reads an fsm file + a phenotype table and runs, per feature, a
  **presence/absence** test (binarize count>0) and a **count** test (use the dosage). `--substrate kmer|graph`
  selects which sample-level file to test; `--mode continuous|binary` selects the trait coding. It writes
  `*.assoc.tsv` with `pa_q` and `count_q`, plus per-feature `pos` (reference), `node_min` (graph order),
  `nodes`, and the traced `variant`.
- **`scripts/plot_gwas.R`** draws the two-panel Manhattan and the QQ (with λ), on either axis.

## 5. Example — real LPA graph, literature-based phenotype

`tests/gwas/make_lpa_phenotype.py` reads each **real** LPA haplotype's KIV-2 copy number from `panphorte`
(`copies.tsv`; the demo finds the KIV-2 unit = 5,547 bp, copies **6–32** across 465 haplotypes), pairs real
haplotypes into 200 synthetic diploid individuals, and assigns `Lp(a) = base − slope·(CN_A+CN_B) + noise`
(synthetic values, literature-plausible **inverse** effect; case/control by median split). One driver runs
the whole thing across **both substrates** and **both trait codings**:

```bash
bash tests/gwas/run_lpa_real.sh build/panvar out/gwas_lpa python3 Rscript
```

It runs `bubble → panphorte → call → describe --samples → gwas_demo.py → plot_gwas.R` and a self-check.
Observed (continuous trait):

```
gwas_kmer_continuous.assoc.tsv    top variant=bubble7_DUP_4790(DUP) count_q=1.4e-29 pa_q=1.0e+00  present-in-all & count-only-sig=921  -> PASS
gwas_graph_continuous.assoc.tsv   top variant=bubble7_DUP_4790(DUP) count_q=1.5e-26 pa_q=1.0e+00  present-in-all & count-only-sig=16   -> PASS
```

Both substrates land on the same locus — the called KIV-2 **`DUP`** (`bubble7_DUP_4790`) — through the
**count** test, while the **presence/absence** test finds nothing there (`pa_q = 1.0`):

- **k-mers** spread the signal across hundreds of repeat-unit markers (≈ 920 present-in-all k-mers that are
  count-significant but presence/absence-blind);
- **node/edge dosage** localizes it to a handful of features — the repeat node `4790` and its self-loop edge
  `4790+>4790+` — the most direct graph read of copy number.

The binary trait behaves the same way (e.g. k-mer `count_q ≈ 1.3e-12`, `pa_q = 1.0`). The Manhattan
(`out/gwas_lpa/plot_kmer_continuous.manhattan.png`) shows the KIV-2 peak only in the **count** panel; `--x
nodes` places the same peak by graph order.

## 6. Reading the outputs

- **`*.assoc.tsv`** — one row per tested feature: `kmer, bubble_id, pos, node_min, nodes, variant,
  n_carriers, max_count, pa_p, count_p, pa_bonf, count_bonf, pa_q, count_q` (the first column holds the
  k-mer sequence or, for `--substrate graph`, the node id / edge key). Sort by `count_q` for the count
  GWAS, `pa_q` for presence/absence.
- **Traceback** — the `variant` column already maps each significant feature back through its `nodes` to the
  `call` variant (here the KIV-2 `DUP`), so a hit reads as a biological variant, not just a sequence.
- **Manhattan / QQ** — `*.manhattan.{png,pdf}` (P/A vs count panels) and `*.qq.{png,pdf}` (λ per method).

## 7. Caveats (honest limits)

- The phenotypes here are **simulated** from the genotype (a literature-plausible inverse KIV-2→Lp(a)
  effect, Gaussian noise, binary trait by median split). This demonstrates *recovery of a planted signal*,
  **not** a real Lp(a) study.
- The simulation has **no covariates** and the harness uses **fixed-effect, single-predictor** tests with
  **no covariate or population-structure correction**. Real studies add covariates (age, sex, …) and a
  kinship/PC or mixed-model correction (watch the QQ λ). `gwas_demo.py` is a transparent illustration, not a
  production GWAS engine.
- Markers are **bubble-local**, and that is all the association needs: the per-sample dosage is read from
  the known graph traversal (`cosigt` assigns the haplotypes; provenance is recorded), so genome-wide
  uniqueness is irrelevant here — it would only matter if you instead screened raw sequencing reads with
  these features, which this pipeline does not do.
- Counts are the *graph-derived* per-haplotype multiplicity; with real reads you would carry read-depth
  uncertainty from `cosigt` into the dosage.

## See also
- [describe](modules/describe.md) — the sample-level files and `--samples` / `--variant-nodes` / `--variant-flank-bp` options.
- [call](modules/call.md) — `node_track.tsv` / `variant_nodes.tsv` used for position and traceback.
