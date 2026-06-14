# Pangenome k-mer GWAS with panvar (LPA copy-number example)

This page is a self-contained worked example of running a **genome-wide association study (GWAS)** on
panvar's k-mer output, and — importantly — of *why* the **k-mer multiplicity** (copy count) panvar emits
matters. It is written for someone who is **not** a GWAS specialist, so it starts with the concepts.

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

- the **unit of analysis is the sample** (an individual), not a haplotype;
- the **genotype is a dosage**: 0/1/2 copies of an allele for a diploid SNP, or a **count** for a copy-number
  variant (CNV).

The output is a list of p-values, one per variant, usually:

- corrected for multiple testing (we report **Benjamini-Hochberg q-values** and Bonferroni), and
- drawn as a **Manhattan plot** (x = genomic position, y = −log₁₀ q; peaks = associated loci) plus a
  **QQ plot** (observed vs expected p-values; the genomic-inflation factor λ should be ≈ 1 if the test is
  well-calibrated and there is no hidden population structure).

## 2. How this maps onto a pangenome

In a pangenome there is no single reference coordinate system; variation lives in the **graph**. panvar
turns the graph into association features as follows:

1. **`cosigt` genotypes a sequenced sample** → it assigns the sample a **pair of haplotype paths** from the
   pangenome (the two alleles it best matches at the locus).
2. **`describe`** counts canonical k-mers per **haplotype path** (per bubble). A k-mer is a proxy for a
   *variant*: it localizes to the graph nodes it comes from, which trace back to a called variant.
3. The **sample genotype** at a k-mer is the **aggregate over its two assigned haplotypes**:
   `sample_count = count(hapA) + count(hapB)`. For KIV-2 this is exactly the **diploid copy number**
   (CN_A + CN_B). `describe --samples <cosigt.tsv>` does this aggregation and writes a sample-level file.
4. **GWAS** = phenotype ~ sample genotype, one test per k-mer, across samples.

Three consequences worth internalizing:

- **It is variant-level, not whole-haplotype.** Every k-mer points back to a bubble/variant; we never test
  "the whole haplotype" as one blob.
- **Testing is reference-free.** k-mer counts exist for every haplotype whether or not the reference
  traverses that sequence. A sample whose haplotypes carry a bubble **not spanned by the reference** is
  still genotyped and tested — its novel k-mers get counts like any other.
- **The reference is only for plotting.** A Manhattan needs an x-axis. panvar offers two
  (`plot_gwas.R --x ref|nodes`): reference bp (familiar; novel insertions plot at their bubble anchor) or
  **graph/node order** (pangenome-native; places *every* variant, including reference-disjoint ones).

## 3. Why multiplicity (counts), not just presence/absence

Most k-mer GWAS tools (**pyseer**, kmersGWAS, DBGWAS) are **presence/absence**: the genotype is "does this
sample contain the k-mer (yes/no)". That works for SNP- and indel-like variants. It **fails for CNVs** like
KIV-2: the repeat-unit k-mer is present in *everyone* (CN ≥ 1), so presence/absence has **no contrast**
(allele frequency = 1) and the variant is filtered out before testing. The information is entirely in the
**count**.

The count-aware approach is what **HAWK** (Rahman et al. 2018, k-mer abundance GWAS) formalizes. panvar's
`fsm_kmers*.txt.gz` carries the true per-strain count, so a count-based test recovers KIV-2 while
presence/absence cannot. The example below demonstrates exactly this contrast.

> **Tools used here.** For the count (multiplicity) side we use a small, transparent **scipy** harness
> (`scripts/gwas_demo.py`) that runs *both* a presence/absence test and a count test, so the comparison is
> apples-to-apples and reproducible. **pyseer** is the reference presence/absence tool; the demo cross-checks
> against it where it is installable (it currently has no build for Python 3.13 / this architecture — the
> exact command is given in §7 to run in a compatible environment). **HAWK** is the published count-based
> tool but ships Linux-only binaries.

## 4. The pipeline

```
graph.gfa ──bubble──▶ sites ──panphorte──▶ (KIV-2 copies) 
                          │                      │
                          ├──call──▶ DUP variant + node_track.tsv + variant_nodes.tsv
                          │
                          └──describe --samples cosigt.tsv──▶ fsm_kmers.samples.txt.gz  (per-SAMPLE dosage)
                                                                      │
                          scripts/gwas_demo.py ◀───────────────────────┘  (P/A vs count, BH q, traceback)
                                   │
                          scripts/plot_gwas.R --x ref|nodes  ──▶ Manhattan + QQ
```

- **`describe --samples <cosigt.tsv>`** (one sample per line: `sample <tab> hap1 <tab> hap2 …`; haplotype
  names must match graph path names) writes, **in addition** to the per-haplotype `fsm_kmers.txt.gz`, a
  **sample-level `fsm_kmers.samples.txt.gz`** whose value is the summed dosage over the sample's assigned
  haplotypes (a haplotype listed twice → counted twice = homozygous). See [describe](modules/describe.md).
- **`scripts/gwas_demo.py`** reads an fsm file + a phenotype table and runs, per k-mer, a **presence/absence**
  test (binarize count>0) and a **count** test (use the dosage). It writes `*.assoc.tsv` with `pa_q` and
  `count_q`, plus per-k-mer `pos` (reference), `node_min` (graph order), `nodes`, and the traced `variant`.
- **`scripts/plot_gwas.R`** draws the two-panel Manhattan and the QQ (with λ), on either axis.

## 5. Example A — synthetic cohort (known ground truth)

`tests/gwas/make_lpa_cnv.py` builds a **haplotype panel** (40 graph haplotypes spanning KIV-2 CN 1–6) and a
**diploid cohort** (150 individuals, each two haplotypes; `samples.tsv`), with three loci of *known* truth:

| locus | model | linked to phenotype? | expected result |
|---|---|---|---|
| **KIV-2 VNTR** | multiplicity (present in all) | yes (causal) | **count-significant, P/A-blind** |
| decoy A | presence/absence | no | non-significant in both |
| decoy B | presence/absence | yes | significant in **both** |

Run it:

```bash
bash tests/gwas/run_synthetic.sh build/panvar out/gwas_syn python3 Rscript
```

This regenerates the cohort, runs `bubble → call → describe --samples → gwas_demo.py` (continuous Lp(a) and
binary high-Lp(a)) → `plot_gwas.R` (both axes) and a **self-check** against `truth.tsv`. Observed:

```
KIV-2 present-in-all unit k-mers: 7 | min count_q=1.8e-83 | max pa_q=1.0e+00   <- count finds it, P/A cannot
decoy_A: min count_q=0.67 | min pa_q=1.0e+00                                   <- null in both (clean)
decoy_B: min count_q=5.0e-49 | min pa_q=3.4e-28                                <- both find it
SELF-CHECK: PASS
```

The Manhattan (`out/gwas_syn/plot_continuous_ref.manhattan.png`) shows the KIV-2 peak **only in the count
panel**; `--x nodes` places the same peak by graph order.

## 6. Example B — real LPA graph, literature-based phenotype

`tests/gwas/make_lpa_phenotype.py` reads each **real** LPA haplotype's KIV-2 copy number from `panphorte`
(`copies.tsv`; the demo finds the KIV-2 unit = 5,547 bp, copies **6–32** across 465 haplotypes), pairs real
haplotypes into 200 synthetic diploid individuals, and assigns `Lp(a) = base − slope·(CN_A+CN_B) + noise`
(synthetic values, literature-plausible **inverse** effect). Run it:

```bash
bash tests/gwas/run_lpa_real.sh build/panvar out/gwas_lpa python3 Rscript
```

Observed (continuous):

```
top count hit: bubble=7 variant=bubble7_INS_4605  count_q=7.2e-29  pa_q=1.0e+00
count q<0.05 concentrate on bubble 7 (1218 k-mers) = the KIV-2 locus
present-in-all-200 k-mers significant by COUNT but not P/A: 1050
SANITY: PASS
```

The real-data Manhattan (`out/gwas_lpa/plot_continuous_ref.manhattan.png`) shows a clean KIV-2 peak at the
real *LPA* locus (~chr6:161.86 Mb in CHM13) in the **count** panel and essentially nothing in the
presence/absence panel — the multiplicity advantage on genuine topology.

> On the **un-normalized** graph the KIV-2 tandem surfaces as an `INS` of repeat units (the k-mers trace to
> `bubble7_INS_...`). Run the intended `panphorte → odgi sort → vg snarls → bubble → call` pipeline first to
> get a clean **`DUP`** with explicit copy number; the k-mer counts (and the association) are identical
> either way because describe counts each traversal of the repeat unit.

## 7. presence/absence with real pyseer (optional cross-check)

The harness P/A test mirrors what pyseer does. To validate against the real tool (in a Python ≤3.11 env):

```bash
pip install pyseer
# phenotype file: header "samples<TAB>binary", sample ids == fsm strain ids
pyseer --phenotypes pheno.pyseer.tsv --kmers out/gwas_syn/desc/fsm_kmers.samples.txt.gz \
       --no-distances --min-af 0 --max-af 1 > pyseer.assoc.tsv
```

Expected: pyseer agrees with the harness P/A column, **cannot flag the KIV-2 unit k-mers** (allele
frequency ≈ 1 → pre-filtered), and *does* flag the presence/absence decoy. `run_synthetic.sh` runs this
automatically when `pyseer` is on `PATH`.

## 8. Reading the outputs

- **`*.assoc.tsv`** — one row per tested k-mer: `kmer, bubble_id, pos, node_min, nodes, variant,
  n_carriers, max_count, pa_p, count_p, pa_bonf, count_bonf, pa_q, count_q`. Sort by `count_q` for the
  count GWAS, `pa_q` for presence/absence.
- **Traceback** — the `variant` column already maps each significant k-mer back through its `nodes` to the
  `call` variant (e.g. the KIV-2 DUP/INS), so a hit reads as a biological variant, not just a sequence.
- **Manhattan / QQ** — `*.manhattan.{png,pdf}` (P/A vs count panels) and `*.qq.{png,pdf}` (λ per method).

## 9. Caveats (honest limits)

- The phenotypes here are **simulated** from the genotype (with a literature-plausible effect for LPA); this
  demonstrates *recovery of a planted signal*, not a real Lp(a) study.
- The harness uses **fixed-effect** tests with **no population-structure correction**. Real studies add a
  kinship/PC correction (pyseer's LMM, or a mixed model); watch the QQ λ.
- k-mer **specificity is within-bubble** (genome-wide uniqueness is a known gap); the traceback to nodes
  mitigates interpretation.
- Counts are the *graph-derived* per-haplotype multiplicity; with real reads you would carry read-depth
  uncertainty from `cosigt` into the dosage.

## See also
- [describe](modules/describe.md) — the `fsm_kmers*.txt.gz` and `--samples` / `--variant-nodes` options.
- [call](modules/call.md) — `node_track.tsv` / `variant_nodes.tsv` used for position and traceback.
- [presentation](presentation.md) — slide-ready schematics of the variants behind the hits.
