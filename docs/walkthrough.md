# Walkthrough — the LPA locus, end to end

A single worked run of the whole pipeline on one real locus: LPA (the `LPA` gene's KIV-2 VNTR — a variable-number tandem repeat whose copy number is a known cardiovascular risk factor). It shows the usage of all the modules in `panvar` with the exception of `associate`, which is covered by the dedicated [GWAS example](gwas/example.md); here it is only pointed at.

Every command runs against the committed test graph `tests/real_data/lpa.gfa.gz` and writes under `results/real_data/lpa/`. The whole run is reproduced by `scripts/regen_results.sh lpa`. The reference path is `GRCh38#0#chr6:160509252-160734894`.

```bash
GFA=tests/real_data/lpa.gfa.gz
REF="GRCh38#0#chr6:160509252-160734894"
OUT=results/real_data/lpa
```

---

## 1. `bubble` — find the variant sites

Sorts/flips the graph along the reference and decomposes it into bubbles (snarl sites where haplotypes diverge). On LPA this yields 12 bubbles; the KIV-2 tandem is bubble 7, spanning ~12–184 kb across the 466 haplotypes (its size is the copy-number variation).

```bash
./build/panvar bubble \
-i "$GFA" \
-o "$OUT/bubble/bubble" \
-r "$REF" \
--gtf tests/real_data/Homo_sapiens.GRCh38.116.gtf.gz
```

Bandage plot, from `$OUT/bubble/bubble.sorted.gfa` + `bubble.bandage_nodes.csv`):

![LPA bubble graph in Bandage](img/lpa_bubble_bandage.png)

The linear reference backbone with the KIV-2 site ballooning into a dense tangle of parallel node paths — one strand per copy-number allele. That tangle is bubble 7.

---

## 2. `inspect` — cluster the haplotypes at the site

A sanity check on one bubble: it pulls every haplotype's `source → sink` walk through bubble 7 and clusters them by walk similarity.

```bash
./build/panvar inspect \
  -i "$OUT/bubble/bubble.sorted.gfa" \
  --bubble-prefix-in "$OUT/bubble/bubble" --bubble-id 7 \
  --cluster \
  --cluster-similarity 0.95 \
  -o "$OUT/inspect/inspect"
```

Then, run `scripts/plot_node_coverage_heatmap.R`:

```bash
Rscript scripts/plot_node_coverage_heatmap.R \
  --table "$OUT/inspect/inspect.bubble_7.node_counts.tsv" \
  --node-lengths "$OUT/inspect/inspect.bubble_7.node_lengths.tsv" \
  --cluster-by "$OUT/inspect/inspect.bubble_7.clusters.tsv" \
  --clusters "$OUT/inspect/inspect.bubble_7.clusters.tsv" \
  --out "$OUT/plots/lpa_node_heatmap"
```

![LPA bubble graph in Bandage](img/lpa_inspect_node_heatmap.png)

Rows are (representative) haplotypes, columns are the repeat-unit nodes. A haplotype with more KIV-2 copies traverses the unit nodes more times, so its row is warmer across the repeat block — the heatmap is a direct picture of the copy-number ladder. Drop `--clusters` to show all the haplotypes.

---

## 3. `panphorte` — fold the tandem into a countable unit

The KIV-2 array is spelled out copy-by-copy in the graph, which would mis-type as a pile of insertions. `panphorte` collapses each detected tandem array onto one repeat-unit (`REP`) node with a self-loop, so copy number becomes the loop-traversal count. On LPA bubble 7 it folds a 5,547 bp unit carried by all 466 haplotypes, collapsing 3,484 nodes, with copies ranging 1–32.

```bash
./build/panvar panphorte \
  -i "$OUT/bubble/bubble.sorted.gfa" \
  --bubble-prefix-in "$OUT/bubble/bubble" \
  -o "$OUT/panphorte/panphorte" \
  --reference-path "$REF" \
  --min-similarity 0.90 \
  --gtf tests/real_data/Homo_sapiens.GRCh38.116.gtf.gz
```

![LPA panphorte-normalized graph in Bandage](img/lpa_panphorte_bandage.png)

The KIV-2 tangle has collapsed to a single node with a self-loop back onto itself — the whole array is now one countable unit. Contrast with the step-1 Bandage view of the same site.

---

## 4. `inspect` again — on the folded graph

The same cluster view, now on the `panphorte` output, confirms the fold preserved the per-haplotype signal: the copy-number ladder is still there, just read off self-loop traversals instead of a long node run.

```bash
./build/panvar inspect \
  -i "$OUT/panphorte/panphorte.normalized.sorted.gfa" \
  --bubble-prefix-in "$OUT/panphorte/panphorte" --bubble-id 7 \
  -o "$OUT/inspect/inspect_panphorte"
```

Then, run `scripts/plot_node_coverage_heatmap.R`. Let's use the clusters coming from the previous `inspect` run to have a clear visual comparison with the previous heatmap:

```bash
Rscript scripts/plot_node_coverage_heatmap.R \
  --table "$OUT/inspect/inspect_panphorte.bubble_7.node_counts.tsv" \
  --node-lengths "$OUT/inspect/inspect_panphorte.bubble_7.node_lengths.tsv" \
  --cluster-by "$OUT/inspect/inspect.bubble_7.clusters.tsv" \
  --clusters "$OUT/inspect/inspect.bubble_7.clusters.tsv" \
  --out "$OUT/plots/lpa_node_heatmap_panphorte"
```


![LPA inspect on the panphorte graph](img/lpa_inspect_node_heatmap_panphorte.png)

---

## 5. `call` — type the variants and read copy number

Diffs every haplotype against the reference and emits a VCF. On LPA: 7 DEL, 13 INS, 3 DUP. The headline record is the KIV-2 DUP on bubble 7: `REF_CN=6`, `RU_LEN=5547`, with a per-sample `CN` spanning 1–32 copies (and `CNBP` giving the actual bp each haplotype gains/loses).

```bash
./build/panvar call \
  -i "$OUT/panphorte/panphorte.normalized.sorted.gfa" \
  --bubble-prefix-in "$OUT/panphorte/panphorte" --reference-path "$REF" \
  -o "$OUT/call/call" \
  --cn \
  --gtf tests/real_data/Homo_sapiens.GRCh38.116.gtf.gz
```

`scripts/plot_vcf_map.R` plot the oncoprint-style variant map), subsetting to 50 haplotypes:

```bash
Rscript scripts/plot_vcf_map.R \
  --vcf "$OUT/call/call.region.vcf" \
  --reference-path "$REF" \
  --out "$OUT/plots/lpa_vcf_map" \
  --max-paths 50
```

![LPA call variant map](img/lpa_vcf_map.png)


Rows are haplotypes, columns are called variants grouped by bubble. The KIV-2 column is the wide blue band shaded by `FORMAT:CN` — dark = high copy number. DEL are red, INS green/purple, INV orange. The KIV-2 blue gradient is the copy-number spectrum the whole pipeline was after.

---

## 6. `describe` — genotype features for association

Turns the called bubbles into per-haplotype dosage features (k-mer, node/edge graph, and per-variant), exported as BIMBAM matrices — the input to `associate`. Restricting to the called variant nodes keeps only features that carry signal.

```bash
./build/panvar describe \
  -i "$OUT/panphorte/panphorte.normalized.sorted.gfa" \
  --bubble-prefix-in "$OUT/panphorte/panphorte" \
  --out-dir "$OUT/describe" \
  --variant-nodes "$OUT/call/call.variant_nodes.tsv"
```
---

## 7. `associate` — the GWAS

Association (`phenotype ~ genotype`, multiple-testing correction, conditional independence) is its own worked example on this exact LPA output: see [GWAS example](gwas/example.md). In short, testing the KIV-2 features against a Lp(a) phenotype recovers the KIV-2 copy-number signal as the region's top hit, and structure correction keeps the genomic-inflation λ ≈ 1.


---

## 8. `benchmark` — round-trip QV of the calls

How faithfully do the calls reconstruct the haplotypes? For each called bubble and each haplotype, `benchmark` rebuilds the haplotype from the reference walk with*only the called events applied (which nodes each call explains comes from `variant_nodes.tsv`), aligns that to the haplotype's true graph walk, and scores a cosigt-style `QV = -10·log10(max(0.5, δ)/S)`. Uncalled variation (SNPs, sub-50 bp indels) is exactly what keeps it short of perfect. 

```bash
./build/panvar benchmark \
  -i "$OUT/panphorte/panphorte.normalized.sorted.gfa" \
  --bubble-prefix-in "$OUT/panphorte/panphorte" \
  --reference-path "$REF" \
  --variant-nodes "$OUT/call/call.variant_nodes.tsv" \
  -o "$OUT/call/benchmark"
```

Because a perfect reconstruction of a short region can't reach a high absolute QV, the length-fair headline is per-haplotype reconstruction identity (`1 − Σδ/ΣS`). `scripts/plot_benchmark_qv.R --category dist` draws it per locus, on a "nines" axis:

![Round-trip reconstruction identity](img/benchmark_qv_dist.png)

Every haplotype across every locus reconstructs at >99.9% identity. On LPA, the copy number lands exactly (the DUP module scores δ=0), so the residual is small sub-threshold variation — e.g. a ~32 bp insertion at one bubble that sits just under `--min-sv-bp=50`. CYP2D6 is near-perfect; GSTM1 sits lowest (its paralog stack is dense with small paralogous sequence variants) but is still >99.9%. This confirms the calls round-trip the haplotypes.

---

## Not just LPA — the other duplicated loci

The same pipeline is validated on C4 (C4A/C4B), GSTM1 (GSTM1/2/5), and CYP2D6 (CYP2D6/2D7/2D8P). `scripts/regen_results.sh` runs all of them and checks copy number against the pangene ground truth:

```bash
scripts/regen_results.sh              # all loci + validation + plots
Rscript scripts/plot_cn_correlation.R \
  --table results/cn_table.tsv \
  --out results/cn_correlation        # the two panels below
```

Total copy number per locus (called vs ground truth):

![Locus total-CN correlation](img/cn_correlation_loci.png)

Per-gene split for the paralog clusters:

![Per-gene CN correlation](img/cn_correlation_genes.png)

Every point on the diagonal is an exactly-correct call. Locus totals are exact (LPA 465/465, C4 131/131, GSTM1 159/159); the per-gene splits land C4A/C4B at 96.2%, CYP2D6 99.2%, CYP2D7 100%, and the whole GSTM cluster (GSTM1, GSTM2, GSTM4, GSTM5) at 466/466 each — the off-diagonal C4/CYP handful are gene-conversion mosaics and unannotated pseudogene/hybrid modules. GSTM3 is stable single-copy and sits outside the folded module, so it is not resolved.
