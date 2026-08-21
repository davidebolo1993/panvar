# Walkthrough — one locus, end to end

A single worked run of the whole pipeline on one real locus: LPA, whose KIV-2 tandem array varies in copy number between haplotypes. It shows every module except `associate`, which has its own [GWAS example](gwas.md).

Every command runs against the committed test graph `tests/real_data/lpa.gfa.gz` and writes under `results/real_data/lpa/`. The numerical results and generated plots are reproduced by `scripts/regen_results.sh lpa`; the Bandage and interactive-viewer screenshots are illustrative views of those committed outputs.

```bash
GFA=tests/real_data/lpa.gfa.gz
REF="GRCh38#0#chr6:160509252-160734894"
OUT=results/real_data/lpa
```

A note on bubble ids before starting. `panphorte` and `refine` both re-decompose the graph after rewriting it, and each decomposition numbers the sites it finds. An id therefore identifies a site within one stage, not across stages: the array below is bubble 7 as `bubble` finds it and bubble 6 by the time `call` runs. Each step here says which stage its numbers come from.

---

## 1. `bubble` — find the variant sites

Sorts and flips the graph along the reference, then decomposes it into snarls: boundary-node pairs whose removal isolates an interior subgraph. On this locus that yields 12 sites across 466 haplotypes.

```bash
./build/panvar bubble \
  -i "$GFA" \
  -o "$OUT/bubble/bubble" \
  -r "$REF" \
  --gtf tests/real_data/Homo_sapiens.GRCh38.116.gtf.gz
```

The tandem array is site 7 at this stage. Its interior span runs from 12,270 bp to 184,333 bp across the panel, and it carries 459 distinct alleles among 466 haplotypes — the size range is the copy-number variation, and the allele count is what makes it hard.

Bandage view of `$OUT/bubble/bubble.sorted.gfa` coloured by `bubble.bandage_nodes.csv`:

![The locus graph in Bandage](img/lpa_bubble_bandage.png)

A linear reference backbone with the array ballooning into a dense tangle of parallel paths, one strand per copy-number allele.

---

## 2. `inspect` — cluster the haplotypes at a site

`inspect` is read-only and can be run at any stage. Before normalization, it exposes the 459 distinct array walks reported above and can cluster paths for closer inspection.

```bash
./build/panvar inspect \
  -i "$OUT/bubble/bubble.sorted.gfa" \
  --bubble-prefix-in "$OUT/bubble/bubble" --bubble-id 7 \
  --cluster --cluster-similarity 0.97 \
  -o "$OUT/inspect/inspect_before"
```

The matrix is large because every physical copy is still spelled out. The compact, countable view appears after `panphorte`.

---

## 3. `panphorte` — fold the tandem into a countable unit

The array is spelled out copy by copy in the graph, which would type as a pile of insertions. `panphorte` collapses it onto one repeat-unit node with a self-loop, so copy number becomes the number of traversals.

```bash
./build/panvar panphorte \
  -i "$OUT/bubble/bubble.sorted.gfa" \
  --bubble-prefix-in "$OUT/bubble/bubble" \
  -o "$OUT/panphorte/panphorte" \
  --reference-path "$REF" \
  --min-similarity 0.95 \
  --gtf tests/real_data/Homo_sapiens.GRCh38.116.gtf.gz
```

At site 7 it folds a 5,547 bp unit, rewriting all 466 traversing paths and collapsing 3,485 nodes, with copy counts spanning 1 to 32. Of those paths, 465 meet the array-carrier definition used for the prevalence calculation; the remaining path carries one copy. The other 11 sites are reported `no_seed` or `below_prevalence` and are left alone — the prevalence gate is what keeps a private duplication out of the fold.

![The normalized graph in Bandage](img/lpa_panphorte_bandage.png)

The tangle has become a single node with a self-loop onto itself. Compare with the step-1 view of the same site.

---

## 4. `inspect` again — on the folded graph

The same cluster view on the normalized graph confirms the fold kept the per-haplotype signal: the ladder is still there, now read off self-loop traversals instead of a long node run.

```bash
./build/panvar inspect \
  -i "$OUT/panphorte/panphorte.normalized.sorted.gfa" \
  --bubble-prefix-in "$OUT/panphorte/panphorte" --bubble-id 7 \
  --cluster --cluster-similarity 0.97 \
  -o "$OUT/inspect/inspect_panphorte"
```

```bash
Rscript scripts/plot_node_coverage_heatmap.R \
  --table "$OUT/inspect/inspect_panphorte.bubble_7.node_counts.tsv" \
  --node-lengths "$OUT/inspect/inspect_panphorte.bubble_7.node_lengths.tsv" \
  --cluster-by "$OUT/inspect/inspect_panphorte.bubble_7.clusters.tsv" \
  --out "$OUT/plots/lpa_node_heatmap_panphorte"
```

![The same site after folding](img/lpa_inspect_node_heatmap_panphorte.png)

Rows are representative haplotypes and columns are nodes. More copies mean more traversals of the repeat-unit node, so the copy-number ladder remains visible after folding.

---

## 5. `refine` — remove graph-builder artifacts

Where the graph builder split one true event across adjacent snarls, `call` would report the split as separate records. `refine` re-aligns the actual per-haplotype interior sequences and collapses such artifacts, holding folded repeat-unit nodes fixed so copy number is untouched.

```bash
./build/panvar refine \
  -i "$OUT/panphorte/panphorte.normalized.sorted.gfa" \
  --bubble-prefix-in "$OUT/panphorte/panphorte" \
  --reference-path "$REF" \
  -o "$OUT/refine/refine" \
  --gtf tests/real_data/Homo_sapiens.GRCh38.116.gtf.gz
```

Here 7 regions are rebuilt and 5 skipped, and the re-decomposition leaves 11 sites where there were 12. The array's region is among those skipped, by design: it carries the copy-number signal, and re-aligning it would linearize the copies. Losslessness is checked rather than assumed, so every haplotype still spells exactly what it spelled on the way in.

From here on, `call` and `describe` read the refine prefix.

---

## 6. `call` — type the variants and read copy number

Diffs every haplotype's walk against the reference's and types each difference.

```bash
./build/panvar call \
  -i "$OUT/refine/refine.normalized.sorted.gfa" \
  --bubble-prefix-in "$OUT/refine/refine" --reference-path "$REF" \
  -o "$OUT/call/call" \
  --cn --allele-vcf \
  --gtf tests/real_data/Homo_sapiens.GRCh38.116.gtf.gz
```

On this locus: 5 deletions, 9 insertions and 3 duplications over 466 haplotype columns.

The headline record is the array, `bubble6_DUP_5100`, on site 6 of the refined graph. It carries `REF_CN=6` and `RU_LEN=5547` with per-haplotype `CN` spanning 1 to 32 and 464 non-reference carriers. `FORMAT:CNBP` is the signed difference in module-walk bases relative to the reference; `CNRESID` records any difference from the ideal `(CN-REF_CN)×RU_LEN`. Its `CN_METHOD` is `REP`: the reported CN is the exact number of traversals of the folded unit in this graph. Because this run allowed 95% sequence similarity during folding, biological copy-number validity is established separately below rather than assumed from the topology.

The other two duplications happen to exercise the other two routes, which is worth knowing because they carry different guarantees:

| record | `CN_METHOD` | what a copy means |
|---|---|---|
| `bubble6_DUP_5100` | `REP` | exact traversal count of the represented repeat unit |
| `bubble5_DUP_3916` | `MODULE_BP` | bases across a collapsed module, divided by a reference-calibrated unit |
| `bubble8_DUP_7468` | `PEAK` | the highest multiplicity any interior node reaches — heuristic, and labelled `CN_CONFIDENCE=HEURISTIC` |

Draw the calls as an oncoprint-style map:

```bash
Rscript scripts/plot_vcf_map.R \
  --vcf "$OUT/call/call.region.vcf" \
  --reference-path "$REF" \
  --out "$OUT/plots/lpa_vcf_map" \
  --max-paths 50
```

![The called variants across haplotypes](img/lpa_vcf_map.png)

Rows are haplotypes, columns are calls grouped by site. The array is the wide band shaded by copy number; deletions, insertions and inversions take their own colours. That gradient is the copy-number spectrum the pipeline exists to recover.

For a per-node view of how each haplotype traverses the graph under every call, build the viewer bundle and launch the app:

```bash
Rscript scripts/build_variant_node_data.R \
  --gfa "$OUT/refine/refine.normalized.sorted.gfa" \
  --vcf "$OUT/call/call.region.vcf" \
  --variant-nodes "$OUT/call/call.variant_nodes.tsv" \
  --out "$OUT/plots/lpa_variant_nodes.rds"

VN_RDS="$OUT/plots/lpa_variant_nodes.rds" Rscript scripts/variant_node_heatmap_app.R
```

![The interactive node-coverage viewer](img/lpa_node_coverage.png)

---

## 7. `describe` — features for association

Turns the called sites into per-haplotype feature matrices on three substrates: k-mers and node/edge dosages spelled from the graph, and variant-level dosages read from the VCF.

```bash
./build/panvar describe \
  -i "$OUT/refine/refine.normalized.sorted.gfa" \
  --bubble-prefix-in "$OUT/refine/refine" \
  --out-dir "$OUT/describe" \
  --variant-nodes "$OUT/call/call.variant_nodes.tsv" \
  --no-wide-matrix
```

A copy-number locus is first-class here: the variant substrate carries the duplication's `CN` as its dosage rather than a presence flag, so the association tests the copy number itself. Pass `--samples` with a sample-to-haplotype table to get diploid per-sample matrices as well.

The association step is covered in the [GWAS example](gwas.md).

---

## 8. `benchmark` — how much of the sequence comes back

Reconstructs each haplotype from the caller's own output and compares it against the truth walk.

```bash
./build/panvar benchmark \
  -i "$OUT/refine/refine.normalized.sorted.gfa" \
  --bubble-prefix-in "$OUT/refine/refine" \
  --reference-path "$REF" \
  --variant-nodes "$OUT/call/call.variant_nodes.tsv" \
  --vcf "$OUT/call/call.region.vcf" \
  -o "$OUT/benchmark/benchmark"
```

Four levels are reported and they are not interchangeable. The first three implant the haplotype's own true block and therefore bound what the graph and the retained records could achieve; only the fourth reconstructs what the VCF actually says.

![Reconstruction anatomy across the reference loci](img/benchmark_qv.png)

- Graph ceiling: near 100% at every locus, because it substitutes the true block wherever a call overlaps the truth event's node space. This is a discovery/representation ceiling, not the region VCF reconstruction.
- From the VCF alone: 85.3% identity at this locus, 99.7% at the least varied of the six.
- Where the loss lives: the residual partitioned into five terms that sum to it exactly.
- Variation found: the residual split into sub-threshold and missed. Missed is zero everywhere.

The third panel is the one to read first. `Not found` is 0.00% at all six loci: every above-threshold truth event overlaps at least one emitted call. That is a discovery statement, not a guarantee that every carrier genotype is correct. Carrier false negatives remain at ACOT and CYP2D6, while most residual bases come from how the compact region VCF represents overlapping or merged calls.

That is also why the allele VCF closes the gap completely. `call --allele-vcf` spells every allele out instead of describing it, and reconstructing from it gives 0 bp residual on all six loci. Both figures live in `results/reconstruction.tsv`; neither replaces the other.

Two metrics are reported and they answer different questions. Identity is a per-base fraction; gap closed is the fraction of the reference-to-truth distance the VCF removes. At this locus they read 85.3% and 75.5% for the same run, so neither should be quoted without its name.

---

## Not just one locus

`scripts/regen_results.sh` runs the same pipeline on six loci and compares copy number with committed assembly-derived truth tables. These checks are independent of the GWAS phenotype simulation, but they are validation on the six supplied loci rather than a universal accuracy guarantee.

![Called copy number against truth, by locus](img/cn_correlation.loci.png)

At LPA every one of the 465 haplotypes with a truth value is called exactly right. This validates the 95%-similarity fold on this dataset; the same conclusion should not be inferred for a new locus without suitable truth. Where paralogs carry usable private k-mers, the module total can also be split per gene:

![Called copy number against truth, by gene](img/cn_correlation.genes.png)

Points off the diagonal are disagreements between the private-marker estimate and the annotation-derived truth. Gene-conversion mosaics, hybrid modules, ambiguous markers and annotation differences can all contribute. Genes outside the represented module are not split out.

---

## `rebuild` — re-inducing a tangled locus

A minority of graphs are too fragmented to decompose sanely: transitive closure over all-pairs alignments can collapse sequence that recurs across a locus onto shared, very high-degree nodes. `rebuild` re-induces such a locus by progressive graph generation before `bubble` sees it.

```bash
./build/panvar rebuild -i tangled.gfa -o rebuilt.gfa
```

It gates first and passes healthy graphs through untouched, so it is not part of the run above. On a locus that does trip the gate:

![A tangled locus before re-induction](img/myom2.original.png)
![The same locus after re-induction](img/myom2.rebuild.png)

The rebuild is accepted only if every haplotype comes back within the identity and coverage bounds; otherwise the original graph is written unchanged and the audit says which haplotype failed and why.
