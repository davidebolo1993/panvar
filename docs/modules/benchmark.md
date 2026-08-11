# Module `benchmark`

CLI: `panvar benchmark`

## What it does

Scores the caller's own output with a round-trip reconstruction, on two levels. Both cover the same bubbles, compare against the same truth and use the same QV scale; they differ only in what the reconstruction is allowed to use.

| level | reconstruction built from | answers |
|-------|---------------------------|---------|
| `graph` | reference walk + **the haplotype's own steps** at blocks the calls explain | can the graph represent this haplotype, and did the caller flag the divergent blocks |
| `genotype` | reference sequence + **only the edits this haplotype's `GT` names**, from the VCF alone | can a consumer reconstruct this sample from the VCF |

`--vcf` accepts either VCF `call` produces. Give it `<prefix>.region.vcf` to score the merged, interpreted output; give it `<prefix>.alleles.vcf` (from `call --allele-vcf`) to score the lossless one. Records whose `ALT` is literal sequence are applied as a plain `REF` to `ALT` substitution; symbolic ones (`<DEL>`, `<INS>`, `<DUP>`) are reconstructed from `SVLEN` / `INSSEQ` / `CNBP`. Running both and comparing is the point: the difference is exactly what the merged representation costs.

`graph` reads no genotype, so it cannot be wrong about *which* haplotype carries what. It is an upper bound on reconstruction, not a genotyping score — a locus can score near-perfect there while every genotype is on the wrong haplotype. `genotype` (enabled by `--vcf`) is the genotyping score: a missed carrier keeps reference and a spurious one edits sequence that was already correct, so both error directions cost bases.

A third **baseline** reconstruction — plain reference, no edits applied — is always computed alongside `genotype`, and the headline is how much of the distance between that baseline and the `graph` bound the VCF actually closes:

```
gap_closed = (baseline_delta - genotype_delta) / (baseline_delta - graph_delta)
```

`1.0` means the VCF recovers everything the graph could offer, `0.0` means it is worth nothing over doing nothing, and a negative value means the genotypes are placing edits where they do not belong. The count of haplotypes scoring *worse than baseline* is reported for the same reason. This also makes the metric self-checking: a haplotype genotyped as carrying nothing must score byte-identically to the baseline, so a coordinate or orientation error in the projection shows up immediately instead of hiding inside a plausible-looking QV.

### The per-level reconstruction

For each called bubble and each haplotype that traverses it, `graph`:
- reconstructs the haplotype on the passed graph by taking the reference walk and substituting only the divergent blocks the calls explain — node-align the haplotype's walk against the reference walk and, between shared anchors, keep the haplotype's block when its nodes are in `variant_nodes.tsv`, else revert to reference. Uncalled variation (SNPs, sub-threshold indels) stays at reference;
- aligns the reconstruction to the haplotype's true walk over the bubble to get the edit distance `δ` and alignment length `S`.

**Reconstruction identity `= 1 − Σδ/ΣS`** (length-weighted over the bubbles a haplotype traverses). `δ` stays above zero wherever a haplotype carries variation `panvar` did not call, so identity is `1.0` only when the called events fully reconstruct every bubble, and drops in linear proportion to the uncalled sequence. Report it per haplotype and aggregated per gene.

**Residual anatomy:** The reconstruction-vs-truth alignment path is walked and the residual `δ` is split by contiguous non-match block size against `--min-sv-bp`: a block shorter than the threshold is variation that could not have been called (sub-threshold SNPs/indels), a block `≥` the threshold is a callable-size event missed or mis-represented. This separates a real caller problem from a metric artifact — e.g. a locus whose only variation is a small indel scores a low raw identity over its tiny denominator, but the split shows it as 100% sub-threshold, not a miss. Aggregated per gene, `Reconstructed (= identity) + Not-callable (sub-threshold) + Mis-called (≥ threshold)` sum to 100% of the aligned sequence. Results can be plotted via `scripts/plot_benchmark.R` — two side-by-side plots: left stacks Reconstructed + Residual as % of the aligned sequence; right splits that residual into Not-callable vs Mis-called, also as % of aligned sequence but with the y-axis auto-scaled to the largest residual.

A logarithmic `QV = -10·log10(max(0.5, δ)/S)` with its ratio `qv_ratio = QV / QV_max` (`QV_max = 10·log10(2S)`) binned into quintiles is also emitted.

Algorithm and worked trace: [algorithms/benchmark.md](../algorithms/benchmark.md).

## Required inputs

- `-i, --gfa <graph.gfa>` — the same passed graph the calls were made on.
- one of `-b, --bubble-prefix-in <prefix>` or `-c, --bubbles-csv-in <path>`.
- `--variant-nodes <path>` — `call`'s `<prefix>.variant_nodes.tsv`.
- `-r, --reference-path <name>` — the diff baseline (full name or unique case-insensitive substring).
- `-o, --out-prefix <prefix>`.

Optional: `--vcf <path>` — `call`'s `<prefix>.region.vcf`, which turns on the `genotype` level. VCF sample columns are joined to graph paths by exact name; each column is one haplotype, so genotypes are haploid and nothing has to be phased. Without it only `graph` is scored.

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `--vcf <path>` | add the genotype-aware level and the baseline it is measured against | off (`graph` only) |
| `--all-bubbles` | score every bubble in the CSV, not just called ones — a haplotype diverging from reference at an uncalled bubble then surfaces as an identity drop / over-threshold miss | off (called bubbles only) |
| `--min-sv-bp <N>` | threshold for the residual split — match the `call` run. Residual blocks `< N` bp are sub-threshold (uncallable) variation, `≥ N` bp are callable-size misses | `50` |
| `--threads <N>` | worker threads over haplotypes | `0` (auto) |
| `-q, --quiet` | disable progress logs | — |

## Outputs

| file | contents |
|------|----------|
| `<prefix>.qv.tsv` | per (haplotype, bubble): `svtypes`, `is_carrier` (a called block was applied), `δ`, `aln_len`, the residual split `sub_threshold_bp` / `over_threshold_bp`, and `qv`. With `--vcf`, also `gt_delta` / `gt_aln_len` / `gt_qv`, the baseline `ref_delta`, and `gt_called_carrier` / `gt_true_carrier` |
| `<prefix>.qv_by_haplotype.tsv` | per-haplotype `identity` and the residual split `sub_threshold_bp` / `over_threshold_bp` (the headline columns `plot_benchmark.R` reads), plus the comparability columns `qv`, cosigt `band`, `qv_max`, `qv_ratio`, `quintile`. With `--vcf`, the `gt_*` mirror of those plus `ref_sum_delta`, `ref_qv` and `gap_closed` |
| `<prefix>.qv_summary.tsv` | the `residual` scope splitting total un-reconstructed bp into `sub_threshold` (`< min-sv-bp`) vs `over_threshold` (`≥ min-sv-bp`); plus, for comparability, `% haplotypes` per `qv_ratio` quintile / cosigt band (`scope` = `quintile`/`haplotype`) and `% observations` per band per SV class (`scope` = `svclass`). With `--vcf`: `gt_quintile` / `gt_haplotype` (the same distributions for the genotype level), `gt_gap` (pooled `baseline_delta` / `genotype_delta` / `graph_delta`, `gap_closed_pooled`, `gap_closed_mean`, `worse_than_baseline`), `gt_carrier` (TP/FP/TN/FN of the carrier call per SV class, with precision and recall), and `gt_records` (how many VCF records were `applied` / `unplaceable` / `clamped` / `unhandled`, so a genotype QV is never read without knowing what it was computed from) |