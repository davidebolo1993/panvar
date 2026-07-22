# Module `benchmark`

CLI: `panvar benchmark`

## What it does

Scores the caller's own output with a round-trip reconstruction. For each called bubble and each haplotype that traverses it:
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

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `--all-bubbles` | score every bubble in the CSV, not just called ones — a haplotype diverging from reference at an uncalled bubble then surfaces as an identity drop / over-threshold miss | off (called bubbles only) |
| `--min-sv-bp <N>` | threshold for the residual split — match the `call` run. Residual blocks `< N` bp are sub-threshold (uncallable) variation, `≥ N` bp are callable-size misses | `50` |
| `--threads <N>` | worker threads over haplotypes | `0` (auto) |
| `-q, --quiet` | disable progress logs | — |

## Outputs

| file | contents |
|------|----------|
| `<prefix>.qv.tsv` | per (haplotype, bubble): `svtypes`, `is_carrier` (a called block was applied), `δ`, `aln_len`, the residual split `sub_threshold_bp` / `over_threshold_bp`, and `qv` |
| `<prefix>.qv_by_haplotype.tsv` | per-haplotype `identity` and the residual split `sub_threshold_bp` / `over_threshold_bp` (the headline columns `plot_benchmark.R` reads), plus the comparability columns `qv`, cosigt `band`, `qv_max`, `qv_ratio`, `quintile` |
| `<prefix>.qv_summary.tsv` | the `residual` scope splitting total un-reconstructed bp into `sub_threshold` (`< min-sv-bp`) vs `over_threshold` (`≥ min-sv-bp`); plus, for comparability, `% haplotypes` per `qv_ratio` quintile / cosigt band (`scope` = `quintile`/`haplotype`) and `% observations` per band per SV class (`scope` = `svclass`) |