# Module `benchmark`

CLI: `panvar benchmark`

## What it does

Scores the caller's own output with a round-trip QV. For each called bubble and each haplotype that traverses it:
- reconstructs the haplotype on the passed graph by taking the reference walk and substituting only the divergent blocks the calls explain — node-align the haplotype's walk against the reference walk and, between shared anchors, keep the haplotype's block when its nodes are in `variant_nodes.tsv`, else revert to reference. Uncalled variation (SNPs, sub-threshold indels) stays at reference. This is coordinate-free, so bubble orientation never matters;
- aligns the reconstruction to the haplotype's true walk over the bubble (edlib global/NW);
- scores `QV = -10·log10(max(0.5, δ)/S)`, where `δ` is the edit distance and `S` the alignment length.

Per-haplotype QV is the length-weighted combine over the bubbles it traverses (`Σδ / ΣS`, then the formula). `δ` stays above zero wherever a haplotype carries variation we did not call (so an accurate call set still scores high but finite); `δ` reaches zero only when the called events fully reconstruct the bubble. No VCF is needed — `variant_nodes.tsv` says which nodes each call explains, and the haplotype walks come from the graph.

Algorithm and worked trace: [algorithms/benchmark.md](../algorithms/benchmark.md).

**Length-fair headline (quintiles).** The `max(0.5, δ)` floor makes the achievable ceiling depend on region length: a perfect reconstruction of total length `S` maxes out at `QV_max = 10·log10(2S)`, so a haplotype spanning only small bubbles can never reach the high cosigt bands however good its calls are. So we also report `qv_ratio = QV / QV_max` — length-fair, `1.0` for a perfect reconstruction of *any* size — and bin it into quintiles (`0.0–0.2` … `0.8–1.0`, top is best). An accurate call set piles into the top quintile regardless of locus size; a haplotype drops only in proportion to the uncalled sequence it carries *relative to its region*. Because QV is logarithmic, `qv_ratio` is sensitive: a handful of uncalled SNPs over tens of kb still moves it down a quintile. The fixed cosigt bands (`<17`, `17–23`, `23–33`, `>33`) are kept alongside for external comparability.

## Required inputs

- `-i, --gfa <graph.gfa>` — the same passed graph the calls were made on.
- one of `-b, --bubble-prefix-in <prefix>` or `-c, --bubbles-csv-in <path>`.
- `--variant-nodes <path>` — `call`'s `<prefix>.variant_nodes.tsv`.
- `-r, --reference-path <name>` — the diff baseline (full name or unique case-insensitive substring).
- `-o, --out-prefix <prefix>`.

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `--all-bubbles` | score every bubble in the CSV, not just called ones — a haplotype diverging from reference at an *uncalled* bubble then surfaces as a low-QV miss (turns the metric into a recall/sensitivity view) | off (called bubbles only) |
| `--threads <N>` | worker threads over haplotypes | `0` (auto) |
| `-q, --quiet` | disable progress logs | — |

## Outputs

| file | contents |
|------|----------|
| `<prefix>.qv.tsv` | per (haplotype, bubble): `svtypes`, `is_carrier` (a called block was applied), `δ`, `aln_len`, `qv` |
| `<prefix>.qv_by_haplotype.tsv` | per-haplotype length-weighted `qv`, cosigt `band`, `qv_max`, `qv_ratio`, `quintile` |
| `<prefix>.qv_summary.tsv` | `% haplotypes` per `qv_ratio` quintile and per cosigt band (`scope` = `quintile`/`haplotype`), plus `% observations` per band per SV class (`scope` = `svclass`) |

## Notes

- Read the `qv_ratio` quintile distribution as the length-fair headline and the per-bubble `qv.tsv` as a diagnostic (a perfectly reconstructed but short bubble scores a modest *per-bubble* QV because of the `max(0.5, δ)` floor, even though it is a perfect reconstruction).
- Sensitivity: rerun `call --min-sv-bp <huge>` to drop records and re-benchmark with `--all-bubbles` — the dropped events are then scored against reference and QV collapses, confirming the metric tracks missed calls.
