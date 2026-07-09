# Module `benchmark` - algorithm

Mechanism for the `benchmark` module. For usage/flags see [modules/benchmark.md](../modules/benchmark.md); References in [references.md](../references.md).

## Terms

- **reconstruction** — the reference walk over a bubble with only the called divergences applied: where the haplotype diverges from reference and that divergence is a call, the haplotype's nodes are spliced in; everything else stays reference.
- **anchor** — a `source→sink` step (node id and orientation) shared by the reference and the haplotype walk, used to pin the two walks together (a longest common subsequence over the step tokens).
- **divergent block** — a maximal stretch between two anchors where reference and haplotype differ (a reference-only run for a deletion, a haplotype-only run for an insertion, or both for a substitution/inversion).
- **called block** — a divergent block whose reference-side or haplotype-side nodes appear in the bubble's `variant_nodes.tsv` set — i.e. a divergence one of our calls explains.
- **δ / S** — edit distance and alignment length from the global (Needleman–Wunsch) alignment of the reconstruction against the true haplotype sequence.
- **QV** — `-10·log10(max(0.5, δ)/S)`. **QV_max** — `10·log10(2S)`, the ceiling for a perfect (`δ≤0.5`) reconstruction of length `S`. **qv_ratio** — `QV/QV_max` ∈ (0,1], length-fair. **identity** — `1 − δ/S`, the linear counterpart.

## Algorithm

Score only bubbles that carry ≥1 call (or every bubble in the CSV with `--all-bubbles`). For each such bubble, get the reference `source→sink` walk; a bubble the reference does not traverse is dropped (no baseline). Then for each haplotype that traverses the bubble:

1. Spell the truth: the haplotype's own canonical `source→sink` walk, orientation-aware. This is the sequence we are trying to reproduce;
2. Node-align the haplotype walk against the reference walk on shared step tokens (LCS). Between consecutive anchors, emit the haplotype's block if it is a called block (any of its nodes ∈ the bubble's `variant_nodes` set), otherwise revert to the reference's block. Uncalled divergence (SNPs, sub-threshold indels) therefore stays at reference. Spell the result — again orientation-aware;
3. Align reconstruction vs truth → `δ`, `S`; `QV = -10·log10(max(0.5, δ)/S)`.

Aggregate per haplotype length-weighted (`Σδ`, `ΣS` over its bubbles, then the QV formula), and derive `qv_ratio`/`quintile`/`identity`. 

The reconstruction is entirely coordinate-free — it works on step tokens and node-set membership — so a reverse-oriented bubble (see [call's coordinate handling](call.md)) needs no special treatment.

## Worked trace

A locus with two bubbles, each 5 nodes of 10 bp; the reference spells nodes `1,2,3,4,5` then `6,7,8,9,10`. Take one haplotype **H**:

| bubble | reference walk | H's true walk | our call at this bubble |
|--------|----------------|---------------|-------------------------|
| A | `1,2,3,4,5` | `1,3,4,5` (node 2 deleted) | **DEL of node 2** (node 2 ∈ `variant_nodes`) |
| B | `6,7,8,9,10` | `6,7,8′,9,10` (node 8 → alt node 8′, a SNP) | **none** (8′ ∉ `variant_nodes`) |

**Bubble A** — anchors `1,3,4,5`. The one divergent block is between anchors `1` and `3`: reference side `[2]`, haplotype side `[]`. Node 2 is in `variant_nodes`, so it is a called block and is dropped. Reconstruction `1,3,4,5` = truth `1,3,4,5`. δ=0, S=40.

**Bubble B** — anchors `6,7,9,10`. The divergent block between anchors `7` and `9`: reference side `[8]`, haplotype side `[8′]`. Neither node is in `variant_nodes` (we did not call the SNP), so it is not a called block. Keep the reference side `[8]`. Reconstruction `6,7,8,9,10` (reference) vs truth `6,7,8′,9,10`: they differ only at node 8/8′, therefore δ=1, S=50. 

Then aggregate: `Σδ = 1`, `ΣS = 90`. `QV = -10·log10(1/90) = 19.5`; `QV_max = 10·log10(180) = 22.6`; `qv_ratio = 0.87` (top quintile) `0.8-1.0`; `identity = 1 − 1/90 = 0.989`. So H reconstructs at ~98.9% — the one uncalled SNP is all that separates it from perfect, and the length-fair views (`qv_ratio`, `identity`) both place it near the ceiling, whereas the raw QV (19.5) looks middling only because the region is short.
