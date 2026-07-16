# Module `benchmark` - algorithm

Mechanism for the `benchmark` module. For usage/flags see [modules/benchmark.md](../modules/benchmark.md); References in [references.md](../references.md).

## Terms

- **reconstruction** — the reference walk over a bubble with only the called divergences applied: where the haplotype diverges from reference and that divergence is a call, the haplotype's nodes are spliced in; everything else stays reference.
- **anchor** — a `source→sink` step (node id and orientation) shared by the reference and the haplotype walk, used to pin the two walks together (a longest common subsequence over the step tokens).
- **divergent block** — a maximal stretch between two anchors where reference and haplotype differ (a reference-only run for a deletion, a haplotype-only run for an insertion, or both for a substitution/inversion).
- **called block** — a divergent block whose reference-side or haplotype-side nodes appear in the bubble's `variant_nodes.tsv` set — i.e. a divergence one of our calls explains.
- **δ / S** — edit distance and alignment length from the global (Needleman–Wunsch) alignment of the reconstruction against the true haplotype sequence.
- **identity** — `1 − δ/S`, the headline metric: the fraction of aligned sequence the calls reconstruct.
- **QV** (comparability) — `-10·log10(max(0.5, δ)/S)`; **QV_max** — `10·log10(2S)`, the ceiling for a perfect (`δ≤0.5`) reconstruction of length `S`; **qv_ratio** — `QV/QV_max` ∈ (0,1].

## Algorithm

Score only bubbles that carry ≥1 call (or every bubble in the CSV with `--all-bubbles`). For each such bubble, get the reference `source→sink` walk; a bubble the reference does not traverse is dropped (no baseline). Then for each haplotype that traverses the bubble:

1. Spell the truth: the haplotype's own canonical `source→sink` walk, orientation-aware. This is the sequence we are trying to reproduce;
2. Node-align the haplotype walk against the reference walk on shared step tokens (LCS). Between consecutive anchors, emit the haplotype's block if it is a called block (any of its nodes ∈ the bubble's `variant_nodes` set), otherwise revert to the reference's block. Uncalled divergence (SNPs, sub-threshold indels) therefore stays at reference. Spell the result — again orientation-aware;
3. Align reconstruction vs truth → `δ`, `S`, contributing to identity (`1 − δ/S`). Walking the edit path also splits `δ` by contiguous non-match block size against `--min-sv-bp` into `sub_threshold_bp` (blocks `< N`, uncallable variation) and `over_threshold_bp` (blocks `≥ N`, callable-size misses). A comparability `QV = -10·log10(max(0.5, δ)/S)` is emitted alongside.

Aggregate per haplotype length-weighted (`Σδ`, `ΣS` over its bubbles) to give identity `= 1 − Σδ/ΣS` and the summed residual split (plus the comparability `QV`/`qv_ratio`/`quintile`). Per gene, `identity + Σsub/ΣS + Σover/ΣS = 1` — the anatomy `scripts/plot_benchmark.R` renders (left: reconstructed + residual; right: the residual's not-callable vs mis-called split).

The reconstruction is entirely coordinate-free — it works on step tokens and node-set membership — so a reverse-oriented bubble (see [call's coordinate handling](call.md)) needs no special treatment.

## Worked trace

A locus with two bubbles, each 5 nodes of 10 bp; the reference spells nodes `1,2,3,4,5` then `6,7,8,9,10`. Take one haplotype H:

| bubble | reference walk | H's true walk |  call at this bubble |
|--------|----------------|---------------|-------------------------|
| A | `1,2,3,4,5` | `1,3,4,5` (node 2 deleted) | DEL of node 2 (node 2 ∈ `variant_nodes`) |
| B | `6,7,8,9,10` | `6,7,8′,9,10` (node 8′is a SNP) | none (8′ ∉ `variant_nodes`) |

**Bubble A** — anchors `1,3,4,5`. The one divergent block is between anchors `1` and `3`: reference side `[2]`, haplotype side `[]`. Node 2 is in `variant_nodes`, so it is a called block and is dropped. Reconstruction `1,3,4,5` = truth `1,3,4,5`. δ=0, S=40.

**Bubble B** — anchors `6,7,9,10`. The divergent block between anchors `7` and `9`: reference side `[8]`, haplotype side `[8′]`. Neither node is in `variant_nodes` (we did not call the SNP), so it is not a called block. Keep the reference side `[8]`. Reconstruction `6,7,8,9,10` (reference) vs truth `6,7,8′,9,10`: they differ only at node 8/8′, therefore δ=1, S=50. 

Then aggregate: `Σδ = 1`, `ΣS = 90`. identity `= 1 − 1/90 = 0.989` — H reconstructs at ~98.9%, the one uncalled SNP all that separates it from perfect. That 1 bp is a SNP below `--min-sv-bp`, so it is `sub_threshold` with zero `over_threshold`: the anatomy reads 98.9% Reconstructed / 1.1% Not-callable / 0% Mis-called. (Comparability: `QV = -10·log10(1/90) = 19.5`, `QV_max = 10·log10(180) = 22.6`, `qv_ratio = 0.87` — the raw QV looks middling only because the region is short.)
