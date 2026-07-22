# Module `benchmark` - algorithm

Mechanism for the `benchmark` module. For usage/flags see [modules/benchmark.md](../modules/benchmark.md); References in [references.md](../references.md).

`benchmark` measures how much of each haplotype the caller's own output reconstructs. For a bubble it takes the reference `source → sink` walk and applies only the called divergences — where a haplotype diverges from the reference and that divergence is one of our calls, the haplotype's nodes are spliced in; everything else stays reference — then aligns that reconstruction against the haplotype's true sequence. The identity, `1 − δ/S` (edit distance `δ` over alignment length `S`), is the headline: the fraction of aligned sequence the calls reproduce. Only bubbles carrying at least one call are scored (or every bubble with `--all-bubbles`); a bubble the reference does not traverse is dropped, having no baseline.

## How it works

For each scored bubble, for each haplotype that traverses it:

### 1. Spell the truth

The haplotype's own canonical `source → sink` walk, orientation-aware — the sequence we are trying to reproduce.

### 2. Reconstruct from the calls

Node-align the haplotype walk against the reference walk on shared step tokens — the anchors, a longest common subsequence over the `(node id, orientation)` tokens. Between two consecutive anchors lies a divergent block, a maximal stretch where the two walks differ. Emit the haplotype's side of the block if it is a called block (any of its nodes are in the bubble's `variant_nodes.tsv` set — a divergence one of our calls explains), otherwise revert to the reference's side, so uncalled divergence (SNPs, sub-threshold indels) stays at reference.

### 3. Align and split the residual

Align the reconstruction against the truth with a global (Needleman–Wunsch) alignment, giving `δ` and `S` and contributing to identity `1 − δ/S`. Walking the edit path splits `δ` by contiguous non-match block size against `--min-sv-bp`: `sub_threshold_bp` (blocks below N, uncallable variation) and `over_threshold_bp` (blocks at or above N, callable-size misses). A comparability `QV = -10·log10(max(0.5, δ)/S)` is emitted alongside, with a ceiling `QV_max = 10·log10(2S)` for a perfect reconstruction and `qv_ratio = QV/QV_max` in (0, 1].

Aggregate per haplotype length-weighted (`Σδ`, `ΣS` over its bubbles) to give identity `= 1 − Σδ/ΣS` and the summed residual split (plus the comparability `QV`/`qv_ratio`/`quintile`). Per gene, `identity + Σsub/ΣS + Σover/ΣS = 1` — the anatomy `scripts/plot_benchmark.R` renders (left: reconstructed + residual; right: the residual's not-callable vs mis-called split).


## Worked trace

A locus with two bubbles, each 5 nodes of 10 bp; the reference spells nodes `1,2,3,4,5` then `6,7,8,9,10`. Take one haplotype H:

| bubble | reference walk | H's true walk |  call at this bubble |
|--------|----------------|---------------|-------------------------|
| A | `1,2,3,4,5` | `1,3,4,5` (node 2 deleted) | DEL of node 2 (node 2 ∈ `variant_nodes`) |
| B | `6,7,8,9,10` | `6,7,8′,9,10` (node 8′is a SNP) | none (8′ ∉ `variant_nodes`) |

**Bubble A** — anchors `1,3,4,5`. The one divergent block is between anchors `1` and `3`: reference side `[2]`, haplotype side `[]`. Node 2 is in `variant_nodes`, so it is a called block and is dropped. Reconstruction `1,3,4,5` = truth `1,3,4,5`. δ=0, S=40.

**Bubble B** — anchors `6,7,9,10`. The divergent block between anchors `7` and `9`: reference side `[8]`, haplotype side `[8′]`. Neither node is in `variant_nodes` (we did not call the SNP), so it is not a called block. Keep the reference side `[8]`. Reconstruction `6,7,8,9,10` (reference) vs truth `6,7,8′,9,10`: they differ only at node 8/8′, therefore δ=1, S=50. 

Then aggregate: `Σδ = 1`, `ΣS = 90`. identity `= 1 − 1/90 = 0.989` — H reconstructs at ~98.9%, the one uncalled SNP all that separates it from perfect. That 1 bp is a SNP below `--min-sv-bp`, so it is `sub_threshold` with zero `over_threshold`: the anatomy reads 98.9% Reconstructed / 1.1% Not-callable / 0% Mis-called. (Comparability: `QV = -10·log10(1/90) = 19.5`, `QV_max = 10·log10(180) = 22.6`, `qv_ratio = 0.87` — the raw QV looks middling only because the region is short.)
