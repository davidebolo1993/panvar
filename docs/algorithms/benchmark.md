# Module `benchmark` - algorithm

Mechanism for the `benchmark` module. For usage/flags see [modules/benchmark.md](../modules/benchmark.md); References in [references.md](../references.md).

`benchmark` measures how much of each haplotype the caller's own output reconstructs, at two levels that share every other ingredient — same bubbles, same truth, same alignment, same QV scale — and differ only in what the reconstruction may draw on. The `graph` level below is the original one; the `genotype` level is described in its own section at the end.

For a bubble it takes the reference `source → sink` walk and applies only the called divergences — where a haplotype diverges from the reference and that divergence is one of our calls, the haplotype's nodes are spliced in; everything else stays reference — then aligns that reconstruction against the haplotype's true sequence. The identity, `1 − δ/S` (edit distance `δ` over alignment length `S`), is the headline: the fraction of aligned sequence the calls reproduce. Only bubbles carrying at least one call are scored (or every bubble with `--all-bubbles`); a bubble the reference does not traverse is dropped, having no baseline.

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

## The `genotype` level

`graph` splices in **the haplotype's own steps** wherever a call explains the divergence. That is deliberate — it isolates "can the graph hold this haplotype, and did the caller find the divergent blocks" from any question of who carries what. But it also means no genotype is ever read, so the level is structurally incapable of being wrong about assignment: a call placed on entirely the wrong haplotypes still scores perfectly. It is an upper bound.

With `--vcf`, a second reconstruction uses **only what the VCF says about this haplotype**, in sequence space, which is what a downstream consumer actually has.

### 1. Join samples to haplotypes

Each VCF sample column is one haplotype path, matched to the graph by exact name, so genotypes are haploid and nothing is phased or imputed.

### 2. Place each bubble on the reference

VCF `POS` lives in the reference path's own forward coordinates, so the bubble's reference span is taken as a substring there rather than from its canonical `source → sink` walk. A bubble the reference crosses `sink → source` has its reconstruction reverse-complemented before scoring, so it meets the canonical truth in the same orientation.

### 3. Apply the genotyped edits, right to left

For each record at the bubble with `GT ≥ 1`, applied in descending `POS` so earlier edits do not shift later coordinates. `POS` is the anchor base and the event occupies `POS+1` onwards:

| type | reconstruction |
|------|----------------|
| `DEL` | erase `\|SVLEN\|` bases after the anchor |
| `INS` | insert `INSSEQ` after the anchor |
| `DUP` | a copy-number record: adjust the span by the sample's own `CNBP`, repeating the span when it grows |

Records that cannot be laid down are counted rather than silently dropped — `unplaceable` (`POS` outside the bubble's reference span), `clamped` (the edit ran past the span end), `unhandled` (no rule for the type, or the record carries no usable payload) — and reported in `gt_records`, so a genotype QV is never read without knowing what produced it.

### 4. Score against the same truth, and against doing nothing

The result is aligned to the haplotype's true walk exactly as in step 3 of the `graph` level. A third reconstruction — the plain reference span with no edits — is scored alongside as the **baseline**, giving

```
gap_closed = (baseline_delta - genotype_delta) / (baseline_delta - graph_delta)
```

the fraction of the available distance the VCF actually closes, with `worse_than_baseline` counting haplotypes the genotypes actively harm.

The baseline is also the correctness check on the whole projection: a haplotype genotyped as carrying nothing at a bubble applies no edits, so its genotype reconstruction must be **byte-identical** to the baseline. Any coordinate, orientation or span error breaks that equality, and it is asserted on real data rather than assumed.

### 5. Carrier confusion

Separately from sequence, each (haplotype, bubble) contributes to a TP/FP/FN/TN table: called carrier is `GT ≥ 1` on any record at the bubble; true carrier is whether the baseline-versus-truth alignment leaves a residual block of at least `--min-sv-bp`. The size threshold matters — plain walk inequality would score a single uncalled SNP as a missed carrier, which is not something the caller set out to emit — and it also means precision and recall here move with `--min-sv-bp` and should be read against it.
