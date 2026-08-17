# Module `benchmark` - algorithm

Mechanism for the `benchmark` module. For usage/flags see [modules/benchmark.md](../modules/benchmark.md); References in [references.md](../references.md).

`benchmark` answers two separable questions about the caller's own output: **what was there to be found** (the truth event ledger) and **how much of the sequence comes back** (three reconstructions plus a do-nothing baseline). Everything is derived from one decomposition of the two walks, so the ledger and the reconstructions cannot disagree about what a divergent block is.

Every reference-traversed bubble is scored by default. A bubble the reference does not traverse has no baseline and is dropped — counted in the `excluded` scope, never silently.

## One decomposition

For each scored bubble, for each haplotype that traverses it:

**Spell the truth.** The haplotype's own canonical `source → sink` walk, orientation-aware.

**Node-align.** Longest common subsequence over `(node id, orientation)` step tokens between the reference walk and the haplotype's walk. The matched positions are anchors; each maximal stretch between consecutive anchors where the two walks differ is one **divergent block**, holding `ref_steps[ri0,ri1)` against `hap_steps[hi0,hi1)`. On walks too large for the `O(nm)` table the whole bubble becomes one block — a coarser decomposition, not a wrong one, and counted as `bubbles_decomposed_coarsely`.

**Size it, and attribute it.** `ref_bp` and `hap_bp` sum the node lengths on each side; the event size is `max(ref_bp, hap_bp)`, which is how `call` gates a linked DEL/INS replacement — on the larger arm. The block is attributed to the first record (by variant id, so the choice cannot depend on file order) whose node set it actually contains. A block that merely sits in the same bubble as a record is not attributed to it.

**Event size is a property of the two walks alone.** It does not move with the alignment, the threshold, the calls or the reconstruction, which is what makes it usable as truth.

## The truth ledger

Each block with a non-zero size is one truth event:

| class | condition |
|-------|-----------|
| `below_threshold` | `size_bp < --min-sv-bp` |
| `called` | `size_bp ≥ --min-sv-bp` and a specific record's node is in the block |
| `missed` | `size_bp ≥ --min-sv-bp` and no record covers it |

`called` says a record covering this block was **emitted**. It does not say this haplotype was genotyped as carrying it — that is the carrier table, and the two genuinely differ: at ACOT bubble 7, 145 haplotypes have an above-threshold event attributed to a record and are genotyped `0` at every record of the bubble.

### Why not classify from the alignment

The residual used to be split by contiguous non-match run length in the reconstruction-vs-truth alignment: runs shorter than `--min-sv-bp` were called "variation that could not have been called", runs at or above it "callable-size events missed". That measures the edit path, not the event.

A clean 60 bp deletion of ordinary non-repetitive sequence, aligned globally by edlib, comes back as **fourteen runs of 1–10 bases** — the co-optimal edit path distributes the gap over chance matches, and every distribution costs the same 60 edits. Consequently `over_threshold_bp` read **0 at all six reference loci**: the bucket meant to hold missed callable-size events was structurally empty, and so was the carrier truth flag built on it. The columns survive as `resid_run_lt_bp` / `resid_run_ge_bp`, describing the residual's shape and nothing more.

## Three reconstructions

Each walks the reference and substitutes the haplotype's own steps at the blocks its rule accepts.

| level | substitutes a block when |
|-------|--------------------------|
| `graph` | any of its nodes is in the **union** of called nodes at that bubble |
| `called` | it is attributed to a **specific** record and `size_bp ≥ --min-sv-bp` |
| `genotype` | not block-based at all — see below |

`graph` is the optimistic upper bound and is what every QV figure previously quoted for this project measured. Sharing a node with some call is not matching one, so a called deletion containing a shared reference node can authorise copying a different, uncalled allele; and no genotype is read, so a call placed on entirely the wrong haplotypes still scores perfectly. `called` is the strict version: the reference with exactly the retained calls implanted and nothing else.

Each reconstruction is globally aligned (Needleman–Wunsch, edlib) to the truth, giving `δ` and `S`. Identity is `1 − Σδ/ΣS` length-weighted over a haplotype's bubbles; `QV = -10·log10(max(0.5, δ)/S)` with ceiling `QV_max = 10·log10(2S)` and `qv_ratio = QV/QV_max` are emitted for comparability.

## Worked trace

Two bubbles, each 5 nodes of 10 bp; the reference spells `1,2,3,4,5` then `6,7,8,9,10`. Haplotype H:

| bubble | reference walk | H's true walk | call at this bubble |
|--------|----------------|---------------|---------------------|
| A | `1,2,3,4,5` | `1,3,4,5` (node 2 deleted) | DEL of node 2 (node 2 ∈ `variant_nodes`) |
| B | `6,7,8,9,10` | `6,7,8′,9,10` (8′ is a SNP) | none (8′ ∉ `variant_nodes`) |

**Bubble A** — anchors `1,3,4,5`. One divergent block between anchors `1` and `3`: reference side `[2]` (10 bp), haplotype side `[]` (0 bp). Size `max(10,0) = 10`. At `--min-sv-bp 50` that is `below_threshold`; at `--min-sv-bp 5` it is `called`, since node 2 is a record's node. Both reconstructions that accept it drop node 2, giving `1,3,4,5` = truth. δ=0, S=40.

**Bubble B** — anchors `6,7,9,10`. Block between `7` and `9`: reference `[8]`, haplotype `[8′]`, both 10 bp, size 10, attributed to no record. `below_threshold` at 50 bp, `missed` at 5 bp. Neither reconstruction substitutes it, so the reconstruction is the reference `6,7,8,9,10` against truth `6,7,8′,9,10`: δ=1 (one substituted base — the nodes differ by a SNP), S=50.

Aggregate: `Σδ = 1`, `ΣS = 90`, identity `= 1 − 1/90 = 0.989`. The ledger says: at 50 bp, two below-threshold events and nothing missed; at 5 bp, one called and one missed. The identity is the same number in both cases — which is the point of keeping the two apart.

## The `genotype` level

With `--vcf`, a reconstruction that uses **only what the VCF says about this haplotype**, in sequence space, which is what a downstream consumer actually has.

**Join samples to haplotypes.** Each VCF sample column is one haplotype path, matched by exact name; genotypes are haploid and nothing is phased or imputed. Duplicate sample columns are refused — which haplotype a genotype belongs to would otherwise depend on column order. A partial join is accepted and reported prominently, with both directions counted, because a QV over a subset is not the QV of the run.

**Place each bubble on the reference.** VCF `POS` lives in the reference path's own forward coordinates, so the bubble's reference span is taken as a substring there rather than from its canonical walk. A bubble the reference crosses `sink → source` has its reconstruction reverse-complemented before scoring, so it meets the canonical truth in the same orientation.

**Apply the genotyped edits, right to left** — descending `POS`, so earlier edits do not shift later coordinates. `POS` is the anchor base and the event occupies `POS+1` onwards.

| type | reconstruction | exact? |
|------|----------------|--------|
| explicit `ALT` | replace the record's `REF` span with the allele the `GT` names | yes |
| `DEL` | erase `\|SVLEN\|` bases after the anchor | yes |
| `INS` | insert `INSSEQ` after the anchor | yes |
| `INV` | reverse-complement `POS+1..END` in place | yes — an inversion carries no new sequence, so the reference supplies all of it |
| `DUP` | tile an inferred reference span to the per-sample `CNBP` delta, or lay down `CN` copies of `RU_LEN` in place of `REF_CN` under `--dup-model cn` | **no** — counted `heuristic`. `CNBP` supplies the length change, not the inserted sequence, so what is reconstructed is the right length of approximately right sequence |

Records that cannot be laid down are counted rather than dropped — `unplaceable` (`POS` outside the bubble's span), `clamped` (the edit ran past the span end), `unhandled` (no rule, or no usable payload), `ref_mismatch` (the record's `REF` is not the reference sequence at `POS`, so the placement is wrong and every score downstream would be fiction) — and reported in `gt_records`.

**Score against the same truth, and against doing nothing.** The plain reference span with no edits is scored alongside as the **baseline**:

```
gap_closed          = (baseline_delta - genotype_delta) / (baseline_delta - graph_delta)
variation_recovered = (baseline_delta - <level>_delta)  /  baseline_delta
```

`gap_closed` measures against the achievable bound, `variation_recovered` against all the variation there was; reported for all three levels, the pair says whether representation or the graph is the limit. **`gap_closed` is `NA` when `baseline_delta == graph_delta`.** Those are equal exactly when nothing was called at the bubble, so the old convention of returning `1.0` on a zero denominator reported a total miss as having closed 100% of the gap.

The baseline is also the correctness check on the whole projection: a haplotype genotyped as carrying nothing applies no edits, so its genotype reconstruction must be byte-identical to the baseline. Any coordinate, orientation or span error breaks that equality, and it is asserted on real data rather than assumed.

**Carrier confusion.** Each (haplotype, bubble) contributes to a TP/FP/FN/TN table: called carrier is `GT ≥ 1` on any record at the bubble; true carrier is whether the ledger holds an event of at least `--min-sv-bp` there. Truth therefore comes from the walks, not from an alignment. Reported for `ALL` only — a bubble-level judgement has no single SV type, and fanning it out over every type in the bubble counted one observation several times; the per-type partition is the ledger, where each event maps to at most one record.
