# Module `benchmark` - algorithm

Mechanism for the `benchmark` module. For usage/flags see [modules/benchmark.md](../modules/benchmark.md); References in [references.md](../references.md).

### The one idea

Every number in this module comes from a single comparison: line up a haplotype's path through a bubble against the reference's path, and see where they part company. Each place they part company is a divergent stretch. That one decomposition feeds both the ledger and every reconstruction, so the two can never disagree about what counts as an event.

From there the module answers two separable questions — what was there to be found (the truth event ledger) and how much of the sequence comes back (four reconstructions plus a do-nothing baseline).

Every reference-traversed bubble is scored by default. A bubble the reference does not traverse has no baseline and is dropped — counted in the `excluded` scope, never silently.

## How it works

### 1. Find where each haplotype differs from the reference

For each scored bubble, for each haplotype that traverses it:

Spell the truth. The haplotype's own `source → sink` walk, orientation-aware. This is what everything is scored against.

Line the two walks up. Both walks are sequences of steps, each step being a node together with the direction it is read in. A longest-common-subsequence match over those steps finds the steps the two walks share; those shared steps are anchors. Between two consecutive anchors, wherever the two walks are not identical, is one divergent stretch — a piece of the reference walk on one side and a piece of the haplotype walk on the other, either of which may be empty.

> On walks too large for the quadratic alignment table, the whole bubble is treated as one stretch. That is a coarser decomposition, not a wrong one, and it is counted as `bubbles_decomposed_coarsely` rather than passed off as a fine one.

Size it. Add up the node lengths on each side, and take the event's size to be the larger of the two. That matches how `call` gates a deletion-plus-insertion replacement: on the larger arm.

Attribute it. A stretch belongs to a record when it actually contains one of that record's nodes. If several qualify, the lowest variant id wins, so the choice can never depend on file order. A record that merely sits in the same bubble is not attributed anything.

The size of an event is a property of the two walks and nothing else. It does not move when the alignment, the threshold, the calls or the reconstruction change — which is exactly what makes it usable as truth.

### 2. Classify the truth events

Each divergent stretch with a non-zero size is one truth event:

| class | condition |
|-------|-----------|
| `below_threshold` | `size_bp < --min-sv-bp` |
| `called` | `size_bp ≥ --min-sv-bp` and a specific record's node is in the stretch |
| `missed` | `size_bp ≥ --min-sv-bp` and no record covers it |

`called` says a specific emitted record shares at least one node with this stretch. It does not say the record spans the stretch, represents it correctly, or that this haplotype was genotyped as carrying it. The two negatives therefore differ in strength: `missed` is firm — nothing emitted touches the event — while `called` is an upper bound on discovery, not a genotyping statement.

#### Why not classify from the alignment

Contiguous mismatch-run lengths from a sequence alignment cannot define event size: co-optimal alignments can distribute one deletion over several short runs through chance matches. The ledger therefore sizes events from the two graph walks. `resid_run_lt_bp` and `resid_run_ge_bp` remain descriptive summaries of residual shape only.

### 3. Reconstruct at each level

Three of the four levels share one mechanism: walk the reference, and at some of the divergent stretches, paste in the haplotype's own true sequence. They differ in nothing but which stretches qualify.

| level | pastes a stretch when | what that makes it |
|-------|-----------------------|--------------------|
| `graph` | any of its nodes is one that some call at that bubble also touches | loosest — mere overlap with the union of called nodes |
| `called` | it is attributed to a specific record, and reaches `--min-sv-bp` | the same idea made strict |
| `carrier` | as `called`, and this haplotype's `GT` names an overlapping record | …with genotypes finally taken into account |
| `region_vcf` | — does not work this way at all; see step 4 | the only level that reads the records' own sequence |

Because all three paste in the truth, none of them tests whether a record encodes its event correctly. They measure where the caller pointed. A record can be flawlessly placed and describe nonsense, and all three still score it as recovered. That is deliberate: it separates finding an event from writing it down properly, and step 4 is where the second question gets asked.

The three are progressively stricter, and each strictness costs something specific:

- `graph` does not read genotypes and does not require a match. Touching a node that a call also touches is not the same as agreeing with that call, so a correctly called deletion sharing one reference node can license pasting an entirely different, uncalled allele nearby. And because no `GT` is consulted, a call placed on completely the wrong samples still scores perfectly here. Read it as a discovery ceiling: could the graph hold this haplotype, and did the caller flag the right neighbourhood?
- `called` requires the stretch itself to map to a record, and to clear the size threshold. Genotypes still play no part.
- `carrier` adds the one thing `called` ignores: whether the VCF says this haplotype carries anything there. That information exists nowhere else — `variant_nodes.tsv` has one row per record holding the union of nodes over every merged event and every carrier, with no haplotype column — so carrier status is knowable only from the VCF, joined by record id.

That join exists only for the region VCF, and the mode is reported as `vcf_mode`. Region record ids must equal the `variant_nodes` id set and agree on `BUBBLE_ID`. An allele VCF instead has one `bubbleN_ALLELES` record per bubble and no per-call id, so `carrier` and the per-call `loss_bp` terms are reported `NA` / `not_applicable`. Anything matching neither contract is refused.

Within region mode, the stretch is pasted when any overlapping record has `GT ≥ 1`: one divergent region can be described by several records and different carriers take different ones, so testing only the primary attribution would be testing which record sorted first, not whether the haplotype was called.

A fifth reconstruction, the in-scope floor, pastes every stretch reaching `--min-sv-bp` regardless of any call. Its residual is purely sub-threshold variation, so it separates out-of-scope sequence from discovery failure.

The consecutive differences along `truth → in-scope floor → called → carrier → region_vcf` are the `loss_bp` partition, and the run asserts that they sum exactly to the region-VCF residual (`sum_check`).

This is a partition, not a descending chain. Only `in-scope floor ≥ called ≥ carrier` are nested ceilings — each substitutes a subset of the previous level's stretches, all implanting true sequence. `region_vcf` sits outside that nesting because it applies every record the haplotype carries rather than only attributed eligible stretches, so it can beat `carrier` locally; the last two terms are signed for exactly that reason.

Two rules make the partition mean what it says:

One population. Every comparative total is over the common set — the `(haplotype, bubble)` observations all levels could score. `graph`/`called` cover every graph haplotype, `carrier` only VCF-joined ones, `region_vcf` only those whose bubble also placed. Totalling each over its own population reported `carrier` above `called` and a negative loss when one haplotype was left out of the VCF.

A false positive is not a representation failure. Where the haplotype has no eligible truth event and the VCF still edits it, `carrier` has no true stretch to paste and leaves the reference alone while `region_vcf` applies the erroneous edit. That difference is `false_positive_damage`, split out per observation from `representation`. The two are signed, because `region_vcf` applies every record the haplotype carries while `carrier` only substitutes attributed eligible stretches — so a record can improve a region `carrier` left alone.

Each reconstruction is globally aligned (Needleman–Wunsch, edlib) to the truth, giving `δ` and `S`. Identity is `1 − Σδ/ΣS` length-weighted over a haplotype's bubbles; `QV = -10·log10(max(0.5, δ)/S)` with ceiling `QV_max = 10·log10(2S)` and `qv_ratio = QV/QV_max` are emitted for comparability.


### 4. Reconstruct from the VCF alone

#### The `region_vcf` level

With `--vcf`, a reconstruction that uses only what the VCF says about this haplotype, in sequence space, which is what a downstream consumer actually has.

Join samples to haplotypes. Each VCF sample column is one haplotype path, matched by exact name; genotypes are haploid and nothing is phased or imputed. Duplicate sample columns are refused — which haplotype a genotype belongs to would otherwise depend on column order. A partial join is accepted and reported prominently, with both directions counted, because a QV over a subset is not the QV of the run.

Place each bubble on the reference. VCF `POS` lives in the reference path's own forward coordinates, so the bubble's reference span is taken as a substring there rather than from its canonical walk. A bubble the reference crosses `sink → source` has its reconstruction reverse-complemented before scoring, so it meets the canonical truth in the same orientation.

Apply the genotyped edits, right to left — descending `POS`, so earlier edits do not shift later coordinates. `POS` is the anchor base and the event occupies `POS+1` onwards.

| type | reconstruction | exact? |
|------|----------------|--------|
| explicit `ALT` | replace the record's `REF` span with the allele the `GT` names | yes |
| `DEL` | erase `\|SVLEN\|` bases after the anchor | yes |
| `INS` | insert `INSSEQ` after the anchor | yes |
| `INV` | reverse-complement `POS+1..END` in place | yes — an inversion carries no new sequence, so the reference supplies all of it |
| `DUP` | tile an inferred reference span to the per-sample `CNBP` delta, or lay down `CN` copies of `RU_LEN` in place of `REF_CN` under `--dup-model cn` | no — counted `heuristic`. `CNBP` supplies the length change, not the inserted sequence, so what is reconstructed is the right length of approximately right sequence |

Records that cannot be laid down are counted rather than dropped — `unplaceable` (`POS` outside the bubble's span), `clamped` (the edit ran past the span end), `unhandled` (no rule or usable payload), and `ref_mismatch` (the record's `REF` does not match the reference at `POS`) — and reported in `gt_records`.

Score against the same truth, and against doing nothing. The plain reference span with no edits is scored alongside as the baseline:

```
gap_closed          = (baseline_delta - region_vcf_delta) / (baseline_delta - graph_delta)
variation_recovered = (baseline_delta - <level>_delta)  /  baseline_delta
```

`gap_closed` measures against the achievable bound, while `variation_recovered` measures against all variation. `gap_closed` is `NA` when `baseline_delta == graph_delta`, because its denominator is then zero.

The baseline is also the correctness check on the whole projection: a haplotype genotyped as carrying nothing applies no edits, so its genotype reconstruction must be byte-identical to the baseline. Any coordinate, orientation or span error breaks that equality, and it is asserted on real data rather than assumed.

Carrier confusion. Each (haplotype, bubble) contributes to a TP/FP/FN/TN table: called carrier is `GT ≥ 1` on any record at the bubble; true carrier is whether the ledger holds an event of at least `--min-sv-bp` there. Truth therefore comes from the walks, not from an alignment. Reported for `ALL` only — a bubble-level judgement has no single SV type, and fanning it out over every type in the bubble counted one observation several times; the per-type partition is the ledger, where each event maps to at most one record.

### 5. Partition the residual and compare against a baseline

The genotype residual is split into five consecutive terms that sum to it exactly: variation below the size threshold, eligible events no record covers, events covered but not genotyped onto this haplotype, sequence a record covers and carries but does not reproduce, and edits applied where there is no eligible truth event. The last two are signed, since the genotype level applies every record the haplotype carries rather than only attributed eligible stretches and can therefore beat the stretch-based ceiling locally.

Every comparative figure is taken over the common set, the haplotype-and-bubble observations all levels could score, so two totals are never compared over different populations.

A do-nothing baseline, the plain reference with no edits, is the denominator of `gap_closed` and `variation_recovered`. It is also the metric's own correctness check: a haplotype genotyped as carrying nothing must reconstruct byte-identically to it, so a coordinate or orientation error shows up immediately rather than hiding inside a plausible score.

## Worked trace

The steps below follow the five above, one for one. The first bubble contains a 60 bp node; the second contains a 20 bp replacement whose two sequences differ by one base. Haplotype H:

| bubble | reference walk | H's true walk | call at this bubble |
|--------|----------------|---------------|---------------------|
| A | `1,2,3,4,5` | `1,3,4,5`, 60 bp node 2 deleted | a deletion of node 2, whose node is in the call's node set |
| B | `6,7,8,9,10` | `6,7,8',9,10`, a 20 bp substitution | none, since node 8' is in no call's node set |

1. Decompose each haplotype against the reference. At bubble A the shared anchors are `1,3,4,5`, leaving one divergent stretch between anchors `1` and `3`: the reference side holds node 2 and the haplotype side nothing. At bubble B the anchors are `6,7,9,10`, leaving a stretch between `7` and `9` holding one node on each side.

2. Classify the truth events. Under a 50 bp threshold, bubble A is a `called` 60 bp truth event and bubble B is a 20 bp `below_threshold` event. Nothing eligible was missed.

3. Reconstruct at each level. Bubble A's stretch shares a node with an attributed record and reaches the threshold, so graph, called and carrier substitute its true steps. Bubble B's stretch is below threshold and shares a node with no record, so every ceiling leaves the reference sequence there.

4. Reconstruct from the VCF alone. The genotype level applies only the edits H's genotype names. The deletion at bubble A is emitted and carried, so it is applied; nothing is emitted at bubble B, so that substitution is left as reference.

5. Partition the residual and compare against a baseline. The only remaining mismatch is the one changed base in the 20 bp substitution. It belongs to `out_of_scope`, not discovery, because the truth event itself is below the reporting threshold.
