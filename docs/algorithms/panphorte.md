# Module `panphorte` - algorithm

Mechanism for the `panphorte` module. For usage/flags see [modules/panphorte.md](../modules/panphorte.md); References in [references.md](../references.md#panphorte).

panphorte finds tandem-repeat bubbles and rewrites each as one repeat unit looped a copy-number of times, so downstream `call` reads copy number off a loop count rather than as a stack of insertions. The unit (its period) is the segment the array is built from; the collapsed form is a `REP` node carrying a self-loop, and a haplotype's copy number (`CN`) is the number of times its path traverses that self-loop. An interruption — a non-unit step interspersed in an array — is tolerated up to `--max-interruption-frac` and kept as a literal step between copies.

## How it works

### 1. Detect tandem arrays

Detection runs on the node-walk, tokenizing each `source → sink` step by the sequence it spells (orientation-aware), so identical copies match even when they are distinct node ids. For each crossing path, the unit period is established from a clean adjacent identical pair, then extended in both directions, collecting further copies and bridging short interruptions. A run is accepted when the unit length is at least `--min-unit-bp`, the copy count is at least `--min-copies`, and interruptions stay within `--max-interruption-frac` of the array span.

### 2. Gate on cohort prevalence

A bubble folds only if a `≥ --min-copies` array is carried by at least `--min-array-prevalence` of the haplotypes that traverse it. This separates a genuine population tandem repeat (carried by a large fraction of the cohort) from a rare or private duplication (carried by only a few haplotypes, and left for `call`). The default `0.5` asks for a majority.

### 3. Collapse to a REP node

Each accepted array collapses to one `REP` node — deduplicated by canonical unit sequence — carrying a self-loop, and every crossing path is rerouted through it once per copy, interruptions kept as literal steps. `call --cn` then reads copy number straight off the self-loop count.

Rerouting a path is not enough on its own. What the paths no longer walk has to leave the graph too, or the replaced route survives beside the normalized one: the re-snarl reads links, not paths, so it picks the dead branch back up as a second allele, and everything downstream — bubbles, calls, tangle detection, gene annotation — sees a site that is half normalized. Both the nodes and the oriented links a replaced span used are therefore marked, and each is dropped unless something still uses it. A node another site still walks is kept, but its arcs *here* are not: the link is what keeps the old route reachable. On LPA that removes 1508 nodes and 51 links; the fragments the rewrite adds in their place are few — 8 at LPA against 8974 folded copies, because a copy boundary that lands on a node edge needs no fragment at all.

A link is one arc read from either end — `a+ → b+` *is* `b- → a-`, and a GFA stores whichever of the two it happens to store — so every key is canonicalized over the pair. Keying on the stored direction alone had two consequences: an arc a path walked backwards was never recognised as obsolete and the whole replaced branch stayed reachable, and `ensure` re-added links the graph already held in the other direction, so a run that folded *nothing* still doubled the link set (C4 6532 → 12374, CYP2D6 8171 → 15606).

Some links cannot be dropped: another path still walks them, and a GFA link is global. At a tandem array that is ordinary rather than exceptional. LPA keeps 2009, and the reason is the fold itself — an arc inside the repeat unit is crossed once per copy, so one arc there is traversed 9436 times across the panel before folding and 466 times after, one crossing per haplotype falling outside the folded span. Judging any single link is therefore meaningless.

The question that is not meaningless: **after the rewrite, can a rewritten path still be walked along its original route?** If every node and link of a replaced span survives for other reasons, the site exists twice, folded and unfolded, and the re-snarl reports two alleles for one haplotype. That is checked directly on the delivered graph and is a hard error — an internal invariant, not a choice about how to fold, so there is no override. At LPA 0 of 466 original walks remain walkable, because 2.9% of each one's steps name nodes that are gone.

## Exact and approximate collapse

Exact collapse (`--min-similarity 1.0`, the default) is sequence-preserving: every path spells exactly the same sequence through the normalized graph (enforced internally). It only folds byte-identical copies, so it misses repeats whose copies differ by SNVs or indels.

Approximate collapse (`< 1.0`) folds divergent copies, lossily, via an in-process banded aligner ([edlib](https://github.com/Martinsos/edlib), bit-parallel edit distance):

- Seed one representative unit per bubble from the exact detector — it must come from a clean adjacent identical pair somewhere in the cohort, so the true period is recovered, and the most-supported unit wins.
- Find copies per path by aligning the unit (and its reverse complement) into the path's spelled sequence: k-mer anchors propose starts, and a banded fit alignment decides each copy's extent and orientation. The edit budget is `(1 − f)·|unit|`, uncapped, so a copy with a large internal indel still folds when `f` is low. Copies need not be adjacent; sequence between them is kept as literal nodes.
- Collapse the copies to one `REP` node traversed once per copy — a self-loop when adjacent, edges through the flanking literals otherwise. Within-copy SNVs and indels are discarded; per-copy orientation is preserved (`REP +,−,+`). A single copy of the seeded unit folds too.

## Seeding when node boundaries do not follow the repeat

Both routes above measure the period in node steps: the unit is a run of consecutive steps that repeats. That works when the graph splits nodes at repeat-unit boundaries, which is what a builder producing one node per distinct observed segment does. It fails when the boundaries fall elsewhere — a graph built by progressive alignment splits an array wherever chaining happened to break it, so a four-copy array may be spread over an arbitrary number of nodes with no step period at all, even though the base-level period is unambiguous. For a bubble where the node-step detector finds nothing anywhere, the unit is instead seeded at base level from the spelled sequence: candidate lags are proposed by k-mer occurrence, and each is accepted or rejected by comparing adjacent copies at that lag. Once a unit is seeded, the approximate route takes over unchanged. **The base-level fallback is reached only when `--min-similarity < 1.0`.** Exact mode measures the period in node steps and nowhere else, so a two-copy repeat stored as one node for the first copy and two nodes for the second is not folded at the default `--min-similarity 1.0`: at node-step level those two traversals are not identical, which is what exact mode means by a copy. A unit must occur as at least two adjacent identical copies in at least one haplotype to be seeded by the node-step route. A duplication with no adjacent identical pair anywhere (paralogs separated by other sequence, say) is left intact unless the base-level fallback seeds it.

## Worked trace

`U` = a 60 bp unit, `x` = an 8 bp interruption, `L`/`R` = flanks. Input path interior:

```text
step :  0    1    2    3    4    5    6
token:  L    U    U    x    U    U    R
bp   :  ..   60   60   8    60   60   ..
```

1. Establish the unit period from a clean adjacent identical pair. Anchoring at step 1 (`U`), the smallest period with an adjacent identical block is 1 — step 1 equals step 2, both `U` — so the repeat unit is `U` and the unit length `= 60`.
2. Extend the array by that period, tolerating short interruptions. Rightward: step 2 is another `U`; step 3 is `x`, a one-step/8 bp gap that is bridged because it stays within the interruption limits (`≤ unit_bp`); step 4 is `U` again, step 5 is `U`, and step 6 (`R`) ends the array. The copies are steps 1, 2, 4, 5 — the interruption is kept but not counted as a copy.
3. Accept the array against the thresholds: 4 copies ≥ `--min-copies` (2), `unit_bp` 60 ≥ `--min-unit-bp` (50), and the 8 bp interruption ≤ `--max-interruption-frac` of the array span (`0.25 · 248 = 62`). All pass, so this array is eligible to fold (subject to the cohort-prevalence gate above).
4. Collapse the copies onto one `REP` node with a self-loop and reroute the path through it once per copy, keeping the interruption as a literal step: the interior becomes `L REP REP x REP REP R` (`REP` traversed 4× total). Downstream `call --cn` reads the copy number straight off the `REP` self-loop count.

Resulting `panphorte.report.tsv` row:

```text
bubble_id  normalized  unit_bp  paths_normalized  min_copies  max_copies  interruptions_bp
1          yes         60       1                 4           4           8
```

The trace above is exact — copies must be byte-identical to seed. Real repeat arrays may differ between copies, so the exact test fails and nothing seeds; at `--min-similarity 0.90` the copies are instead found by alignment (e.g. identities 0.98, 0.95, 0.93, 0.97 give 4 copies and one `REP = U`), which can also bridge a large internal indel between divergent copies when the edit budget is wide enough. The approximate collapse is lossy: the per-copy small differences are discarded in favour of the single representative.

Two limits of the approximate detector, stated because `--min-similarity` reads like a sensitivity dial and is not one:

- **A copy is found only if it carries an exact 16-mer shared with the seed.** Identity is checked after anchoring, not before it, so a copy at 95% identity whose substitutions are spread such that every 16-mer window contains one is invisible at any `--min-similarity`. The threshold governs how much divergence is *accepted* once a copy is anchored, not how much can be *found*.
- **A copy boundary inside a node keeps its flanks.** Copy boundaries are base coordinates and rarely land on a node edge. Rounding to the nearest edge and replacing whole steps deletes the bases between the node edge and the copy edge — sequence *outside* the copy, which this mode may not discard. Declining such a copy instead is worse: at a tandem array a declined copy is a whole repeat unit missing from the copy number, while the bases at stake are a handful. So the containing step range is replaced by `fragment · REP · fragment`, where the fragments carry those outside bases verbatim as per-occurrence nodes. Nothing is lost and nothing is undercounted, and no global node split is needed — a fragment is local to the one path that needed it.

  Several copies can share one step range — two copies inside a single node have no step boundary between them — and are emitted as one block, `fragment · REP · fragment · REP · fragment`, so both are folded and the bases between them survive. A copy is declined only when its containing range collapses to nothing, which the detector's own bound on a copy (at least half a unit consumed) makes unreachable: `copies_declined_partial_boundary` reads 0 on every measured input, LPA KIV-2 included, where all 466 haplotypes fold 8974 copies in total, median 20 each. The guard behind it stays because the condition it prevents is severe — a haplotype left with no foldable copy reaches the site literally while its neighbours reach it through the REP, and copy number is counted from REP occurrences, so it reads CN 0 while carrying copies. Then the site is refused (`status=partial_boundary`) and `--allow-partial-boundary` overrides.

- **Copies are grouped into arrays by `--max-interruption-frac`.** The interrupting bases between successive copies are measured against the span of the group so far; past the fraction the group ends and a new one begins.

  The contract precisely: `--max-interruption-frac` decides whether copies form a qualifying array and so contribute to prevalence. Once the site is confirmed, every safely representable occurrence of the motif becomes a REP step — including copies in groups that did not qualify, because leaving them literal would recreate the mixed representation above. Copy number is therefore the number of REP occurrences, not the number of adjacent self-loop transitions, and `copies.tsv` is one row per (haplotype, site motif) rather than one row per array. `interruptions_bp` reports the bases accepted *inside* arrays.

## Differences from the original panphorte

panvar re-implements the [GenoGra/Panphorte](https://github.com/GenoGra/Panphorte) idea, deliberately:

| | original panphorte | this re-implementation |
|---|---|---|
| repeat detection | full base sequence (`O(bp²)`) | node-walk, compared by spelled sequence (`O(node²)`, ~1000× speed-up on big bubbles) |
| aligner | external | in-process banded aligner (edlib, approximate mode) |
| threading | — | per-haplotype parallel seed scan and copy detection (`--threads`) |
| modes | single collapse | explicit exact (sequence-preserving) vs approximate (lossy) via `--min-similarity` |
\