# Module `panphorte` - algorithm

Mechanism for the `panphorte` module. For usage/flags see [modules/panphorte.md](../modules/panphorte.md); References in [references.md](../references.md#panphorte).

## Terms

- **unit/period** — the tandem repeat (TR) segment an array is built from.
- **interruption** — non-unit steps interspersed in an array, tolerated up to `--max-interruption-frac` and kept as nested literal steps.
- **`REP` node** — the single self-looping node an array collapses to; a haplotype's copy number (`CN`) is the number of self-loop traversals.

## How it works

Detection runs on the node-walk, tokenizing each step by the sequence it spells (orientation-aware), so identical copies match even when they are distinct node ids. For each bubble, for each crossing path:

1. tokenize the `source→sink` walk by spelled-sequence hash;
2. find non-overlapping tandem arrays: establish the unit period from a clean adjacent identical pair, then extend bidirectionally collecting further copies, tolerating short interruptions (kept as nested literal steps). Accept a run when the TR unit length  `≥ --min-unit-bp`, the number of copies `≥ --min-copies`, and interruptions stay within `--max-interruption-frac`;
3. apply a cohort-prevalence gate: fold the bubble only if a `≥ --min-copies` array is carried by `≥ --min-array-prevalence` of the bubble-traversing haplotypes. This separates a genuine population TR (carried by a large fraction of the cohort) from a rare/private duplication (carried by only a few haplotypes). Default `0.5` (a real TR is carried by a majority);
4. collapse each array to one `REP` node (deduplicated by canonical unit sequence) with a self-loop; reroute the path through it `copies` times.

Exact (`--min-similarity 1.0`, default) is sequence-preserving — every path spells exactly the same sequence through the normalized graph (enforced internally), but misses divergent repeats (copies differing by SNVs/indels). 
Approximate mode (`< 1.0`) does single-block, lossy collapse via an in-process banded aligner ([edlib](https://github.com/Martinsos/edlib), bit-parallel edit distance):

- Seed one representative unit per bubble from the exact detector (must come from a clean adjacent identical pair somewhere in the cohort, so the true period is recovered; most-supported unit wins).
- Find copies per path by aligning the unit (and its reverse complement) into the path's spelled sequence: k-mer anchors propose starts, a banded fit alignment decides each copy, its extent, and orientation. Edit budget `(1 − f)·|unit|` (uncapped): a copy folds while its edits (mismatches + indels) stay within it, so a copy with a large internal indel still aligns when `f` is low. Copies need not be adjacent; sequence between them is kept as literal nodes.
- Collapse copies to one `REP` node traversed once per copy (self-loop when adjacent, edges through flanking literals otherwise). Lossy: within-copy SNVs/indels are discarded. Per-copy orientation is preserved (`REP +,−,+`). Single copies of the seeded unit fold too.

A limitation worth considering is that the unit must occur as `≥ 2` adjacent identical copies in `≥ 1` haplotype to be seeded. A duplication with no adjacent identical pair anywhere (e.g. paralogs separated by other sequence) is left intact.

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

## Differences from the original panphorte

panvar re-implements the [GenoGra/Panphorte](https://github.com/GenoGra/Panphorte) idea, deliberately:

| | original panphorte | this re-implementation |
|---|---|---|
| repeat detection | full base sequence (`O(bp²)`) | node-walk, compared by spelled sequence (`O(node²)`, ~1000× speed-up on big bubbles) |
| aligner | external | in-process banded aligner (edlib, approximate mode) |
| threading | — | per-haplotype parallel seed scan and copy detection (`--threads`) |
| modes | single collapse | explicit exact (sequence-preserving) vs approximate (lossy) via `--min-similarity` |
