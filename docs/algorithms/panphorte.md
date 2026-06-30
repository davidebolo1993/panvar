# panphorte — algorithm & worked example

Mechanism and a hand-traced example for **Module 2**. For usage/flags see
[modules/panphorte.md](../modules/panphorte.md); citations: [references.md](../references.md#panphorte).

## How it works

Detection runs on the **node-walk**, tokenizing each step by the sequence it spells (orientation-aware), so
identical copies match even when they are distinct node ids. For each bubble, for each crossing path:

1. tokenize the source→sink walk by spelled-sequence hash;
2. find non-overlapping tandem arrays: establish the **unit period** from a clean adjacent identical pair,
   then extend bidirectionally collecting further copies, tolerating short interruptions (kept as nested
   literal steps). Accept a run when `unit_bp ≥ --min-unit-bp`, `copies ≥ --min-copies`, and interruptions
   stay within `--max-interruption-frac`;
2b. **cohort-prevalence gate** — only fold the bubble if a `≥ --min-copies` array is carried by
   `≥ --min-array-prevalence` of the bubble-traversing haplotypes. This is what separates a genuine
   population VNTR (folded; e.g. LPA KIV-2, carried by ~all haplotypes) from a **rare private duplication
   of a gene/segmental module** (left untouched). A single haplotype with a `CYP2D6×2` tandem, or the
   C4/RCCX module present in a minority, no longer collapses the whole bubble — those are gene-level
   events the `call` module resolves per-gene/total downstream, and folding them here would fuse paralogs
   (e.g. CYP2D6 with CYP2D7/2D8P). Default `0.5` (a real VNTR is carried by a majority);
3. collapse each array to one **REP node** (deduplicated by canonical unit sequence) with a self-loop;
   reroute the path through it `copies` times; drop now-unused duplicate nodes.

**Exact (`--min-similarity 1.0`, default)** is sequence-preserving — every path spells exactly the same
sequence through the normalized graph (enforced internally). **Approximate (`< 1.0`)** is intentionally
lossy (below).

## Approximate collapse (`--min-similarity < 1.0`)

Exact mode misses divergent repeats (copies differing by SNVs/indels). Approximate mode does single-block,
lossy collapse via a built-in banded aligner:

1. **Seed** one representative unit per bubble from the exact detector (must come from a clean adjacent
   identical pair somewhere in the cohort, so the true period is recovered; most-supported unit wins).
2. **Find copies** per path by aligning the unit (and its reverse complement) into the path's spelled
   sequence: k-mer anchors propose starts, a banded global alignment decides each copy, its extent, and
   orientation. Band width `(1 − f)·|unit|` (uncapped), so a copy with a large internal indel still aligns
   when `f` is low. Copies need not be adjacent; sequence between them is kept as literal nodes.
3. **Collapse** copies to one REP node traversed once per copy (self-loop when adjacent, edges through
   flanking literals otherwise). Lossy: within-copy SNVs/indels are discarded. Per-copy orientation is
   preserved (`REP +,−,+`). Single copies of the seeded unit fold too (`copies = 1`).
4. **Limitation:** the unit must occur as ≥ 2 adjacent identical copies in ≥ 1 haplotype to be seeded. A
   duplication with no adjacent identical pair anywhere (e.g. paralogs separated by other sequence) is left
   intact — recovered downstream by `call` from node-traversal multiplicity. See
   [modules/call.md](../modules/call.md).

## Worked trace — tandem-repeat collapse

`U` = a 60 bp unit, `x` = an 8 bp interruption, `L`/`R` = flanks. Input path interior:

```text
step :  0    1    2    3    4    5    6
token:  L    U    U    x    U    U    R
bp   :  ..   60   60   8    60   60   ..
```

1. **Unit period.** Anchor `i=1` (`U`); smallest period `p` with an adjacent identical block is `p=1`
   (`block_equal(1,2,1)`: U==U). So `unit_bp = 60`.
2. **Extend, tolerating interruptions.** `copy_starts = {1}` walking right by `p`:
   ```text
   1→2 : U==U ✓                                   {1,2}
   2→3 : token[3]=x ≠ U ✗ → gap g=1: candidate 4, gap_bp=8 ≤ 60 ✓ and U==U ✓   {1,2,4}  (x = interruption)
   4→5 : U==U ✓                                   {1,2,4,5}
   5→6 : R ≠ U → stop
   ```
   A gap is bridged only while ≤ `kMaxGapSteps` steps **and** cumulative `gap_bp ≤ unit_bp`.
3. **Accept test.** `copies=4 ≥ 2` ✓; `unit_bp=60 ≥ 50` ✓; interruption `8 ≤ 0.25·(4·60+8)=62` ✓ → collapse.
4. **Collapse.** One REP node carries `U` with a self-loop; the path becomes `L REP REP x REP REP R` (REP
   traversed 4× total), `x` kept as a literal step. Report row: `unit_bp=60 copies=4 interruption_bp=8`.
   Downstream `call` reads CN straight off the REP self-loop multiplicity (the always-on self-loop path, no
   flag; `--cn-from-multiplicity` is for folded bubbles that have *no* self-loop).

**Exact vs approximate.** The trace above is exact (`block_equal` requires byte-identical copies). Real
repeats (KIV-2) differ by SNVs, so the exact pair test fails and nothing seeds. At `--min-similarity 0.90`
copies are found by alignment (`identity 0.98, 0.95, 0.93, 0.97 → 4 copies → one REP=U`); because alignment
bridges a copy, a large internal indel between divergent copies is also bridged when `f` is low enough
(this lets a C4 short module fold onto the long one at `--min-similarity 0.70`). The collapse is lossy: the
per-copy small events are discarded in favour of the single representative.

## Differences from the original panphorte

panvar re-implements the [GenoGra/Panphorte](https://github.com/GenoGra/Panphorte) idea, deliberately:

| | original panphorte | this re-implementation |
|---|---|---|
| repeat detection | full base sequence (`O(bp²)`) | node-walk, compared by spelled sequence (`O(node²)`, ~1000× less on big bubbles) |
| aligner | external | built-in banded aligner (approximate mode) |
| threading | — | per-haplotype parallel seed scan + copy detection (`--threads`) |
| modes | single collapse | explicit exact (sequence-preserving) vs approximate (lossy) via `--min-similarity` |
