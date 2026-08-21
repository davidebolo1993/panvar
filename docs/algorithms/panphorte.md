# Module `panphorte` - algorithm

Mechanism for the `panphorte` module. For usage/flags see [modules/panphorte.md](../modules/panphorte.md); references in [references.md](../references.md#panphorte).

A tandem array spelled out copy by copy in a graph looks like a long insertion. Every extra copy is sequence one haplotype has and another does not, so a caller comparing walks reports a pile of insertions of unrelated sizes, and the quantity that actually varies — how many copies each haplotype carries — is nowhere in the output. Collapsing the array onto one repeat-unit node with a self-loop makes that quantity explicit: a haplotype's copy number becomes the number of times it traverses the loop, and a downstream caller reads it off directly.

## How it works

### 1. Detect the repeat unit

Within a bubble, each haplotype's interior is taken as a sequence of node steps. The unit is established from the smallest period at which two adjacent blocks of steps spell the same sequence; that block is the repeat unit and its spelled length is the unit length.

Where no step period exists, because the graph's node boundaries were set by alignment rather than by the repeat, the unit is seeded from the spelled bases instead. Both routes end with the same thing: a unit sequence and the positions of its copies.

### 2. Extend the array and accept it

From the seed, copies are extended in both directions at the detected period. Bases lying between successive copies are interruptions: they are kept but not counted as copies, and are tolerated while they stay within `--max-interruption-frac` of the span covered so far. Past that fraction the array ends and a new one may begin.

An array qualifies when it carries at least `--min-copies` copies and its unit spans at least `--min-unit-bp`.

### 3. Gate on how many haplotypes carry it

Folding rewrites the site for every haplotype, so it should only happen where the array is a property of the locus rather than of one sample. The fraction of the bubble's haplotypes carrying a qualifying array must reach `--min-array-prevalence`. Below that the site is left alone and reported as `below_prevalence`, and a rare private duplication is typed by `call` as an ordinary event instead.

### 4. Collapse onto a repeat-unit node

The copies are replaced by a single node carrying a self-loop, and each haplotype's walk is rerouted through it once per copy. Interruptions stay in the walk as literal steps, so the sequence between copies is not lost.

Copy boundaries are base coordinates and rarely land on a node edge. Rounding to the nearest edge would delete the bases between the edge and the copy boundary, which lie outside the copy. Instead the containing step range is replaced by the unit node flanked by fragment nodes carrying those bases verbatim, local to the haplotype that needed them. Several copies inside one node share a range and are emitted as one such block, so all of them fold and the bases between them survive.

At `--min-similarity 1.0` only byte-identical copies fold, and the graph still spells every haplotype exactly. Below that, copies are found by alignment and folded onto one representative, which discards the differences between them.

### 5. Re-snarl and emit

With `--reference-path`, the normalized graph is sorted along the reference and decomposed again, so the output is ready for `call` without an external tool. That re-snarl applies its own `--resnarl-min-variant-bp`, a separate decision from the one that produced the input sites. Bubble ids are reassigned by it.

Alongside the graph, the report gives one row per input site and the reason any site was not folded, and the provenance table records which site and motif each created unit node stands for, so a copy number can be traced back to where it came from.

## Worked trace

The steps below follow the five above, one for one. `U` is a 60 bp unit, `x` an 8 bp interruption, `L` and `R` the flanks. One haplotype's bubble interior:

```text
step :  0    1    2    3    4    5    6
token:  L    U    U    x    U    U    R
bp   :  ..   60   60   8    60   60   ..
```

1. Detect the repeat unit. Anchoring at step 1, the smallest period with an adjacent identical block is 1, since step 1 and step 2 spell the same sequence. The unit is `U` and the unit length is 60.

2. Extend the array and accept it. Rightward from step 2: step 3 is `x`, an 8 bp interruption, bridged because it stays within the tolerance; step 4 and step 5 are `U` again; step 6 is the flank and ends the array. The copies are steps 1, 2, 4 and 5, and the interruption is kept without being counted. Four copies clears `--min-copies` of 2, the 60 bp unit clears `--min-unit-bp` of 50, and the 8 bp interruption is well inside a quarter of the 248 bp the array spans.

3. Gate on how many haplotypes carry it. Suppose most haplotypes crossing this bubble carry a similar array, so prevalence clears `--min-array-prevalence` and the site is folded. Had only one haplotype carried it, the site would be left alone with `status` reading `below_prevalence`.

4. Collapse onto a repeat-unit node. The interior becomes the unit node traversed four times with the interruption still between the second and third, and the array's original nodes are removed. Copy number for this haplotype is now 4, readable as the number of traversals.

5. Re-snarl and emit. With a reference given, the normalized graph is sorted and decomposed again, and the site's row is written to the report:

```text
bubble_id  normalized  unit_bp  paths_normalized  min_copies  max_copies  interruptions_bp  status
1          yes         60       1                 4           4           8                 normalized
```

At `--min-similarity 1.0` this trace requires the copies to be byte-identical. Real arrays usually differ between copies, so nothing seeds and the site is reported `no_seed`; below 1.0 the copies are found by alignment instead and folded onto one representative, at the cost of the differences between them.

## Relationship to the original method

This is a re-implementation of the [Panphorte](https://github.com/GenoGra/Panphorte) idea, with detection over node walks rather than raw base sequence, an in-process aligner for the approximate mode, per-haplotype parallelism, and an explicit split between a sequence-preserving exact collapse and a lossy approximate one.
