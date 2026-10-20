# M10-D: COSIGT-style sequence-similarity evaluation with the existing exact edit distance (frozen before execution)

**Scope (Davide).** Structural similarity means **sequence similarity**. Small substitutions and indels count
proportionately. The evaluator does not identify or classify events. The structural-alphabet and event-comparator
investigation is closed; a clarification is appended to the M10-B, B.1 and C records. M10-A remains the computational
checkpoint.
- **Not done:** no predictions, GQ, confidence, filters, priors or caller code are changed. No cutoff, cluster or posterior
  pooling is defined.
- **Branch:** `m10-d-sequence-similarity`. **Results:** `tests/results/m10d_sequence_similarity_20261020/`.

## 1. Contract

For predicted haplotypes A, B and truth C, D:
- `D1 = d(A,C) + d(B,D)`, `D2 = d(A,D) + d(B,C)`, `D = min(D1, D2)`.
- **One assignment** for the whole evaluated region, chosen by summed edits (never by QV, never per block).
- **Tie** (D1 = D2): assignment 1 (A–C, B–D), recorded.

**The distance:**
- `d` is the exact global (NW) edit distance, `panvar::nw_edit_distance` (edlib). Every base is accounted for: no local
  clipping.
- L1 and L2 are its alignment-column counts (the benchmark's S) for the selected assignment, from edlib's deterministic
  path. Co-optimal paths can differ in column count; that is recorded, not resolved.

**The metrics:**
- `error_rate = D/(L1+L2)`: unfloored, exactly 0 for identical sequences, NA if L1 + L2 = 0;
- `identity = 1 − error_rate`;
- `QV = −10 log10(max(0.5, D)/(L1+L2))`, or 60 when L1 + L2 = 0.

QV is the benchmark's `qv_from` (`benchmark_command.cpp`), with its length-weighted haplotype combination (Σδ/ΣS).
`qv_from` is in an anonymous namespace, so the expression is reused verbatim and checked on edge cases.

**Reported:** both totals, the selected assignment and tie, per-haplotype edits and alignment lengths, aggregate metrics, the
scope, and orientation.

## 2. Scope, endpoints and orientation

**Truth** comes through the existing evaluation path: `project_path_onto_locus` via `genotype_probe block-sequences
--truth-gfa --donor`, all blocks, reused from M10-C (`m10c.truth`).
- It does not require truth to be panel-spellable.
- The built model's fingerprint must equal the M10-B.1 dump's.
- It expresses truth in the model's **canonical chain orientation**, the frame in which the predictions are spelled. So
  whole-haplotype orientation is matched by construction.
- **Verified, not inferred:** per haplotype, a whole-sequence check counts truth 31-mers found in the prediction as given
  against in its whole reverse complement. Internal pieces are never flipped.

**Scope:**
- **Evaluated:** the **chain interior** (block 1 .. block n−2), which is the region with matched endpoints. It is named as
  such, never "the complete locus".
- **Not scored:** the flank-inclusive sequences (each path's own assembly extent). Their lengths are reported with status
  UNSUPPORTED, unmatched endpoints.

## 3. Contract gate (`m10d.py gate`)

**Fixtures:** M10-C's frozen sequences, reused, plus the missing checks:
- swapped truth homologues;
- whole-sequence reverse-complemented truth;
- an empty truth pair, all-empty, and empty haplotypes matched by the assignment.

**Checks:**
- identical DNA gives D = 0, error exactly 0, and the floored QV;
- swapped homologues give the same D under the other assignment;
- reverse complement gives a much larger D, and the orientation check prefers the given strand;
- small differences: 0 < D ≤ 9;
- unequal lengths give exactly the length-difference bound: F03 55, F04 600, F05 1,200 at equal total length, F11 55, F14
  60;
- order change and inversions: the exact distance is reported, not named events;
- tie recorded; empty handling;
- alignment columns lie within [max(len), len_a + len_b] and are ≥ edits (no clipping);
- **independent textbook DP oracle:** 400 tiny pairs, including inversions, order changes and empty sequences, and 60 tiny
  diploids against the brute-force assignment minimum;
- QV normalization edge cases.

## 4. Frozen development cases (the M10-C set; saved predictions; verified truth inputs)

| case | predictions | truth (donor) |
|---|---|---|
| C4 allocation failure | M1-B reported 3+1, saved 2+2 | `c4.gfa` NA19650 |
| C4 control | M10-A `D.c4_5` | `c4.gfa` HG00096 |
| CYP2D6 dosage | M1-B 5× 2+2, 30× 3+2, saved 2+3 | `cyp2d6.gfa` NA19650 |
| ACOT L4 / L2 | M9-H JP on and off | `acot.gfa` HG02083 / HG02074 |

- The M1-B 30× call is **also** the CYP2D6 synonymous-labelling case (M9-I).
- Truth GFA hashes and the fingerprint match are recorded.
- Unsupported comparisons stay in the table with their reason.
- **Earlier per-block metrics** (M1-B excess, M9-H interval distance) are a different evaluation unit. They are shown in a
  separate column, **never subtracted**, and no new floor search is run.

## 5. Execution

**Units (manifest P):** `build`, `gate.contract`, `truth.<5>`, `cases`.
- **Admission:** spent + declared timeout ≤ 900 s, including 10 s of disclosed pre-freeze development.
- **Watchdog (frozen) per unit:** its declared timeout, and a **6 GiB envelope** on the summed resident size of its process
  tree, polled every 0.2 s. A breach kills the tree and is recorded.
- A failed gate stops before any case.
- There is no approximate fallback. A unit that fails stays failed and visible.

## Amendment 1 (after `cases` exited 1; the truth units had succeeded)

**Failure.** `cases` raised `TypeError: dict() got multiple values for keyword argument 'case'` before writing any row. The
provenance file written by the truth units has its own `case` key, and my harness also passed `case=`.

**Fix.** The column is renamed to `evaluated_case`. This is a harness defect: no metric, case, tie rule or scope changes.

**Rerun:** `cases.a1`, on the same truth units. `cases` stays recorded.

## Amendment 2 (after `cases.a1`; no metric or rule changed)

**Observation.** All three CYP2D6 interior comparisons are UNSUPPORTED. `project_path_onto_locus` returns no frame for
**sv:17 of truth haplotype 2** (`NA19650#2#haplotype2-0000078`): that path's walk is shorter than the bytes the divergence
trim moves out of the block. The bytes are **not** reconstructed. The interior rows stay UNSUPPORTED, with the block named.

sv:9 truth is available (39,609 + 27,459 bp), unlike M1-B's single-block run.

**Addition.** For a case whose interior is unsupported, one further row reports the **largest framed prefix**: blocks
1..k, where blocks 1..k+1 are projected for both truth haplotypes, so the moved bytes at the end boundary are framed.
- The row is labelled as that region, **not** the chain interior, and the contract is otherwise unchanged.
- For CYP2D6 it is sv:1..sv:15.
- It is a narrower named scope, not an approximation of the interior.

**Rerun:** `cases.a2`. `cases.a1` stays recorded.
