# M10-D result: COSIGT-style sequence-similarity evaluation, implemented on the existing exact edit distance; contract verified; seven development cases evaluated (CYP2D6 on a named framed prefix)

**Run.** Branch `m10-d-sequence-similarity`.
- **Commits:** freeze `40bc3d0`; Amendment 1 `0ff47ef` (harness column name); Amendment 2 `eba1e76` (a named
  framed-prefix scope); result after.
- **Results:** `tests/results/m10d_sequence_similarity_20261020/`, with `runs/cases.a2/{cases,provenance}.tsv` and
  `runs/gate.contract/`.
- **Budget:** 44.4 of 900 s, including 10 s of disclosed development.
- **Watchdog:** no unit hit its timeout or the 6 GiB envelope. The peak process-tree RSS was 1.47 GiB, in the ACOT truth
  extraction.
- **Failed unit, recorded:** `cases`, my harness bug.
- **Closed:** the structural-alphabet and event-comparator investigation. Davide's scope clarification is appended to the
  M10-B, B.1 and C results. M10-A remains the computational checkpoint.

## 1. Contract and implementation

**Contract:**
- `D = min(d(A,C)+d(B,D), d(A,D)+d(B,C))`, with one whole-region assignment chosen by summed edits. A tie selects
  assignment 1 and is recorded.
- `error_rate = D/(L1+L2)` (unfloored), `identity = 1 − error_rate`, `QV = −10 log10(max(0.5,D)/(L1+L2))`, and 60 for an
  empty alignment.
- This is the benchmark's `qv_from` and its length-weighted combination. L1 and L2 are the alignment columns of edlib's
  deterministic path; co-optimal paths can differ in column count.

**Reuse:**
- the distance is `panvar::nw_edit_distance` (`m10c_align_eval --distance`, one call per pair);
- truth comes from the existing `project_path_onto_locus` (`block-sequences`), which needs no panel spellability.

No new aligner.

**Gate: 21/21.**
- identical DNA gives D = 0, error exactly 0 and the floored QV;
- swapped homologues give the same D under the other assignment;
- a whole-sequence reverse complement gives D 2,909 against 9, and the orientation check prefers the given strand;
- small differences: 0 < D ≤ 9;
- unequal lengths give exactly the length-difference bound (55, 600, 60), and 1,200 for 3+1 vs 2+2 at equal total length;
- order and inversion: exact distances reported (637; 315 and 35);
- tie handling; empty handling (all-empty is EMPTY with QV 60);
- alignment columns within [max(len), len_a+len_b], so there is no clipping;
- **the textbook DP oracle agrees on 400 tiny pairs** (inversions, order changes and empty sequences included) and on 60
  tiny diploids against the brute-force assignment minimum.

**Orientation.** Truth and predictions are both in the model's canonical chain orientation by construction. The
whole-sequence 31-mer check returns AS_GIVEN for every evaluated haplotype (for example C4 123,491 shared k-mers as given
against 910 in the reverse complement). Nothing was flipped.

## 2. Development cases (inspected, not held-out accuracy or calibration)

**Scope.**
- **Evaluated:** the **chain interior**, the region with matched endpoints. C4 is blocks 1..9 of 11; CYP2D6 and ACOT are
  blocks 1..17 of 19.
- **Not scored:** each flank-inclusive comparison is listed as UNSUPPORTED, because the flank blocks hold each path's own
  assembly extent.

| case | prediction | scope | D1 | D2 | assignment | h1 / h2 edits | L1 + L2 | **D** | error rate | identity | QV | earlier per-block metric (different unit; not subtracted) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| C4 NA19650 | reported 3+1 | interior | 65,617 | 65,588 | 2 | 32,747 / 32,841 | 189,660 + 156,925 | **65,588** | 1.892e-01 | 0.8108 | 7.23 | M1-B sv:7 excess 65,444 |
| C4 NA19650 | saved 2+2 | interior | 244 | 168 | 2 | 11 / 157 | 156,922 + 156,927 | **168** | 5.353e-04 | 0.99947 | 32.71 | M1-B sv:7 excess 86 |
| C4 HG00096 control | call | interior | 0 | 66,280 | 1 | 0 / 0 | 150,343 + 117,856 | **0** | **0** (exact) | 1 | 57.29 (floored) | — |
| CYP2D6 NA19650 | 5× 2+2 | interior | — | — | — | — | — | **UNSUPPORTED** | truth2 sv:17 not projectable | | | M1-B sv:9 excess 12,159 |
| | | framed prefix sv:1..sv:15 | 12,275 | 13,179 | 1 | 12,174 / 101 | 200,905 + 188,936 | **12,275** | 3.149e-02 | 0.9685 | 15.02 | |
| CYP2D6 NA19650 | 30× 3+2 (also the synonymous case) | interior | — | — | — | — | — | **UNSUPPORTED** | (same) | | | M1-B sv:9 excess 12; M9-I compound target exact |
| | | framed prefix | 25,218 | 120 | 2 | 95 / 25 | 188,929 + 200,911 | **120** | 3.078e-04 | 0.99969 | 35.12 | |
| CYP2D6 NA19650 | saved 2+3 | interior | — | — | — | — | — | **UNSUPPORTED** | (same) | | | M1-B sv:9 excess 281 |
| | | framed prefix | 133 | 25,220 | 1 | 27 / 106 | 200,911 + 188,939 | **133** | 3.412e-04 | 0.99966 | 34.67 | |
| ACOT L4 HG02083 | JP on | interior | 637 | 1,011 | 1 | 556 / 81 | 162,188 + 161,687 | **637** | 1.967e-03 | 0.99803 | 27.06 | M9-H interval distance 0 |
| ACOT L4 HG02083 | JP off | interior | 1,025 | 651 | 2 | 82 / 569 | 161,687 + 162,201 | **651** | 2.010e-03 | 0.99799 | 26.97 | M9-H interval distance 14 |
| ACOT L2 HG02074 | JP on | interior | 485 | 34 | 2 | 21 / 13 | 161,656 + 161,697 | **34** | 1.051e-04 | 0.99989 | 39.78 | M9-H interval distance 3 |
| ACOT L2 HG02074 | JP off | interior | 42 | 478 | 1 | 18 / 24 | 161,696 + 161,654 | **42** | 1.299e-04 | 0.99987 | 38.86 | M9-H interval distance 5 |

**Unsupported, with reasons (kept in the table):**
- CYP2D6 interior ×3: `project_path_onto_locus` gives no frame for **sv:17 of truth haplotype 2**, whose walk is shorter
  than the bytes the trim moves out. The bytes are not reconstructed. The CYP2D6 rows above are therefore on the named
  framed prefix, which is **not** the interior.
- The flank-inclusive scope ×7: unmatched path extents.

**Provenance.** The truth GFAs are SHA-256-recorded (`provenance.tsv`). Every truth extraction's model fingerprint equals
the M10-B.1 dump's.

## 3. What changes under the sequence-similarity interpretation

- **Labelling no longer matters.** The CYP2D6 30× call, whose block labels were wrong by 118 + 118 bp but whose compound
  target was exact, is scored only by sequence: D = 120 on the framed prefix. ACOT's 3 + 3 bp label artefacts likewise
  contribute nothing.
- **Allocation is captured without classifying events.** The whole-haplotype assignment makes C4's 3+1 call cost
  65,588 edits, about two RCCX units, against 168 for the saved 2+2, at equal total dosage. CYP2D6's missing unit at 5×
  costs 12,275, against 120 for the 3+2 twin. No unit or event detector is involved.
- **Small differences count in proportion.** ACOT's residuals are a few hundred edits over about 323 kb of interior
  (identity ≈ 0.998–0.9999). JP1.1 lowers D by 14 at L4 (651 → 637, matching M9-H's interval change of 14 → 0) and by 8
  at L2 (42 → 34). These are **whole-interior** numbers: they include differences outside sv:5..ctx:6.
- **Units differ.** The earlier per-block availability-floor metrics are another unit. They are shown beside D and never
  subtracted. No floor search was run.
- **What the metric does not distinguish.** A large structural difference and many scattered small differences of equal
  total edits score the same. That is the accepted consequence of the clarified notion.

**Not defined or changed:** no categorical "structurally correct" cutoff; no clusters; no posterior pooling; no change to
GQ, confidence, filters, predictions, priors or caller code. No reads, families, parser, aligner, merge, push, cleanup or
promotion. WIP, stash, indexes and earlier packages are preserved.
