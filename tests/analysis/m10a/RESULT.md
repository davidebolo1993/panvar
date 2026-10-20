# M10-A result: `genotype-dev-m10a-v1` is computationally validated on the six development cells, within the stated gaps; JP1.1 changed only one block's conditional GQ in four cells

**Run.** Branch `m10-a-consolidation`.
- **Commits:** plan `585f99a`; Amendment 1 (preflight) `69f014d`; Amendment 2 `c9cc5fc`; Amendment 3 text `50804fc`, code
  `4af69fc`; Amendment 4 `e853ce0`.
- **Results:** `tests/results/m10a_dev_configuration_20261020/`, containing `ledger.tsv`, `runs/summary/{runs,comparisons,index}.tsv`
  and the index files under `index/`.

**Budget: 1,150.5 of 2,400 s.** That is 1,100.5 s of ledger wall time plus 50 s of disclosed development (Amendment 1: 30
s; Amendments 2–4: 20 s). No caps were relaxed, no defaults changed, no truth was used and no donors were added.

**Scope.** Each locus is **one donor at two coverages** (C4 HG00096, CYP2D6 HG02769, ACOT HG00097). This is not
validation across biological samples.

## Recorded failures and amendments

Every failure was in comparator classification. **No scientific file ever differed in a computational comparison.**

| amendment | trigger | cause | fix |
|---|---|---|---|
| 1 (prospective preflight) | Sol's reproducer: `jp` overrode a settings failure with JP_ATTRIBUTABLE | wrapper defect | J1–J5 invariants that eligibility never waives; 18 controls, including a real JP off/on pair and the named negatives; JP failures stop caller runs; comparison reserve before caller admission |
| 2 | `cmp.route.c4_5` and `cmp.jp.c4_5` failed; all four arms' scientific files byte-identical | unclassified status keys: `compound_seconds` (timer), footprint-checkpoint measurements in `jp`, the status echo of the JP setting | keys classified (checkpoints also added to `prepared`, prospectively); the stop had recorded the CYP2D6 5× units NOT RUN, so they ran as `.a2` |
| 3 | `cmp.jp.acot_5` failed, only on `time_limit_seconds` 130 vs 160 | **my plan** gave D larger timeouts than D0 | the time guard is a declared operational difference in `jp` (both runs must complete; values reported). Commit `50804fc` holds only the text; the code is `4af69fc` (no unit ran in between) |
| 4 | `cmp.prepared.c4_5` / `c4_30` failed, only on `junction_read_table_seconds` (0.0211 vs 0.0208 s) | unclassified JP timer | both JP timers classified as operational; `write.cyp` (recorded NOT RUN) ran as `.a4` |

Controls were rerun after each comparator change: 18/18 every time. Each negative failed with exactly its intended reason.
Sol's reproducer also lacks the JP setting in B, so it fails on J1 in addition to J3.

## Complete matrix

**Key:** P = PASS (identity); J = JP_ATTRIBUTABLE (invariants J1–J5 held); F→P = failed, then passed after the amendment
shown; NR = NOT RUN.

| cell | REF / O / D0 / D runs | backend REF–O | route O–D0 | JP D0–D | cold–loaded D–P | repeat |
|---|---|---|---|---|---|---|
| C4 5× | all COMPLETE | P | F→P (A2) | F→P (A2): eligible 0, identity | F→P (A4) | — |
| CYP2D6 5× | first attempt NR (stop); `.a2` all COMPLETE | P | P | J: eligible 1 | P | — |
| ACOT 5× | all COMPLETE | P | P | F→P (A3) J: eligible 1 | P | — |
| C4 30× | O / D0 / D COMPLETE | **not planned** | P | P: eligible 0, identity | F→P (A4) | — |
| CYP2D6 30× | O / D0 / D: completed sv:3/sv:5 tie (`RESCORE_INCOMPLETE` / distinct tie, comparable) | **not planned** | P | J: eligible 1 | P | — |
| ACOT 30× | D0 / D COMPLETE | **gap** | **gap** | J: eligible 1 | P | P (D vs D.r) |

**NOT RUN entries, kept in the ledger:**
- `REF/O/D0/D.cyp_5`, `cmp.{backend,route,jp}.cyp_5` and `write.cyp`: recorded by the stop rule, superseded by `.a2` /
  `.a4`;
- `loaded.acot.single30` and its comparison: the condition was not met, because the sequential batch delivered a completed
  ACOT 30× sample.

## Computational equivalence

- **Backend** (reference against optimised backends; P1 chain plus every sequence and compound file): **identical** on C4,
  CYP2D6 and ACOT at 5×.
- **Retained route** (ordinary against blocked union + round-1 admission + early incumbent): **identical** on all three
  5× cells and on C4 and CYP2D6 30×.
- **Cold against persisted index:** **identical** on all six cells, including the compound annotation file.
- **Determinism:** ACOT 30× D repeated gives **identical** results.
- **The annotation** (`sequence_compound.tsv`) was byte-identical in every computational pair.
- **Gap kept explicit: ACOT 30× has no reference or ordinary-route comparison.** The ordinary route refuses at 500 M
  physical origins under 6 GiB (S2-L), and REF at 30× was not planned. Its evidence remains S2-G (earlier source,
  validation-only cap extensions) plus this campaign's cold–loaded identity and repeat. Backend equivalence was tested at 5×
  only.

## JP1.1-caused differences (no accuracy claim; no truth used)

| cell | eligible blocks | pool change | common hypotheses (objective exactly equal) | prediction, representative haplotypes, sequence status | compound target / bounds | changed outputs |
|---|---|---|---|---|---|---|
| C4 5×, 30× | 0 | none | — | identical | identical | none: full identity |
| CYP2D6 5× | 1 | 8,781 → 8,780 sequence candidates | 8,969 | unchanged | sv:3..sv:5 target unchanged; error upper 5.09e-7 (last digits only) | sv:11 `conditional_GQ` 162.4 → 125.8 |
| CYP2D6 30× | 1 | 8,767 → 8,766 | 8,964 | unchanged (completed tie) | unchanged | equivalence set, candidates and scores only |
| ACOT 5× | 1 | 109,887 → 109,886 | 109,824 | unchanged | sv:5..ctx:6 target and bounds unchanged | sv:5 `conditional_GQ` 503.5 → 275.9; dosage blocks |
| ACOT 30× | 1 | 109,911 → 109,912 | 109,848 | unchanged | unchanged | sv:5 `conditional_GQ` 1,858.8 → 1,444.0; dosage blocks |

JP1.1 cost: one invocation per run, 0.007–0.015 s.

## Persisted index: three separate questions

1. **Cold against indexed scientific equivalence:** identical for all six samples.
2. **Sequential reuse (5×, then 30×, one process, release between samples, 6 GiB ceiling):** **completed at all three
   loci**: 2/2 samples each, `RESOURCE_COMPLIANT`.
3. **Fresh-process fallback:** not needed and not run. Nothing here rests on a fallback.

| locus | index size | write (prepare + encode) | verify (includes a full load) | load in batch | 5× sample | 30× sample | current after 5× → after 30× | process lifetime peak | cold peaks 5× / 30× |
|---|---|---|---|---|---|---|---|---|---|
| C4 | 45.1 MB | 2.0 s (1.80 + 0.19) | 0.87 s | 0.47 s | 7.2 s | 26.3 s | 0.40 → 0.31 GB | 1.33 GB | 0.81 / 1.69 GB |
| CYP2D6 | 39.8 MB | 2.2 s (2.03 + 0.17) | 0.78 s | 0.41 s | 7.8 s | 28.1 s | 0.45 → 0.35 GB | 1.17 GB | 0.83 / 1.60 GB |
| ACOT | 183.0 MB | 7.3 s (6.43 + 0.79) | 3.61 s | 1.92 s | 27.4 s | 103.1 s | 1.17 → 1.92 GB | 4.31 GB (30× sample) | 2.82 / 4.40 GB |

- **Memory carried between samples:** ACOT holds 1.17 GB after the 5× sample. The 30× sample's peak (4.31 GB) is still
  below its cold run's (4.40 GB) and within the 6 GiB ceiling.
- **Timings are descriptive only.** Runs were sequential, but machine load was not controlled.

## Is the unified development configuration validated?

**Computationally, yes, on these six cells, within the stated gaps.** Specifically, all of the following hold for
`genotype-dev-m10a-v1` (performance backends, retained route, JP1.1, sequence_first_v2, compound annotation, per-model
declarations, per-cell resource layer):
- it reproduces the reference backends (5×) and the ordinary route (5×; C4/CYP2D6 30×) exactly;
- it runs identically from a persisted index, including sequential 5×→30× reuse in one process;
- it is deterministic at ACOT 30×;
- JP1.1 alters only what candidate generation legitimately can, under invariants that held exactly.

**Not validated:**
- ACOT 30× against the reference or ordinary route, and backend equivalence at 30×;
- other donors or loci (each locus is one donor at two coverages);
- concurrent in-process callers;
- accuracy or calibration.

No default was changed.
