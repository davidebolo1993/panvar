# M9-J.1 result: decision-set resolution corrected (no saved annotation changes), validation closed in 126 of 600 s, development checkpoint frozen, production quality target proposed

**Run.** Branch `m9-j1-resolution-fix`.
- **Commits:** fix and freeze `e61c09c`; Amendment 1 `59c9919`.
- **Results:** `tests/results/m9j1_closure_20261020/`.

**Budget.** 126.0 of 600 s, including 25 s of disclosed pre-freeze development. No genotyping run.

## 1. The correction

`compound_assess` now decides resolution **before and independently of** the bounds, over decision-set groups only:

| condition | `target_resolution` | `mass_reference` |
|---|---|---|
| no group flagged | `NO_DECISION_SET` | `…representative_projection_only_no_decision_set` |
| any decision-set group without a recorded labelling | `INCOMPLETE_DECISION_SET` (new) | `…representative_projection_only_decision_set_incompletely_projected` |
| complete, projections differ | `UNRESOLVED_DECISION_SET` | unchanged |
| complete, one projection | `RESOLVED` | `the_resolved_target` |

- Missing **non-decision** groups do not affect resolution (fixture R4). They still widen both upper bounds.
- Bounds, status, competitor, counts and the prediction are unchanged.
- **Declared schema change:** `sequence_compound.tsv` gains `decision_set_groups` and `decision_set_missing_groups` after
  `target_resolution`.

## 2. Validation

| unit | result |
|---|---|
| build | exit 0, 10.5 s |
| `compound_annotation_test` | **26/26** (21 earlier + R1–R5) |
| JP1.1, sequence-report, e2e `sequence_first` | pass (0.7, 0.8, 42.9 s) |
| `sol_repro` | **failed at link time** (my harness linked too few libraries); kept recorded |
| `sol_repro.a1` (Amendment 1: the test target's link line) | both of Sol's cases → `ASSESSED_WITH_MISSING_PROJECTIONS / INCOMPLETE_DECISION_SET / representative_projection_only_decision_set_incompletely_projected / missing=1`; no `the_resolved_target` |
| `replay_compare` | M9-J replay against Python **34/34** (all `RESOLVED`); **every saved column byte-identical** to `m9j…/gate.replay.a2/cpp.tsv`; every case has 0 missing decision-set groups |
| `saved_smoke` | ACOT 4-thread, ACOT 1-thread and CYP2D6 saved caller annotations re-assessed from their own saved tables: `RESOLVED` as saved, 2 decision-set groups, 0 missing. The C4 control still has no rows |

**No saved annotation changes.** The M9-J records stand, and a correction note is appended to M9-J's RESULT. The defect was
unreachable in every saved run, because every decision-set group had recorded labellings.

## 3. Development checkpoint

See `CHECKPOINT.md` for the source revision, the effective settings per validated configuration, declaration artifacts
with their sha256, required inputs, and limitations. It records:
- **Concurrency:** concurrent caller invocations in one process are outside the validated contract. Helper concurrency and
  within-sample threads are not multi-sample validation.
- **Timings:** M9-J timings overlapped storage compaction, so the off/on wall-time differences are not an annotation speed
  effect.
- **No combined profile:** none combining JP1.1, the M1-B overlay and the annotation has been validated.

Not merged, pushed or promoted.

## 4. Production quality target

See `QUALITY_TARGET.md`. It is a specification only.
- **Target:** the structural diplotype per structural interval. It treats synonymous labellings identically and tolerates
  declared small differences. It preserves copy number, homologue allocation (C4 3+1 ≠ 2+2), order, orientation and
  junction class.
- **Phase:** reported and decided separately.
- **Interface:** one prediction (always kept), one `structure_quality` decision with one reason code, and one phase row.
- **Evaluation:** representability against the active panel, with confident calls on unavailable content counted as errors.
- **Approval:** eight choices need Davide's approval; none is set here.
