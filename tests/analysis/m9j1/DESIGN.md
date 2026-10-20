# M9-J.1: decision-set resolution correction, validation closure, development checkpoint

**Results:** `tests/results/m9j1_closure_20261020/`. Branch `m9-j1-resolution-fix` from `m9-j-compound-integration` (`3481bca`).

## 1. The defect (Sol's review of M9-J)

`compound_assess` resolved the target from the **projections it could collect** from decision-set groups. A decision-set
group with no recorded labelling contributed nothing:
- one known decision-set projection plus one missing decision-set group gave `RESOLVED` / `the_resolved_target`;
- a decision set whose only group was missing gave `NO_DECISION_SET`, yet still `the_resolved_target`.

Reproducer: `/private/tmp/panvar-m9j-review.lm2ynx/repro.cpp` (copied verbatim into the `sol_repro` unit).

## 2. The correction (frozen)

Resolution is decided **independently of the support / error calculation**, over the decision-set groups only:

| target_resolution | condition | mass_reference |
|---|---|---|
| `NO_DECISION_SET` | no group flagged as a decision-set member | `representative_projection_only_no_decision_set` |
| `INCOMPLETE_DECISION_SET` | at least one decision-set group has no recorded labelling | `representative_projection_only_decision_set_incompletely_projected` |
| `UNRESOLVED_DECISION_SET` | every decision-set group projected; projections differ | `representative_projection_only_not_a_resolved_target` (unchanged) |
| `RESOLVED` | every decision-set group projected; one projection | `the_resolved_target` (unchanged) |

- An incompletely projected decision set is never `RESOLVED`, even when its known projections differ (an unknown projection
  is not evidence of disagreement either, so it is reported as `INCOMPLETE`, not `UNRESOLVED`).
- Missing **non-decision** groups do not affect resolution. They still enter both upper bounds, as before.
- Bounds, status, competitor, counts and the prediction are unchanged.
- **Schema (declared):** two columns are added to `sequence_compound.tsv`, after `target_resolution`:
  `decision_set_groups` and `decision_set_missing_groups`. The replay harness output gains the same two columns.

## 3. Fixtures (added to `compound_annotation_test`, 21 → 26 checks)

| fixture | expectation |
|---|---|
| R1 known + missing decision-set group (Sol case 1) | `INCOMPLETE_DECISION_SET`, representative-only wording, `decision_set_missing_groups` 1, support [0.6, 1.0] unchanged |
| R2 entirely missing decision set (Sol case 2) | `INCOMPLETE_DECISION_SET`, never `the_resolved_target` |
| R3 no decision set flagged | `NO_DECISION_SET`, not `the_resolved_target` |
| R4 complete decision set + missing non-decision group | `RESOLVED`, `the_resolved_target`, status `ASSESSED_WITH_MISSING_PROJECTIONS` |
| R5 differing projections + a missing decision-set group | `INCOMPLETE_DECISION_SET` |

## 4. Validation units (no genotyping)

| unit | check |
|---|---|
| `build` | every touched target |
| `gate.compound_unit` | 26/26 |
| `gate.junction_unit`, `gate.sequence_report`, `gate.e2e_sequence_first` | unchanged behaviour, publication and failure handling |
| `sol_repro` | the reviewer's file compiled against the fixed library: both cases `INCOMPLETE_DECISION_SET`, no `the_resolved_target` |
| `replay_compare` | the unchanged M9-J 34-case replay against the M9-I.1 Python reference (same declared tolerances: 34/34), **and** every column of the saved M9-J replay (`gate.replay.a2/cpp.tsv`) byte-identical, with `decision_set_missing_groups` 0 and `decision_set_groups` > 0 in all 34 |
| `saved_smoke` | the four saved M9-J caller annotation files (ACOT 4 and 1 threads, CYP2D6, C4 control) re-assessed by the fixed helper from each run's **own saved tables** (manifests saved by `smoke_check.a4`): resolution equal to the saved one, decision set complete, `the_resolved_target` iff `RESOLVED`; the C4 control still has no rows |

**Expected:** no saved annotation changes (every saved case has complete decision-set projections). If one changes, it is
reported as a correction of that historical record, appended, never edited. **No caller rerun** unless a reproducible
regression requires one (it would be a declared amendment).

## 5. Budget

**600 s** including builds and gates; 25 s of disclosed pre-freeze development (one compile of the unit test, one run).
Admission: spent + declared timeout of the next unit. Failed units stay recorded.

## 6. Then (documents only)

- `CHECKPOINT.md`: the development checkpoint (revision, effective experimental settings, declaration artifacts, required
  inputs and indexes, limitations). No merge, push or default promotion.
- `QUALITY_TARGET.md`: the production quality target proposal, specification only; no threshold, no calibration.
