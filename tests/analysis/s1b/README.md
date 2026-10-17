# S1-B: opt-in sequence-first reporting (`output_schema sequence_first_v1`)

**What it is.** The implemented reporting change from S1-A (contract section 7 of `tests/genotype_sequence_label_reporting_design.md`, and `tests/analysis/s1a/S1_PROPOSAL.md`), with the corrections listed below. It is built on the validated performance baseline: branch `s1-b-sequence-first`, from `p1-consolidated-baseline` `49e6dc3`.

**Code:**
- `src/genotype_sequence_report.cpp` and `include/panvar/genotype_sequence_report.hpp` hold the reporting layer;
- `src/genotype_command.cpp` holds the integration.

**Untouched:** scoring, priors, candidate generation, adaptive search, caps, the performance backends, every legacy output, and dosage and phase decisions.

## 1. Switch and provenance

| setting | effect |
|---|---|
| unset (default) | nothing changes: outputs, status, settings and exit codes are identical to before |
| `output_schema legacy`, given explicitly (profile or `--output-schema legacy`) | exactly two rows: settings `output_schema legacy <source>` and status `output_schema legacy`. Nothing else |
| `output_schema sequence_first_v1` | the same two rows (value `sequence_first_v1`), plus four status rows (`sequence_prediction`, `sequence_status`, `sequence_completion`, `sequence_exit_code`) and five new files (section 2). The dosage stage copies the caller status to `call.original.mosaic_status.tsv`, so that copy carries the same rows |

**Refused under `sequence_first_v1`:**
- runs without the adaptive repair (`repair_rounds 0` / `--adaptive-repair off`), and `diagnostic_union_only`, with an explicit error;
- a legacy call on a route the report does not cover. It would be written as `UNSUPPORTED_ROUTE`, never a silent legacy fallback.

**Computation cost.** Nothing extra is computed unless the schema is `sequence_first_v1`.

## 2. Outputs (`sequence_first_v1` only)

| file | content |
|---|---|
| `.sequence_locus.tsv` | Key/value. `prediction` (PREDICTION / PROVISIONAL / NO_PREDICTION); `sequence_status`; `provisional`; `completion` and its note; `decision_scope`; the meaning of RESOLVED; the representative group, its SHA-256 and lengths, one explored labelling, the tie-break and its scope; decision-set, pool and labelling counts; block-status counts; `representation_intervals`; representative and decision-set probabilities (candidate-conditional, uncalibrated); `structure_confidence NA` with reason `UNCALIBRATED`; `projection_coverage EXPLORED` and its scope; legacy GQ scope; `legacy_phase_support`; `phase_assertion`; `dosage_scope`; `legacy_run_status`, `legacy_call_class` and `legacy_no_call_reason` (diagnostics); `exit_code` |
| `.sequence_blocks.tsv` | `block_status` (IDENTIFIED / REPRESENTATION_AMBIGUOUS / GENOTYPE_UNRESOLVED / PROVISIONAL_LABEL); `genotype` (identified or provisional only); `observed_genotypes` with the decision-set groups carrying each; `error_mass_lower` and `error_mass_upper` (identified only); `projection_coverage`; `representation_interval`; `legacy_representative_GQ` and `legacy_numerically_separated` (legacy diagnostics) |
| `.sequence_haplotypes.fa` | The representative sequence pair, ordered by SHA-256. Each header carries `prediction`, `sequence_status`, `representative_of_set`, `provisional`, the intervals, `phase=not_asserted` and `structure_confidence=NA`. Written for PREDICTION and PROVISIONAL; removed for NO_PREDICTION |
| `.representation_intervals.tsv` | Per decision-set group: interval, first and last block, the blocks, `conservative_enclosing_not_minimal`, and the observed tuples (`hapA/hapB` labellings over the interval) |
| `.sequence_alternatives.tsv` | The decision set, then the next five groups by rank. Score components to 12 significant digits, probability, explored labellings and projection ambiguity, internal digests and lengths (no sequence copied), and the blocks whose explored genotypes differ from the representative's |

**Exit codes under `sequence_first_v1`:**

| code | meaning |
|---|---|
| 0 | PREDICTION: a completed ranking, even when uncertain |
| 3 | PROVISIONAL: an incomplete search |
| 1 | NO_PREDICTION |
| 2 | errors, unchanged |

Under legacy the codes are unchanged: 0 for a call, 1 for a no-call.

## 3. Semantics (corrected contract)

**Representative.** The canonical member of the decision set (`PooledDecision::top_groups`):
- with one group, that group (`highest_score`);
- with several, every group is numerically tied with the best under the caller's own certified bounds, so the tie is broken by the sorted pair of SHA-256 digests of the spelled sequences (`canonical_sha256`).

Bitwise-equal scores are counted as `exact_score_ties` for information only, because sample order changes a prior in its last bits (amendment 1). The representative never depends on discovery or input order.

**Sequence status:**
- `RESOLVED`: one group in the decision set. This is resolution within the declared decision set, not confidence and not correctness.
- `GENOTYPE_RESOLVED_SEQUENCE_SET`: several groups, every block identified.
- `SEQUENCE_TIED`: several groups that disagree somewhere. Identified blocks stay identified.

**Block status, over EXPLORED labellings:**
- `IDENTIFIED`: one unordered genotype across every explored labelling of every decision-set group.
- `REPRESENTATION_AMBIGUOUS`: one group whose labellings disagree here.
- `GENOTYPE_UNRESOLVED`: several groups that disagree here.

**Error-mass bounds (identified blocks only).** Each sequence group counts once, with its candidate-set probability:
- lower = groups whose every explored labelling carries another genotype;
- upper = groups with any explored labelling that carries another genotype.

They are conditional on the explored pool and labellings. The old representative-based GQ is kept as `legacy_representative_GQ`.

**Representation intervals.** For each decision-set group, the spans between the first and last differing block of every pair of explored labellings, merged when they overlap. Coupled alternatives therefore include the unchanged blocks between them. The interval is a conservative enclosing range, not claimed minimal, and its tuples are whole observed labellings, never per-block Cartesian combinations.

**PROVISIONAL.** A `round_cap` or `cumulative_cap` stop after at least one committed round, reported from the last committed, exactly scored incumbent:
- `sequence_status` is NOT_DETERMINED and `decision_scope` is NOT_DETERMINED;
- blocks are PROVISIONAL_LABEL only;
- there are no bounds, intervals or alternatives;
- exit code 3.

**NO_PREDICTION** covers:
- watchdog stops, which write no sequence file at all;
- refusals before a ranking;
- a `neighbourhood_bound` refusal, which happens before any search;
- improving-move ties and tie-branching failures, which are not declared caps.

**Downstream:**
- no PASS threshold and no calibrated probability;
- the representative's labels never drive dosage repair or phase;
- strict-call dosage and phase behaviour is unchanged;
- the persisted prepared-locus index contract is unchanged, so its version is unchanged; F2 and the e2e prepared route show identical sequence outputs.

## 4. Deviations from `S1_PROPOSAL.md`, declared

- **Annotations go in new `.sequence_*` files, not new columns in `locus_genotype.tsv` and `block_genotypes.tsv`.** The dosage stage rewrites those tables (repair mode) and downstream scripts parse them, so leaving every legacy file byte-identical in both schemas is the safer contract. The legacy classifications are repeated in `.sequence_locus.tsv` as diagnostics.
- **Intervals are coupled enclosing spans, not maximal runs of adjacent ambiguous blocks** (the P1.1 correction).
- **The representative tie-break uses SHA-256 of the spelled sequences** over the numerically tied decision set (amendment 1).
- **The caller-level tiny fixture for coupling across an unchanged context was not constructible.** The block builder merges every small periodic layout tried into one block: four layouts, including context SNV alleles. The regression is therefore at the reporting level (`genotype_sequence_report`, `coupled`). The CYP2D6 30× sentinel shows it at real scale: `sv:3..sv:5` encloses an unchanged `ctx:4`.

## 5. Reproducing the validation

`python3 tests/analysis/s1b/s1b_run.py next` runs the frozen manifest plus amendments 1 and 2.

**Results:** `tests/results/s1b_sequence_first_20261018/README.md`.

**Fixture gates** (`s1b_fixtures.py`) compare against `s1b_oracle.py`, whose projections are complete, using exhaustive enumeration on tiny fixtures only.

## 6. Remaining limitations

- **EXPLORED coverage.** Uniqueness and bounds are conditional on the explored labellings (fixture F9). Complete projection sets from model-level intervals are later work.
- **Intervals** are conservative and not minimal.
- **No calibration.** There is no PASS threshold and no structural-tolerance model.
- **The unrepaired view (`repair_rounds 0`) and union-only diagnostics** are refused, not reported. Panel-only runs with the adaptive repair are supported (F7, F9).
- **Neighbourhood-bound and improving-move-tie stops** give no prediction. Neither is a declared-cap stop reached after a committed round.
