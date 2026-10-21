# M15-D.1 result: declared-library preflight and termination statuses corrected; CYP2D6 audit limitation reported honestly; every unit and the full suite pass

**Run.**
- **Branch and commits:** `m15-d1-correctness`, from `6cbd717`. Freeze `7f2def1`; result after.
- **Results:** `tests/results/m15d1_correctness_20261021/`.
- **Budget:** 1,405.0 of 1,800 s, including 240 s of disclosed development.
- **Units:** 10 of 10 passed first time, with no amendment and no watchdog or envelope event.
- **Resources:** 4 threads, one process at a time, peak 1.5 GiB (the build).
- **Not done:** no merge, push, tag, cleanup, model tuning or reserved rerun.

**Final source and binaries:**

| item | value |
|---|---|
| source revision | `7f2def16af7ed36ff2c9d12b51fe21cd989b767c` |
| `panvar` SHA-256 | `457b03cc7d26d38b4451f4f7edff7ed23989da70658b4bb8ae4dd6d44204b3ba` |
| `genotype_command_e2e_test` SHA-256 | `1279bfe360a93773ae104b5c7f628b5d342cb9ce8a434c41afcb8fc04e92bbac` |

## 1. Declared-library preflight (M15-D F2)

- **The fix:** a `--library FILE` is read and validated in full before the run clears its prefix. The check sits at the
  start of the run body, so the CLI, the direct API and the prepared route all get it.
- **e2e gate** (`test_declared_library_preflight`, CLI and API):
  - seven cases are refused with exit 2 and the problem named: missing, unreadable, repeated key, unknown key, missing
    key, negative sd and non-finite value;
  - all six sentinel files at the prefix stay byte-identical, and nothing is added.
- **Binary** (`ex.library_preflight`, C4, five sentinel files): missing, unreadable, repeated-key and unknown-key files are
  refused (8 of 8 checks), with the prefix byte-identical.
- **The documented guarantee** (guide §6, module page): only argument, profile and declared-library errors are refused
  before the prefix is touched. Errors found while running (reads, graph or model refusals, resource stops, failed writes)
  happen after the run owns the prefix, and are **not** transactional.

## 2. Termination statuses (M15-D F3)

The process exit, `.sequence_locus.tsv` `exit_code`, the status's `sequence_exit_code` and (for stops and failures) the
status's new `exit_code` row now agree. The fixtures are deterministic: test hooks and write blockers, not wall-clock
races.

| outcome | exit | e2e (tiny locus) | binary (C4) | before M15-D.1 |
|---|---|---|---|---|
| time stop before a checkpoint | 124 | pass | pass: `TIME_LIMIT`, NO_PREDICTION, no sequences (4/4) | record and status said 1 |
| footprint stop before a checkpoint (new hook) | 125 | pass | pass: `RESOURCE_REFUSAL (footprint_limit)`, NO_PREDICTION (4/4) | record and status said 1; no deterministic fixture existed |
| stop after a committed checkpoint | 3 | pass | pass: PROVISIONAL, `budget_stop`, sequences kept (4/4) | the kept record had no `exit_code` |
| ordinary no prediction | 1 | pass | pass: `neighbourhood_bound`, early incumbent off (3/3) | unchanged |
| completed prediction | 0 | pass | — | unchanged |
| checkpoint publication failure (existing contract) | 2 | pass | — | unchanged |
| locus table cannot be written after the sequence report was published | 2 | pass | pass: every exit field 2, `locus_output FAILED: cannot write …` in report and status (3/3) | record and status said 0 |

**Also:**
- **A v2 run without the early incumbent** now writes its NO_PREDICTION sequence record and status rows on a stop or
  failure. This route has no deterministic test hook (the hooks need the early incumbent), so it is covered by code
  reading only.
- **The existing S2-A stop assertion** (`sequence_exit_code == "1"` for an exit-124 stop) encoded the inconsistency, and
  is corrected to 124.

## 3. CYP2D6 audit limitation (M15-D F1, not fixed)

`ex.cyp_audit` (HG02769, declared library, `--dosage-check`, `--quiet`), 7 of 7 checks:
- **The prediction is unchanged:** identical to the archived M11-D run and to M15-D's run without the audit.
- **The audit reports FAILED with its concrete reason:** "unit library: alleles disagree on their flanks or anchor trims
  in sv:9".
- **Dosage is never reported as consistent:** locus DOSAGE_FILTER is NOT_EVALUATED (`audit_FAILED` at every block), with
  no CONSISTENT reason, and no block reports a dosage verdict.
- **The failure is visible to the user who asked for it**, on stderr even under `--quiet`: `WARNING: the requested dosage
  audit did not complete (FAILED: …); the prediction is unchanged and DOSAGE_FILTER is NOT_EVALUATED`.
- **The record is kept as it ran:** the qualification is appended to M15-D's RESULT/README (its "audit completes"
  expectation was incorrect for this known case), and the original failed unit `a2.ex.dosage_absent_cyp2d6` stays
  recorded.

## 4. Recorded regression run

- **Focused reproductions:** `e2e.sequence_first` passed (48.1 s). It covers:
  - termination;
  - preflight;
  - the corrected S2-A stop;
  - **the restored API dosage-off provenance** (`sf21_off`: `dosage_mode off`, source `caller_options`, no audit
    artefact).
- **Then the full enabled suite on the final revision:** **63 of 63 passed** (1,006.3 s, serial, 900 s per test).
- **Disabled, reported separately:** `panvar_real_regressions` and `panvar_rebuild_stats`.
- **What is development, not campaign units:** the three earlier subset runs and the hook smoke run (240 s). M15-D's
  ledgers are not rewritten.

## 5. Remaining limitations

- **Preflight covers declared-library errors only.** Unreadable reads, graph or model refusals, resource stops and failed
  writes are found after the run has cleared its prefix.
- **A locus-table failure leaves the published sequence report in place.** Its sequences are complete, and it is marked
  `exit_code 2` and `locus_output FAILED`; it is not withdrawn.
- **The v2-without-early-incumbent stop record** is untested by a deterministic fixture.
- **The CYP2D6 unit-library limitation stands:** the audit cannot assess this cell. Other loci or cells may hit the same
  refusal.
- **The legacy schema's completed status has no `exit_code` row.** Agreement there is between the process and its
  `status` / `call_class` only.
- **Still open from M15-D:** the filter review and the advanced help.

**Stopped for review.**
