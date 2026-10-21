# M15-D.1: basic-command correctness (frozen before execution)

**Scope.**
- **Unchanged:** the basic interface from M15-D. Only M15-D's findings F2 (preflight) and F3 (exit agreement) are fixed,
  and F1 (the CYP2D6 audit) is verified and reported honestly.
- **Not done:** no index CLI, filter redesign, advanced-help reorganisation, model tuning, reserved-family rerun,
  panphorte integration, graph regeneration, worktree cleanup, merge, push or tag.

**Logistics:**
- **Branch:** `m15-d1-correctness`, from `m15-d-basic-command` `6cbd717`.
- **Results:** `tests/results/m15d1_correctness_20261021/`.
- **Budget:** 1,800 s, including 240 s of disclosed development (see `m15d1_run.py`).
- **Resources:** 4 threads, one process at a time, 6 GiB.

## 1. Changes

**Declared-library preflight.**
- **What:** a `--library FILE` is read and validated in full (`read_library`: existence, readability, keys, values,
  format) **before** the run clears its prefix. The check sits at the start of the run body, so the CLI, the direct
  command API and the prepared-locus route all get it.
- **Effect:** a missing, unreadable or invalid file refuses the run (exit 2, with the file and the problem named), and
  every file already at the prefix is left byte-identical.
- **The guarantee, documented as such:** only argument, profile and declared-library errors are refused before the prefix
  is touched. Errors found while running (reads, graph, model, resources, writing) are not transactional.

**Termination statuses.** The process exit, `.sequence_locus.tsv` `exit_code` and the status (`sequence_exit_code`,
plus a new `exit_code` row for stops and failures) state the same value:

| outcome | exit | what changed |
|---|---|---|
| time stop before a checkpoint | 124 | the sequence record and status said 1 |
| footprint stop before a checkpoint | 125 | the sequence record and status said 1 |
| stop after a committed checkpoint | 3 | the kept PROVISIONAL record now carries `exit_code 3` |
| an exception during the run | 2 | the status gains `exit_code 2` |
| a v2 run without the early incumbent | — | a stop or failure now also writes the NO_PREDICTION sequence record and status rows |
| locus table cannot be written after the sequence report was published | 2 | the record and status said 0. Now both say 2, and `locus_output FAILED: <reason>` names the failure. The published sequences stay |
| checkpoint publication failure | 2 | unchanged (the existing contract) |
| ordinary no prediction | 1 | unchanged |
| completed prediction | 0 | unchanged |

**A deterministic footprint stop:** a new test hook, `test_interrupt_point footprint_before_first_evaluation` (a
RESOURCE_REFUSAL / 125 stop at the same point as the existing time hook).

**A requested dosage audit that does not complete** is reported on stderr, even under `--quiet`, as a warning. The
prediction and every filter are unchanged by this.

**Docs:** the guide's and module page's outcome tables state the actual preflight guarantee and the exit agreement.

## 2. Tests (frozen; development controls only)

| unit | what | frozen expectations |
|---|---|---|
| `build` | every target | build passes; hashes of `panvar` and the e2e test binary recorded |
| `e2e.sequence_first` | the e2e sequence-first gates (fixtures as the full run writes them) | pass, including: **`test_termination_consistency`** (124 and 125 before a checkpoint, 3 PROVISIONAL, 1 no prediction, 0 prediction, 2 for a checkpoint publication failure and for a locus-table failure; every exit field checked, not only the return code); **`test_declared_library_preflight`** (missing, unreadable, repeated key, unknown key, missing key, negative sd, non-finite; CLI and API; sentinels byte-identical, nothing added); the **restored API dosage-off provenance** (`sf21_off`: `caller_options`, no audit artefact); the updated S2-A stop (`sequence_exit_code` and `exit_code` 124) |
| `ex.library_preflight` | C4 binary: missing, unreadable, repeated-key and unknown-key library files, with five sentinel files at the prefix | exit 2, the message names the problem, every sentinel byte-identical, nothing added |
| `ex.term_time` / `ex.term_footprint` | C4 binary with a one-key test profile (`before_first_evaluation` / `footprint_before_first_evaluation`) | 124 / 125; `NO_PREDICTION`; `TIME_LIMIT` / `RESOURCE_REFUSAL (footprint_limit)`; every exit field equal; no sequences |
| `ex.term_checkpoint` | C4 binary, `after_checkpoint` | 3; PROVISIONAL (`budget_stop`); every exit field 3; the checkpoint's sequences kept |
| `ex.term_nocall` | C4 binary, profile `early_incumbent off` + `--repair-max-neighbours 1` | 1; NO_PREDICTION; report and status exit 1 |
| `ex.term_locus_failure` | C4 binary, a directory blocking the locus table's staging file | 2; the report and status exit fields 2; `locus_output FAILED: …` in both |
| `ex.cyp_audit` | CYP2D6 HG02769, declared library, `--dosage-check` | exit 0; prediction identical to the archived M11-D run **and** to M15-D's run without the audit; `dosage_audit_status FAILED` with "alleles disagree on their flanks or anchor trims"; locus DOSAGE_FILTER NOT_EVALUATED (`audit_FAILED`) with no CONSISTENT reason; no block reports a dosage verdict; the stderr warning present (the run uses `--quiet`) |
| `ctest.full` | the full enabled suite on the final revision, serial, 900 s per test | every enabled test passes; disabled tests reported separately |

**Failure handling:**
- A failed check is a result. It is recorded and never retuned, and an amendment is declared before any rerun.
- Development checks are disclosed as development, never as campaign units.
- Earlier ledgers are not rewritten.

## 3. Deliverable

- the final source revision and binary hashes;
- the unit results;
- the remaining limitations.

The batch then stops for review.
