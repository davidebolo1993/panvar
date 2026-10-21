# M15-D result: the basic short-read command runs without a profile; the accelerated configuration, the automatic library and the opt-in dosage check behave as frozen, with three findings

**Run.**
- **Branch and commits:** `m15-d-basic-command`, from `c0e2dea`. Freeze `0c22bd1`; amendments A1/A2 `f4513f8` (declared
  before rerun); result after.
- **Results:** `tests/results/m15d_basic_command_20261021/`.
- **Budget:** 2,170.7 of 2,400 s. That includes 660 s of disclosed development (300 s before the freeze, 360 s for the
  amendments).
- **Resources:** 4 threads, one process at a time, peak footprint ≤ 2.6 GiB.
- **Not done:** no reserved family, filter redesign, advanced-help reorganisation, model tuning, merge, push or tag.

## 1. The basic command

```
panvar genotype --graph locus.gfa --reads1 R1.fq.gz --reads2 R2.fq.gz --threads 4 --out sample
```

`--reference-path` is needed only when the graph names no single reference path. The development graphs carry no `RS:Z`,
so the examples pass it. The refusal lists the candidate paths.

**Before / after defaults** (no profile; a flag or profile still overrides every row; full table in `DESIGN.md` §1):

| | before (`c0e2dea`) | after |
|---|---|---|
| exact backends and routes (indexed fine union, fine adaptive discovery, fused seed, fixed-geometry gap DP, interned census, blocked union, reuse-only round-1 admission) | reference routes | the release candidate's, recorded `preset short-read (default)` |
| cache budgets and consumption capacities; `max_state_classes` | mostly 0; 2e8 | the release candidate's; 6e8 |
| `footprint_limit_bytes` / `time_limit_seconds` | none / none | **6 GiB** / none |
| `repair_max_neighbours` | 670 (CYP2D6 and ACOT needed profiles) | **`model`**: the model's own structural maximum (C4 670, CYP2D6 958, ACOT 998), recorded |
| `repair_junction_proposals` (JP1.1) | off | **on**: candidate generation, recorded as such |
| `early_incumbent` | off | on: completed runs unchanged; a budget stop keeps a PROVISIONAL checkpoint |
| report | legacy | `sequence_first_v2` + `compound_annotation contract_v2` |
| dosage | off (legacy), **audit implied by v2** | **off unless `--dosage-check`**; reported NOT_EVALUATED (`dosage_mode_off`) |
| `--dosage-mode repair` | legacy only | unchanged: `--output-schema legacy --dosage-mode repair` |
| library | omitted = auto | unchanged, and now tested identical to `--library auto` |
| profile | needed | optional expert override |

**Preset features whose requirements a flag removes** (`--adaptive-repair off`, `--repair-rounds 0`, `--output-schema
legacy`) are switched off and recorded `(inactive: <reason>)`. The same features set explicitly keep their refusals.

## 2. Tested examples (development controls: C4 NA19705, CYP2D6 HG02769, ACOT HG03209, 5×)

| test | result |
|---|---|
| **no profile, no `--library`, no dosage flag** | **3 of 3 PREDICTION, exit 0** (C4 12/12, CYP2D6 12/12, ACOT 12/12 checks). Every accelerated key comes from the preset; no time limit; nothing from a profile; `model` bound recorded; DOSAGE_FILTER NOT_EVALUATED (`dosage_mode_off`). FASTA rule: C4 and CYP2D6 off-panel → one FASTA with both records (SHA-256 match); ACOT panel pair → no FASTA, extraction references verified |
| omitted vs `--library auto` (C4) | identical prediction, estimate, sensitivity tables, sequence and block tables, and every setting. Library status differs only in timing and the provenance row `library_invocation` (12/12, A2) |
| **declared library + `--dosage-check` vs the archived M11-D accelerated runs** | **identical prediction at all three cells** (rank-1 score, block paths, haplotype SHA-256 and origin, status, fingerprint). Sequence tables are identical; the FASTA is identical apart from prefix-derived record names. Settings differ only in the declared set (watchdog, diagnostics, neighbour cap, absent built-in `repair_joint_search off`); dosage status is equal apart from timing (8/8 each; C4 and CYP2D6 via A2) |
| dosage absent vs `--dosage-check` | identical prediction, sequence tables, FASTA, scores and candidates at all three cells. GENOTYPE_FILTER is identical. The audit adds only the eleven `.dosage*` diagnostics and no prediction FASTA. **CYP2D6 12/13: see finding F1** |
| unrepaired view (`--adaptive-repair off`) | runs under legacy (exit 0). The four repair-only preset features are recorded inactive with reasons; the exact backends stay (7/7) |
| automatic library refusal (300 read pairs) | `LIBRARY_REFUSED`, exit 1, no prediction, no FASTA: "only 175 decisively placed concordant pairs (at least 500 are required)" (5/5) |
| binding time limit (one-key expert profile `time_limit_seconds 8`, ACOT) | exit 124 before the first checkpoint, `NO_PREDICTION` (`budget_stop_before_first_evaluated_candidate`); the minimal status names the profile and limit (3/3, A2). **See finding F3** |
| invalid combinations | 9 of 9 refused, exit 2, with the message naming the problem, the prefix untouched for the other eight: `--dosage-check` with `--dosage-mode off` (both orders), with `--dosage-audit off`, and with the legacy repair; `--dosage-mode repair` under v2; explicit v2 with `--adaptive-repair off`; no reference on a graph without `RS:Z`; `--preset long-read`. **For a missing library file, see finding F2** |

**Regression suite:**
- **`ctest.full`:** 62 of 63 passed (981.9 s). The one failure was `command_e2e`, which still expected "v2 implies audit".
- **A1 updated that test.** It found that my resolver edit had dropped S1-D.1's `caller_options` provenance for an API
  caller who declares dosage off; that is restored. The CLI is unaffected, so the passed example units stand.
- **The A1 rerun unit (`a1.ctest`) was not admitted by the budget** (2,170.7 + 240 > 2,400) and is recorded NOT RUN.
  `command_e2e` passed in disclosed development on exactly the committed source (169.7 s). The CLI-contract,
  sequence-report and library-estimate tests also passed there. **This is development evidence, not a recorded unit.**

## 3. Estimated parameters and runtime (observed; nothing tuned)

**Declared library** (the simulation's): insert 350 / 40 bp, intensity 0.0125125.

| cell | auto insert median / robust sd (pairs) | auto intensity (anchor-depth ratio) | sensitivity | auto prediction = declared prediction | wall: auto, no dosage | declared, no dosage | declared + `--dosage-check` | peak GiB (auto) |
|---|---|---|---|---|---|---|---|---|
| C4 | 349 / 41.5 (2,885) | 0.0126563 (1.020) | STABLE, PASS | yes | 27.9 s | 20.9 s | 20.1 s | 1.10 |
| CYP2D6 | 350 / 40.0 (3,238) | 0.0125429 (1.018) | STABLE, PASS | yes | 46.0 s | 40.8 s | 40.2 s | 1.28 |
| ACOT | 348 / 40.8 (2,578) | 0.0123615 (1.017) | STABLE, PASS | yes | 70.9 s | 51.1 s | 51.7 s | 2.08 |

**How to read the timings:**
- **Automatic library and sensitivity audit:** cost 5–20 s per run.
- **Dosage audit:** at C4 and ACOT, 1.2 and 3.4 s of stage time, within run-to-run noise.
- **Declared no-profile runs:** match M15-C's profile-based runs of the same cells on the same day (20.2 / 39.9 / 53.1 s).
  The archived M11-D `seconds_total` (10.9 / 20.6 / 27.7 s) came from an earlier run and is not a comparable baseline.

## 4. Findings

- **F1. The dosage audit fails at the CYP2D6 development cell.** `dosage_audit_status FAILED`: "unit library: alleles
  disagree on their flanks or anchor trims in sv:9".
  - This also happened in the archived M11-D run, and the status is otherwise identical, so it predates M15-D.
  - The prediction is unaffected; DOSAGE_FILTER reads NOT_EVALUATED (`audit_FAILED`).
  - This is an added audit that does not complete. The frozen check "the audit completes" was wrong for this cell, and the
    unit is recorded failed.
- **F2. A missing `--library` file is refused** (exit 2, the file named), but only after the run has started. The prefix
  then holds `.mosaic_status.tsv` and `.sequence_locus.tsv`, and a sentinel `.settings.tsv` there was removed. The
  contradictions are refused before any output; this input error is not.
- **F3. Exit-code mismatch at a pre-checkpoint time stop:** the process exits 124 (documented), but `.sequence_locus.tsv`
  records `exit_code 1`.

## 5. Remaining limitations

- **Accuracy evidence still comes from declared libraries only.** The automatic library gave the declared library's
  predictions at all three development cells (§3). That is descriptive, not an accuracy result.
- **No time limit by default.** A run is bounded by the 6 GiB sampled footprint guard and the work caps; users who need a
  wall-clock bound set one with a one-key profile.
- **The `model` neighbourhood bound** records the model's own structural maximum. It is not a tuned cap: any numeric cap
  must come from the user.
- **JP1.1 is on by default** and can change predictions. It is the configuration M15-B evaluated, and it is recorded as
  candidate generation.
- **Still open:** F1–F3 are reported, not fixed. The filter review and the advanced help are the next steps.

**Stopped for review.**

---

## Qualification appended by M15-D.1 (2026-10-21)

**F1 (CYP2D6 audit).**
- **The original failed unit stays recorded as it ran:** `a2.ex.dosage_absent_cyp2d6`, 12 of 13 checks.
- **The frozen M15-D expectation was incorrect for this known case.** It said that with `--dosage-check` "the audit
  completes and the dosage filter is evaluated from it".
- **What actually happens at the CYP2D6 development cell:**
  - the audit's unit library refuses sv:9 ("alleles disagree on their flanks or anchor trims"), so the audit reports
    FAILED, as it did in the archived M11-D run;
  - the prediction is unchanged;
  - DOSAGE_FILTER reads NOT_EVALUATED (`audit_FAILED`), never consistent.
- **The unit-library model is not changed.** M15-D.1 verifies the honest reporting of this limitation and adds a stderr
  warning when a requested audit does not complete.

**F2 and F3** are corrected in M15-D.1: the declared library is preflighted, and the termination exit fields agree.
