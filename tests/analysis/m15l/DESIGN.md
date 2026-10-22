# M15-L design: primary results beside the prefix, every other report in `<prefix>.diagnostics/` (frozen before implementation)

**Branch:** `m15-l-diagnostics-dir`, from `m15-k-unrefined-v2` (`14b58fc`).

**This is a breaking output-path change** for `panvar genotype`.
- **Not changed:** no file content (apart from file references, §3), inference, assessment or filter rule.
- **Not done:** no automatic deletion of diagnostics, no change to archived outputs, no merge, push or release.

**Budget:** 2,400 s, raised to 3,600 s before the harness freeze. The disclosed development, about 1,500 s, was larger
than usual: rebuilding the 6,000-line e2e test several times; one lost background attempt; e2e reruns that turned out to
have used a stale test binary, because a helper was defined after its first use and the build failed silently in a
pipeline; and trial runs. With each example run on two binaries and the full suite, 2,400 s could not admit the suite.

**Also updated (an internal reader):** `scripts/genotype_report/load.py` (the report tool behind
`plot_genotype_results.py`). It still finds runs by their locus table, now in `P.diagnostics/`, and identifies each run by
its primary prefix `P` (run_id unchanged).

## 1. The layout

For a command run with `--out P`, where `P = <dir>/<name>`:

| location | files |
|---|---|
| beside the prefix (`P.*`) | `.sequence_locus.tsv`, `.sequence_blocks.tsv`, `.sequence_haplotypes.tsv`, and `.sequence_haplotypes.fa` when the v2 FASTA rule writes it |
| `P.diagnostics/` (`P.diagnostics/<name>.*`) | every other report the run writes, under its usual basename: the status, settings, locus and block tables, candidates, scores, equivalence set, alternatives, representation intervals, compound annotation, library estimate and sensitivity, dosage audit files, repair logs, partial files and diagnostics |

**Unchanged:**
- **Basenames:** `<name>.<suffix>`; FASTA record names (`<name>_seq_hapN`) are unchanged too.
- **The primary outcome report** (`P.sequence_locus.tsv`) is written at its normal location whenever it is produced:
  predictions, provisional checkpoints, no-predictions, watchdog stops and reporting failures.
- **Not affected:** `--build-panel-index` (it writes only the index), and other commands.

## 2. Mechanism

**The run's working prefix is the diagnostics prefix.**
- The command (`parse_genotype_arguments` sets `diagnostics_directory`) runs with `out_prefix = D = P.diagnostics/<name>`.
- One layout module (`genotype_output_layout.{hpp,cpp}`) registers `D → P` for the run.
- `genotype_output_path(D, suffix)` returns `P + suffix` for the four primary suffixes (and their `.writing` temporaries),
  and `D + suffix` otherwise.

**What changes in the code:**
- Every writer, reader and remover of a primary file goes through that function:
  - the sequence report's publication and failure handling;
  - the watchdog's stop records;
  - the committed-record check;
  - the failure handler's exit-code marking;
  - the assessment's reads and column rewrites;
  - the locus table's `main_fasta` link.
- Every other path is `D + suffix` as before. So the dosage stage, headline filters, locus table, assessment mirrors,
  settings, status, partial files and run-owned cleanup move with no further change.

**The internal API** (`GenotypeCommandOptions` built directly, as most e2e tests do) keeps the flat layout unless
`diagnostics_directory` is set. The command line always sets it.

## 3. Content and file references

**Contents are unchanged, apart from:**
- the locus table's `main_fasta`, a path relative to the table, now `../<name>.sequence_haplotypes.fa`;
- paths the run records that contain the output prefix. Those now contain `P.diagnostics/<name>` where they name a diagnostic,
  and `P` where they name a primary file. Each is listed per unit.

Every link between diagnostics is relative to its own table, so it stays correct.

## 4. Cleanup, checkpoints and failures

**Before a run** (after every preflight refusal, so a refused run still touches nothing):
- `P.diagnostics/` is created;
- the run-owned reports are removed at their new locations;
- the same run-owned names are also removed at the old flat locations `P.<suffix>`, so an earlier run's flat reports cannot
  sit beside the new layout and be mistaken for this run's. Only this command's own report names are removed; nothing else
  in the directory is touched.

**During and after a run:**
- **Checkpoints:** an early-incumbent checkpoint publishes its PROVISIONAL record at `P.sequence_locus.tsv`.
- **Watchdog stops:** the run-owned files are removed at their locations; the stop record is written at
  `P.sequence_locus.tsv` and its status in `P.diagnostics/`.
- **No automatic deletion:** diagnostics are never removed after a completed run.

## 5. Verification (units)

| unit | checks |
|---|---|
| `build` | all targets |
| `tests` | the CLI, curation and e2e command-route tests (paths updated), the layout test |
| `ex.layout_{c4,cyp2d6,acot}` | the ordinary command with the M15-L binary and with the M15-K binary (`baseline/panvar.m15k`). Exactly the primary files beside `P`; every other file in `P.diagnostics/`; every file's content equal to the M15-K file of the same name (M15-H rule: names, timing, scheduling counters), apart from the declared reference changes (§3); exit codes equal |
| `ex.dosage_check_cyp2d6` | the same with `--dosage-check` and the declared library |
| `ex.indexed_c4` | the same with `--panel-index` (the M15-K index, so both binaries load one file) |
| `ex.provisional` | the same for a provisional run (a time limit with the early incumbent): `P.sequence_locus.tsv` PROVISIONAL, exit 3 |
| `ex.failed` | the same for a watchdog stop (exit 124/125): the stop record at `P.sequence_locus.tsv` where the contract writes one. A refused run (exit 2) touches nothing: no directory is created |
| `ex.stale_flat` | a prefix holding an earlier run's flat reports: those report names are removed; an unrelated file is kept |
| `ctest.full` | the full enabled suite on the final revision |
