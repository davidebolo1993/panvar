# M15-L result (BREAKING output-path change): only the primary sequence reports beside the prefix, every other report in `<prefix>.diagnostics/`; every file's content equal to M15-K's on C4, CYP2D6, ACOT, the dosage check, an indexed run, provisional and failed outcomes; full suite 65/65

**Run.**
- **Branch:** `m15-l-diagnostics-dir`, from `m15-k-unrefined-v2` (`14b58fc`).
- **Commits:** design frozen `5bc55c7`; implementation and harness frozen `bc19abc`; amendment A1 (harness only); result
  after.
- **Results:** `tests/results/m15l_diagnostics_dir_20261022/`, including `baseline/panvar.m15k`, the M15-K binary.
- **Budget:** 2,468.7 of 3,600 s (raised from 2,400 s before the harness freeze; 1,500 s of disclosed development, an upper
  estimate; DESIGN.md).
- **Units:** 9 of 11 passed on the first pass. `ex.provisional` and `ex.failed` failed on harness expectations and passed as
  `a1.*` (§4).
- **Not done:** no content change apart from the declared reference, no inference, assessment or filter change; no
  automatic deletion of diagnostics, no change to archived outputs; no merge, push or release.

| binary | SHA-256 |
|---|---|
| `panvar` (`bc19abc`) | `7dc25db6cd5de2c82d883d0996edb4a98d732c651027b8f2cd6f1b3e131c9d30` |
| baseline `panvar.m15k` | `5c4707a575f23ed4eab3b392143ad9ed1dbef8d9a3feb1e7ee78f12a9cc182e6` |

## 1. The breaking change

`panvar genotype ... --out P`, with `P = <dir>/<name>`:

| before (M15-K and earlier) | now |
|---|---|
| every report beside the prefix: `P.<suffix>` | **beside the prefix:** `P.sequence_locus.tsv`, `P.sequence_blocks.tsv`, `P.sequence_haplotypes.tsv`, and `P.sequence_haplotypes.fa` when the v2 FASTA rule writes it |
| | **`P.diagnostics/<name>.<suffix>`:** every other report: status, settings, locus and block tables, candidates, scores, equivalence set, alternatives, representation intervals, compound annotation, library estimate and sensitivity, dosage audit files, repair logs |

**Contents are unchanged.** The one declared reference change is the locus table's `main_fasta`, now relative to the table:
`../<name>.sequence_haplotypes.fa`. Links between diagnostics are relative to their table and unchanged.

**Consumers that read a diagnostic at `P.<suffix>` must read `P.diagnostics/<name>.<suffix>`.**
- The repository's report tool (`plot_genotype_results.py`, `genotype_report/load.py`) is updated. It still finds runs by
  their locus table and identifies each by its primary prefix `P`, so run ids are unchanged.
- Older campaign scripts in `scripts/` read archived flat outputs and are unchanged.

## 2. Implementation

**The layout module** (`genotype_output_layout.{hpp,cpp}`):
- The command runs under `D = P.diagnostics/<name>`, registered against `P` for the run.
- Every writer, reader and remover of a primary file asks `genotype_output_path`:
  - the sequence report's publication and failure handling;
  - the watchdog's stop records and the exception path;
  - the committed-record check;
  - the failure handler's exit-code marking;
  - the assessment's reads and column rewrites;
  - the locus table link.
- Every other path is `D + suffix` as before. So the dosage stage, headline filters, locus table, assessment mirrors,
  settings, status, partial files and cleanup moved with no other change.
- **The assessment** is generated from the moved legacy tables (unchanged results, §3).

**Lifecycle:**
- **Directory creation:** when the run first touches its prefix, after every preflight refusal, one level only. A refused
  run creates nothing.
- **Earlier runs' reports:** the run removes its own report names at their new locations and at the old flat locations.
  Nothing else is touched.
- **No automatic deletion:** diagnostics are never deleted automatically.

**The internal API** stays flat unless `diagnostics_directory` is set. The command line always sets it.

## 3. Verification

Each example ran the same command with the M15-L and the M15-K binaries. Every report must exist under its mapped name, and
its content must equal M15-K's (M15-H rule: names, timing, scheduling counters), apart from the declared `main_fasta` change.

| unit | reports (beside / in diagnostics) | result |
|---|---|---|
| `ex.layout_c4` | 20 (4 / 16) | 23/23, every file equal to M15-K; exit 0 both |
| `ex.layout_cyp2d6` | 20 (4 / 16) | 23/23 |
| `ex.layout_acot` | 19 (3 / 16), no FASTA (panel pair) | 22/22 |
| `ex.dosage_check_cyp2d6` | 19 (4 / 15), the dosage audit files included | 22/22 |
| `ex.indexed_c4` (the M15-K index on both binaries) | 20 (4 / 16) | 24/24; the status, now in the diagnostics directory, records the loaded index |
| `a1.ex.provisional` | 7 (3 / 4) | 13/13: exit 3; `P.sequence_locus.tsv` PROVISIONAL, `exit_code 3` |
| `a1.ex.failed` | stop: 2 (1 / 1) | 12/12. A watchdog stop: exit 124, the NO_PREDICTION stop record beside the prefix, the status in the diagnostics directory. A reporting failure: exit 2, `locus_output FAILED` in `P.sequence_locus.tsv`, equal to M15-K's record. A refused run (missing declared library): exit 2, **nothing created**, not even the diagnostics directory |
| `ex.stale_flat` | | 4/4: an earlier run's flat reports (`.mosaic_status.tsv`, block and locus tables, settings, legacy FASTA) are removed; unrelated files are kept; the run's reports are in the new layout |
| `tests` | | CLI contract (including the pre-block no-call's layout) and curation contract |
| `ctest.full` | | **65 of 65 passed** (603.0 s); disabled, as before: `panvar_real_regressions`, `panvar_rebuild_stats` |

**Tests updated for the new paths:**
- **e2e command-route tests** read through a layout-aware helper; the blocked-union comparison lists the diagnostics
  directory.
- **CLI contract:** the pre-block no-call's locus table is in `pre.diagnostics/`.
- **Curation test:** uses the internal API's flat layout explicitly, since it sets the internal legacy contract.

## 4. Amendment A1 (harness only; failures stay recorded)

**`ex.provisional` failed 12/13:**
- **Cause:** a provisional record carries the checkpoint's timing rows (`checkpoint_seconds_*`), which differ by run. The
  primary record now gets the same timing-row allowance as every keyed file.
- **Never in question:** layout and outcome had passed.

**`ex.failed` failed 11/12:**
- **Cause:** the record comparison replaced the primary prefix before the longer diagnostics path it begins with. The
  longer path is now replaced first.
- **Before that fix:** the records were already identical once properly normalized (checked by hand).

Both passed on rerun.

## 5. Notes for review

- **This changes paths every downstream consumer of `panvar genotype` reads:** the release notes, guide §6, module page and
  help say so.
- **`--build-panel-index`** writes only its index and is unaffected.

**Stopped for review.**
