# S1-D.1: completing the sequence_first_v2 output contract

**What it is:** two narrow fixes on top of S1-D (`38ef26e`). Scoring, search, priors, caps and the persisted prepared-locus format are unchanged. So are `legacy` and `sequence_first_v1` output.

## 1. The FASTA policy covers the whole prefix

`.sequence_haplotypes.fa` is the only prediction FASTA under v2:
- for an off-panel pair, it holds both haplotypes, each with its origin;
- for a panel pair, no FASTA is written (extraction references only);
- for no prediction, no FASTA.

| file | under v2 |
|---|---|
| `.mosaic_haplotypes.fa` (and any `.original` copy) | removed after its internal consumers have read it: the dosage stage's eligibility check and the locus table (hashes and lengths) |
| empty `.dosage_{original,repaired,repaired_block_alleles}.fa` | removed |
| non-empty dosage FASTA | kept as an audit diagnostic, linked from the locus table (not produced by any fixture here; see limitations) |
| every legacy table and audit table | unchanged |

**`.locus_genotype.tsv` under v2:** `main_fasta` links `.sequence_haplotypes.fa` (or `NA`). `original_fasta` and `hap1_id`/`hap2_id` are `NA`. `dosage_representative_fasta` links only a file that stays.

**`.sequence_locus.tsv`** gains a `legacy_fasta` row stating the policy.

**When it applies:** the policy also holds when the locus output fails. Run start still removes every stale run-owned FASTA.

## 2. One dosage contract at every command entry

`resolve_dosage_contract(options)` is called both by argument resolution and by `execute_genotype_command`, which serves the direct API and the prepared-locus command route.

| case | result | `dosage_mode_source` |
|---|---|---|
| mode not declared | audit | `contract_default_sequence_first_v2` |
| declared `off` (`--dosage-mode off`, `--dosage-audit off`, or `dosage_mode_explicit`) | off | `flag --dosage-mode`, `flag --dosage-audit` or `caller_options` |
| `repair`, declared or set by value | refused (exit 2) before the prefix is touched | — |

**Caller-only APIs.** `run_genotype` and `run_genotype_prepared` call the caller alone, with no dosage stage. Under v2 they record `dosage_mode not_run`.

## 3. Validation tools

- **`s1d1_fixtures.py`:** whole-prefix inventory and references per case; CLI dosage arms; legacy identity with the declared differences; loaded index; a lifecycle on a stale-seeded prefix; cross-library legacy/v1 reruns; `--check-prefix` for the real-data unit.
- **`s1d1_run.py`:** the frozen manifest runner.
- **API and prepared routes:** gated in `tests/genotype_command_e2e.cpp`, `test_sequence_first_v2_contract`.

## 4. Limitations (stated, not validated)

- **Whole-path frame.** Extraction assumes each template path spans the locus graph whole; no path offsets are persisted (S1-D).
- **No audit on a legacy no-call.** The audit is NOT_RUN when the legacy classification is a no-call (F1, and the provisional F7, here).
- **Non-empty dosage FASTA.** No fixture or real-data run here produced one. The retention path is implemented but not exercised.
