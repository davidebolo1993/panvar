# M15-J design: simplify the public interface (frozen before implementation)

**Branch:** `m15-j-interface`, from `m15-i-release` (`31d99d7`; code identical to the M15-H freeze `cd51d37`).

**Scope:** the command interface only. No change to genotype inference: scores, candidates, search, priors, decisions,
assessment rules and output files are unchanged. No donors, reserved data, tuning, merge, push or tag.

**Budget:** 2,400 s for the recorded units (the runner admits a unit only if spent + its timeout fits).

## 1. Public boundary

The **public interface** is the command line and the profile file. The **internal API** is `GenotypeCommandOptions` and
the functions that take it. Historical behaviour stays reachable through the internal API only, where regression tests
and diagnostics need it.

## 2. Changes

### 2.1 One reporting contract

- `--output-schema` is removed: refused with a message naming the fixed contract.
- The profile key `output_schema` is refused ("not a profile setting"), so no profile reaches an older contract.
- Ordinary genotyping keeps the preset's `sequence_first_v2`. The settings record keeps its existing `output_schema` row
  (value and source), which is the reporting-contract version in metadata. No row or file is renamed (item 4).
- The internal API keeps `output_schema` (legacy, `sequence_first_v1`) for regression tests.

### 2.2 The unrefined route: a dependency, reported, not changed

`--adaptive-repair off` / `--repair-rounds 0` today switch the preset's `sequence_first_v2` to `legacy`, recorded as
`output_schema legacy` with source `preset short-read (inactive: the unrepaired view has no pooled decision to report)`.
This batch **does not change** that behaviour. It is reported as a decision (RESULT §decisions), with these findings:

**Can the unrefined route keep v2 without adding scoring? Yes, for the scores; no, for the decision rule.**
- The v2 report is a reporting layer over `decide_pooled_genotype`, a pure function of exactly scored pool candidates
  (allele labels, sequence digests, `DiplotypeScore`, sequence and path priors).
- The refined route seeds its pool from the union candidates with exactly those quantities, all already computed by the
  union pass (`adaptive_completion`, "the pool: the exactly scored union"). A pool of the union candidates alone needs no
  new scoring.
- But the unrefined route decides by a different rule set over the same scores (`rank_set`, then its own tail):
  - **ties:** block consensus (tied candidates agreeing at a block emit that block; disagreeing blocks are no-called)
    against v2's top-set decision (`SEQUENCE_TIED`, `GENOTYPE_RESOLVED_SEQUENCE_SET`, `GENOTYPE_UNRESOLVED` blocks);
  - **phase:** the optional phase audit over a K-prefix (`phase_audit_candidates`, `READS_AUDITED`) against the pooled
    decision's phase support;
  - the full-beam audit on a tied, exhaustive search (the union candidates it uses are the same set a pool would get).
- So v2 for the unrefined route is feasible as a reporting adapter (zero-round pool, no new scoring), but the published
  outcome would follow v2's decision rules, which differ from today's unrefined outcome on ties and phase. Whether that
  counts as changing inference is Davide's decision.

**Options for the decision:**
- **A (recommended):** the unrefined route reports `sequence_first_v2` from the zero-round pool (union candidates,
  `decide_pooled_genotype`); the unrefined legacy tail (block consensus, the K-prefix phase audit) remains internal.
  Scores, candidates and the union ranking are unchanged; tie and phase outcomes follow v2.
- **B:** keep the unrefined route as an explicitly labelled development view with legacy reporting (today's behaviour),
  documented as the one exception.
- **C:** remove the public refinement controls. Not proposed: the brief says not to without approval.

### 2.3 Dosage: audit only

- `--dosage-mode` and `--dosage-audit` are removed: refused with a pointer to `--dosage-check`.
- `--dosage-check` is unchanged: the audit of the published alleles, never replacing the prediction.
- No profile key enables dosage repair: `dosage_mode` is not and stays not a profile setting; it is refused by name
  ("not a profile setting") rather than as an unknown key. `output_schema legacy`, the only contract allowing repair, is
  unreachable from the command line and profiles (2.1).
- The repair machinery stays internal (the `dosage_mode` field; `genotype_dosage_mode_gate.py` and e2e use it through the
  API).

### 2.4 Legacy and development spellings removed

Refused with a clear message, never ignored: `--phase-audit-candidates`, `--dosage-allocation-test`,
`--library-diagnostic`, `--repair-origin-cache-bytes`, `--repair-tie-seed-cache-bytes`.
- Kept as expert profile keys: `repair_origin_cache_bytes`, `repair_tie_seed_cache_bytes` (cache budgets) and
  `phase_audit_candidates` (it affects only the unrefined route's phase audit; it goes with decision 2.2).
- Kept internal: the `dosage_allocation_test` and `library_diagnostic` fields.
- `--max-candidates` stays refused, as before.
- Refusal happens during argument parsing, before the output prefix is touched.

### 2.5 `--profile`

Documented: optional overrides, never required; literal `key<TAB>value`; a required `profile<TAB>name` row; `#` comments
and blank lines; a key set twice is refused; precedence default < preset < profile < flag.

**One behaviour change, declared here:** a repeated `--profile` is refused. Today each file is applied in order and later
files override earlier keys, but the status's `profile`, `profile_name` and `profile_sha256` rows name only the last file,
so file-level provenance is incomplete. Refusing is the smallest consistent rule; merging profiles is left to the user.

### 2.6 Outputs (item 4)

No output file, column, assessment category or filter rule is deleted, renamed or restructured. The assessment's
dependency on compatibility reports is inventoried in §3 with a retirement proposal; nothing is retired here.

## 3. The assessment's dependency on compatibility reports

`src/genotype_assessment.cpp` (`read_facts`, `assess`) reads:

| file | fields read | used for |
|---|---|---|
| `.sequence_locus.tsv` (v2) | `prediction`, `completion`, `sequence_status`, `dosage_mode`, `legacy_call_class`, `legacy_no_call_reason`, `legacy_phase_support` | outcome and context |
| `.sequence_blocks.tsv` (v2) | `block_id`, `genotype`, `block_status` | the published blocks |
| `.mosaic_status.tsv` (legacy status) | `status`, `call_class`, `no_call_reason`, `library_provenance`, `library_sensitivity_status` / `library_sensitivity`, `library_sensitivity_incomplete`, `library_sensitive_blocks`, `library_sensitivity_not_evaluated_blocks`, `phase_support` | verified context, library check, phase |
| `.block_genotypes.tsv` (legacy block table) | `allele1`, `allele2` | verifying that the legacy block equals the v2 genotype |
| | `filter` | representation-guard conflict (`genotype_check`) |
| | `genotype_support` | `read_support` |
| | `dosage_status` and the dosage fields, `dosage_filter_reasons`, `dosage_stage_filter`, `dosage_stage_filter_reason`, `dosage_repair_applied` | `dosage_check` |
| `.dosage_status.tsv`, `.dosage_blocks.tsv` | `dosage_audit_status`, `dosage_audit_note`; per-block dosage rows | `dosage_check` |

It also writes its columns into `.block_genotypes.tsv`, `.genotypes.tsv` and `.locus_genotype.tsv` (mirrors).

**Proposal for retiring the compatibility reports without losing assessments (not executed):**
1. Have the producers that already compute these per-block facts (the pooled decision's block read support, the
   representation guard's verdict, the dosage audit) write them into one internal facts file per run, or into v2 tables
   as new columns, at the same point they write the legacy table. No new scoring: the values are already computed.
2. Point `read_facts` at that source; the "verified context" check (legacy alleles equal the v2 genotype) becomes
   "the facts were computed for the published representative", recorded by the producer.
3. Gate: `assessment_v1` byte-identical on the 606 saved development prefixes (the M15-G replay) and on the C4, CYP2D6 and
   ACOT controls, before and after.
4. Only then retire the writers (`.block_genotypes.tsv`, `.genotypes.tsv`, `.locus_genotype.tsv` mirrors, the headline
   columns) in a separately approved step.

## 4. Verification (units)

| unit | checks |
|---|---|
| `build` | all targets; `panvar` SHA-256 |
| `cli_contract` | the rewritten CLI contract test (refusals, profile rules, help, inventory) |
| `ex.ordinary_{c4,cyp2d6,acot}` | the M15-G `a3.ex.complete_*` command (automatic library): every output identical apart from names, timing and the exact caches' scheduling-dependent counters (the M15-H rule) |
| `ex.dosage_check_cyp2d6` | the M15-G `ex.dosage_cyp2d6` command: outputs identical (as above); prediction identical to the ordinary run |
| `ex.dosage_check_c4` | `--dosage-check` on C4: prediction, sequences and scores identical to `ex.ordinary_c4` |
| `ex.refusals` | every removed option and forbidden profile key (and a repeated profile): exit 2, the expected message, the output prefix untouched (an existing file kept byte-identical, no new file) |
| `ex.index_c4` | the M15-E index build and load: index bytes identical to M15-E's apart from recorded paths; indexed run equal to the cold run |
| `ex.unrefined` | `--adaptive-repair off` on C4: behaviour unchanged (legacy, recorded as inactive v2); evidence for decision 2.2 |
| `ctest.full` | the full enabled suite on the final revision, serially, 900 s per test |

Amendments are declared in `m15j_run.py` before any rerun; failed units stay recorded.
