# assessment_v1 (amended in M15-G): the primary v2 assessment

**Reference implementation:** `assessment_spec.py`. Its strings are normative, and the C++ must reproduce them byte for
byte. This page is the prose.
- **Supersedes:** M15-F's draft (`tests/analysis/m15f/CONTRACT.md`), with the corrections C1–C5 below.
- **Not changed:** the old headline fields keep their meanings.

## 1. Where it is written

| level | file | how |
|---|---|---|
| locus | `.sequence_locus.tsv` | rows appended, in this order: `assessment_rule`, `sequence_resolution`, `genotype_check`, `genotype_check_coverage`, `dosage_check`, `dosage_check_coverage`, `read_support`, `read_support_coverage`, `phase`, `phase_reason`, `findings`, `not_assessed`, `dosage_alternatives`, `assessment_summary` |
| locus mirror | `.locus_genotype.tsv` | the same names as columns, plus `headline_filters_status superseded_by_assessment_v1` |
| block | `.sequence_blocks.tsv` | columns appended: `block_resolution`, `genotype_check`, `genotype_check_reasons`, `dosage_check`, `dosage_check_reason`, `dosage_alternative`, `read_support`, `read_support_reason` |
| block mirror | `.block_genotypes.tsv`, `.genotypes.tsv` (when they exist) | the same columns, matched by `block_id` |

**v2 runs** (`sequence_first_v2`) write the view on every outcome:
- **Completed and capped runs:** after the locus table.
- **Watchdog stops:** in the stop record (NO_PREDICTION) or the kept checkpoint (PROVISIONAL).
- **Failures:** in the failure record.

**Legacy-schema runs** do not write it. Neither does the caller-only API (`run_genotype`, no command layer); that is
documented.

## 2. The rules

**Verified context.** A block's legacy facts (representation guard, library sensitivity, dosage row, `genotype_support`)
are used only when all of these hold:
- the v2 prediction is a completed PREDICTION;
- the block is IDENTIFIED with a known genotype;
- the legacy decision exists (a strict or genotype-equivalent call);
- the legacy reported alleles, sorted, equal the v2 published genotype.

**Otherwise every check of that block is NOT_ASSESSED, with a reason:** `provisional`, `genotype_not_identified`,
`legacy_decision_<no_call_reason>`, `no_legacy_block` or `legacy_block_differs`. No fact is borrowed from another
candidate or from stale output.

**Block checks:**

| check | values |
|---|---|
| `block_resolution` | IDENTIFIED · LABELS_AMBIGUOUS (v2 `REPRESENTATION_AMBIGUOUS`: one sequence pair, different labels) · GENOTYPE_UNRESOLVED · PROVISIONAL_LABEL |
| `genotype_check` | CONFLICT (LOWQUAL_REPRESENTATION_CONFLICT, or library SENSITIVE at the block) · NONE_DETECTED (the guard ran) · NOT_ASSESSED (unverified, or the guard did not run). Dosage is never folded in |
| `dosage_check` | NOT_REQUESTED (no `--dosage-check`) · UNAVAILABLE (`AUDIT_*`, `NOT_VERIFIED:<why>`, UNDERPOWERED, `NOT_EVALUATED:<r>`) · NOT_APPLICABLE · AGREES · SUSPECT · DISAGREES (TRIGGERED, REPAIRED-not-applied, CONTENT_REPAIRED). REPORTED_IS_DOSAGE_REPAIR is AGREES (legacy repair mode only) |
| `dosage_alternative` (diagnostic) | NA · NONE · COMPLETE_REPAIR_AVAILABLE · `COMPLETE_REPAIR_REFUSED:<guard>` · CONTENT_CORRECTION_AVAILABLE · `REPAIR_NOT_FOUND:<reason>` · APPLIED. The audit-mode not-applied record comes from this run's own headline reasons |
| `read_support` | READS · PRIOR_DEPENDENT · NOT_ASSESSED (unverified, or `genotype_support` not READS or PRIOR_ONLY: **never READS**) |

**Locus aggregation.** A finding is reported whenever present, at any coverage:

| check | value | coverage |
|---|---|---|
| `genotype_check` | CONFLICT > NONE_DETECTED > NOT_ASSESSED | `<assessed>/<blocks>` |
| `dosage_check` | DISAGREES > SUSPECT > AGREES; otherwise NOT_APPLICABLE (every block not sensitive) or UNAVAILABLE; NOT_REQUESTED when off | `<assessed>/<applicable>` (NA when not requested) |
| `read_support` | PRIOR_DEPENDENT > READS > NOT_ASSESSED | `<assessed>/<blocks>` |

The **`assessment_summary`** token carries `[a/n]` whenever coverage is partial (0 < a < n), so a positive value is never
shown as complete when it is not.

**Phase** is locus-level, decided in this order:
1. NO_PREDICTION → NA.
2. PROVISIONAL → NOT_ASSESSED.
3. Several sequence pairs → UNRESOLVED.
4. Every block's zygosity known and **fewer than two HET blocks anywhere** → NOT_APPLICABLE.
5. The legacy decision is not a strict call → NOT_ASSESSED (`legacy_decision_<…>`). The legacy classification never
   replaces the v2 resolution.
6. A block unverified → NOT_ASSESSED.
7. Otherwise from `phase_support`: READ_SUPPORTED · PRIOR_ONLY · NOT_ASSESSED.

**The lists:**
- **`findings`:** `CODE@ids` in first-occurrence order. Codes: REPRESENTATION_CONFLICT, LIBRARY_SENSITIVE,
  DOSAGE_DISAGREES, DOSAGE_SUSPECT, PRIOR_DEPENDENT, LABELS_AMBIGUOUS, GENOTYPE_UNRESOLVED, PHASE_PRIOR_ONLY,
  PHASE_UNRESOLVED.
- **`not_assessed`:** `GENOTYPE_CHECK:<why>`, `LIBRARY_SENSITIVITY:<why>`, `DOSAGE:<why>` (audit-level once, at `@locus`,
  with `AUDIT_FAILED(<note>)`), `READ_SUPPORT:<why>`, `PHASE:<why>`.
- **Both lists are always written,** and `@all` abbreviates the full block list.

## 3. The corrections and their fixtures

| | correction | fixture |
|---|---|---|
| C1 | the legacy classification never overrides v2 resolution | G1 (a unique v2 pair beside a legacy genotype-equivalent call: RESOLVED, phase NOT_ASSESSED, not UNRESOLVED) |
| C2 | phase applicability is two or more HET blocks anywhere | G2 (HET · HOM · HET: READ_SUPPORTED; the old headline read NOT_APPLICABLE) |
| C3 | unknown support is never READS | G3 (READS, unknown, unverified: `READS` with coverage `1/3`, `SUPPORT_READS[1/3]`, both missing blocks listed) |
| C4 | provisional never yields a decision set or completed support | G4 (PROVISIONAL with legacy rows carrying PRIOR_ONLY and SUSPECT: every check NOT_ASSESSED, no findings) |
| C5 | findings and missing assessments both visible; coverage counts | G5 (AGREES `1/3` beside UNDERPOWERED and NOT_EVALUATED blocks); S10 (a conflict beside an incomplete library audit) |
| — | verified context | G6 (a legacy block whose alleles differ from the published genotype: its conflict is not reported; NOT_ASSESSED `legacy_block_differs`) |

## 4. Compatibility

**Kept unchanged** for one experimental release: `genotype_filter`, `dosage_filter`, `phase_filter`, their reasons and
counts, `call_summary`, `headline_rule`.
- Same values, same historical meanings.
- Labelled `headline_filters_status superseded_by_assessment_v1` in the locus table.
- Documented as superseded.

**Unchanged component diagnostics:** `filter`, `quality_filter`, `locus_filter`, `genotype_support`, `phase_support`, the
library and dosage tables.
