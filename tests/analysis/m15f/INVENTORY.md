# M15-F inventory: the current filter, support and status fields

**Source:** the code at `m15-e-panel-index` `9a7e750`, read directly. File and line references are to that revision.
- **Abbreviations:** GC `src/genotype_command.cpp`, SR `src/genotype_sequence_report.cpp`, HL `src/genotype_headline.cpp`,
  DI `src/genotype_dosage_integrate.cpp`, RP `src/genotype_dosage_repair.cpp`, SH `src/genotype_dosage_shadow.cpp`,
  AU `src/genotype_dosage.cpp`, GR `src/genotype_repair.cpp`.
- **Column codes:**
  - **Rep / Diag:** whether the field describes the reported prediction or a diagnostic alternative or candidate set.
  - **NE ≠ adverse:** whether "not evaluated" is distinct from adverse evidence.
  - **N / D:** native or derived.

## 1. Run and locus outcome

| field (file) | values and exact trigger | Rep / Diag | scope | NE ≠ adverse | N / D |
|---|---|---|---|---|---|
| `prediction` (`.sequence_locus.tsv`, v2) | PREDICTION: completed pooled ranking (SR:114). PROVISIONAL: a cap or budget stop after a fully evaluated checkpoint. NO_PREDICTION: no report built, or publication failed | Rep | locus | — | N |
| `completion` (v2) | `completed_ranking`; the stop name for PROVISIONAL (`round_cap`, `cumulative_cap`, `neighbourhood_bound`, `budget_stop`, …); `budget_stop_before_first_evaluated_candidate`; `report_publication_failed`; `run_failed`; else the legacy run status | Rep | run | — | N |
| `sequence_status` (v2) | RESOLVED: one top sequence group. GENOTYPE_RESOLVED_SEQUENCE_SET: several groups, every block IDENTIFIED. SEQUENCE_TIED: otherwise. NOT_DETERMINED: PROVISIONAL. NA: no prediction (SR:259–261) | Rep (decision set) | locus | — | D (from `block_status`) |
| `status` / `run_status` (`.mosaic_status.tsv`, `.locus_genotype.tsv`) | COMPLETE (unique decision); CONDITIONAL_GENOTYPE_EQUIVALENT; RESCORE_INCOMPLETE (ties, pool mismatches); REPAIR_INCOMPLETE; MODEL_INCOMPLETE; CANDIDATE_INCOMPLETE; LIBRARY_REFUSED; TIME_LIMIT; RESOURCE_REFUSAL; diagnostic-only statuses | Rep (legacy decision) | run | — | N |
| `call_class` | strict_call / genotype_equivalent_call / no_call, mapped from `status` only (GC:4961) | Rep (legacy) | run | — | D |
| `no_call_reason` | `none`; `distinct_tie`, `divergent_tie`, `neighbourhood_bound`, `round_cap`, `cumulative_cap`, `tie_branching_incomplete`, …; HL synthesises `status_<STATUS>` when absent | Rep (legacy) | run | — | N |
| `legacy_run_status` / `legacy_call_class` / `legacy_no_call_reason` (v2) | copies of the above | — | run | — | copy |

**The legacy and v2 outcomes diverge by design.** The v2 report is built before the legacy tie check (GC:6814, then the
tie check at 6821), so a decision that ties only over labellings reports **PREDICTION** with exit 0 beside
`legacy_call_class no_call` / `distinct_tie`. Saved development outputs show 119 such runs:
- **66** are PREDICTION with sequence_status RESOLVED;
- **12** are PREDICTION with SEQUENCE_TIED;
- **41** are PROVISIONAL.

**No legacy block table is written for them** (`.block_genotypes.tsv`, `.genotypes.tsv`): only the v2 `.sequence_blocks.tsv`
exists.

## 2. Headline filters (`headline_filters_v2`, HL)

| field | values and exact trigger | Rep / Diag | scope | NE ≠ adverse | N / D |
|---|---|---|---|---|---|
| `genotype_filter` (+ `_reasons`) | **NO_CALL** if legacy `call_class` is no_call (`RUN_NO_CALL:<reason>`). Otherwise the most severe of: NOT_EVALUATED (`REPRESENTATION_NOT_EVALUATED`, `LIBRARY_SENSITIVITY_NOT_EVALUATED`); LOW_QUALITY (`REPRESENTATION_CONFLICT`, `LIBRARY_SENSITIVE`, `EXACT_DIPLOTYPE_UNRESOLVED`, the last when the audit ran and the block is CONTENT_REPAIRED, TRIGGERED_UNREPAIRED or REPAIRED-not-applied); else PASS (HL:54–71) | Rep, but `EXACT_DIPLOTYPE_UNRESOLVED` is dosage evidence | block → locus | yes, but NOT_EVALUATED **outranks** LOW_QUALITY | D |
| `dosage_filter` (+ `_reasons`) | NO_CALL (legacy no-call); NOT_EVALUATED (`dosage_mode_off`, `audit_<status>`, `DOSAGE_UNDERPOWERED`, `DOSAGE_NOT_EVALUATED`); PASS (`DOSAGE_CONSISTENT`, `COMPLETE_REPAIR_APPLIED`, **`COMPLETE_REPAIR_AVAILABLE_NOT_APPLIED`**, **`CONTENT_CORRECTED_EXACT_DIPLOTYPE_UNRESOLVED`**); LOW_QUALITY (`COMPLETE_REPAIR_NOT_APPLIED:<guard>`, `DOSAGE_SUSPECT`, `TRIGGERED_UNREPAIRED`, the unreachable `CONTENT_CORRECTION_CONDITIONAL`); NOT_APPLICABLE (`NOT_DOSAGE_SENSITIVE`) (HL:73–100) | **mixed**: the two bold PASS reasons describe the alternative, not the reported alleles | block → locus | yes, but NOT_EVALUATED outranks LOW_QUALITY | D |
| `phase_filter` (+ `_reasons`) | per adjacent-block boundary: NOT_APPLICABLE unless both are HET; NO_CALL `PHASE_UNRESOLVED` (`UNRESOLVED` or genotype-equivalent); PASS (READS_*); LOW_QUALITY `MODEL_ONLY_PHASE` (PRIOR_ONLY); NOT_EVALUATED (`PHASE_ASSESSMENT_INCOMPLETE`, `PHASE_NOT_ASSESSED_AFTER_DOSAGE_REPAIR`) (HL:102–114) | Rep | boundary → locus, but **inherited**: the only per-boundary inputs are zygosity and an applied repair | yes | D |
| `call_summary` | `GENOTYPE_<s>;DOSAGE_<s>;PHASE_<s>` | — | locus | — | D |
| `*_filter_units / applicable / pass / low_quality / not_evaluated / no_call / not_applicable` | per-status unit counts | — | locus | — | D |

**Locus aggregation:** the most severe unit in the order NO_CALL > NOT_EVALUATED > LOW_QUALITY > PASS > NOT_APPLICABLE,
with `REASON@ids` for every reason (HL:116–151). Two consequences:
- **Any missing optional evidence hides a concrete finding in the headline value.** One UNDERPOWERED block makes a locus
  with a SUSPECT block read NOT_EVALUATED.
- **A legacy no-call makes all three NO_CALL,** even beside a v2 PREDICTION.

## 3. Block quality, representation and support

| field (`.block_genotypes.tsv` / `.genotypes.tsv`) | values and exact trigger | Rep / Diag | scope | NE ≠ adverse | N / D |
|---|---|---|---|---|---|
| `filter` (representation guard) | PASS; LOWQUAL_REPRESENTATION_CONFLICT (marker, fragment and exposure evidence each ≥ ln 1000 against the call); NOT_EVALUATED (`competitor_not_scored`, `guard_off`, `guard_not_run_on_this_path`, after a dosage repair) (GC:6853–6936) | Rep | block | yes | N |
| `locus_filter` | LOWQUAL_REPRESENTATION_CONFLICT if any block is; else NOT_EVALUATED if any block is not PASS; else PASS; `;`-appended LOWQUAL_LIBRARY_SENSITIVE / NOT_EVALUATED_LIBRARY_SENSITIVITY (GC:6945, 7582) | Rep | locus | yes | D |
| `quality_filter` (`block_quality_rule_v1`) | PASS iff none of: `representation_conflict`, `representation_not_evaluated` (fails closed), `no_read_support` (genotype_support ≠ READS), `library_sensitive` (GC:2238–2251) | Rep | block | **no:** not-evaluated folds into LOW_QUALITY (the reason shows it) | D |
| `quality_pass_blocks` / `quality_lowqual_blocks` | counts of the above | — | locus | no | D |
| `genotype_support` | READS; **PRIOR_ONLY** when a fragment contender (fragment term within bounds of the winner's) projects to a different unordered genotype at the block (GR:1548–1560); NOCALL (consensus route) | Rep vs explored rivals | block | — | N |
| `conditional_GQ` / `legacy_representative_GQ`, `error_mass_lower/upper`, `numerically_separated`, `separation(_certified)` | conditional on the explored candidate set; uncalibrated | Rep vs alternatives | block / locus | — | N / D |
| `structure_confidence` | constant NA (`UNCALIBRATED`) | — | locus | — | constant |
| `block_status` (`.sequence_blocks.tsv`, v2) | IDENTIFIED (one unordered genotype across the top groups' members); **REPRESENTATION_AMBIGUOUS** (one top group whose explored labellings disagree: the same sequence pair, different block labels); GENOTYPE_UNRESOLVED (several top groups disagree); PROVISIONAL_LABEL (SR:232–250, 330) | Rep (decision set) | block | — | N |
| `representation_interval` (`.representation_intervals.tsv`) | the block span where the representative group's labellings differ | Rep | interval | — | N |
| `phase_support` | READS_WITHIN_CANDIDATES / READS_AUDITED / PRIOR_ONLY / UNRESOLVED. **Computed once per locus** (GR:1562–1566: UNRESOLVED if the decision is not unique, READS_WITHIN_CANDIDATES if no fragment contender, else PRIOR_ONLY) and **copied to every block row** (GC:7647) | Rep | **locus** | UNRESOLVED ≠ PRIOR_ONLY | N |
| `zygosity` | HOM / HET of the reported alleles | Rep | block | — | D |

## 4. Library sensitivity (`--library auto` only)

| field (`.mosaic_status.tsv`) | values and exact trigger | Rep / Diag | scope | NE ≠ adverse | N / D |
|---|---|---|---|---|---|
| `library_sensitivity(_status)` | STABLE / SENSITIVE (a block's 0.9× or 1.1× intensity rank or pipeline projection differs) / NOT_EVALUATED (with `library_sensitivity_incomplete`: audit_off, library_refused, not_adaptive_completion, primary_no_call, call_not_reached). A declared library writes no key; HL treats it as NOT_APPLICABLE | Rep under perturbation | locus | yes | N |
| `library_sensitive_blocks`, `library_sensitivity_not_evaluated_blocks` | block lists | Rep | blocks | yes | N / D |

## 5. Dosage

**The stage (`.dosage_status.tsv`, run; native):**
- `dosage_audit_status` takes COMPLETE, REFUSED_LIBRARY (no insert distribution), NOT_RUN (no rank-1 call or no
  genotype emitted) or FAILED (an exception, with `dosage_audit_note`).
- **All non-COMPLETE values mean "not evaluated", never adverse.**
- The file is not written when dosage is off.

**Per block (`.dosage_blocks.tsv`; native):**

| `dosage_status` | trigger | about the reported alleles | NE ≠ adverse |
|---|---|---|---|
| NOT_SENSITIVE | `!sensitive()` (AU:532) | not applicable | — |
| NOT_EVALUATED | `no_tested_feature`, `no_reference`, `no_reference_hits`, `bootstrap_zero_reference`, `unsupported_library` | not evaluated | yes |
| CONSISTENT | no feature triggered, suspect or underpowered (AU:681–688) | agrees | — |
| SUSPECT | `robust_z ≥ 2.5`, not triggered (wins over UNDERPOWERED) | weak adverse | — |
| UNDERPOWERED | `robust_sd > sd_max`, not triggered or suspect | not evaluated (power) | yes |
| REPAIRED | triggered (`robust_z ≥ zc` and the call is not admissible), and `repair_class` COMPLETE_DIPLOTYPE_REPAIR | **the reported alleles are dosage-inconsistent**; a complete alternative exists | — |
| CONTENT_REPAIRED | triggered, class CONDITIONAL_DOSAGE_CONTENT | **dosage-inconsistent**; corrected content exists, exact diplotype unresolved | — |
| TRIGGERED_UNREPAIRED | triggered, no repair class. **Includes repairs that could not run** (`cap_pairs`, `cap_refused`, `resource_refusal`, `winner_not_reproduced`, `top_genotype_tie`, `no_admissible_candidate`, …; visible only in `dosage_reason`) | dosage-inconsistent | — |

**Diagnostic fields** (all describe the alternative, never the reported alleles):
- `repair_class`, `dosage_content_scope`, `diplotype_scope`, `allocation_equivalence_scope`;
- `serialization_representative`, `repaired_alleles_role`, `allocation_status`, `corrected_content`;
- `.dosage_repaired*.fa`.

`gq_status` reads UNAVAILABLE_REPAIRED for a REPAIRED or CONTENT_REPAIRED block even though, under audit, the reported GQ
is unchanged.

**Repair mode only** (legacy; refused under v2): `dosage_stage_filter` (`PASS` also for NOT_SENSITIVE) and its reason,
`dosage_repair_applied`, `filter_basis`, and the status keys `dosage_repaired_blocks` and following.

**Locus (`.locus_genotype.tsv`):**
- `dosage_mode`, `dosage_audit_status`;
- counts (`dosage_complete_class_blocks`, `dosage_content_corrected_blocks`, `dosage_triggered_unrepaired_blocks`,
  `dosage_suspect_blocks`, `dosage_underpowered_blocks`).
- **Under audit:** `dosage_complete_not_applied_blocks/_detail` read **0 / none even when complete repairs exist**, and
  `locus_repair_scope` can read COMPLETE while nothing was applied. The audit-mode not-applied record reaches the outputs
  only as the headline reason `COMPLETE_REPAIR_AVAILABLE_NOT_APPLIED`.

**In the v2 report:** `dosage_mode`, `dosage_mode_source` and the constants `dosage_scope` / `dosage_contract`.
`dosage_mode audit` does not mean the audit completed.

## 6. Documentation mismatches found

1. **`docs/guides/genotype_short_read.md` §6** lists `distinct_tie` as a NO_PREDICTION example. The code reports
   PREDICTION (RESOLVED or SEQUENCE_TIED) with exit 0.
2. **`docs/algorithms/genotype.md`** describes `phase_support` as "PRIOR_ONLY otherwise". On the repair route, any
   non-unique decision is UNRESOLVED.
