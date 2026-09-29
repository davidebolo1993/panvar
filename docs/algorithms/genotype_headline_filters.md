# Headline filters (`headline_filters_v2`)

## 1. What they are

`panvar genotype` reports three independent headline verdicts, **GENOTYPE_FILTER**, **DOSAGE_FILTER** and
**PHASE_FILTER**, with deterministic reasons, and renders them together as

    CALL_SUMMARY=GENOTYPE_<status>;DOSAGE_<status>;PHASE_<status>

- **There is no global PASS/FAIL.**
- **They are derived summaries.** The representation filter, library sensitivity, `genotype_support`, GQ, dosage
  status and repair class, allocation status and `phase_support` are all still reported. A stricter downstream view,
  for example excluding `genotype_support=PRIOR_ONLY`, is built from those fields.
- **Provenance:** the design was frozen and audited against the sealed C4 replay before implementation
  (`tests/results/c4_headline_filters_26401b4/SPEC_v2.md`). The implementation reproduces that audit's truth-free
  derived tables exactly.

## 2. Values and precedence

**Values:** PASS, LOW_QUALITY, NOT_EVALUATED, NO_CALL and NOT_APPLICABLE.

**Precedence is fail-closed:**

    NO_CALL > NOT_EVALUATED > LOW_QUALITY > PASS > NOT_APPLICABLE

- A unit's value is the most severe applying component. Every applying reason is kept, comma-separated, in rule order.
- A locus takes the most severe of its units: blocks for genotype and dosage, adjacent-block boundaries for phase.
- A locus is NOT_APPLICABLE only when every unit is, so NOT_APPLICABLE never hides an evaluated failure.
- Locus reasons are `REASON@id1,id2;...`, listing the responsible block or boundary IDs (`left>right`).
- **Per-status unit counts** are reported beside every locus value: units, applicable, PASS, LOW_QUALITY,
  NOT_EVALUATED, NO_CALL and NOT_APPLICABLE. A fail-closed locus headline with one non-PASS block is then visibly
  different from a locus whose blocks all failed.

## 3. GENOTYPE_FILTER (per called block)

| component | status | reason |
|---|---|---|
| the run emitted no genotype | NO_CALL | `RUN_NO_CALL:<no_call_reason>` |
| representation filter neither PASS nor LOWQUAL_* | NOT_EVALUATED | `REPRESENTATION_NOT_EVALUATED` |
| library sensitivity not evaluated (auto library), locus-wide or at the block after a dosage repair | NOT_EVALUATED | `LIBRARY_SENSITIVITY_NOT_EVALUATED` |
| representation conflict | LOW_QUALITY | `REPRESENTATION_CONFLICT` |
| library-sensitive block | LOW_QUALITY | `LIBRARY_SENSITIVE` |
| the dosage stage ran and found the reported alleles dosage-inconsistent without resolving the exact diplotype | LOW_QUALITY | `EXACT_DIPLOTYPE_UNRESOLVED` |
| otherwise | PASS | none |

"The run emitted a genotype" means the status `call_class` is strict_call or genotype_equivalent_call. A run with
`--adaptive-repair off` (the development view) writes no `call_class`; the headline then derives it from the run
`status`: COMPLETE is a strict call, CONDITIONAL_GENOTYPE_EQUIVALENT a genotype-equivalent call, and any other status
is a no-call with reason `status_<status>`. The component `call_class` column of locus_genotype.tsv is left as the
run wrote it (NA).

**What counts as dosage-inconsistent and unresolved:**
- a content correction (`CONTENT_REPAIRED`);
- a triggered block left unrepaired;
- a complete repair that was not applied, whether audit mode kept it out or a guard refused it. The reported alleles
  stay the caller's in both cases.

**Rules outside the table:**
- `genotype_support=PRIOR_ONLY` is reported and does not alone fail the filter.
- GQ is conditional on the scored candidate set. It never creates PASS and never overrides a component.

## 4. DOSAGE_FILTER (per called block)

| dosage stage state | status | reason |
|---|---|---|
| `--dosage-mode off` | NOT_EVALUATED | `DOSAGE_NOT_EVALUATED:dosage_mode_off` |
| stage did not complete | NOT_EVALUATED | `DOSAGE_NOT_EVALUATED:audit_<status>` |
| CONSISTENT | PASS | `DOSAGE_CONSISTENT` |
| complete diplotype repair, applied | PASS | `COMPLETE_REPAIR_APPLIED` |
| complete diplotype repair, audit mode, eligible (repair mode would apply it) | PASS | `COMPLETE_REPAIR_AVAILABLE_NOT_APPLIED` |
| complete diplotype repair refused by an eligibility, reproduction or safety guard (repair mode, or repair mode's guards evaluated in audit mode) | LOW_QUALITY | `COMPLETE_REPAIR_NOT_APPLIED:<reason>` |
| content correction, content scope COMPLETE | PASS | `CONTENT_CORRECTED_EXACT_DIPLOTYPE_UNRESOLVED` |
| content correction, content scope conditional | LOW_QUALITY | `CONTENT_CORRECTION_CONDITIONAL` |
| SUSPECT | LOW_QUALITY | `DOSAGE_SUSPECT` |
| TRIGGERED_UNREPAIRED | LOW_QUALITY | `TRIGGERED_UNREPAIRED` |
| UNDERPOWERED | NOT_EVALUATED | `DOSAGE_UNDERPOWERED` |
| NOT_SENSITIVE | NOT_APPLICABLE | `NOT_DOSAGE_SENSITIVE` |
| unrecognised | NOT_EVALUATED | `DOSAGE_NOT_EVALUATED` |
| no-call | NO_CALL | `RUN_NO_CALL` |

**A complete repair audit mode did not apply** has usable dosage content: the stage resolved the whole diplotype and
repair mode's guards (strict call, a single complete repair, panel alleles, reproduced haplotypes) all hold. Only the
selected mode kept it from the reported genotype, which stays the caller's and dosage-inconsistent, so the genotype is
LOW_QUALITY `EXACT_DIPLOTYPE_UNRESOLVED`. A guard's refusal is different: the repair could not be confirmed, so the
dosage is LOW_QUALITY with the guard's reason (`not_a_strict_call`, `several_complete_repairs_not_jointly_scored`,
`repaired_allele_not_a_panel_allele`, `reported_haplotypes_not_reproduced`). A reason that was not recorded reads
`COMPLETE_REPAIR_NOT_APPLIED:not_recorded`, never PASS. Audit mode evaluates the guards as a record only
(genotype_dosage.md §22.2); no candidate, repair decision or reported genotype changes.

**A content-only correction is never an exact allele call.** Its genotype is LOW_QUALITY `EXACT_DIPLOTYPE_UNRESOLVED`
and its dosage is PASS. The reported alleles stay the caller's, and the corrected content is reported only in the dosage
fields (`corrected_content`).

**The default dosage mode is off,** so dosage is honestly NOT_EVALUATED by default. It is never silently PASS.

## 5. PHASE_FILTER (per adjacent-block boundary, then the locus)

Phase support is locus-level evidence. No block-specific read phasing is invented.

| boundary state, first applicable | status | reason |
|---|---|---|
| no-call | NO_CALL | `RUN_NO_CALL` |
| an endpoint called homozygous | NOT_APPLICABLE | `NO_PHASE_INFORMATIVE_BOUNDARY` |
| genotype-equivalent call, or `phase_support` UNRESOLVED | NO_CALL | `PHASE_UNRESOLVED` |
| READS_AUDITED or READS_WITHIN_CANDIDATES | PASS | none |
| PRIOR_ONLY | LOW_QUALITY | `MODEL_ONLY_PHASE` |
| anything else | NOT_EVALUATED | `PHASE_ASSESSMENT_INCOMPLETE` |

An endpoint whose genotype an applied dosage repair replaced adds NOT_EVALUATED
`PHASE_NOT_ASSESSED_AFTER_DOSAGE_REPAIR`.

## 6. Where they are written

**`<prefix>.block_genotypes.tsv` and `<prefix>.genotypes.tsv`:** four trailing columns, `genotype_filter`,
`genotype_filter_reasons`, `dosage_filter` and `dosage_filter_reasons`, after every existing column.

**`<prefix>.locus_genotype.tsv`:** trailing columns.
- `headline_rule`
- `genotype_filter`, `genotype_filter_reasons`, `dosage_filter`, `dosage_filter_reasons`, `phase_filter` and
  `phase_filter_reasons`
- `call_summary`
- For each of `genotype`, `dosage` and `phase`: `<f>_filter_units`, `_applicable`, `_pass`, `_low_quality`,
  `_not_evaluated`, `_no_call` and `_not_applicable`.

The headline is computed from the run's own final outputs, after the dosage stage, in every mode.

**Schema changes:**
- **Renamed:** the dosage stage's own per-block filter columns in repair mode, formerly `dosage_filter` and
  `dosage_filter_reason`, are now `dosage_stage_filter` and `dosage_stage_filter_reason`. Their content is unchanged.
- **Appended:** the new columns are all trailing, so readers that select by name are unaffected. Readers that assume a
  fixed column count, or read the old repair-mode `dosage_filter` column, must be updated.

## 7. Usage policies (documented selections, not output filters)

    EXACT_GENOTYPE_USABLE  = GENOTYPE_FILTER == PASS and DOSAGE_FILTER in {PASS, NOT_APPLICABLE}
    DOSAGE_CONTENT_USABLE  = DOSAGE_FILTER == PASS
    PHASED_GENOTYPE_USABLE = EXACT_GENOTYPE_USABLE and PHASE_FILTER == PASS

Ordinary genotype analyses never require PHASE_FILTER = PASS.

**Measured on the sealed C4 replay** (`tests/results/c4_headline_filters_26401b4/REPORT_v2.md`):
- **EXACT_GENOTYPE_USABLE:** retention 0.889 and precision 0.858 over all experiments; 1.000 on leave-zero-out. It
  keeps 2 of the 15 catastrophic errors that GENOTYPE_FILTER alone keeps.
- **DOSAGE_CONTENT_USABLE:** content precision 0.997.
- **PHASED_GENOTYPE_USABLE:** almost empty on C4, because the caller's phase there is prior-only.

## 8. Known short-read limitation: allocation

**What it is.** At C4 sv:7, which homologue carries each HERV-bearing module (allocation) is not observable from short
reads. At `loo.NA19650` 10× and 30× the call is:
- dosage-CONSISTENT;
- representation PASS and read-supported;
- nevertheless a catastrophic allocation error.

**What that means here:**
- It passes EXACT_GENOTYPE_USABLE, and no existing truth-free field identifies it (`allocation_status` is NOT_TESTED,
  as for 600 sv:7 blocks).
- No truth-tuned exception is made. Allocation should be treated as unverified wherever `allocation_status` is not
  resolved.

## 9. Settings and the advanced CLI

**Scientific settings** are exposed as named flags (`panvar genotype --help-advanced`). Each is exactly its profile
key: the same parser and the same range validation.

**Precedence** is `flag > profile > preset > built-in default`, whatever the argument order, applied once.

**The preset is short-read by default.** An omitted `--preset` means `--preset short-read`, so the plain command runs
the production algorithm (adaptive repair, exact tie branching, the representation guard). An explicit
`--preset short-read` is equivalent; `--preset long-read` is refused as unimplemented. A profile still overrides the
preset, key by key. The unrepaired development view needs an explicit `--adaptive-repair off`.
`--phase-audit-candidates` is a development-only control: adaptive repair refuses it, because its pooled rerank
replaces the audit.

**Resource caps, caches and time and footprint limits** stay profile keys.

**Provenance:** `<prefix>.settings.tsv` records every public setting's effective value and its source (`default`,
`preset short-read (default)` when the preset was not given, `preset short-read` when it was, `profile <path>` or
`flag --…`). `mosaic_status.tsv` records `preset`, `preset_source` (`default` or `flag`) and `preset_keys`.

**Invalid values** are refused before any work, and no output is written.
