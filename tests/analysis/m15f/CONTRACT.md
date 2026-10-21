# M15-F: proposed filter contract (`assessment_v1`), for approval before implementation

**Status:** a proposal.
- **Not done:** nothing in the production reporting is changed; no model, prediction, confidence or threshold change; no
  campaign.
- **The executable specification:** `contract.py`. Its `proposed()` is the contract below, and its `current()` is a port
  of today's `headline_filters_v2`, verified against every recorded headline (§6).
- **The field inventory:** `INVENTORY.md`.

## 1. What the current view gets wrong (with evidence)

| # | problem | evidence |
|---|---|---|
| P1 | **The legacy call class overrides the v2 prediction.** A decision that ties only over labellings is a v2 PREDICTION with exit 0, but the headline reads `GENOTYPE_NO_CALL;DOSAGE_NO_CALL;PHASE_NO_CALL`. No legacy block table is written for it | 119 saved development runs: 78 PREDICTION (66 with a RESOLVED sequence pair) and 41 PROVISIONAL, all NO_CALL in all three filters, none with a block table |
| P2 | **Missing evidence outranks concrete findings.** The precedence NO_CALL > NOT_EVALUATED > LOW_QUALITY lets one underpowered block hide a disagreement | ACOT 5× (S2-J): `DOSAGE_SUSPECT@sv:13` under a locus DOSAGE_FILTER of NOT_EVALUATED. C4 NA19434: a dosage-inconsistent sv:7 under NOT_EVALUATED (because sv:1 is underpowered) |
| P3 | **DOSAGE_FILTER PASS describes the alternative, not the reported alleles.** `COMPLETE_REPAIR_AVAILABLE_NOT_APPLIED` and `CONTENT_CORRECTED_EXACT_DIPLOTYPE_UNRESOLVED` are PASS, while the same fact appears as LOW_QUALITY `EXACT_DIPLOTYPE_UNRESOLVED` in GENOTYPE_FILTER: one fact, two places, opposite signs | C4 NA19434 and NA24385 (M11-D): sv:7 dosage PASS and genotype LOW_QUALITY |
| P4 | **Prior dependence is invisible in the headline.** `genotype_support PRIOR_ONLY` reaches only `quality_filter` | PRIOR_ONLY blocks in 7 of the 8 saved predictions in §6 that have block tables; one of the two largest reserved ACOT errors (M15-B, NA19468 sv:9) was PRIOR_ONLY with GENOTYPE_FILTER PASS |
| P5 | **Phase looks local but is inherited.** `phase_support` is computed once per locus and copied to every boundary; only zygosity and an applied repair vary per boundary | GR:1562–1566, GC:7647, HL:102–114 |
| P6 | **"Not evaluated" is folded into adverse or into PASS in component filters** | `quality_filter` fails closed on `representation_not_evaluated`; `dosage_stage_filter` maps NOT_SENSITIVE to PASS and an unrun repair to TRIGGERED_UNREPAIRED |
| P7 | **Dosage diagnostics under-report under audit** | `dosage_complete_not_applied_blocks` reads 0 when complete repairs exist; `locus_repair_scope` can read COMPLETE with nothing applied; TRIGGERED_UNREPAIRED hides "repair could not run" (cap, resource, reproduction) except in `dosage_reason` |
| P8 | **The docs say otherwise** | guide §6 lists `distinct_tie` as NO_PREDICTION; the algorithm page says `phase_support` is "PRIOR_ONLY otherwise" |

## 2. The proposed primary view

**Principles:**
1. **The primary outcome is the v2 sequence report's.**
2. **Every assessment concerns the reported prediction.** Corrections and alternatives are listed separately.
3. **"Not assessed" is never a quality failure,** and a finding is never hidden by a missing check: each assessment states
   what the checks that ran found, and every check that could not run is listed.
4. **There is no overall PASS and no single severity.**

### 2.1 Locus level

In `.sequence_locus.tsv`, which every v2 run writes, mirrored at the front of `.locus_genotype.tsv`.

| column | values | rule |
|---|---|---|
| `assessment_rule` | `assessment_v1` | |
| `prediction` | PREDICTION · PROVISIONAL · NO_PREDICTION | v2, unchanged. The legacy `call_class` never overrides it |
| `prediction_detail` | `completed_ranking`, the stop name, or the no-prediction reason | v2 `completion` |
| `sequence_resolution` | RESOLVED · GENOTYPE_RESOLVED_SEQUENCE_SET · SEQUENCE_TIED · NOT_DETERMINED · NA | v2 `sequence_status`, renamed. Sequence ambiguity, distinct from labels |
| `genotype_check` | NONE_DETECTED · CONFLICT · NOT_ASSESSED · NA | CONFLICT if any block has one; else NONE_DETECTED if any block's representation guard ran; else NOT_ASSESSED |
| `dosage_check` | AGREES · SUSPECT · DISAGREES · UNAVAILABLE · NOT_REQUESTED · NOT_APPLICABLE · NA | the first present in DISAGREES > SUSPECT > AGREES > UNAVAILABLE > NOT_REQUESTED > NOT_APPLICABLE. **Concrete evidence first;** unavailability never masks it |
| `read_support` | READS · PRIOR_DEPENDENT · NOT_ASSESSED · NA | PRIOR_DEPENDENT if any block is PRIOR_ONLY |
| `phase` (+ `phase_reason`) | READ_SUPPORTED · PRIOR_ONLY · UNRESOLVED · NOT_ASSESSED · NOT_APPLICABLE · NA | **locus-level only** (phase support is a property of the decision, §3.3) |
| `findings` | `CODE@ids;…` (every concrete finding, nothing collapsed) | codes: REPRESENTATION_CONFLICT, LIBRARY_SENSITIVE, DOSAGE_DISAGREES, DOSAGE_SUSPECT, PRIOR_DEPENDENT, LABELS_AMBIGUOUS, GENOTYPE_UNRESOLVED, PHASE_PRIOR_ONLY, PHASE_UNRESOLVED |
| `not_assessed` | `CHECK:why@ids;…` | REPRESENTATION_GUARD:<why>, LIBRARY_SENSITIVITY:<why>, DOSAGE:UNDERPOWERED, DOSAGE:NOT_EVALUATED:<why>, DOSAGE:AUDIT_FAILED:<note>, DOSAGE:AUDIT_NOT_RUN, DOSAGE:AUDIT_UNSUPPORTED_LIBRARY, BLOCK_CHECKS:<why>, PHASE:<why> |
| `dosage_alternatives` | `ALT@ids;…` or `none` (**diagnostic**) | COMPLETE_REPAIR_AVAILABLE, COMPLETE_REPAIR_REFUSED:<guard>, CONTENT_CORRECTION_AVAILABLE, REPAIR_NOT_FOUND:<reason> (this includes caps and resource refusals), APPLIED (legacy repair mode) |
| `summary` | `<prediction>;SEQUENCE_<v>;GENOTYPE_<v>;DOSAGE_<v>;SUPPORT_<v>;PHASE_<v>` | replaces `call_summary` as the scan line |

### 2.2 Block level

In `.sequence_blocks.tsv`, written for every v2 prediction, and mirrored into `.block_genotypes.tsv` / `.genotypes.tsv`
when they exist.

| column | values | rule |
|---|---|---|
| `block_resolution` | IDENTIFIED · LABELS_AMBIGUOUS · GENOTYPE_UNRESOLVED · PROVISIONAL_LABEL | v2 `block_status`. `REPRESENTATION_AMBIGUOUS` is renamed LABELS_AMBIGUOUS: the same sequence pair with different block labels |
| `genotype_check` (+ `_reasons`) | NONE_DETECTED · CONFLICT · NOT_ASSESSED | CONFLICT: LOWQUAL_REPRESENTATION_CONFLICT, or library SENSITIVE at this block. NOT_ASSESSED only when the representation guard did not run. A missing library audit is listed in `not_assessed`, never a value. **Dosage is never folded in** |
| `dosage_check` (+ `dosage_check_reason`) | as at the locus | §3.1 |
| `dosage_alternative` | as at the locus (one value) | diagnostic |
| `read_support` | READS · PRIOR_DEPENDENT · NOT_ASSESSED | `genotype_support` |

## 3. Rules

### 3.1 Dosage: the reported alleles first, the alternative separately

| facts | `dosage_check` | reason | `dosage_alternative` |
|---|---|---|---|
| no prediction | NA | | NA |
| `--dosage-check` not given | NOT_REQUESTED | dosage_check_not_requested | NA |
| audit FAILED / NOT_RUN / REFUSED_LIBRARY | UNAVAILABLE | AUDIT_FAILED:<note> / AUDIT_NOT_RUN / AUDIT_UNSUPPORTED_LIBRARY (recorded once, at the locus) | NA |
| NOT_SENSITIVE | NOT_APPLICABLE | NOT_DOSAGE_SENSITIVE | NONE |
| NOT_EVALUATED (block) | UNAVAILABLE | NOT_EVALUATED:<dosage_reason> | NONE |
| UNDERPOWERED | UNAVAILABLE | UNDERPOWERED | NONE |
| CONSISTENT | AGREES | CONSISTENT | NONE |
| SUSPECT | SUSPECT | SUSPECT | NONE |
| TRIGGERED_UNREPAIRED (any reason, including caps and resources) | DISAGREES | TRIGGERED | REPAIR_NOT_FOUND:<dosage_reason> |
| REPAIRED, not applied, eligible (audit) | DISAGREES | TRIGGERED | COMPLETE_REPAIR_AVAILABLE |
| REPAIRED, not applied, refused by a guard | DISAGREES | TRIGGERED | COMPLETE_REPAIR_REFUSED:<guard> |
| CONTENT_REPAIRED | DISAGREES | TRIGGERED | CONTENT_CORRECTION_AVAILABLE |
| REPAIRED and applied (legacy `--dosage-mode repair` only) | AGREES | REPORTED_IS_DOSAGE_REPAIR | APPLIED |

**Why REPAIRED-not-applied is DISAGREES:** a triggered block is one whose reported alleles are dosage-inconsistent. That
an alternative exists does not change what the check found about the reported prediction.

### 3.2 Genotype conflicts, read support, ambiguity

- **Kept apart:**
  - detected conflicts: representation, library sensitivity;
  - prior dependence: a characteristic of the evidence, not a conflict;
  - ambiguity: `sequence_resolution` and `block_resolution`.

  None of these implies another.
- **The library sensitivity audit runs only with `--library auto`.** With a declared library it is not applicable, and it
  is not listed.
- **Prior dependence** is reported as a finding (PRIOR_DEPENDENT@ids), not as a quality failure.

### 3.3 Phase is locus-level, by construction

`phase_support` is decided once per locus, from whether a fragment contender projects differently. So `phase` is reported
**only at the locus**, labelled as such. Per-boundary phase assessments are not offered until they are computed per
boundary.

| facts | `phase` | `phase_reason` |
|---|---|---|
| no prediction | NA | no_prediction |
| PROVISIONAL | NOT_ASSESSED | provisional_checkpoint |
| several sequence pairs (GENOTYPE_RESOLVED_SEQUENCE_SET or SEQUENCE_TIED, or genotype-equivalent) | UNRESOLVED | several_sequence_pairs |
| no adjacent HET/HET blocks | NOT_APPLICABLE | no_phase_informative_boundary |
| READS_WITHIN_CANDIDATES / READS_AUDITED | READ_SUPPORTED | |
| PRIOR_ONLY | PRIOR_ONLY | prior_only |
| RESOLVED sequence pair, but the legacy decision tied over its labellings (`phase_support` UNRESOLVED) | NOT_ASSESSED | label_tie_support_not_computed |
| a dosage repair applied (legacy repair mode) | NOT_ASSESSED | after_dosage_repair |

**Why the RESOLVED pair is NOT_ASSESSED rather than UNRESOLVED:** one sequence pair fixes phase, because the haplotypes
are whole-locus paths. What is missing is its read support, not the phase itself.

### 3.4 The decision table (the brief's cases; `synthetic.tsv`)

| case | current `call_summary` | proposed `summary` | proposed findings · not assessed · alternatives |
|---|---|---|---|
| no dosage requested | GENOTYPE_PASS;DOSAGE_NOT_EVALUATED;PHASE_PASS | PREDICTION;…;DOSAGE_NOT_REQUESTED;… | none · none · none |
| requested audit failed | …DOSAGE_NOT_EVALUATED… | …DOSAGE_UNAVAILABLE… | none · DOSAGE:AUDIT_FAILED:<note>@locus · none |
| underpowered at one block | …DOSAGE_NOT_EVALUATED… | …DOSAGE_AGREES… | none · DOSAGE:UNDERPOWERED@sv:3 · none |
| dosage consistent | …DOSAGE_PASS… | …DOSAGE_AGREES… | none |
| complete correction available, not applied | GENOTYPE_LOW_QUALITY;DOSAGE_PASS | GENOTYPE_NONE_DETECTED;DOSAGE_DISAGREES | DOSAGE_DISAGREES@sv:3 · none · COMPLETE_REPAIR_AVAILABLE@sv:3 |
| correction refused by a guard | GENOTYPE_LOW_QUALITY;DOSAGE_LOW_QUALITY | GENOTYPE_NONE_DETECTED;DOSAGE_DISAGREES | DOSAGE_DISAGREES@sv:3 · COMPLETE_REPAIR_REFUSED:<guard>@sv:3 |
| conditional content correction | GENOTYPE_LOW_QUALITY;DOSAGE_PASS | GENOTYPE_NONE_DETECTED;DOSAGE_DISAGREES | DOSAGE_DISAGREES@sv:3 · CONTENT_CORRECTION_AVAILABLE@sv:3 |
| triggered, repair not run (cap) | GENOTYPE_LOW_QUALITY;DOSAGE_LOW_QUALITY | …DOSAGE_DISAGREES | REPAIR_NOT_FOUND:cap_refused@sv:3 |
| concrete conflict + library audit incomplete | **GENOTYPE_NOT_EVALUATED** | GENOTYPE_CONFLICT | REPRESENTATION_CONFLICT@sv:1 · LIBRARY_SENSITIVITY:primary_no_call@sv:1,sv:3 |
| concrete conflict + dosage audit failed | GENOTYPE_LOW_QUALITY;DOSAGE_NOT_EVALUATED | GENOTYPE_CONFLICT;DOSAGE_UNAVAILABLE | REPRESENTATION_CONFLICT@sv:1 · DOSAGE:AUDIT_FAILED@locus |
| prior-dependent, otherwise unflagged | GENOTYPE_PASS;DOSAGE_PASS;PHASE_LOW_QUALITY | …SUPPORT_PRIOR_DEPENDENT;PHASE_PRIOR_ONLY | PRIOR_DEPENDENT@sv:3;PHASE_PRIOR_ONLY@locus |
| one sequence pair, ambiguous labels (legacy tie) | **NO_CALL ×3** | PREDICTION;SEQUENCE_RESOLVED;GENOTYPE_NOT_ASSESSED;DOSAGE_UNAVAILABLE;SUPPORT_NOT_ASSESSED;PHASE_NOT_ASSESSED | LABELS_AMBIGUOUS@sv:3 · BLOCK_CHECKS:legacy_decision_distinct_tie; DOSAGE:AUDIT_NOT_RUN; PHASE:label_tie_support_not_computed |
| genotype available, phase unresolved | …PHASE_NO_CALL | PREDICTION;SEQUENCE_GENOTYPE_RESOLVED_SEQUENCE_SET;…;PHASE_UNRESOLVED | PHASE_UNRESOLVED@locus |
| provisional / resource-stopped after a checkpoint | **NO_CALL ×3** | PROVISIONAL;SEQUENCE_NOT_DETERMINED;GENOTYPE_NOT_ASSESSED;…;PHASE_NOT_ASSESSED | BLOCK_CHECKS:provisional_checkpoint; PHASE:provisional_checkpoint |
| stop before a checkpoint | NO_CALL ×3 | NO_PREDICTION;SEQUENCE_NA;GENOTYPE_NA;DOSAGE_NA;SUPPORT_NA;PHASE_NA | none |
| one problematic block among ten | GENOTYPE_PASS;DOSAGE_LOW_QUALITY | …DOSAGE_SUSPECT… | DOSAGE_SUSPECT@b8 (the locus value plus the id list; nothing about the other nine is lost) |
| library sensitive at one block | GENOTYPE_LOW_QUALITY | GENOTYPE_CONFLICT | LIBRARY_SENSITIVE@sv:1 |
| legacy repair applied | …PHASE_NOT_EVALUATED | DOSAGE_AGREES;PHASE_NOT_ASSESSED | APPLIED@sv:3 · PHASE:after_dosage_repair |

The **bold** entries are where the current headline is wrong or masks a finding.

## 4. Compatibility treatment for the old fields

| field | treatment |
|---|---|
| `genotype_filter`, `dosage_filter`, `phase_filter` (+ reasons and counts), `call_summary`, `headline_rule` | **kept unchanged** for one release as `headline_filters_v2` compatibility columns, documented as superseded by `assessment_v1`, with the P1–P3 caveats stated next to them; then removed from the default output (kept behind a flag) |
| `filter`, `quality_filter`, `locus_filter`, `quality_*_blocks`, `genotype_support`, `phase_support`, `library_sensitivity*` | kept as **component diagnostics** (the inputs of the primary view), unchanged. `quality_filter`'s fail-closed fold is documented (P6) |
| `status`, `call_class`, `no_call_reason` | kept as the **legacy decision**, documented as never overriding `prediction` |
| dosage tables and columns | unchanged as diagnostics. `dosage_complete_not_applied_*` and `locus_repair_scope` gain the audit-mode record (P7); `gq_status` is documented as describing the alternative |
| `REPRESENTATION_AMBIGUOUS` (`block_status`) | kept. `block_resolution` carries the clearer name LABELS_AMBIGUOUS |
| docs (P8) | the guide's §6 example and the algorithm page's `phase_support` sentence corrected with the implementation |

## 5. Example user-facing output

`--library FILE --dosage-check`, C4 NA19434 (development; M11-D's outputs read, unchanged):

```
prediction           PREDICTION (completed_ranking)
sequence_resolution  RESOLVED
genotype_check       NONE_DETECTED
dosage_check         DISAGREES          DOSAGE_DISAGREES@sv:7
                                        alternative: CONTENT_CORRECTION_AVAILABLE@sv:7 (diagnostic; not applied)
read_support         PRIOR_DEPENDENT    PRIOR_DEPENDENT@ctx:0,sv:1,ctx:2,sv:3,ctx:4,sv:5,sv:7,ctx:8,sv:9
phase                PRIOR_ONLY         (locus-level)
not_assessed         DOSAGE:UNDERPOWERED@sv:1
summary              PREDICTION;SEQUENCE_RESOLVED;GENOTYPE_NONE_DETECTED;DOSAGE_DISAGREES;SUPPORT_PRIOR_DEPENDENT;PHASE_PRIOR_ONLY
```

**Today the same run reads** `GENOTYPE_LOW_QUALITY;DOSAGE_NOT_EVALUATED;PHASE_LOW_QUALITY`: the disagreement is hidden by
the underpowered sv:1, and sv:7's dosage reads PASS.

CYP2D6 HG02769 with `--dosage-check`, the known audit failure (M15-D.1):

```
prediction           PREDICTION (completed_ranking)
genotype_check       NONE_DETECTED
dosage_check         UNAVAILABLE        DOSAGE:AUDIT_FAILED:unit library: alleles disagree on their flanks or anchor trims in sv:9
read_support         PRIOR_DEPENDENT    PRIOR_DEPENDENT@ctx:0,ctx:2,sv:3,ctx:4,ctx:14,ctx:16
phase                PRIOR_ONLY         (locus-level)
```

A CYP2D6 30× run whose legacy decision tied over labellings (S2-L control):

```
prediction           PREDICTION (completed_ranking)
sequence_resolution  RESOLVED           LABELS_AMBIGUOUS@sv:3,sv:5 (one sequence pair; block labels differ)
genotype_check       NOT_ASSESSED       BLOCK_CHECKS:legacy_decision_distinct_tie
dosage_check         UNAVAILABLE        DOSAGE:AUDIT_NOT_RUN
phase                NOT_ASSESSED       label_tie_support_not_computed
```

Today this run reads `GENOTYPE_NO_CALL;DOSAGE_NO_CALL;PHASE_NO_CALL`.

## 6. Demonstration (`tests/results/m15f_filter_contract_20261021/`)

- **`fidelity` (with A1):** the Python port of `headline_filters_v2` reproduces **1,382 of 1,382** recorded headlines (all
  seven locus fields).
  - One run was recorded under an earlier rule (`a1e3ffd`, before `02f846a` added the `:<guard>` suffix).
  - 909 pre-headline runs have no headline to compare.
  - **The first run's 155 differences were only the `@ids` list of no-calls:** production passes the model's block ids in
    memory, and no file lists them. A1 takes them from the recorded list, as an input fact.
- **`synthetic`:** the 20-case decision table, current against proposed (`synthetic.tsv`).
- **`saved`:** 11 representative development runs, current against proposed (`saved.tsv`), including the CYP2D6 audit
  failure, two unapplied content corrections, a v2 prediction beside a legacy tie, PROVISIONAL, a representation conflict,
  a genotype-equivalent call and a suspect block. **250 of 250 source files byte-identical** before and after.
- **Budget:** about 76 of 600 s: 58.8 s recorded in the ledger (including 30 s of disclosed development), plus about 17 s of
  verification scans after the run (the 119-run block-table check and the guide line), disclosed here.

## 7. Decisions for Davide

1. **Approve `assessment_v1`** as the primary view (§2–§3), or amend it.
2. **v2 predictions whose legacy decision tied (P1)** get `NOT_ASSESSED` block checks under this contract. Computing the
   representation guard, read support and dosage on the v2 representative is separate implementation work. Proceed now,
   or later?
3. **The compatibility period** for the `headline_filters_v2` columns: one release, then behind a flag?
4. **Names:** `genotype_check` / `dosage_check` / `read_support` / `phase`, and the `findings` / `not_assessed` lists. Are
   these acceptable?

**Not proposed:** `--min-GQ`, a "reliable call" class, automatic dropping, any truth-dependent filter, or any calibrated
confidence claim.
