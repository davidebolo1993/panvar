# M15-G result: the corrected assessment_v1 is the primary v2 view; the C++ agrees byte for byte with the executable contract on 26 fixtures and 606 saved prefixes; predictions, scores, sequences, phase inference and the old headline fields are unchanged

**Run.**
- **Branch:** `m15-g-assessment`, from `5ebe2e9`.
- **Commits:**
  - contract and specification frozen `642b960`;
  - implementation and harness frozen `6c3b45f`;
  - amendments A1 `99cfa44`, A2 `ee80159` and A3 `961a7d8` (harness only);
  - docs `e2463ed`;
  - result after.
- **Results:** `tests/results/m15g_assessment_20261022/`.
- **Budget:** 1,761.3 of 2,400 s, including 240 s of disclosed development.
- **Units:**
  - **First pass:** 11 of 12 passed. `ex.complete_c4` failed on a harness name-normalisation defect; the outputs were
    identical (§4).
  - **Final:** passed after the declared reruns.
- **Not done:** no model, threshold, dropping, calibrated claim, index change, advanced-help reorganisation, reserved
  rerun, merge, push or tag.

**Binaries** (source at `6c3b45f`; A1–A3 and the docs changed no C++ and no registered test):

| binary | SHA-256 |
|---|---|
| `panvar` | `07ddd909907ad7c8f320d79e2d87f3c32ce1582ee9335708b1ac96e2b8b0b7d6` |
| `genotype_assessment_test` | `ef6edae4ed4dd5aeefe39e82baed9384c3415c1431dda555438cdfa379589ffc` |

## 1. The amended contract (`CONTRACT.md`, normative `assessment_spec.py`; frozen first)

| correction | fixture | the specification now gives |
|---|---|---|
| C1: the legacy classification never overrides v2 resolution | G1: a unique v2 pair beside a legacy genotype-equivalent call | RESOLVED; phase NOT_ASSESSED (`legacy_decision_genotype_equivalent`), not UNRESOLVED |
| C2: phase applicability is two or more HET blocks anywhere | G2: HET · HOM · HET | phase READ_SUPPORTED (the old headline: NOT_APPLICABLE) |
| C3: unknown support is never READS | G3: READS, unknown, unverified | `read_support` READS with coverage `1/3`, `SUPPORT_READS[1/3]`, both gaps listed |
| C4: provisional yields no decision set or completed support | G4: PROVISIONAL with legacy rows carrying PRIOR_ONLY and SUSPECT | every check NOT_ASSESSED, no findings |
| C5: findings and missing assessments both visible; coverage | G5: AGREES beside UNDERPOWERED and NOT_EVALUATED | `DOSAGE_AGREES[1/3]` and both listed |
| verified context | G6: a legacy block whose alleles differ from the published genotype | its conflict is not reported; NOT_ASSESSED `legacy_block_differs` |

**Legacy block facts are used only in verified context:**
- a completed v2 PREDICTION;
- an IDENTIFIED block;
- an existing legacy decision;
- sorted legacy alleles equal to the v2 genotype.

**Otherwise NOT_ASSESSED, with a reason.** No fact is borrowed from another candidate or from stale output, and nothing is
scored or recomputed.

## 2. Implementation (reporting only)

**The module.** `src/genotype_assessment.cpp` is a line-for-line port of the specification. It writes:
- **locus:** rows in `.sequence_locus.tsv`, mirrored as columns in `.locus_genotype.tsv` with
  `headline_filters_status superseded_by_assessment_v1`;
- **block:** columns in `.sequence_blocks.tsv`, mirrored in `.block_genotypes.tsv` and `.genotypes.tsv` when they exist.

**When it is written:**
- **The command route:** after the locus table. On a reporting failure it is still written, best effort, beside the exit-2
  record.
- **Watchdog stops** (PROVISIONAL checkpoint or NO_PREDICTION).
- **The failure handler.**

**Old headline fields:** unchanged, with their historical meanings, for one experimental release. They are documented as
superseded (guide §7.0, module page).

## 3. Verification

| unit | result |
|---|---|
| `spec_equivalence` (ctest `panvar_genotype_assessment`) | **26 of 26 fixtures byte-identical** (C++ = specification), and 9 of 9 counterexample expectations |
| `ex.replay_saved` | **606 of 606** saved development v2 prefixes byte-identical (438 PREDICTION, 120 PROVISIONAL, 48 NO_PREDICTION); **3,636 source files unchanged** (read-only) |
| `a3.ex.complete_{c4,cyp2d6,acot}` | **25/25 each.** The written rows and columns equal the specification. **Unchanged from M15-E's run of the same command:** prediction key, sequences, candidate scores and candidates, equivalence set, intervals, alternatives, compound, FASTA, every pre-existing block column (old headline columns included), the locus table, status / call class / **phase support** / library fields, and the sequence record. The old fields are labelled superseded |
| `ex.dosage_cyp2d6` | 6/6: dosage UNAVAILABLE, `DOSAGE:AUDIT_FAILED(unit library: alleles disagree on their flanks or anchor trims in sv:9)@locus`; prediction identical to M15-D.1's |
| `ex.legacy_tie` | 7/7: v2 PREDICTION beside legacy `distinct_tie`, **no legacy block table**, primary locus rows and block columns present; resolution RESOLVED, LABELS_AMBIGUOUS@sv:1,sv:2; the unverifiable checks NOT_ASSESSED with reasons |
| `ex.provisional` | 7/7: exit 3; the stop path writes PROVISIONAL rows; NOT_DETERMINED, no findings, every block PROVISIONAL_LABEL with support NOT_ASSESSED |
| `ex.absent` | 6/6: exit 124; NO_PREDICTION rows, all NA |
| `ex.reporting_failure` | 5/5: exit 2, `locus_output FAILED: …`; the primary view of the published prediction still written and equal to the specification |
| `ctest.full` (final C++ revision) | **64 of 64 passed** (1,009.6 s), including the new `panvar_genotype_assessment`; disabled, reported separately: `panvar_real_regressions`, `panvar_rebuild_stats` |

**Test change, declared in the freeze:** the S1-D e2e "audit vs off" comparison now excludes exactly the assessment's
dosage columns and rows. It also requires the off run to record NOT_REQUESTED, and every non-dosage assessment row to
match.

## 4. Amendments (harness only; every failed unit stays recorded)

**A1.**
- **The defect:** the completed-run comparison normalised the prefix name as a bare substring, and the C4 prefix `c4` is a
  hexadecimal substring of SHA-256 values.
- **What was actually true:** checked outside the harness, `.sequence_haplotypes.tsv` is byte-identical to M15-E's, the
  sequence record has no differing row, and the locus table differs only in its `prefix` (name) column.
- **The fix:** names are replaced only as names.

**A2.** The A1 units never ran: a runner dispatch error, 0.2 s each, recorded. It was fixed.

**A3.** A1 had left one call (the locus-table check) passing the bare basename. It was fixed, and all three units passed
25/25.

## 5. Example output (C4 NA19705, declared library, `--dosage-check`; `.sequence_locus.tsv`)

```
assessment_rule          assessment_v1
sequence_resolution      RESOLVED
genotype_check           NONE_DETECTED
genotype_check_coverage  11/11
dosage_check             AGREES
dosage_check_coverage    1/2
read_support             PRIOR_DEPENDENT
read_support_coverage    11/11
phase                    PRIOR_ONLY
phase_reason             prior_only
findings                 PRIOR_DEPENDENT@sv:1,ctx:2,ctx:4,sv:5,sv:7;PHASE_PRIOR_ONLY@locus
not_assessed             DOSAGE:UNDERPOWERED@sv:1
dosage_alternatives      none
assessment_summary       PREDICTION;SEQUENCE_RESOLVED;GENOTYPE_NONE_DETECTED;DOSAGE_AGREES[1/2];SUPPORT_PRIOR_DEPENDENT;PHASE_PRIOR_ONLY
```

**Block columns (`.sequence_blocks.tsv`):**

```
block_id  block_status  genotype  block_resolution  genotype_check  genotype_check_reasons  dosage_check    dosage_check_reason   dosage_alternative  read_support     read_support_reason
ctx:0     IDENTIFIED    21/48     IDENTIFIED        NONE_DETECTED   none                    NOT_APPLICABLE  NOT_DOSAGE_SENSITIVE  NONE                READS            reads
sv:1      IDENTIFIED    51/60     IDENTIFIED        NONE_DETECTED   none                    UNAVAILABLE     UNDERPOWERED          NONE                PRIOR_DEPENDENT  prior_only
```

**A legacy-tied v2 prediction** (S1-B F1): `PREDICTION;SEQUENCE_RESOLVED;GENOTYPE_NOT_ASSESSED;DOSAGE_NOT_REQUESTED;SUPPORT_NOT_ASSESSED;PHASE_NOT_ASSESSED`.
- `findings`: `LABELS_AMBIGUOUS@sv:1,sv:2`.
- `not_assessed`: `GENOTYPE_CHECK:legacy_decision_distinct_tie@ctx:0,ctx:3;…;PHASE:legacy_decision_distinct_tie@locus`.
- The old headline read `GENOTYPE_NO_CALL;DOSAGE_NO_CALL;PHASE_NO_CALL`.

## 6. Remaining limitations

- **Legacy-tied v2 predictions** report their block checks as NOT_ASSESSED. Computing them on the v2 representative needs
  new scoring and dosage work, which was not authorised.
- **Phase support is locus-level** by construction. Per-boundary phase is not offered.
- **The caller-only API** (`run_genotype` without the command layer) does not write the view. The command and
  prepared-locus routes do.
- **The old headline fields remain in the outputs** for this experimental release, with their known misreadings
  documented. They are not removed and not changed.
- **The basic `--help` text** still describes the headline filters; the docs now lead with `assessment_v1`. The help is
  untouched in this batch because its C++ string was frozen with the tested revision.

**Stopped for review.**
