# M15-G: assessment_v1 implementation (frozen before the recorded run)

**Logistics:**
- **Branch:** `m15-g-assessment`, from `m15-f-filter-contract` `5ebe2e9`.
- **Contract and specification frozen:** `642b960` (`CONTRACT.md`, `assessment_spec.py`, fixtures S01–S20 and G1–G6).
- **Results:** `tests/results/m15g_assessment_20261022/`.
- **Budget:** 2,400 s, including 240 s of disclosed development.
- **Resources:** 4 threads, one process at a time.

## 1. Implementation (reporting only)

**The module.** `src/genotype_assessment.cpp` (+ header) is a line-for-line port of `assessment_spec.py`.
- **What it reads:** only the run prefix's own outputs; legacy block facts in verified context only.
- **What it writes** (idempotent):
  - locus rows in `.sequence_locus.tsv`;
  - block columns in `.sequence_blocks.tsv`;
  - mirrors in `.locus_genotype.tsv` (with `headline_filters_status superseded_by_assessment_v1`) and the legacy block
    tables when they exist.
- **What it never does:** score, recompute dosage, or change any existing column.

**Where it is written** (sequence_first_v2 runs only):
- **The command route:** after the locus table, inside the same try. A failure there is the existing exit-2 reporting
  failure, and the view is still attempted best effort.
- **Watchdog stops:** after the minimal status (NO_PREDICTION record or kept checkpoint).
- **The failure handler:** after its record.

**Tests:**
- **`panvar_genotype_assessment` (ctest):** C++ against the specification, byte for byte, on all 26 fixtures, plus the 9
  counterexample expectations.
- **The S1-D audit-vs-off e2e comparison:** now excludes exactly the assessment's dosage columns and rows. It requires the
  off run to record NOT_REQUESTED and every non-dosage assessment row to match.

## 2. Verification units (frozen expectations; `examples.py`)

| unit | expectation |
|---|---|
| `spec_equivalence` | 26 of 26 fixtures identical (C++ = specification), and 9 of 9 expectations |
| `ex.replay_saved` | every saved development v2 prefix (reserved excluded): C++ `show` = specification `render`; every source file read is byte-identical before and after |
| `ex.complete_{c4,cyp2d6,acot}` | exit 0; the written rows and columns equal the specification on the run's final files. **Unchanged from M15-E's run of the same command:** prediction key, sequences, candidate scores and candidates, equivalence set, intervals, alternatives, compound, FASTA, every pre-existing block column (old headline columns included), the locus table's pre-existing columns, the status's status / call class / phase support / no-call / library fields, and the sequence record's pre-existing rows (timing aside). Headline columns labelled superseded |
| `ex.dosage_cyp2d6` | declared library + `--dosage-check`: dosage UNAVAILABLE with `AUDIT_FAILED(<reason>)@locus`; prediction identical to M15-D.1's |
| `ex.legacy_tie` | the S1-B F1 fixture: v2 PREDICTION beside legacy `distinct_tie`, **no legacy block table**, primary locus rows and block columns present, unverifiable checks NOT_ASSESSED, equal to the specification |
| `ex.provisional` | checkpoint stop (exit 3): PROVISIONAL rows written by the stop path; no decision set and no completed support |
| `ex.absent` | stop before a checkpoint (exit 124): NO_PREDICTION rows, all NA |
| `ex.reporting_failure` | locus table blocked (exit 2, `locus_output FAILED`): the primary view of the published prediction still written, equal to the specification |
| `ctest.full` | the full enabled suite on the final revision passes; disabled tests reported separately |

**Failure handling:** a failure is recorded and never retuned. An amendment is declared before any rerun.
