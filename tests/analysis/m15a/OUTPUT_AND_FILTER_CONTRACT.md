# Output and filter contract (proposed; for approval)

**Scope:** the release-candidate configuration (`RELEASE_CANDIDATE.md`) with `output_schema sequence_first_v2`. The
user-facing definitions remain in `docs/modules/genotype.md`. This contract fixes **which fields answer which question**, and
proposes the filter procedure.

**Nothing here is calibrated confidence.** GQ, error mass and the conditional posterior are conditional on the explored
candidates, and are uncalibrated. QVfrac is a descriptive evaluation quantity, never a reported score.

## 1. Six questions, six places

| question | where it is answered | values |
|---|---|---|
| **Prediction** | `.sequence_locus.tsv` `prediction`; `.sequence_haplotypes.tsv` (SHA-256, length, `origin` PANEL/OFF_PANEL, every sequence-equivalent panel path, extraction reference, GFA SHA-256, model fingerprint, block path); `.sequence_haplotypes.fa` only for an off-panel pair; legacy `.locus_genotype.tsv` call class (strict / genotype-equivalent / no-call) | `PREDICTION` · `PROVISIONAL` · `NO_PREDICTION` |
| **Completion / resource status** | `run_status`, `completion`, exit code, `.mosaic_status.tsv` caps and limits | completed ranking; `PROVISIONAL` (an early-incumbent checkpoint or declared cap; exit 3); no-call (exit 1); time/footprint limit (exit 124/125, minimal status); errors (exit 2) |
| **Sequence and genotype ambiguity** | `sequence_status`; per-block `block_status`; `.representation_intervals.tsv`; `.sequence_alternatives.tsv`; legacy `CONDITIONAL_GENOTYPE_EQUIVALENT` | `RESOLVED` · `GENOTYPE_RESOLVED_SEQUENCE_SET` · `SEQUENCE_TIED`; blocks `IDENTIFIED` · `REPRESENTATION_AMBIGUOUS` · `GENOTYPE_UNRESOLVED` (over **explored** labellings only) |
| **Phase** | `phase_support`; PHASE_FILTER; FASTA header `phase_support` / `representation=` | read-supported · `MODEL_ONLY_PHASE` (prior-only) · unresolved · not applicable |
| **Dosage-audit disagreement** | `.dosage_status.tsv` (audit COMPLETE / FAILED / NOT_RUN with note); `.dosage_blocks.tsv` `dosage_status`; DOSAGE_FILTER | CONSISTENT · CONTENT_REPAIRED · TRIGGERED_UNREPAIRED · SUSPECT · UNDERPOWERED · NOT_SENSITIVE · NOT_EVALUATED |
| **Pool-conditional support** | per block `genotype_support`; `error_mass_lower/upper`; `legacy_representative_GQ`; `candidate_scope CONDITIONAL`; `structure_confidence NA` | READS / PRIOR_ONLY: PRIOR_ONLY means the reads alone do not separate the call from at least one **explored** rival (pool- and labelling-conditional; M11-C) |

**Contract rules:**
- **The dosage audit never replaces or reconstructs haplotypes.**
  - Under v2, `--dosage-mode repair` is refused before the output prefix is touched.
  - Audit representatives go only to dosage diagnostics, labelled `REPRESENTATIVE_ONLY_NOT_RESOLVED`.
  - A FAILED audit (for example, the CYP2D6 sv:9 unit library, or ACOT `KMAX > 3`) is reported as **not evaluated**,
    never as consistent.
- **Panel predictions** carry identifiers sufficient to recover the sequence from the graph:
  - the path name, interval and orientation, plus the GFA SHA-256 and the model fingerprint;
  - checked by `tests/analysis/s1d/extract_panel_haplotype.py --check`;
  - **no new FASTA is required for evaluation.**
- **Off-panel mosaic outputs** follow the v2 FASTA rule: one FASTA holding both records, annotated with their origin.
- **Predictions stay accessible when quality is uncertain.** An uncertain completed prediction is still `PREDICTION`
  (exit 0), with annotations.
- **Incomplete or invalid computations never masquerade as completed predictions:**
  - provisional results are `PROVISIONAL`;
  - watchdog stops before a checkpoint, refusals and invalid input give `NO_PREDICTION`;
  - the run-owned reports are removed on a later run or a watchdog stop.

## 2. Proposed filter procedure (no new classifier, no new score)

### A. Computational validity checks (they decide whether a completed prediction exists)

A run yields a **completed prediction** only if all of the following hold:
1. exit 0 and `prediction = PREDICTION`;
2. `completion = completed_ranking` and not provisional;
3. the representation guard is complete, and `sequence_publication` did not fail;
4. the model fingerprint is recorded, and every panel reference passes the extraction check.

**Anything else is reported by its own status** — PROVISIONAL, NO_PREDICTION with its reason, or a resource stop. It is
never filtered into or out of a quality class.

### B. Uncertainty and disagreement annotations (annotate; never hide a prediction)

| annotation | trigger | development basis |
|---|---|---|
| `SEQUENCE_AMBIGUITY` | `sequence_status` ≠ RESOLVED; any block not IDENTIFIED | definitional |
| GENOTYPE / DOSAGE / PHASE headline filters | the existing v2 filters (LOW_QUALITY / NOT_EVALUATED / NO_CALL) | existing contract (selected on the C4 campaign; **development**) |
| `PRIOR_ONLY` block support | `genotype_support = PRIOR_ONLY` | M11-B/C/D: associated with non-top outcomes (22 of 27 non-top against 228 of 427 top), **not** a calibrated risk |
| `DOSAGE_AUDIT_DISAGREES` / `DOSAGE_AUDIT_UNAVAILABLE` | block `CONTENT_REPAIRED` / `TRIGGERED_UNREPAIRED` / `SUSPECT`; audit FAILED | M11-E/M12-A: where the audit ran it moved toward the truth's direction in 2 of 2 catastrophic cases; it FAILED in 3 of 5 |
| **proposed:** `REPEATED_UNIT_STRUCTURE` | census route `repeated_unit` for the block | M11-D: **all 5** catastrophic blocks were `repeated_unit` (5 of 29 such blocks; 0 of 425 others). **Inspected development evidence; an annotation only, pending approval** |

### C. Quality threshold: **none is proposed**

**Existing evidence cannot support a reliable-call category:**
- **M11-D:** per-block error mass (C1), and error mass plus PRIOR_ONLY (C2), did **not** beat a constant baseline out of
  family. The pooled Δ log-loss intervals include 0.
- **The catastrophic blocks carried high confidence:** GQ up to 443, with predicted top probability ≥ 0.91 out of family.
- **M12/M13/M14** produced no reliable within-state or variation-tolerant discriminator.

**Proposal: an experimental release with explicit uncalibrated support.**
- Every completed prediction is reported with the B annotations and the uncalibrated support fields.
- There is no RELIABLE / PASS-quality class, and no threshold is invented.

## 3. Decisions requiring Davide's approval

1. **Acceptable error:** the per-block excess and normalized excess that count as acceptable for reporting. Evaluation
   only, not a filter.
2. **Catastrophic-error handling** (excess ≥ 1,000 bp): accept as a reported limitation; annotate repeated-unit blocks
   (proposal B); or withhold some category (no evidence supports a safe withholding rule).
3. **Acceptable error among "reliable" calls:** needed only if a reliable class is created. The proposal is that none is.
4. **Whether `REPEATED_UNIT_STRUCTURE` and `DOSAGE_AUDIT_*` become formal annotations,** and whether either should lower
   GENOTYPE_FILTER, which today does so only for some dosage states.
5. **Reporting of all-empty and representation-dependent blocks:** M11-E eligibility (counted / interval-level /
   unresolved).
6. **RC packaging:** a named profile against default promotion (`RELEASE_CANDIDATE.md` §4).
