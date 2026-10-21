# M15-B frozen output and filter contract (corrects M15-A's proposal; frozen before any reserved genotyping)

**Recorded decisions** (Sol's recommendations, approved for this batch):
- **Packaging:** a **named experimental short-read profile**, `genotype-dev-m10a-v1` (M10-A arm D). Defaults are not
  promoted.
- **Reliable-call category:** **none.** Existing `PASS` values in the headline filters (`GENOTYPE_FILTER`, `DOSAGE_FILTER`,
  `PHASE_FILTER`), `quality_filter` and `filter` are computational and representation checks. **They do not establish
  calibrated reliability.**
- **Repeated-unit and dosage annotations:** descriptive warnings, never rejection or withholding rules.
- **Catastrophe:** excess ≥ 1,000 bp, a **descriptive category only**. Every case is listed, and the predictions stay
  accessible.
- **Eligibility:**
  - all-empty blocks are visible but **excluded from informative accuracy summaries**;
  - representation-dependent regions are assessed at interval level where possible, and are otherwise UNRESOLVED.
- **Acceptable-error threshold:** **none is invented.** This campaign is a **descriptive reserved evaluation**: continuous
  excess distributions, and no accuracy PASS. A production acceptance decision is separate.
- **Neighbour limits:** the locus-specific `repair_max_neighbours` (CYP2D6 958, ACOT 1004; C4 default) are retained
  unchanged and documented as the supported per-locus profiles. No general rule is introduced.

## Correction to M15-A §2B: two separate ambiguity annotations

M15-A's `SEQUENCE_AMBIGUITY` conflated two different things. A representation-ambiguous block does not mean an ambiguous
sequence. The annotations are now:

| annotation | trigger (from existing outputs; nothing is changed) | meaning |
|---|---|---|
| `SEQUENCE_AMBIGUOUS` | `.sequence_locus.tsv` `sequence_status` ≠ `RESOLVED` (`GENOTYPE_RESOLVED_SEQUENCE_SET` or `SEQUENCE_TIED`) | more than one sequence pair is in the decision set |
| `BLOCK_REPRESENTATION_AMBIGUOUS` | a block's `block_status` = `REPRESENTATION_AMBIGUOUS` (explored labellings of the same sequence pair differ at that block) | the block labels are ambiguous; **the sequence may be unique** |
| `BLOCK_GENOTYPE_UNRESOLVED` | a block's `block_status` = `GENOTYPE_UNRESOLVED` | the decision set disagrees on that block's genotype |

The underlying outputs (`sequence_status`, `block_status`, `.representation_intervals.tsv`) are preserved unchanged.

## Descriptive warnings (formalised; none withholds a prediction)

| annotation | values | rule |
|---|---|---|
| `REPEATED_UNIT_STRUCTURE` | `YES` / `NO` / **`UNKNOWN`** | from `block_diagnostics` `census_route`: `repeated_unit` → YES; `not_sensitive` or `pair` → NO; any `NOT_EVALUATED:*` or `UNAVAILABLE:*` → **UNKNOWN**. An unavailable census route is never evidence of a non-repeated block |
| `DOSAGE_AUDIT_DISAGREES` | block `dosage_status` ∈ {CONTENT_REPAIRED, TRIGGERED_UNREPAIRED, SUSPECT} | the audit's count evidence contradicts the call's dosage |
| `DOSAGE_AUDIT_UNAVAILABLE` | audit FAILED or NOT_RUN, or block `dosage_status` ∈ {UNDERPOWERED, NOT_EVALUATED} | no audit evidence |
| `PRIOR_ONLY_SUPPORT` | `genotype_support = PRIOR_ONLY` | the reads alone do not separate the call from an explored rival (pool-conditional) |

The block-diagnostics output (`block_diagnostics diagnostic_v1`, M11-C) is **enabled for this campaign** as a diagnostics
profile, because it is the only source of the census route. Its off/on scientific identity was verified in M11-C.

**Computational validity** (from M15-A §2A, unchanged): a completed prediction requires exit 0, `prediction = PREDICTION`,
`completion = completed_ranking`, not provisional, and a complete representation guard. Every other run is reported by its
own status. **Resource refusals and no-calls are outcomes, not exclusions.**
