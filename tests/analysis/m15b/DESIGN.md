# M15-B: final preflight and reserved ACOT evaluation (frozen before any reserved genotyping)

**What this is:** a descriptive reserved evaluation of the frozen release candidate (RC), run as the **named experimental
profile** `genotype-dev-m10a-v1`.

**What this is not:** a default promotion, a reliable-call category, an accuracy PASS, a confidence fit, or a model or
prediction change.
- **Branch:** `m15-b-reserved-acot`, from `m15-a-release-candidate` `9756d24`.
- **Results:** `tests/results/m15b_reserved_acot_20261021/`.
- **Budget:** 7,200 s, including 90 s of disclosed development (inventory and curation checks, archived-hash lookups).
- **Resources:** 4 threads; one process at a time; 6 GiB.
- **Contract:** `CONTRACT.md` (the corrected M15-A contract and the recorded decisions).

## 1. Regression preflight

- **`build`:** every target is built, and `panvar` must equal the RC SHA-256 (`47042a5c…`), or the unit fails.
- **`ctest.full`:** the **full enabled suite**, run serially, including the three heavy tests omitted in M15-A.
  - Per-test timeout 900 s; unit watchdog 2,400 s.
  - Disabled tests are reported separately (`Not Run (Disabled)`).
  - **A failure stops reserved execution:** `prepare` needs `ctest.full`.
  - `genotype_command_e2e`'s timing is recorded as observed. Stale scheduling metadata neither proves nor excludes a
    performance regression.

## 2. Provenance (frozen in `prepare` and `MANIFEST.tsv` before any call)

**Families:** the 18 reserved ACOT families, with the donor from `m11a_confidence_design_20261020/family_exposure.tsv`.
`prepare` asserts that each has **exactly two graph paths, both in the LOO curation**, so they are outside the active panel.

**Disclosed:**
- **prior exposure:** their assemblies were used in M8-B masked prior scoring;
- **the full-panel graph frame:** their paths, like every path, help define the blocks, because curation never changes
  blocks.

**They are not pristine independent validation.**

**Graph, curation and binaries:** `acot.gfa`, `curation.LOO.txt` (M9-H), and the binaries, all hashed.

**Profiles** (hashed):
- the base `acot_5` (`repair_max_neighbours` 1004, retained as the supported per-locus profile);
- the performance overlay;
- the M10-A overlay;
- `block_diagnostics diagnostic_v1` (diagnostics only; needed for the census route).

**Reads** (the M8-C recipe):
- `simulate`: depth 30, 100 bp, seed 20260916, seed-salt = donor;
- then `subsample` to 5×, name-salt = donor.
- **Lineage gate:** regenerating M9-H's HG02723 reads must reproduce its three recorded SHA-256 values byte for byte,
  otherwise `prepare` fails.

**Limits:** caller guard 280 s; watchdog 300 s; 6 GiB.

## 3. Execution

- **Calls:** one per donor, in the manifest order (sorted donor), with no outcome-driven change and no selective rerun.
  Resource refusals and no-calls are **outcomes**.
- **Admission:** a call runs only if spent + 300 + 3,600 (the evaluation reserve) ≤ 7,200.
- **`freeze`:** hashes every prediction output **before** any truth unit can run (`eval.*` and `truth.*` need `freeze`).
- **Evaluation:** `genotype_probe evaluate --scheme loo`, plus a block-sequences truth projection, used only for the
  representation interval.

## 4. Analysis (`m15b.py analyse`, frozen)

**Per block:** raw distance, availability floor, excess, normalized excess (excess / max(truth bp, called bp)) and floor
attainment.

**Category** (mutually exclusive, in this order):

| category | rule |
|---|---|
| `UNDEFINED_TRUTH` | |
| `NO_PREDICTION` | missing prediction |
| `ALL_EMPTY` | visible; excluded from informative summaries |
| representation-dependent | ACOT sv:5..ctx:6, when summed block distance − concatenated interval distance > 0: `REPRESENTATION_INTERVAL_EXACT` if the concatenated interval distance is 0, else `REPRESENTATION_UNRESOLVED` |
| `BLOCK_OUTCOME` | informative |

**Denominators:** every category is counted. Unavailable truth or failed evaluation never silently shrinks a denominator.
Runs without a completed prediction are listed with their status.

**Per family:** informative blocks; excess sum, median and max; floor-attaining blocks; catastrophes.

**Family bootstrap** (2,000 resamples, seed 11): mean excess per informative block, mean normalized excess, catastrophic
block rate, floor-attainment rate, and the fraction of families with at least one catastrophe.

**Every catastrophic block (≥ 1,000 bp)** is listed individually, with its support fields and annotations.

**Annotations** (descriptive associations, with block and family counts):
- `REPEATED_UNIT_STRUCTURE` (YES / NO / UNKNOWN);
- `DOSAGE_AUDIT_DISAGREES` / `DOSAGE_AUDIT_UNAVAILABLE`;
- `PRIOR_ONLY_SUPPORT`;
- `BLOCK_REPRESENTATION_AMBIGUOUS`, `BLOCK_GENOTYPE_UNRESOLVED`;
- the headline genotype and dosage filters.

**Per run:** completion, sequence ambiguity (`SEQUENCE_AMBIGUOUS`), no-call reason, dosage audit status, wall time, peak
memory and refusals.

**Not computed:** QVfrac (descriptive only and not needed for this report). No threshold or acceptability judgement is
applied.

## 5. Deliverable

Four conclusions, kept separate:
1. computational and regression status;
2. reserved ACOT simulated-read performance (descriptive);
3. remaining limitations;
4. whether the evidence supports publishing an **experimental** release.

**Not extrapolated:** independent C4/CYP2D6 accuracy, real reads, other coverages, or calibrated confidence. The batch stops
for review before any release operation.
