# M15-B result: regression preflight passed; reserved ACOT evaluation complete (18/18 predictions, 0 catastrophes, median excess 0); the evidence supports an experimental release with stated limits

**Run.** Branch `m15-b-reserved-acot`; freeze `13b5f13`; result after.
- **Results:** `tests/results/m15b_reserved_acot_20261021/`.
- **Budget:** 4,121.9 of 7,200 s, including 90 s of disclosed development.
- **Units:** 58 of 58 passed first time, with no amendment and no watchdog or envelope event.

**Not done:**
- no default promotion, model change, prediction correction, threshold, confidence fit or selective rerun;
- no merge, push or tag.

**This is a descriptive reserved evaluation.** It is not an accuracy PASS and not calibrated confidence.

## 1. Computational and regression status

| check | result |
|---|---|
| build | every target; **`panvar` = the RC SHA-256 `47042a5c…`** (asserted) |
| full enabled suite (serial, 900 s per test) | **63 of 63 passed**, 1,004.8 s; 2 disabled in the suite (`panvar_real_regressions`, `panvar_rebuild_stats`) |
| the three heavy tests omitted in M15-A | passed: `genotype_report` 197.7 s, `genotype_branch_lifetime` 189.2 s, `genotype_block_delta` 144.8 s |
| simulator lineage | M9-H's HG02723 30×→5× reads regenerated **byte for byte** (all three SHA-256) |
| provenance | 18 families, each with exactly two graph paths, **both in the LOO curation** (outside the active panel); graph, curation, profiles, binaries and reads hashed before any call |
| prediction freeze | every call output hashed (`runs/freeze`) **before** any truth unit |

**Timings recorded as observed:**
- `genotype_command_e2e` took 174.4 s, and the heavy tests ran longer than their stored cost data (9–105 s).
- Stale scheduling metadata neither proves nor excludes a performance regression; this is recorded only.

## 2. Reserved ACOT simulated-read performance (descriptive; `runs/analyse/`)

**Runs:**
- **18 of 18** produced completed predictions: 14 strict and 4 genotype-equivalent (`SEQUENCE_AMBIGUOUS`: HG01150, HG03688,
  NA18943, NA21110);
- **0** no-calls, provisional results or resource refusals;
- dosage audit COMPLETE in 18 of 18;
- call wall time 80–198 s (median 117.7, sum 2,145 s); peak footprint at most 2.56 GiB.

**Blocks:** 342 in total, with no block silently dropped:

| category | blocks |
|---|---|
| BLOCK_OUTCOME (informative) | **296**, across all 18 families |
| ALL_EMPTY (visible; excluded from informative summaries) | 36 |
| REPRESENTATION_UNRESOLVED (sv:5..ctx:6; artefact 1–2 bp; concatenated interval distance 2–10 bp, so not exact) | 10 (5 families) |
| REPRESENTATION_INTERVAL_EXACT | 0 |
| UNDEFINED_TRUTH / NO_PREDICTION | 0 / 0 |

**Excess over the availability floor**, on the 296 informative blocks:

| measure | value | family-bootstrap 95% (2,000 resamples) |
|---|---|---|
| excess: median / q90 / q99 / max | 0 / 8 / 98 / **478** bp | |
| mean excess per block | 6.26 bp | [2.30, 11.36] |
| normalized excess: median / q90 / q99 / max | 0 / 3.2e-4 / 0.31 / 0.36 | |
| mean normalized excess | 0.0062 | [0.0019, 0.0114] |
| floor attainment | 66.6% (197 / 296) | [60.6%, 72.0%] |
| catastrophic blocks (≥ 1,000 bp) | **0** (0 families) | |

**On zero catastrophes:** with 0 in 18 families, a crude one-sided 95% upper bound on the family-level catastrophe rate is
about 15% (rule of three: 3/18 = 17%). The reserved set does not show the rate is zero. For comparison, the development
ACOT set (12 families) contained one catastrophe (M11-D NA19909, 27 kb).

**Per family** (`families.tsv`): excess sums range 4–613 bp. Two families dominate:
- **GB21 (HG02572) and NA19468:** sv:9 called as 902 + 424 bp against a truth of 424 + 424. That is one extra ~478 bp
  segment, with excess 478 and 475, a floor of 0, and normalized excess **0.36 and 0.35**.
- **Annotations:** both blocks are GENOTYPE_FILTER **PASS**, with `REPEATED_UNIT_STRUCTURE` NO (census route not
  repeated-unit) and no dosage disagreement. Support is READS for HG02572 and PRIOR_ONLY for NA19468.
- **These are the largest relative errors.** They are below the 1,000 bp catastrophe definition, and none of the
  annotations flags them.

**Next largest errors:**
- NA21110 sv:13 (110 bp; 1.3e-3 normalized);
- PK35 sv:13 (98 bp);
- NA20346 ctx:16 (78 bp);
- all other blocks ≤ 43 bp.

**Annotations** (descriptive associations; informative blocks, with block and family counts):

| annotation | value: blocks (families), mean excess bp | | |
|---|---|---|---|
| `REPEATED_UNIT_STRUCTURE` | YES: 18 (18), 13.9 | NO: 278 (18), 5.8 | UNKNOWN: 0 |
| `DOSAGE_AUDIT_UNAVAILABLE` (block underpowered / not evaluated) | True: 61 (18), 19.6 | False: 235, 2.8 | |
| `DOSAGE_AUDIT_DISAGREES` | True: 5 (5), 2.0 | False: 291, 6.3 | |
| `PRIOR_ONLY_SUPPORT` | True: 123 (17), 6.5 | False: 173, 6.1 | |
| block representation ambiguous / genotype unresolved | 0 / 0 | | |
| GENOTYPE_FILTER | PASS on all 296 | | |
| DOSAGE_FILTER | NOT_EVALUATED 61 (mean 19.6); NOT_APPLICABLE 224; LOW_QUALITY 5; PASS 6 | | |

No catastrophe exists, so no annotation can be assessed against catastrophes here. The differences above are small-sample
descriptions, not risk estimates.

## 3. Remaining release limitations

- **Scope of the evidence:**
  - ACOT only, at 5×, on **simulated reads** from real assemblies;
  - the families are **prior-exposed** (M8-B masked prior scoring) and inside the **full-panel graph frame**, so this is not
    pristine independent validation;
  - 18 families, so the rate bounds are wide.
- **C4 and CYP2D6:** development evidence only. All five development catastrophes were at these loci or development ACOT;
  no independent holdout exists.
- **No calibrated confidence and no reliable-call category.**
  - GENOTYPE_FILTER PASS covered both near-500 bp sv:9 errors.
  - GQ and error mass are uncalibrated.
- **The annotations missed the two largest reserved errors.** Dosage evidence was unavailable for 61 of 296 informative
  blocks.
- **Known open items** (from `RELEASE_BLOCKERS.md`):
  - the ACOT 30× route-equivalence gap;
  - the packaging decision (named profile);
  - the locus-specific neighbour limits, retained as per-locus profiles;
  - release notes.

## 4. Does the evidence support publishing an experimental release?

**Yes, as an experimental short-read release** with the named profile `genotype-dev-m10a-v1`. The evidence:
- reproducible computation (frozen source and binary, full regression suite, byte-reproduced inputs);
- complete, refusal-free execution on 18 reserved families;
- near-floor accuracy on most blocks (median excess 0; 99% of informative blocks ≤ 98 bp);
- no catastrophic error in the reserved set.

**Not supported:**
- a validated or high-confidence production caller;
- any accuracy claim for C4 or CYP2D6;
- real-read performance, other coverages or calibrated confidence;
- any reading of PASS as reliability.

**The release decision itself remains Davide's.** The release notes must carry §3's limitations, including the two missed
sv:9 length errors.

**Stopped for review before any release operation.**
