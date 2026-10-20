# M11-B result: offline development audit of per-block confidence against the floor-adjusted block outcome (descriptive only)

**Run.** Branch `m11-a-confidence-design`.
- **Commits:** freeze `25a0a36`; Amendment 1 `4140a02`; result after.
- **Results:** `tests/results/m11a_confidence_design_20261020/m11b/`, containing `ledger.tsv`, `runs/manifest.a1/manifest.tsv`
  and `runs/audit/{block_observations, binned_summary, extremes, rounding_and_artefacts, locus_level}.tsv`.
- **Budget:** 19.2 of 300 s, including 15 s of disclosed development.
- **Not done:** no reads, genotyping, fitting, thresholds, calibration or prediction changes.

**Amendment 1** (declared before the audit).
- **Duplicates:** the first manifest marked 85 C4 coverage replicates as duplicates. Their evaluator tables are identical,
  but their reads differ, and the table has no coverage column. The duplicate key now includes the reads' SHA-256.
- **No-call runs:** these were "excluded"; they are now recorded as NO_CALL and kept in the run-level table.
- `manifest` stays recorded.

## 1. Units, quantities and outcomes

**Outcome** (the established contract, reused verbatim from `scripts/genotype_report/load.py` and `aggregate.py`):
- `excess` = the evaluator's recorded unordered distance minus its availability floor, in bp;
- `normalized_excess` = excess / max(truth₁+truth₂, called₁+called₂);
- `floor_attained` = (excess = 0), defined only where both truth alleles exist;
- `catastrophic` = excess ≥ 1,000 bp (the report's existing default), **descriptive only**.

**Per-block truth-free quantities:**
- `error_mass_upper`: projection-aware and pool-conditional, `sequence_first_v2` runs only;
- `conditional_GQ` (posterior), in every run.

**Tabulated, not used:**
- **Locus-level support** (`decision_set_probability`; legacy `winner_candidate_set_probability`) is tabulated **separately**
  (`locus_level.tsv`) and never broadcast onto blocks.
- **Completion status** is tabulated, not used as evidence. A completed bounded search does not prove that every better
  candidate was explored.

**Strata are configurations, never pooled:**

| stratum | runs (no-call) | donors | truth-defined blocks | configuration provenance |
|---|---|---|---|---|
| v2:M1-B | 4 | 1 (NA19650: C4 10×/30×, CYP2D6 5×/30×) | 58 | fresh combined-route caller, M1-B profiles |
| v2:M9-C | 7 (+1 NO_CALL: LOO.d5 RESCORE_INCOMPLETE) | 5 | 133 | ACOT v2, r = 0.01, JP off |
| v2:M9-G | 5 JP-on (the 5 JP-off runs are duplicates of M9-C) | 4 | 95 | ACOT JP1.1 on |
| v2:M9-H | 12 (JP on and off, L1–L4, Z1–Z2) | 4 | 228 | ACOT JP off / on |
| legacy:C4-production-11fe681 | 483 (+1 NO_CALL) | 57 (incl. mosaic subjects) | 5,309 | C4 production caller 11fe681 |
| legacy:CYP2D6-pilot-a5 | 27 (+9 NO_CALL) | 7 | 511 | CYP2D6 pilot, amendment-5 profile |

Blocks from one donor are not independent. Donor counts are given beside every block count, and no test is made.

## 2. What the audit shows

**Ordering information is weak to moderate.** The rank separation is P(a non-floor block has more pool error mass, or lower
GQ, than an at-floor block):

| stratum | `error_mass_upper` | `conditional_GQ` | floor attained |
|---|---|---|---|
| M1-B | 0.559 | 0.547 | 39 / 58 |
| M9-C | 0.574 | 0.574 | 106 / 133 |
| M9-G | 0.761 | 0.765 | 75 / 95 |
| M9-H | 0.741 | 0.739 | 168 / 228 |
| C4 production | — | 0.792 | 4,599 / 5,309 |
| CYP2D6 pilot | — | 0.613 | 421 / 511 |

- **The two v2 quantities are nearly interchangeable.** The block error mass carries essentially the same ordering as the
  legacy GQ.
- **Floor-attainment rates rise with support**, for example C4 GQ bins: 0–10 → 0.34; 60–70 → 0.91; ≥100 → 0.95. The bins
  are reporting decades, not thresholds.
- **The relation is not monotone enough to trust.**

**High-support errors exist and are structural** (`extremes.tsv`, HIGH_SUPPORT_ERROR):
- **C4 production:** 44 catastrophic blocks in 14 donors. For an illustrative split only (not a threshold), 31 of them, in
  11 donors, carry conditional GQ ≥ 50, **all at sv:7**. Examples: HG04217 sv:7 GQ 4,454, excess 32,709 bp; NA20847 GQ 1,008,
  32,708 bp.
- **M1-B NA19650:** C4 sv:7 error mass 7.7e-75, excess 65,444 bp (the known 3+1 allocation error).
- **M9-C / M9-G NA19909:** ACOT sv:13, 1.1e-38 and 5.4e-45, excess 27,270 bp.
- **CYP2D6 pilot:** NA19650 sv:9 (the known dosage error), catastrophic at high GQ.
- **Small high-support errors at ACOT sv:5:** M9-C HG00639, 1.3e-68 (6 bp); M9-H NA18971, 2.8e-63 (4 bp), a synonymous
  artefact.

These are the cases earlier batches attributed to unrepresentable or unexplored alternatives (M1-B, M9-D/E). The pool is
confident because the competitor was never in it.

**Low-support predictions at the floor:**
- context blocks with PRIOR_ONLY support, error mass 0.10–0.22, exact (M1-B ctx:0, M9-C ctx:0/2/16);
- M1-B C4 sv:7 at 6.5e-3, exact.

**Rounding:**
- `error_mass_upper` printed as **0** in 2 blocks. One of them is an error: M9-H HG02083 sv:5, excess 13 bp. A printed 0 is
  **rounding or underflow, not certainty**.
- `conditional_genotype_posterior` printed as **1** in 1,740 C4, 205 CYP2D6 and 104 v2 blocks.
- GQ itself stays finite (no `inf` and none missing), so ranks are preserved in GQ but not in the printed posterior.

**Synonymous-boundary artefacts**, identified from existing records (block distance summed over the model's degeneracy
interval, minus M9-I's compound distance):
- **Quantified (YES):** 20 v2 blocks (M1-B 4, M9-C 2, M9-G 2, M9-H 12), all in ACOT sv:5..ctx:6 or CYP2D6 sv:3..sv:5.
- **Unquantified:** 28 CYP2D6-pilot blocks with excess in sv:3..sv:5, where no compound distance exists.

They are kept in every count and flagged in `block_observations.tsv`. Several high-support "errors" at ACOT sv:5 are among
them.

## 3. Reading (development only; not calibration or validation)

- **The per-block pool-conditional quantity carries ordering information but is not, by itself, a reliable confidence.**
  Its most consequential failures are confident structural errors whose true competitor was absent from the pool.
- **Calibrating this scalar cannot detect those failures.** Any eventual assessment must report them as a class. Whether an
  additional truth-free input could flag them is open: genotype support, block kind, search-completeness facts. **Nothing is
  chosen here.**
- **No scalar predictor is committed to.** If the validation pilot shows the pool error mass is uninformative for a block
  class, the assessment for that class stays unassessed.

## 4. Revised prospective pilot (for approval): M11-C

**Primary outcome:** the floor-adjusted **block** outcome (excess in bp and normalized excess, as defined in §1). Floor
attainment and catastrophic counts are descriptive.

**Confidence input:** the per-block pool error mass (`error_mass_upper`). Locus-level support is reported separately, and
completion and resource status only gate assessability.

**Configuration:** the unified M10-A configuration (D), leave-one-out panels. Every run's model, profile and rule versions
are recorded. Nothing is pooled across configurations.

| stage | families | data | truth |
|---|---|---|---|
| **Development** | exposed families with existing reads: the 6 CYP2D6 pilot donors; 12 C4 campaign families chosen by a frozen hash; the 12 ACOT families of M8-C/D, M9-C..G and M9-H | **existing reads** (no simulation), 5× | opened: descriptive relation per locus and block kind, by donor; the calibration **procedure** (not values) frozen at the end |
| **Validation** | ACOT's 18 never-genotyped outer-evaluation families. Disclosed: **prior-score-exposed** in M8-B, not pristine | reads simulated with the existing `genotype_probe simulate`, 5× | **sealed:** genotyped, but not evaluated until the procedure is frozen |
| C4 / CYP2D6 validation | **none available** (every family genotyped or used for selection) | — | stated as a limitation |

**Runtime estimates** (from measured unit costs):
- **Development:** CYP2D6 about 6 × 30 s; C4 about 12 × 30 s; ACOT about 12 × 75 s; evaluation about 20–40 s per run.
  About **45 min** in all.
- **Validation genotyping:** about 18 × 80 s, about **25 min**.
- **Envelope:** 6 GiB; watchdog per run.

**Decisions Davide must make before M11-C is frozen** (block-level, not QV bands):
1. **The acceptable block outcome.** Floor attainment (excess = 0), or an excess allowance in bp and/or normalized excess,
   possibly different for context and SV blocks.
2. **The catastrophic definition.** Reuse the report's 1,000 bp default, or another value.
3. **The acceptable confident-error rate.** Among blocks reported confident, the maximum rate of unacceptable outcomes, and
   whether any confident **catastrophic** error is tolerable.
4. **The independent unit.** Donors as the unit (family-level resampling) for every rate.
5. **Synonymous-boundary artefacts.** Counted as block errors, reported beside them, or assessed at the declared compound
   interval.
6. **The validation design.** Accept ACOT's 18 prior-score-exposed families as validation, with C4 / CYP2D6 development-only;
   permit read simulation; 5× first.
7. **The quantity to calibrate.** Fix the input (`error_mass_upper`), and any stratification (block kind, genotype support)
   to declare **before** the development data is examined.

---
## Qualifications (appended after review; the audit above is unchanged)

**1. The audit does not establish that every confident error comes from an absent candidate.** The documented attributions
differ:

| confident error (development) | documented cause | source |
|---|---|---|
| ACOT sv:5, M9-C d1 / d3 / d6 (E1, E3, E6) | **search omission**: a one-block replacement scores higher but was never evaluated (proposal emissions tied at sv:5) | M9-D, M9-E |
| ACOT sv:5 d2 (E2), sv:9 d3 (E4) | **copying-prior effect**: the likelihood favours truth, and the LOO-context prior favours the call | M9-D |
| ACOT sv:5 d4 (E5) | **weak likelihood preference** for the call (about 1 nat) | M9-D |
| C4 NA19650 sv:7 (65 kb) | **mixed**: the truth is unspellable from the LOO panel (availability), **and** the saved 2+2 class is present but scores lower (likelihood/model; residual content, M6-B) | M1-B, M1-C, M6-B |
| CYP2D6 NA19650 sv:9 at 5× (12 kb) | **likelihood/model limitation**: the true 2+3 class is in the pool (best rank 4) but loses by about 31 nats; a count model fixes this LOO dosage (M6-A) | M1-B, M6-A |
| ACOT NA19909 sv:13 (27 kb, M9-C d1) | **unattributed** | — |
| ACOT sv:9 ≈ 475 bp errors (M9-H L1, L4; M9-C d2) | **unattributed** (PRIOR_ONLY blocks; resembles E4, not shown) | — |
| C4 campaign sv:7 catastrophic errors in 10 further donors (e.g. HG04217, NA20847) | **unattributed** (only NA19650 was studied) | — |

**2. Scope of calibration.** A scalar calibration **can** correct average over-confidence for a class of predictions. It
**cannot** distinguish good from bad predictions whose inputs are indistinguishable. This does **not** make calibration
categorically impossible. Whether the documented failure classes have distinguishable truth-free inputs is an empirical
question, addressed below.

---
## Correction to Q1's meaning (appended by M11-C; the associations above are unchanged)

**What `PRIOR_ONLY` is not.** It does **not** mean "no fragment informs the block" (§2 and CALIBRATION_DESIGN Q1 said so).

**What it is** (`block_genotype_read_supported`, `genotype_repair.cpp` / `genotype_command.cpp`). A block is `PRIOR_ONLY`
when at least one **fragment contender** has a different unordered genotype at that block. A fragment contender is a ranked
rival whose read-only fragment term is not certifiably below the winner's, within numeric bounds.

**What follows:**
- Both `READS` and `PRIOR_ONLY` are conditional on the explored, ranked pool and on representative labellings.
- The lower floor attainment of PRIOR_ONLY blocks therefore means: where the reads did not separate the call from an
  explored rival, the call was correct less often.
