# Final evaluation protocol (frozen draft; for approval before any reserved outcome is opened)

**Status:** nothing in this document has been run. **No reserved family has been genotyped, and no reserved truth outcome
has been opened.**

**Order:**
1. The protocol, `OUTPUT_AND_FILTER_CONTRACT.md` and the RC are approved.
2. Inputs are prepared.
3. Every reserved run executes.
4. The evaluator runs.
5. The analysis runs.

No rule, setting, annotation or threshold changes after step 3 begins.

## 1. Evidence classes

| class | members | what it can support |
|---|---|---|
| **Development** (inspected; used for design and selection) | all C4 families with reads (52; headline filters v2 and dosage policy v2 were selected on the C4 campaign); the 6 CYP2D6 pilot donors (amended profile selected there); the 12 genotyped ACOT families (M8-C/D, M9-C..H); every M-series cell | descriptive development performance only |
| **Reserved, prior-exposed** | the **18 ACOT outer-evaluation families never genotyped**: ST016, NA20809, NA21110, NA20806, NA18943, PR04, CLM31, GBR002, PK35, VN046, CLM12, BB46, GB23, GB31, ST232, GB21, 2490, NA19468 (donors as in `m11a_confidence_design_20261020/family_exposure.tsv`) | the RC's ACOT performance on families **never genotyped or inspected**. **Disclosed:** their assemblies were used in M8-B masked prior scoring. They are excluded from the calling panel by the M9-H LOO curation (every outer-evaluation template) |
| **Independent** | none available for C4 or CYP2D6 (no clean holdout; M11-A) | — |

## 2. Manifest and inputs (prepared in the approved batch, not before)

- **Families:** the 18 above. One donor per family: the donor listed in `family_exposure.tsv`; each family has 2 templates,
  which form the truth diploid.
- **Reads:** the M8-C recipe as used by M9-H `prepare.py`, unchanged:
  - `genotype_probe simulate`: depth 30, read length 100, insert 350 ± 40 (4σ), error 0.002, outlier 0.001, seed 20260916,
    seed-salt = donor;
  - then `subsample` to 5×, name-salt = donor.
  - The input SHA-256 values are recorded before any call.
- **Model and panel:** `acot.gfa` with `curation.LOO.txt` (M9-H; every outer-evaluation template excluded). The truth graph
  is `acot.gfa` itself.
- **Configuration:** the RC exactly — `panvar` sha256 `47042a5c…`, the composed arm-D profile with `acot_5` base, 4
  threads, dosage audit.
  - **Watchdog** 300 s per call, time guard 280 s, 6 GiB envelope.
  - **Admission:** the evaluation reserve is set before any call is admitted.

## 3. Outcomes and analysis rules

**Per block, from the evaluator (`genotype_probe evaluate --scheme loo`):**
- **excess** = unordered distance − availability floor (the active model's floor, as in M11-D);
- **normalized excess** = excess / max(truth bp, called bp);
- floor attainment;
- **catastrophic error:** excess ≥ 1,000 bp.

**Block categories** (M11-E eligibility, frozen):

| category | rule |
|---|---|
| `ALL_EMPTY_COUNTED` | counted in coverage; not informative |
| `UNDEFINED_TRUTH` / `NO_PREDICTION` | reported apart |
| representation-dependent | blocks in an evaluated degeneracy interval with artefact bp > 0 (ACOT sv:5..ctx:6). **Interval-exact** when the concatenated predicted interval equals the truth; otherwise **unresolved**. Never relabelled correct |
| `BLOCK_OUTCOME` | everything else |

**Coverage and retained calls:**
- **completed-prediction rate,** by contract §A;
- **provisional / no-prediction / resource-stop counts,** each with its reason;
- **retained-call rate** under any approved annotation policy. With no reliable class, this is just the completed-prediction
  rate.

**Annotation association** (descriptive): non-top and catastrophic rates by B annotation — PRIOR_ONLY, repeated-unit,
dosage-audit status.

**Donor-level uncertainty:**
- every summary is reported per family;
- intervals come from a **family bootstrap** (2,000 resamples, seed 11);
- block counts are always given beside family counts.

**Every catastrophic error is listed individually,** with its annotations and support fields.

**Descriptive only:**
- **QVfrac:** never a pass criterion.
- **The M10-D chain-interior locus distance:** complementary.

**Resources:**
- **per run:** wall time, peak footprint, stage seconds, and resource refusals or caps;
- **per batch:** the totals.

## 4. What the reserved evaluation can and cannot establish

**It can establish:**
- the RC's ACOT block accuracy and catastrophic-error rate at 5× on **simulated reads from real assemblies**;
- the completion and resource behaviour;
- all of it on 18 never-genotyped families, with family-level uncertainty, conditional on the disclosed prior exposure.

**It cannot establish:**
- **C4 or CYP2D6 performance:** development evidence only, with no clean holdout in this repository;
- **real-read performance:** library artefacts, mapping biases, coverage non-uniformity;
- **other coverages;**
- **calibrated confidence,** or a reliable-call category;
- **generalization** beyond the panel's ancestry composition.

**The pass/fail reading needs Davide's approved acceptable-excess and catastrophic-error criteria** (contract §3). The
analysis reports against them; it does not choose them.
