# M11-D: bounded unified-configuration development pilot (frozen before execution)

M11-C is closed. M10-A plus M11-C's gated diagnostic output is the computational caller checkpoint. No default is promoted.
- **Branch:** `m11-d-dev-pilot`, from `m11-c-confidence-inputs` `ed10d53`.
- **Results:** `tests/results/m11d_dev_pilot_20261020/`.
- **Budget:** 3,600 s, including 60 s of disclosed development (inventory reconciliation, no caller run).
- **Every estimate here is a development estimate.** There is no validation, cross-locus reliability or calibrated-confidence
  claim.

## 1. Outcomes (frozen)

**Primary continuous:** per-block **excess** and **normalized excess** over the existing availability floor (evaluator fields;
`scripts/genotype_report` normalization).

**Also reported:**
- QVfrac and all five QVfrac quintiles (M11-C definition: evaluator floor pairs, minimum-edit pairing; values marginally
  above 1 preserved and flagged, never clipped, never read as a probability);
- floor attainment;
- catastrophic excess ≥ 1,000 bp.

**Binary label for the confidence comparison:** `top = QVfrac ≥ 0.8`.
- **Source:** the historical COSIGT top availability-relative category.
- **Scope:** an **evaluation label** only, not a production threshold, and not chosen from M11-C's data.

## 2. Development inputs (frozen; inventory reconciled before freezing; no simulation)

| locus | cells (5×, existing reads) | panel |
|---|---|---|
| CYP2D6 | the 6 pilot donors: HG02769, HG02953, HG03456, NA19650, NA19705, NA19983 | `loo.<donor>.gfa` + the CYP2D6 curation |
| C4 | the first 12 **eligible** development families by sha256("M11-D:C4:" + family), with the donor by sha256("M11-D:C4:donor:" + donor): NA19705, NA19317, HG03065, HG02018, HG00513, NA12329, HG00358, NA19331, HG01457, HG03732, NA24385, NA19434 | `<donor>.loo.gfa` |
| ACOT | the 12 previously genotyped families: HG03209, HG01975, NA19909, NA18747, HG02178, HG00639, NA21106, HG03831, HG02723, HG02074, NA18971, HG02083 | `acot.gfa` + the M9 LOO curation (every outer-evaluation template excluded) |

**C4 eligibility:** 5× reads plus a leave-one-out graph. 51 of the 52 development families are eligible. **HG00096** is
excluded because it is the reference donor and has no LOO graph. That is a stated reason, not an accuracy-based choice.

**Families** are the grouping unit, and a family stays together across loci: NA19705 is family 2368 in both C4 and CYP2D6.

**Configuration:**
- the unified configuration (`compose_m10a` arm D, the per-locus resource profile, no declarations), plus
  `block_diagnostics diagnostic_v1`;
- `--preset short-read --dosage-mode audit --threads 4`;
- time guard = watchdog − 20 s, with watchdogs C4 150 s, CYP2D6 150 s and ACOT 300 s;
- unchanged scientific caps; 6 GiB envelope.

**Evaluation:** `genotype_probe evaluate`, LOO scheme, against the full graph. QVfrac uses `block-sequences` truth with the
same arguments; recomputed distances must equal the evaluator's.

**Recorded:** executable SHA-256, input reads SHA-256 and model paths in `MANIFEST.tsv`. The 18 prospective ACOT validation
families are **untouched**: no run, no evaluation, no truth.

## 3. Confidence candidates (frozen before outcomes are opened)

- **C1:** x = −10·log10(`error_mass_upper`). Where it prints as 0, `legacy_representative_GQ` is used, **only** when
  lower = upper (the correspondence is checked on every block where both are defined), otherwise NA.
- **C2:** C1 plus PRIOR_ONLY (the M11-C-corrected read-support quantity: a fragment contender disagrees at the block).
- **B0:** a constant-probability baseline.
- **Q2 / Q3** (census route, JP tie) are reported as **diagnostics only**.

**Model:** logistic regression, Newton iterations, slopes standardized, L2 = 1e-4 on slopes only (for numerical stability
under separation).

**Fitted on:** assessable, truth-defined blocks. **Assessable** means a completed ranking, not provisional, and the
representation guard complete.

**Comparison:**
- **leave-one-family-out** over the pooled development set; each family is predicted by a model that never saw it;
- metrics per locus and pooled: log loss, Brier and AUC for B0, C1 and C2; Δlog loss (C1 − B0, C2 − B0, C2 − C1) with a
  **family bootstrap** (2,000 resamples, seed 11) 95% interval.

A high overall top rate alone cannot make an uninformative model look useful, because everything is judged against B0.

## 4. Reporting

**Counts and outcomes:**
- donor and family counts; completeness;
- prediction coverage (PREDICTION / PROVISIONAL / NO_PREDICTION);
- completed ambiguous rankings, kept distinct from resource-incomplete searches;
- continuous excess; the categorical outcomes; development out-of-family performance by locus.

**Catastrophic failures** are listed individually, with their out-of-fold C1 and C2 probabilities and their inputs.

**Synonymous-boundary artefacts:**
- **Intervals:** the model degeneracy intervals (ACOT sv:5..ctx:6 and CYP2D6 sv:3..sv:5 from the M9-I inventories of the
  same LOO models; C4 from this pilot's model-only inventory).
- **Measure:** artefact bp = the interval's summed block distances minus its concatenated-interval distance.
- **Handling:** flagged blocks stay visible, and are **not** treated as proven sequence failures.

**Not done:** no new scientific no-calls, no production cutoff, no additional predictors or filters.

## 5. Units (`m11d_run.py`)

`build`, `gate.selftest` (logistic recovers simulated slopes; the baseline equals the prevalence; LOFO covers every row;
separation stays finite), `call / eval / truth` per cell, `inventory`, `outcomes`, `fit`.
- **Expected:** about 1,800 s.
- **Budget rule:** a unit is admitted only if spent + its watchdog ≤ 3,600; otherwise it is recorded NOT RUN. Watchdogs are
  never shortened.
