# M11-A: confidence-calibration design (planning and inventory only)

**Status.**
- **Checkpoints:** M10-D is accepted as the sequence-similarity evaluation checkpoint. M10-A remains the computational
  caller checkpoint.
- **Closed:** no further structural-alphabet, event-comparator or aligner work.
- **Not done in this batch:** no genotyping, simulation, threshold fitting or caller change.
- **Execution:** under 1 s, for the inventory script.
- **Where it lives:** branch `m11-a-confidence-design`. Inventory: `tests/results/m11a_confidence_design_20261020/family_exposure.tsv`.

## 1. Truth-free quantities the caller already emits (unified configuration, `sequence_first_v2`)

None of these is a calibrated accuracy. Every probability below is **conditional on the explored candidate pool**. The
caller says so itself: `probability_scope = candidate_conditional_over_the_explored_pool;not_calibrated`,
`structure_confidence = NA (UNCALIBRATED)`.

| group | quantity (file) | what it measures |
|---|---|---|
| **A. Pool separation / support** | `decision_set_probability`, `representative_probability`, `decision_set_groups`, `exact_score_ties`, `pool_groups` (sequence_locus) | mass of the reported sequence pair, or its numerically inseparable set, among explored candidates |
| | `error_mass_lower` / `upper`, `block_status`, `legacy_representative_GQ` (sequence_blocks) | per-block pool mass against the call, projection-aware bounds |
| | `support` / `error` bounds, `target_resolution`, `competitor` (sequence_compound, when enabled) | the same over declared compound intervals |
| | `winner_candidate_set_probability`, `winner_set_size`, `score_gap_to_runner_up`, `best_other_genotype_score`, `winner_numeric_error_bound`, `candidate_set_log_weight` (mosaic_status); `separation`, `separation_certified` (locus_genotype) | score margin to the best rival; numeric error of the winner |
| | `sequence_status`, `representation_ambiguous_blocks`, `genotype_unresolved_blocks`, `representation_intervals` | whether one sequence pair is identified within the decision set (`resolution_meaning`: not a confidence claim) |
| | `phase_support` / `legacy_phase_support` (for example PRIOR_ONLY) | phase support, separate from genotype |
| **B. Search completeness** | `prediction` (PREDICTION / PROVISIONAL / NO_PREDICTION), `completion`, `provisional`, `projection_coverage` | whether the ranking is a completed exact ranking of the explored pool |
| | `legacy_run_status` (COMPLETE, CONDITIONAL_GENOTYPE_EQUIVALENT, RESCORE_INCOMPLETE, REPAIR_INCOMPLETE …), `no_call_reason`, `repair_rounds_run`, `tie_branched`, `tie_branches_run`, `early_incumbent_status`, `checkpoint_last_stage`, `representation_guard_complete` | how the search ended, and whether the guard completed |
| | `junction_eligible_blocks` (JP1.1) | where the search used the junction proposal heuristic |
| **C. Resource status** | `time_limit_seconds`, `footprint_limit_bytes`, peak footprint, configured caps (`max_*`), cap refusals (for example "physical-origin cap", "neighbourhood_bound") | whether a guard or cap stopped or limited the run; never evidence about correctness |
| **D. Existing production filters (do not add more)** | `locus_filter`, `quality_rule`, `quality_pass` / `lowqual_blocks`, `call_class` (headline filters v2); dosage audit statuses (CONSISTENT / SUSPECT / UNDERPOWERED / NOT_SENSITIVE) | rule-based filters selected on the C4 campaign; the dosage audit is diagnostic only |

## 2. Simplest route to one user-facing confidence assessment

**The prediction is always emitted:** PREDICTION or the PROVISIONAL representative. NO_PREDICTION only when the caller has
none.

**Exactly one assessment field, in two steps.**

**Step 1, assessability (completeness and resource only; no score).** The assessment is `NOT_ASSESSED`, with the single
governing reason, whenever:
- `completion` is not a completed ranking;
- the prediction is PROVISIONAL;
- a cap, time or footprint guard limited the run;
- the representation guard did not complete.

These are facts about the search, not about correctness.

**Step 2, one calibrated mapping of one pool quantity.** For assessable runs, a single monotone calibration maps
`e = 1 − decision_set_probability` (locus-level pool error mass) to the approved outcome scale of §3. That scale is either
an expected error rate or QV, or the probability of a band; Davide decides which.
- **Fitted later**, on development families only.
- **Checked** on the family-disjoint validation set.
- **No other quantity enters the user-facing field.** The rest of §1 stays as diagnostics.

**If calibration shows that `e` does not separate errors,** the field stays `NOT_ASSESSED (uncalibrated)`. This is
plausible: M9-C found ACOT sv:5 errors at extreme pool confidence. There is **no fallback combination** of filters, and the
existing headline filters are not extended.

**Rationale.**
- **One input and one map** avoid stacking independently tuned filters.
- **Separating Step 1 from Step 2** keeps "the search did not finish" distinct from "the search finished and the pool is
  confident".
- **The input is locus-level** because the M10-D outcome is a whole-region sequence distance.

## 3. Evaluation outcomes (kept separate; never mixed)

| | **O1: sequence distance (M10-D contract)** | **O2: per-block excess over the availability floor (existing evaluator)** |
|---|---|---|
| unit | one sample × locus × coverage × **scope** (the chain interior; or a named framed prefix) | one block |
| quantities | D1, D2, D (whole-region assignment), selected assignment and tie, L1, L2, `error_rate = D/(L1+L2)` (unfloored), identity, QV (benchmark `qv_from`) | unordered distance, availability floor, `excess = distance − floor`, truth availability |
| denominator | **every run attempted**, by status: EVALUATED; EMPTY; UNSUPPORTED (named reason: truth projection missing at block X, no prediction, unmatched endpoints); incomplete runs (PROVISIONAL, cap or time stop) are evaluated on the emitted prediction **and** carried as their own stratum | blocks with defined truth; undefined-truth blocks counted separately |
| never | subtract a floor from D; average across scopes | add to or compare numerically with D |

**Historical COSIGT/Locityper thresholds in the code, reusable explicitly rather than refit** (`benchmark_command.cpp`):
- **QV bands:** `<17`, `17–23`, `23–33`, `>33`.
- **Length-fair quintiles:** QV / QV_max with QV_max = 10·log10(2S), binned 0.0–0.2, …, 0.8–1.0.

**No correctness cutoff is chosen here.**

## 4. Family-disjoint development / validation inventory

**Source.** `family_exposure.tsv`: one row per locus × family, built from existing records only. Families come from M8-A's
grouping, roles from the M8-B splits.
- **Genotyping exposure:** the C4 production campaign (484 jobs, 66 subjects), the CYP2D6 pilot (6 donors), and ACOT's
  M8-C/D, M9-C..G, M9-H and the development cell.
- **Also recorded:** failure inspection and parameter selection.

| locus | DEVELOPMENT (genotyped or inspected) | validation candidate | panel-exposed, never genotyped | excluded (reference genome) |
|---|---|---|---|---|
| ACOT | 13 | **18** | 201 | 3 |
| C4 | 52 | **0** | 7 | 2 |
| CYP2D6 | 6 | **0** | 53 | 2 |

**What counts as exposure:**
- **C4.** Every family with reads was genotyped in the production campaign (M8-A's "uninspected" HG02018, HG03520,
  NA19705, NA19836 and NA20355 included), and the headline filters v2 and dosage policy were **selected on that campaign**.
  HG04036 has no reads; the other six never-genotyped families are in M8-A's development set. **No C4 family is a clean
  holdout.**
- **CYP2D6.** Six donors were genotyped in the pilot (HG02769, HG02953, HG03456, NA19650, NA19705, NA19983), and the
  amended profile was selected there. Every other family served as **training panel** in M8-B prior scoring and in every
  LOO panel. **None is a clean holdout.**
- **ACOT.** 18 of M8-B's 30 outer-evaluation families have never been genotyped: ST016, NA20809, NA21110, NA20806, NA18943,
  PR04, CLM31, GBR002, PK35, VN046, CLM12, BB46, GB23, GB31, ST232, GB21, 2490, NA19468.
  - **Disclosed:** their assemblies were used in M8-B masked prior scoring.
  - The 12 others were genotyped in M8-C/D, M9-C..G and M9-H, and are development.

**Already-inspected cases** (C4 NA19650 and HG00096, CYP2D6 NA19650, ACOT L2/L4 and every other M-series cell) **are
development, never "held out".**

**Not frozen.** No campaign is frozen until review.

## 5. Proposed next experiment (one): M11-B, ACOT family-disjoint calibration data, 5×

**Design:**
- **Families:** the 18 ACOT validation candidates.
- **Hash split:** fixed sha256("M11-B:" + family), into **9 development and 9 sealed validation**.
- **Donor per family:** a deterministic hash rule among complete diploid donors.
- **Configuration:** the unified M10-A configuration (D), LOO panel (every outer-evaluation template excluded, as M9-C/H),
  5×, existing read simulation and subsampling.

**Measured per run:**
- §1 groups A–C;
- O1 on the chain interior, or a named framed prefix;
- O2 per block.

**Outputs:**
- **For the 9 development families only:** a descriptive association between `e` and O1 (rank order, the distribution of
  O1 by `e`, the assessability rate). **No map is fitted** until Davide approves the outcome scale.
- **The 9 validation families:** genotyped in the same run, but **O1 and O2 are not computed** (the truth comparison stays
  sealed) until the calibration rule is frozen.

**Runtime** (from measured M9-H costs on the same LOO panel):
- simulation and subsampling: under 1 s per family;
- LOO 5× call: 41–127 s (mean about 75 s);
- truth projection: about 4.5 s;
- O1: about 2 s;
- O2: about 21 s.

Total about **25–30 min** for 18 families. Declare 2,400 s, with a 300 s per-call watchdog and a 6 GiB envelope.

**Minimum missing inputs:**
- simulated reads for the 18 families (existing `genotype_probe simulate`; none exist yet);
- no other input is missing (`acot.gfa` and the LOO curation exist).

**Decisions needed from Davide before freezing M11-B:**
1. **The calibration outcome scale.** Expected error rate or QV, or the probability of a QV band. Reuse the COSIGT bands
   (`<17 / 17–23 / 23–33 / >33`) explicitly, or the QV/QV_max quintiles?
2. **The single truth-free input.** Proposed: `e = 1 − decision_set_probability`.
3. **The evaluated locus.** The chain interior (flanks excluded), with the named framed-prefix rule when truth is not
   projectable.
4. **Whether O2** (per-block excess) remains a reported secondary outcome.
5. **The validation set.** Accept ACOT's 18 families, panel-exposed through M8-B, as the family-disjoint validation, and
   accept **development-only** evidence for C4 and CYP2D6, where no clean holdout exists.
6. **Permission to simulate reads** for those families; coverage 5× only in the first pass.

## Housekeeping: corrected M10-D wording (no rerun)

M10-D scored **10 prediction rows over 5 cases**: **8 caller predictions** and **2 saved alternatives**.
- **Caller predictions:** C4 NA19650 reported; C4 HG00096 call; CYP2D6 5× and 30×; ACOT L4 and L2, JP on and off.
- **Saved alternatives:** C4 saved 2+2 and CYP2D6 saved 2+3, from M1-B.
- **Scored on the chain interior:** 7 rows (C4 ×3, ACOT ×4).
- **Scored on the named framed prefix sv:1..sv:15 only:** 3 rows (CYP2D6), whose interior is UNSUPPORTED.

"Seven development cases" and "all seven development predictions" were inaccurate. The correction is appended to M10-D's
RESULT.
