# M11-C result: confidence inputs prepared; Q1 semantics corrected; default-off block diagnostics gated; offline QVfrac computed; calibration design simplified

**Run.** Branch `m11-c-confidence-inputs`; freeze commit `edcb0e6`; result after.
- **Results:** `tests/results/m11c_confidence_inputs_20261020/`.
- **Budget:** 410.7 of 900 s, including 150 s of disclosed development.
- **Units:** 19 of 19 passed first time.
- **Watchdog:** no timeout or envelope event. The peak was 3.08 GiB, in the ACOT smoke run.
- **Not done:** no simulation, calibration, fitting, threshold or prediction change.

## 1. Q1, corrected

`genotype_support` is `READS` iff no **fragment contender** carries a different unordered genotype at the block.
- **Fragment contender:** a **ranked rival** whose fragment term (the read likelihood only, without exposure or prior) is not
  certifiably below the winner's, i.e. `rival.fragment + bound ≥ best.fragment − bound`.
- **The genotype compared:** read from the contender group's representative labelling.

**So PRIOR_ONLY means the reads alone do not separate the call from at least one explored rival with a different block
genotype.** It does **not** mean "no fragment informs the block".
- **Both values are conditional on the explored pool,** and depend on representative labellings.
- **M11-B's association results stand;** its wording is corrected (appended to `RESULT_M11B.md`).

## 2. Diagnostic outputs (default off; `block_diagnostics diagnostic_v1`)

`<prefix>.block_diagnostics.tsv` has one row per model block: `census_route` and `jp_emission_tie`, with explicit
NOT_EVALUATED / UNAVAILABLE values.

**Code (output only):**
- `JunctionProposalReport::block_eligible` copies the per-block eligibility the search already computed, merged by OR
  across searches;
- `DosageShadowOptions::census_routes_out` receives the routes immediately after `run_census`;
- the table is written after the dosage stage;
- the file is run-owned;
- the settings row appears only when the key is declared.

**Gates:**

| gate | result |
|---|---|
| build; JP1.1, compound and sequence-report units; e2e `sequence_first` | pass |
| off/on identity, C4 5× (unified D configuration) | **science identical**; only the diagnostic file and the declared setting differ |
| off/on identity, ACOT 5× (D, JP1.1 on) | **science identical** |
| off/on identity, CYP2D6 NA19650 5× (M1-B command, JP off) | **science identical** |
| block association | one row per model block, in model order, in all four on-runs |
| census retained past audit failure | CYP2D6: the audit **FAILED** ("unit library: alleles disagree on their flanks … in sv:9"), and every block keeps its route; sv:9 = `repeated_unit` |
| census not run | `--dosage-mode off`: every block `NOT_EVALUATED:dosage_audit_off`; no route inferred |
| JP condition | ACOT: exactly sv:5 `TIED_ELIGIBLE`, equal to the status count (1). C4: 0 eligible, equal to the status. CYP2D6 (JP off): `NOT_EVALUATED:junction_proposals_off` |
| determinism | the C4 table is byte-identical on a repeat run |
| C4 route | sv:7 = `repeated_unit` (as in the saved dosage tables) |

## 3. Offline QVfrac (no caller change)

**Scope:** 514 truth-defined blocks of the current-schema development runs (M1-B, M9-C, M9-G JP-on, M9-H), plus 2 with
undefined truth.

**Truth and alleles:** from `genotype_probe block-sequences` with each evaluation's own model, curation, truth graph and
donor. 15 truth units; every fingerprint matches.

**Consistency:**
- the recomputed predicted distance equals the evaluator's unordered distance in **514 / 514** blocks;
- the best-pair edits equal the evaluator's closest-allele distances in **514 / 514**.

So QVfrac uses exactly the evaluator's floor pairs.

| measure | value |
|---|---|
| QVfrac at the floor (388 blocks) | 1.0000 (2 blocks marginally > 1: equal edits, different alignment lengths) |
| QVfrac above the floor (126 blocks) | min 0.105, median 0.939, max 0.999 |
| availability-relative QVfrac quintiles (0–0.2, …, 0.8–1.0) | 2 / 3 / 8 / 14 / **487**. They barely resolve these data |
| rank agreement with normalized excess (above-floor pairs) | 0.846 |
| summed-QV pairing ≠ minimum-edit pairing | 6 blocks (M1-B NA19650 ctx:2; M9-C/G HG03831 sv:5 and ctx:16; M9-H HG02723 sv:9 ×2). Reported, **not** substituted |

**What QVfrac is not:**
- not QV / perfect-match ceiling;
- QV-maximisation among equal-edit best alleles is **not** evaluated, because the evaluator's tie rule fixes the best allele
  and anything else would need a new floor search.

**Legacy strata** (C4 production 483 runs, CYP2D6 pilot 27): QVfrac was **not** computed. **What it needs:**
- per run, a truth projection with the evaluation's own model, curation and truth graph (`block-sequences`);
- six alignments per truth-defined block.

The evaluator's closest-allele identities are already saved, so **no new floor search** is needed.

## 4. Simplified calibration design (supersedes CALIBRATION_DESIGN.md §2's stratified scheme)

**Two pre-declared candidates:**
- **C1:** `error_mass_upper` alone;
- **C2:** C1 + Q1 (`genotype_support`, with the corrected meaning).

Q2 (census route) and Q3 (emission tie) stay **diagnostic only** at first. There is no three-way stratification or merge
scheme.

**Outcomes:**
- **primary:** per-block excess and normalized excess over the existing floor;
- **secondary / descriptive:** QVfrac and floor attainment.

**Recorded historical interpretation (Davide):**
- COSIGT QV categories `<17 / 17–23 / 23–33 / >33`;
- availability-relative QVfrac quintiles (not QV / QV_max).

Any categorical acceptance rule is a **proposal pending approval**.

**Not done or claimed:**
- no model is fitted, and calibrated confidence is not claimed;
- the small, ACOT-only validation cannot establish cross-locus reliability or a very low confident-error rate.

**Decisions still needed:**
- the acceptable block outcome (excess and/or normalized excess, or QVfrac class);
- the acceptable confident-error rate;
- validation acceptance and read simulation;
- whether the two diagnostic outputs may be enabled in the development runs.
