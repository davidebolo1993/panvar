# M4-B: default-off distance-prior integration and bounded fresh-search validation (frozen before any unit)

**Ancestry.** Branch `m4-b-distance-prior` from `9ed5a9e` (M4-A). The M4-A uniqueness qualification is appended to
`tests/analysis/m4a/CONTRACT.md`.

**Baseline binary.** `a728eb71…` (the M3/M4-A caller), copied to `results/bin/panvar_baseline` before any library change.

## Implementation (default off; constant behaviour unchanged)

### Adapter (`genotype_locus`)

- **`LocusGeometry`** is computed inside `build_locus_input` against the **frame**:
  - the frame is every complete-frame template's allele after the divergence trim, curated templates included, so it
    is fixed before curation;
  - the reference path's allele per block gives the reference axis;
  - cores are trimmed against the frame's alleles;
  - positions are exact doubled integers.
- **Neutral states** (a frame with one distinct sequence) have an empty core at the block midpoint.
- **Unavailable geometry** carries an explicit reason: no reference path, not a graph path, no consistent frame,
  antiparallel, or partial frame.
- **Digests:** a frame digest and a geometry digest.
- **`compute_locus_geometry`** is the rule, exposed for gates.

### Library: `copying_transitions`

- **One builder,** `build_copying_transitions(blocks, geometry, spec)`, with two modes.
  - `constant_v1`: the legacy vector, bit for bit.
  - `distance_ref_core_v1_experimental`: r_b = −expm1(−d_b/λ), clamped to [0, 1]; d = 0 gives r = 0, and saturation
    gives r = 1, i.e. T = U.
- **Validation:** every invalid input throws `UnsupportedCopyingTransition` (λ ≤ 0 or non-finite, constant r outside
  (0, 0.5], missing or unavailable geometry, block-count mismatch, unknown mode).
- **Recorded:**
  - SHA-256 vector digest;
  - zero-distance inventory, total and informative–informative;
  - `hazard_matched_copying_length`.

### Caller

- **Profile keys:** `copying_transition` (`constant_v1` default | `distance_ref_core_v1_experimental`) and
  `copying_length_bp`. They are recorded in settings and status only when set, so the default output is unchanged.
- **One vector per run** (`g_transitions`), built right after the model, feeds every supported stage:
  - early incumbent, proposal beam, union ranking (sequence priors) and path priors;
  - refinement (`RepairSearchOptions::switch_vector`), seeds, tie branches and the start check;
  - the representation guard;
  - dosage shadow and integrate (`switch_vector`).
- **Constant mode** passes the scalar exactly as before, so the repair setup text and cache identities are unchanged.
  Distance mode appends `;transitions=<digest>`.
- **Refused, never reverting to constant transitions:**
  - distance mode with `switch_probability` set by a profile or flag (the preset's default is not a conflict);
  - distance mode without `copying_length_bp`, or `copying_length_bp` under `constant_v1`;
  - distance mode with joint search, the refinement replay, diagnostic tie branching, or a prepared locus (in-process
    or loaded index: no geometry).
- **Unavailable geometry** writes MODEL_INCOMPLETE with the reason `unsupported copying transition: …`.
- **Status rows (only when set):**
  - `copying_transition`, `copying_length_bp(_hex)`, `copying_redraw_measure`;
  - `locus_geometry_{status,reference,digest,frame_digest,distances_bp}`;
  - `transition_vector_sha256`, `transition_r_{min,max}`, `transition_zero_distance_*`.

### Probes

- `score-set` takes `--copying-transition` / `--copying-length-bp` through the same builder. It refuses
  `--switch-probability` under distance mode, and refuses when the run's recorded mode or vector digest differs.
- `prior-geometry` reports the adapter geometry, its digests, the zero-distance inventory and the hazard-matched length.

## Gates before real runs (manifest G)

| unit | what it checks |
|---|---|
| `gate.lib` | `genotype_copying_transitions_test`: M4-A G1–G10 re-run against the shared rule and builder with the brute-force oracle; N1 validation and refusal; N2 constant bit identity; N3 saturation and small-distance accuracy; N4 zero-distance exclusion (two informative states at one reference point exclude a recombinant: −inf under distance mode, finite under constant r; no floor); N5 digests; N6 hazard-matched formula |
| `gate.e2e` | `genotype_copying_transitions_e2e_test`: explicit `constant_v1` = default apart from record rows; distance completes and records; 1/4-thread agreement; all refusals; same-graph curation keeps geometry while curated alleles go inactive; adapter geometry = independent recomputation; a reverse-written graph keeps the distances; insertion/deletion locus |
| `gate.m4a_prototype` | the M4-A prototype, still passing |
| `gate.regress` | existing ctests for prior, beam, repair, adapter, CLI, curation, dosage, report, prepared-locus, v6 profile and source tree, all in constant mode |
| `gate.cli.*` | CLI refusals: distance + `--switch-probability`; distance + profile `switch_probability`; `copying_length_bp` alone |

**Real-data gates** (manifest C/E; listed beside the cells below):
- constant-mode identity, new binary vs baseline, on HG00171 5× and `three_switches_1` 2×;
- caller vs score-set agreement on both arms of HG00171 5×;
- incompatible specifications refused: constant spec on a distance run, λ off by 1e-6 relative, distance plus scalar.

## λ freeze (manifest L), before any matched run

- `geometry.<cell>` runs `prior-geometry` on each cell's own model (graph, reference, curation).
- `lambda.freeze` writes `LAMBDA.tsv`.
  - **Formula:** λ_HM = Σd / ((B−1)·(−ln 0.99)), from the **implemented** adapter geometry, at r = 0.01 (every cell's
    constant switch probability).
  - **Binding:** each value is bound to the model fingerprint and the geometry and frame digests.
  - **Format:** 17 significant digits, round-trip checked against the hex value.
- `LAMBDA.tsv` is committed before manifest C. It is never fitted to accuracy and never taken from rounded M4-A values.
- The comparator checks that each distance run's recorded geometry digest, model fingerprint and λ equal the frozen ones.

## Matched comparisons (manifests C then E); only the prior specification differs

**Cells:**

| # | cell | validated settings | evaluation |
|---|---|---|---|
| 1 | C4 LZO HG00171 5× | v6 campaign profile | lzo |
| 2 | C4 `panel_control_1` 2× | v6 + mosaic curation | mosaic |
| 3 | CYP2D6 LZO HG02769 5× | the P1 "opt" composition of `perf/inputs/profiles/cyp_5.tsv` (its own guards: 900 s, 6 GiB), CYP2D6 curation, `--dosage-mode repair` (as P1) | lzo |
| 4 | ACOT HG00097 5× | the same route with `acot_5.tsv` | lzo |
| 5 | C4 `both_homologues_1` 5× | v6 + mosaic curation | mosaic |
| 6 | C4 `three_switches_1` 2× (the error case) | v6 + mosaic curation | mosaic |

**Arms.**
- **const:** the cell's validated profile, unchanged.
- **dist:** the same profile with any `switch_probability` line removed (only the C4 v6 profile has one), plus
  `copying_transition distance_ref_core_v1_experimental` and `copying_length_bp <frozen λ_HM>`.

**Order and stop rule.**
1. Cells 1–5 are the controls. Each runs const then dist.
2. `control.check` freezes the control predictions (digest) before any evaluation.
3. Controls are evaluated; then `control.gate`. It reuses M3-F's stop rule: a control that becomes HARMFUL or MIXED,
   refuses only in the distance arm, or is downgraded in status stops every further real run.
4. Then cell 6.
5. `compare` freezes all predictions before evaluation.
6. Cell 6 is evaluated, then `accuracy` runs.

**Reported per cell, kept separate:**
- **Predictions:** selected unordered pair and block projections; winner F/E/prior/score.
- **Search:** repair trajectory, ties and pool sizes.
- **Accuracy:** sequence excess and SV and context exactness; phase.
- **Coverage:** candidate coverage (pool membership).
- **Outcome:** refusals and no-calls, kept in every denominator; whole-caller wall time and RSS.

**Changed calls (manifest D, truth-free).** `decomp.<cell>` scores both winners under both specifications. It separates
the change into likelihood, prior under each specification, and pool membership: PRIOR_RESCORING vs POOL_CHANGE.
Neither class is read as better or worse.

## Budget

- 1,200 s counting every attempt, 4 threads, one process at a time.
- 6 GiB execution envelope. Each profile's own guards are unchanged: C4 v6 4 GiB / 3,600 s; CYP2D6 and ACOT base
  profiles 6 GiB / 900 s. These guards are not the envelope.
- A unit is refused when the time spent plus its declared timeout exceeds the budget.
- Builds are not counted, as in earlier batches. Inputs are hashed in INPUTS.tsv.

## Amendment 1 (declared after `gate.lib` and before any further unit; 0.5 s spent)

**What failed.** `gate.lib` exited 1 on one row, `G4_shared_affix carved_mono_block:distances_compose`, which printed
1.5 against an expected 0. Every other row passed, including the same case's sequence-level prior comparison
(`carved_mono_block:max_abs_prior_difference`, ≤ 1e-12).

**Cause: an assertion defect.** The carved neutral block precedes b0. The correct property is that the original blocks'
distances are unchanged: tc[1] = t[0] and tc[2] = t[1]. The assertion instead required tc[1] = t[0] − tc[0]. The neutral
state's own boundary enters no prior.

**Fix.** The assertion is corrected as row `A1:carved_mono_block:original_distances_unchanged`. Nothing else changes in
the test or the library.

**Units.** Manifest A1, `gate.lib.a1` (60 s). `gate.lib` stays recorded as failed.
