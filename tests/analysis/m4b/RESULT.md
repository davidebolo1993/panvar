# M4-B result: default-off distance-prior integration and six matched fresh-search comparisons

**Run.** Branch `m4-b-distance-prior`:
- freeze `320ffea`;
- Amendment 1 (G4 carved-case assertion) and Amendment 2 (E6 reverse-graph expectation);
- λ freeze `04c365d`.

Results are in `tests/results/m4b_distance_prior_20261019/`.

**Budget.** 506.1 of 1,200 s, every attempt counted, at 4 threads, one process at a time.

**Recorded failures.** Four units failed and stay recorded:
- `gate.lib` and `gate.e2e`: test defects, fixed by the amendments;
- the two ACOT evaluations: they exceeded their declared 60 s timeout.

## Integration

**Adapter.** `LocusGeometry` is computed inside `build_locus_input`:
- against the **frame**, meaning every complete-frame allele after the divergence trim, curated templates included,
  before curation;
- with the reference path's allele per block;
- with exact doubled positions.

Unavailability is explicit, with reasons. The record holds the reference identity, a frame digest and a geometry digest.

**Library.** `build_copying_transitions` is the one builder:
- `constant_v1` gives the legacy vector bit for bit;
- `distance_ref_core_v1_experimental` gives r = −expm1(−d/λ), clamped to [0, 1];
- invalid inputs raise `UnsupportedCopyingTransition`;
- it records a SHA-256 vector digest and the zero-distance inventory, and provides `hazard_matched_copying_length`.

**Caller.** One vector per run feeds every supported stage: early incumbent, proposal beam, union sequence priors, path
priors, refinement, seeds and tie branches, the start check, the representation guard, and dosage shadow and integrate.
- **Constant mode** keeps its scalar path and cache identity.
- **Distance mode** adds its digest to the repair setup key.
- **Refused explicitly:** a scalar set by profile or flag; a missing or stray `copying_length_bp`; joint search; the
  refinement replay; diagnostic tie branching; a prepared locus (no geometry).
- **Unavailable geometry** gives MODEL_INCOMPLETE `unsupported copying transition: …`.
- **Status rows** are written only when the key is set.

**Probes.** `score-set` uses the same builder. It refuses on a mode or digest mismatch with the run, and refuses a scalar
alongside the distance mode. `prior-geometry` reports the adapter geometry and λ_HM.

**Prepared-index format: unchanged.** Distance mode refuses any prepared locus, built in-process or loaded from a v1
file; both reach the same code path. Constant-mode prepared runs are unchanged (E4).

## Gates (all pass after the two test-defect amendments)

| gate | result |
|---|---|
| `gate.lib.a1` | Shared rule and builder against the brute-force oracle. M4-A G1–G10 hold on the shared code. **N1:** 10/10 invalid specs refused explicitly. **N2:** constant vector equal to legacy bit for bit. **N3:** saturation r = 1 (T = U) matches the oracle; small-d relative error 8.5e-12. **N4:** two informative states at one reference point exclude the recombinant (−inf under distance mode, finite under constant r; no floor). **N5:** digests. **N6:** λ_HM formula. |
| `gate.e2e.a2` | 24/24 caller checks. Explicit `constant_v1` equals the default apart from record rows. Distance completes; 1 vs 4 threads identical. All refusals hold. Curation keeps geometry, frame digest and distances while 4 curated alleles go inactive. Adapter geometry equals an independent recomputation, on the valley and on an insertion/deletion locus. A reverse-written graph gives mirrored distances. |
| `gate.m4a_prototype` | pass |
| `gate.regress` | 16/16 existing ctests in constant mode: prior, beam, repair, adapter, CLI, curation, dosage, headline, report, prepared-locus file, v6 profile, source tree |
| `gate.cli.*` | 3/3 CLI refusals |
| constant identity | **HG00171 5×:** 12/12 output files identical; **`three_switches_1` 2×:** 11/11. New binary vs the baseline `a728eb71`, volatile rows excluded. |
| caller vs score-set agreement | HG00171 5×: both arms' saved winning totals reproduce under their own specification; the distance digest matches the run's |
| incompatible specifications | 3/3 refused: constant spec on a distance run; λ·(1 + 1e-6) (digest mismatch); distance plus `--switch-probability` |

**Zero-distance inventory on the real models: none.** All six cell models have 0 zero-distance boundaries (C4 10, CYP2D6
18, ACOT 18 boundaries). The exclusion consequence is demonstrated on the N4 fixture only.

## Frozen λ_HM (`LAMBDA.tsv`)

Frozen from the implemented adapter geometry at r = 0.01 before any run, and bound to the model fingerprint and the
geometry and frame digests.

| cells | λ_HM |
|---|---|
| C4, all four (the curated mosaic models equal `c4.gfa`: same fingerprint `d2b8c6c2e0b4079d`, geometry `fbed6033…`) | 1,872,698.61 bp |
| CYP2D6 | 1,130,346.42 bp |
| ACOT | 1,222,908.28 bp |

Every distance run recorded the frozen geometry digest, model fingerprint and λ (comparator checks: all True).

**Resulting r ranges:** C4 9.8e-5 to 0.028; CYP2D6 7.2e-5 to 0.042; ACOT 6.3e-5 to 0.027. The constant prior is 0.01
everywhere.

## Matched comparisons (controls first; predictions frozen before evaluation)

| cell | status (const / dist) | call | excess | blocks / SV exact | phase | winner prior Δ (dist − const) | pool identities | wall (s) | RSS (GiB) |
|---|---|---|---|---|---|---|---|---|---|
| C4 LZO HG00171 5× | COMPLETE / COMPLETE | identical | 0 / 0 | 11/5 both | PRIOR_ONLY | −0.00005 | 9,232 / 9,231 | 18.5 / 18.5 | 1.33 / 1.34 |
| C4 `panel_control_1` 2× | COMPLETE / COMPLETE | identical | 0 / 0 | 11/5 both | PRIOR_ONLY | −0.00003 | 9,248 / 9,252 | 12.3 / 12.4 | 1.02 / 1.03 |
| CYP2D6 LZO HG02769 5× | COMPLETE / COMPLETE | identical | 0 / 0 | 19/9 both | PRIOR_ONLY | −0.0016 | 9,032 / 9,030 | 8.7 / 8.7 | 0.97 / 0.96 |
| ACOT HG00097 5× | COMPLETE / COMPLETE | identical (sequence, rows and projections) | **not assessed** (evaluator timeout) | – | PRIOR_ONLY | −0.0010 | 109,887 / 109,891 | 31.5 / 31.0 | 2.61 / 2.80 |
| C4 `both_homologues_1` 5× | COMPLETE / COMPLETE | identical | 0 / 0 | 11/5 both | PRIOR_ONLY | **−0.485** | 10,392 / 10,396 | 20.6 / 20.9 | 1.31 / 1.30 |
| C4 `three_switches_1` 2× (error case) | GE / GE | identical (same unordered rows) | **8 / 8** | 10/5 both | UNRESOLVED | −1.79 | **12,754 / 10,441** | 17.7 / 13.4 | 0.99 / 0.98 |

**Summary.**
- **Controls:** the stop gate passed. No control changed its call; there are no recoveries and no harmful changes.
- **Changed calls:** 0 of 6, so every decomposition unit reports UNCHANGED.
- **Search, not just rescoring.** In `three_switches_1` 2× the distance prior changed the refinement trunk without
  changing the call.
  - The first two moves (ctx:2, sv:5) were identical in both arms.
  - The constant arm went on to improve at sv:3 (round 3, +3.92) and continued for 7 rounds and 6 moves.
  - Under distance transitions, round 3 found no improving neighbour, so the trunk converged after 3 rounds and 2 moves,
    at −31889.37 instead of −31871.40.
  - Repair-scored candidates fell from 4,476 to 1,920. The reported call is still the union winner, with the same
    unordered genotype and the same ctx:0 error.
- **Costs** are unchanged within noise. ACOT peak RSS is +0.19 GiB, within its 6 GiB profile guard and the 6 GiB
  execution envelope.

## Limitations

- **λ_HM is a diagnostic control** (block-count dependent). Nothing here selects a physical copying length; that
  needs an independent truth-free study.
- **Small sample.** Six cells, one coverage each. Identical calls on these cells are not evidence of equivalence
  elsewhere.
- **ACOT absolute accuracy was not assessed**: the evaluator exceeded its declared 60 s, costing 120.6 s for both arms.
  The two arms' predictions are byte-identical, so the comparison stands.
- **Constant-mode identity was verified against the baseline binary on two C4 cells.** The CYP2D6 and ACOT constant
  arms were not compared against the baseline; they are covered by the gates and the unchanged code path.
- **Prepared-index refusal** is tested through the in-process prepared locus, the same code path a loaded v1 index
  takes, not by loading a file.
- **Still unresolved:** reference-allele length at structural blocks (G-U1). Zero-distance behaviour is not exercised by
  any real model here.

## Recommendation: **retain default-off as an experimental mode; do not promote, and do not run a wider campaign yet**

- **The integration works.** The geometry is computed in the adapter against the frame. One shared vector reaches
  every supported stage and the diagnostic. Constant mode is unchanged, and every unsupported combination is refused.
- **No accuracy effect yet.** At the hazard-matched control, the distance prior changed no call on five controls or on
  the known error. It moved priors (by up to 0.49 nats on controls) and, in one cell, the refinement trajectory.
- **Next step.** An independent, truth-free copying-length study, so that a physically interpretable λ can be tested
  against this control on the same cells. It is not threshold tuning, and not another search mechanism.
- **Recorded separately for later work:**
  - the unrefined-pool-winner issue (the ts1 2× call is again a union winner never refined);
  - the trunk sensitivity seen here.
