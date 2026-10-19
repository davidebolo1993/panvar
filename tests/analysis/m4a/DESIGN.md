# M4-A: distance-aware copying-prior contract and invariance gates (design and diagnostic; frozen before execution)

**Ancestry.** Branch `m4-a-distance-prior` from `af96ddf` (M3-H), plus the M3-H closure (`7f14009`).

**No production change.**
- No caller integration, index format, default, campaign, simulation or tuning.
- New code is diagnostic only:
  - the probe subcommand `prior-geometry` (read-only model geometry);
  - the test `tests/m4a_distance_prior.cpp`, which holds the prototype construction (not in `src/`);
  - scripts in this directory.
- M3-A's grid and fixed-pool comparison are **not** rerun. M3-A's representation test runs once, as an unchanged
  regression (`gates.m3a_regress`).

## The contract under test (fixed before any unit runs)

**Transition.** T_b = (1 − r_b) I + r_b U.
- r_b = 1 − exp(−d_b / λ), where λ is a fixed copying length in bp.
- U is the declared redraw measure: uniform over the n curated template rows, the current one included.

**Geometry: reference-coordinate core midpoints.**
- Each block's copying state sits at the midpoint of the block's **polymorphic core** on the reference axis. The core
  is the reference path's interval for the block along the canonical chain, with the maximal prefix and suffix shared
  by every allele of the **frame** removed.
- d_b = core_mid(b) − core_mid(b−1).
- A monomorphic block has an empty core.

**Scale.** λ is one length. It is never re-derived from the block count. The hazard-matched length is a diagnostic
control only.

## Units (manifest A)

| unit | what it runs | what it checks |
|---|---|---|
| `gates.m4a` | `tests/m4a_distance_prior.cpp`; oracle = mosaic-core `brute_haploid_prior`; λ and U frozen in every comparison | G1 valid distances and transitions; G2 zero distance (T = I); G3 neutral monomorphic insertion at every split point (constant-r control); G4 shared affix written in the block, as invariant sequence, or as a monomorphic block (controls: untrimmed midpoints, M3-A panel-mean); G5 reverse-chain orientation; G6 polymorphic split classified NON_NEUTRAL; G7 sequence-level summation over equivalent spellings; G8 donor multiplicity and order, redraw semantics; G9 donor removal (frame cores vs callable cores vs panel-mean); G10 scale under subdivision (hazard-matched control) |
| `gates.m3a_regress` | M3-A's representation test | unchanged regression |
| `geometry.<model>` | `genotype_probe prior-geometry` on 8 models: C4 full; C4 LOO NA19650; ts1 staged with and without mosaic curation; CYP2D6 with and without curation; CYP2D6 LOO NA19650; ACOT | Gates: the reference projection tiles the chain, and agrees with the reference template when it is one. Output: per block, reference interval and core; per boundary, all five candidate distances. |
| `geometry.compare` | per-model totals, concentration, ratio to M3-A panel-mean span, hazard-matched control length | **Pairs:** curated vs uncurated on the same graph (ts1, CYP2D6); full vs rebuilt LOO graph (C4, CYP2D6). |

**The five candidate distances:**
1. panel-mean span (M3-A);
2. reference start-to-start;
3. reference block midpoints;
4. reference core midpoints;
5. invariant gap only.

**Note on the probe.** It trims cores against the model's **callable** alleles. The frame core is measured as the
uncurated model of the same graph. The gap between the two is what the G9 control predicts.

## Outputs after execution

- `CONTRACT.md`: the contract, geometry and scale recommendation, unresolved choices and references.
- `INTEGRATION.md`: the consumer inventory, the shared construction and provenance proposal, and the validation set.
- `RESULT.md`: gate and geometry results.

## Bounds

- 600 s counting every attempt (development runs counted).
- **Declared timeouts:** gates 60 s each; geometry 60 s each; compare 30 s.
- Builds are not executions and are not counted, as in earlier batches.
- Inputs are hashed in INPUTS.tsv.
