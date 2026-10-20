# M8-C result: r = 0.1 vs r = 0.01 is non-discriminating on two ACOT donors; no call changed, the LOO error is shared, and one donor's LOO timed out in both arms

**Run.** Branch `m8-c-copying-challenge`: freeze `8caf41f`, Amendment 1 (subsample directory), Amendment 2 (missing
outputs in the analysis). Results in `tests/results/m8c_copying_challenge_20261020/`.

**Budget.** 828.5 of 1,200 s.

**Recorded failures:**
- `prepare`: a script defect, re-run as `prepare.a1`;
- `call.HG01975.LOO.r0.01` and `call.HG01975.LOO.r0.1`: **TIMEOUT** at the declared 120-s unit bound, in both arms;
- `checks`: exits 1 only because those two runs have no outputs;
- `eval.HG01975.LOO.*`: exit 1, no call to evaluate;
- `contrast` and `summary`: a missing-file defect, re-run as `.a2`.

**Scope.** A mechanistic development test on two donors, not release validation. No default promotion, fitting, C4/CYP2D6
campaign, distance mode, merge or push.

## Design as executed

- **Donors:** the first two complete diploids from M8-B's 30 ACOT outer-evaluation families, by sha256("M8-C:" + family):
  **HG03209** (family SL17) and **HG01975** (PEL25).
- **Panels:** M8-B's frozen training bank, as curation on the fixed full graph `acot.gfa`.
  - **LOO:** every outer-evaluation template excluded.
  - **LZO:** the bank plus the donor's two haplotypes.
  - No other outer-evaluation donor enters any panel.
  - The frame is derived from the full 467-template panel; this is disclosed.
- **Reads:** the current build's simulator, re-verified by regenerating the archived C4 NA19650 30× source byte for byte;
  30× (seed 20260916, salt = donor) subsampled to 5×. **Byte-identical reads and library across each donor's four runs.**
- **Profiles:** the P1 composition of `acot_5.tsv` (validated combined route, 4 threads, base 6 GiB and 900-s guards).
  The arms differ only in `switch_probability` (0.01 / 0.1), with `copying_transition constant_v1`.
- **Executables:** panvar `98f339e8…`, genotype_probe `af2ed281…`.
- **Arm check:** in every completed run, every recorded `switch_probability` (settings, mosaic status, pre-repair mosaic
  status) equals the arm's r.
  - No transition digest rows are written for constant mode.
  - That the single run vector also feeds dosage and repair rests on M4-B's implementation (source-based), not on a
    recorded field.

## Results (every cell in the denominator)

| cell | call | haplotypes (sequence-based membership) | floor / distance / excess | exact blocks (SV) | phase | seconds / peak GiB |
|---|---|---|---|---|---|---|
| HG03209 LZO r0.01 | strict_call | **truth templates HG03209#1, #2** | 0 / 0 / 0 | 19/19 (9/9) | PRIOR_ONLY | 54.9 / 1.56 |
| HG03209 LZO r0.1 | strict_call | **truth templates HG03209#1, #2** | 0 / 0 / 0 | 19/19 (9/9) | PRIOR_ONLY | 53.8 / 1.71 |
| HG03209 LOO r0.01 | strict_call | panel members HG03225#2 + HG01496#2 | **9 / 14 / 5** | 10/19 (7/9) | PRIOR_ONLY | 53.9 / 1.94 |
| HG03209 LOO r0.1 | strict_call | **the same** HG03225#2 + HG01496#2 | **9 / 14 / 5** | 10/19 (7/9) | PRIOR_ONLY | 54.8 / 1.62 |
| HG01975 LZO r0.01 | strict_call | **truth templates HG01975#1, #2** | 0 / 0 / 0 | 19/19 (9/9) | PRIOR_ONLY | 58.2 / 2.10 |
| HG01975 LZO r0.1 | strict_call | **truth templates HG01975#1, #2** | 0 / 0 / 0 | 19/19 (9/9) | PRIOR_ONLY | 50.9 / 2.01 |
| HG01975 LOO r0.01 | **no call: TIMEOUT** (declared 120 s; the run was still in tie branching, with partial tie-branch and repair memory files) | — | not assessable | — | — | 120.3 / — |
| HG01975 LOO r0.1 | **no call: TIMEOUT** (same) | — | not assessable | — | — | 120.2 / — |

**Run details:**
- **Projection ambiguity:** none flagged.
- **Spellability:** all 19 blocks defined in every evaluated run. Truth is spellable from the panel in LZO; in LOO it is
  not (`truth_candidate_source = absent`), and the floor is not attained in the candidates.
- **No off-panel predictions:** every called haplotype in every arm is either the donor's own truth template (LZO) or a
  single calling-panel template (LOO). Each called haplotype has 0 switches.
- **Search work** is identical between arms: HG03209 LZO 84,533 sequence candidates, 128 mosaic candidates, 178,760,311
  pruned states; HG03209 LOO 83,719 / 128 / 177,353,079; HG01975 LZO 84,524 vs 84,525 / 128 / 178,792,784.

## Objective vs search (§5)

**No reported call changed in any assessable donor × scheme** (`contrast`: three UNCHANGED, one incomplete), so there is
no changed-call contrast to decompose.

**Same call, same fragment and exposure terms; only the prior differs:**

| cell | prior at r = 0.01 | prior at r = 0.1 |
|---|---|---|
| HG03209 LZO | −11.69 | −15.11 |
| HG03209 LOO | −11.69 | −15.11 |
| HG01975 LZO | −11.69 | −15.11 |

A zero-switch pair is penalised by the larger redraw probability.

**The visible effect is on confidence, not the call.** In HG03209 LOO the winner's numeric separation gap fell from **5.06
to 2.63 nats**, and its candidate-set probability from **0.993 to 0.911**.

## Decision

- **Non-discriminating.** Across every assessable cell, neither arm made an error that the other avoided:
  - both LZO controls are exact in both arms;
  - HG03209 LOO carries the same 5-bp excess over its floor in both arms;
  - HG01975 LOO produced no call in either arm, because of the declared 120-s unit timeout, not a caller guard.
- **No LZO deterioration:** r = 0.1 introduced no off-panel or mosaic prediction in either control.
- **No LOO improvement either:** the shared HG03209 error is not touched by the parameter, since the winner and the search
  are identical.
- **Not established:** that r = 0.1 helps or harms genotyping. The masked-score preference from M8-B does not translate
  into a different call here.
- **Limitations:**
  - two donors, one coverage;
  - two of eight runs are resource failures under the declared unit bound;
  - a full-panel frame;
  - all calls were zero-switch single-template haplotypes, a regime where r affects only the prior's level and the
    confidence margin.
