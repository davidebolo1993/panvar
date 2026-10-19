# M3-H result: the ctx:0 comparison shows no inconsistency, only sampling at 2×

**Run.** Branch `m3-h-evidence-decomposition`: freeze `4a36b9c`, M3-G qualifications `7f02018`. No model or caller
change. Results in `tests/results/m3h_evidence_decomposition_20261019/`.

**Budget.** 10.3 of 300 s. Seven units, all exit 0 and first-time, with no amendments.

**Matched pair, frozen in DESIGN.md.** A = R_t0, which carries (14,56) at ctx:0 and is the truth path. B = R_c0, which
carries the called (4,52). Both are in the repair-incumbent context, assignment 0.
- **Identity:** `spell` confirms that R_t0 spells the truth pair, and that the two sequences are identical outside
  ctx:0 on both homologues.
- **Same sequences at both coverages:** the same two sequence pairs are scored against the 2× and 5× reads, each with its
  own frozen library. The libraries differ only in fragment start intensity.
- **Nesting:** verified. All 2,145 2× pairs are in the 5× pool with identical bytes, and 5× adds 3,229.

## Checks (all pass)

- **Identity expectations:** 2× R_c0 equals the saved caller value and R_t0 the M3-G value; 5× R_t0 equals the saved 5×
  equivalence-set member.
- **Per-pair sums:** equal the scorer's fragment terms at both coverages.
- **Canonical-origin reconstruction:** exact (difference 0); `placement_evidence_interpretable` is true.
- **Exposure difference:** 0 at both coverages.
- **Prior gates:** the independent forward/backward over the 131 curated templates reproduces the sequence-level pair
  prior within 7e-15. `prior-pool` re-derives all 12,754 saved pool priors (largest difference 5e-11) and gives the same
  two priors.

## Score differences (A − B)

| | 2× | 5× |
|---|---|---|
| fragment | **−1.332** | **+27.986** |
| exposure | 0 | 0 |
| prior | −10.164 | −10.164 |
| total | −11.496 | **+17.82** |

## Read evidence

| | 2× | 5× |
|---|---|---|
| pairs | 2,145 | 5,374 |
| exactly neutral | 2,134 | 5,341 |
| non-neutral favouring truth / called (abs Δ ≥ 0.01) | 4 / 6 | 18 / 12 |
| placement class carrying Δ | ctx:0 interior: 122 pairs, all of Δ; junction and outside: 0 | interior: 290 pairs, all of Δ |
| concentration | the top 10 pairs carry 99.9% of the summed abs Δ | top 4 = +25.7 of +28.0 |

**What the 2× reads are.** All ten 2× contributors are worth about ln 2 (0.693; one is 0.635).
- Each read places with **0–1 mismatches** at one locus per haplotype, with a proper insert (300–423 bp; library
  350 ± 40).
- It matches **both haplotypes** of one candidate but **only one** of the other. Under the diploid mixture, that is a
  copy-number factor of 2.
- 3 such reads favour the truth and 5 the called pair (plus three small ones), which nets about −2·ln 2 ≈ −1.33.

**Read detail.**
- **Base qualities:** uniformly Q40 (simulated); the scorer uses a fixed substitution rate.
- **Placement ambiguity:** each read has one gapless locus per haplotype, and the scorer's origin counts are consistent.
- **Discriminating positions:** at ctx:0, 61–65 discriminating 31-mers per homologue, in 2–3 clusters.
  - hap1: block offsets about 274–336 and 10,816–10,879.
  - hap2: about 7,585–7,629, 10,982–11,033 and 11,471–11,534.

**Nesting.** At 5×, the shared 2,145 pairs give exactly the 2× delta (−1.332; the largest per-pair change is 0, with no
sign flips). The 3,229 added pairs give **+29.318**.
- Their largest contributors (+8.05, +7.73, +6.62, +3.30) place, with 0–2 mismatches, **only on the truth's hap2
  (allele 56)**, or fit it better by one mismatch.
- These are reads from 56-specific sequence. None of the 2× subsample's 2,145 pairs is such a read.

**No inconsistency demonstrated.** The frozen flag found 0 contributors scored against the candidate they match better.
The ±ln 2 reads are a correct property of the diploid mixture, not errors. The 2× deficit is the balance of about ten
shared-sequence reads. The decisive, candidate-specific reads appear only when more molecules are sampled.

## Prior (sequence-level, summing compatible histories)

The prior uses the current copying model: r = 0.01, a **donor-redraw parameter of the copying model**, not an
established biological recombination rate. No alternative spelling contributes, as the gate shows.

| homologue | truth vs called at ctx:0 | compatible donors | continue into sv:1 without a redraw | expected donor changes | Δlog prior |
|---|---|---|---|---|---|
| hap1 | 14 vs 4 | 2 (HG02769#1, NA19434#1) vs 6 | 0 vs 0 (both need a redraw) | 3.01 vs 3.01 | **−1.099** = log(2/6), donor compatibility only |
| hap2 | 56 vs 52 | 1 (HG02769#2) vs 2 (HG00731#2, HG02011#2) | 0 vs **2** | **3.00 vs 2.00** | **−9.065**: one extra copying transition plus fewer donors |

- **Total:** −10.164.
- **Where the hap2 cost comes from.** The truth's hap2 copies ctx:0 from a donor that does not carry the downstream
  alleles. The copying model charges that as a donor change, which is what a simulated mosaic switch at that boundary
  looks like.
- **Scope.** This is not evidence for or against a particular value of r. No grid or distance rule was evaluated.

## Decision: **NO_INCONSISTENCY. Retain as an unresolved low-depth error.**

At 2× the read evidence between the two candidates is about ten copy-number-type reads, balanced 4:6 by sampling. The
prior adds one copying transition against the truth. The same two sequences scored against the nested 5× reads reverse
the likelihood (+28.0) through candidate-specific reads absent from the 2× subsample. Nothing in the decomposition
indicates a likelihood inconsistency. One coverage pair does not establish calibration.

## Recorded for later (not addressed here)

1. **Unrefined pool winner.** The final pooled call (a union candidate) is never refined. A one-block move inside it
   scored +0.613 above the call in M3-G. Trunk convergence does not establish the local optimality of the reported call.
   It is a separate search limitation; in this cell the move goes to another wrong pair.
2. **Distance-dependent prior and flank questions.** The prior's constant per-boundary redraw parameter is not
   distance-aware, and here it decides about −9 of the −10.2 at the ctx:0/sv:1 boundary. Whether redraw cost should
   depend on block span or flank context remains open (M3-A's distance-matched model was measured, not adopted). No
   grid or rule is proposed or implemented in this batch.

## Closure appended at M4-A (2026-10-19; Sol's review; the record above is not rewritten)

**The ctx:0 outcome stands:** no inconsistency demonstrated, retained as an unresolved low-depth error. Joint search
stays parked and default-off.

Recorded separately and not reopened in M4-A:
1. **The final pooled winner may never have been refined.** A union candidate can be reported without any refinement
   from it.
2. **Conditional data support is not calibrated user confidence.** The support labels (v1x/v2x) are conditional on the
   evaluated candidates.
3. **The C4 sv:7 structural error is a different model/content problem.** The fragment term decides it on content the
   leave-one-out panel cannot spell (M1-B/M1-C). It is unrelated to the copying prior or to this ctx:0 case.
