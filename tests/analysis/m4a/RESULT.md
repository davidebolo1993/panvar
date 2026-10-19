# M4-A result: distance-aware copying-prior contract, invariance gates and model geometry

**Run.** Branch `m4-a-distance-prior`: freeze `add2801`, Amendment 1 `087b254` (plus its appended correction),
Amendment 2. Results in `tests/results/m4a_distance_prior_20261019/`.

**Budget.** 12.5 of 600 s. 13 units: 11 exit 0, and two recorded failures, both explained by the amendments.

**No production change.** Nothing changed in the caller, defaults, index format or campaigns.

**Deliverables:**
- `CONTRACT.md`: contract, geometry, scale, unresolved choices, references;
- `INTEGRATION.md`: inventory, shared construction, default-off proposal, validation set;
- this file.

## Invariance gates (`tests/m4a_distance_prior.cpp`; brute-force history oracle; λ and U frozen)

| gate | result |
|---|---|
| G1 validity | 400 boundaries: 0 violations (d ≥ 0, r in [0, 1), T stochastic and non-negative) |
| G2 zero distance | r = 0, T = I. A parental sequence equals the oracle and a recombinant is −inf. Adjacent single-base cores give d = 1, not 0. |
| G3 neutral monomorphic insertion (3 split points) | distances add exactly; priors equal (≤ 9e-16). **Control:** constant r moves the prior by 0.688 nats. |
| G4 shared affix in the block / invariant / own monomorphic block | distances and priors identical. **Controls:** untrimmed block midpoints move (1 bp); M3-A panel-mean moves (2 bp). |
| G5 reverse-chain orientation | distances mirror exactly; reverse-complement priors equal (≤ 5e-16) |
| G6 polymorphic split | 4 new admissible recombinants: classified **NON_NEUTRAL**, never used as an invariance case |
| G7 equivalent spellings | 2 sequences with more than one spelling; library equals summed oracle (≤ 9e-16) |
| G8 donors | A duplicated row leaves geometry unchanged. U counts rows: equals the oracle, and the reweighting is explicit (0.54). Row order has no effect. Actual change probability is r(n−1)/n. |
| G9 donor removal (Amendment 1 fixture) | frame cores: 0 change. **Controls:** callable-allele cores 1.5 bp; panel-mean 0.5 bp. |
| G10 scale | fixed λ: total hazard unchanged by subdivision. **Control:** λ_HM changes by 497 bp. |
| M3-A regression | `m3a_prior_representation` PASS, unchanged |

**Amendment record.**
- `gates.m4a` failed on the G9 callable-core control because of a fixture defect: the removed donor carried no
  core-defining allele.
- `gates.m4a.a1` passes on a fixed fixture. Every non-G9 row is identical between the two runs.
- The original G9 rows remain in `gates.m4a`'s record. The A1 run emitted only the `A1:` rows, a deviation from the
  declaration recorded in the appended correction.

## Model geometry (`prior-geometry`; truth-, read- and candidate-free)

**Projection gates.** For all 7 buildable models, the reference projection tiles the chain and agrees with the
reference template.

**Refusal.** `cyp_uncurated` is unbuildable: NA18989#1 has a partial frame, which is why the CYP2D6 curation exists.
Amendment 2 leaves it and its pair out.

| model | blocks | panel-mean span (M3-A) | ref start | ref block mid | **ref core mid** | invariant gap | λ_HM, core mid (control) |
|---|---|---|---|---|---|---|---|
| C4 | 11 | 164.8 kb | 162.7 kb (1 zero boundary) | 188.2 kb | **188.2 kb** | 0 at all 10 | 1.87 Mb |
| CYP2D6 (curated) | 19 | 212.1 kb | 211.5 kb (2 zero) | 204.4 kb | **204.5 kb** | 0 at all 18 | 1.13 Mb |
| ACOT | 19 | 189.4 kb | 202.9 kb (2 zero) | 221.2 kb | **221.2 kb** | 0 at all 18 | 1.22 Mb |

**The alternating chain has no inter-block invariant sequence.** So the anchor rule decides where a block's length
goes.
- **Start-based anchors** (reference start, M3-A panel-mean) give near-zero distance after an insertion or tiny variant
  block. At C4, sv:5→ctx:6 is 0 bp (start) or 110 bp (panel-mean), and sv:9→ctx:10 is 1 or 2.6 bp. The next boundary
  then carries the whole context block (48 kb, 63 kb).
- **Core midpoints** split each block symmetrically: sv:5→ctx:6 23.9 kb, sv:9→ctx:10 31.7 kb.
- The ratio of core-midpoint to panel-mean distance per boundary ranges from 0.49 to about 12,000 (C4).

**Core vs block midpoints.** They agree within 0.5% on real models; core trimming matters for representation
invariance (G4), not for magnitude.

**Panel composition and LOO:**
- **Same graph, curated vs uncurated (ts1):** every reference geometry changes by 0 bp; M3-A panel-mean changes by up to
  227 bp per boundary.
- **Rebuilt LOO graphs:** C4 NA19650 moves one divergence trim (sv:3 from 166 to 120 bp), so core-midpoint distances
  change at two boundaries (by 23 bp) with the total unchanged; panel-mean changes everywhere (total −69 bp). CYP2D6
  NA19650 also changes the frame (total panel-mean −87 bp; reference totals unchanged).

**Reference dependence at structural blocks.** C4 sv:7 counts the reference haplotype's 59 kb RCCX allele (panel mean
61 kb). This is the declared unresolved choice G-U1.

## Recommendation (detail in CONTRACT.md)

- **Transition:** T_b = (1 − r_b) I + r_b U, with r_b = 1 − exp(−d_b/λ) and U uniform over the curated template rows,
  current one included and multiplicity explicit.
- **Geometry:** reference-coordinate core midpoints, cores trimmed against the frame (adapter output, digested). Refuse
  when the reference coordinates are unavailable.
- **Scale:**
  - (a) λ_HM, a per-model diagnostic control only;
  - (b) one physical κ = 1/λ, chosen later by an independent, truth-free development study (leave-one-template-out
    copying likelihood) on development loci;
  - (c) biological calibration: unestablished.
- **Integration:** one `build_copying_transitions` feeding every consumer, a digest-checked no-mixing rule, and v1
  indexes refused under distance mode. The validation set is in INTEGRATION.md.

**Unresolved:**
1. G-U1: reference-allele length at structural blocks.
2. G-U2: frame cores need adapter support.
3. The value of κ.
4. Whether U should ever be anything other than uniform over rows. Version 1 says no; any change is declared, never
   implicit.
