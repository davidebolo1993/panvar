# M9-B: output-target mapping and confidence-contract closure

This batch defines what a genotype-confidence value describes. It does not implement or calibrate one.

**Scope:**
- **Read-only:** code inspection, saved outputs, and focused fixtures.
- **Not done:** no genotyping rerun, fitting, clustering, threshold, GQ/PASS change, caller change, default promotion,
  merge, push or cleanup.
- **Results:** `tests/results/m9b_confidence_contract_20261020/`.

## Units

| unit | timeout | what it does |
|---|---|---|
| `fixtures` | 10 s | joint versus marginal agreement through M9-A's frozen core aggregation, executed unchanged from `tests/analysis/m9a/audit.py` |
| `evidence` | 90 s | extracts the three illustrative failures from saved outputs only |

**Fixtures:**
- **J1:** two block marginals of 0.9 with joint agreement 0.8. 0.8 is the Fréchet lower bound for these marginals; the
  product, 0.81, is never used.
- **J2:** the same marginals with joint agreement 0.9 (perfectly dependent errors).
- **J3:** C covers only the reportable blocks, and a block's D is unaffected by the other blocks.

**Evidence:**
- ACOT HG01975 LOO sv:5 (M8-D, r0.01 and r0.1):
  - block posterior and support;
  - the archived evaluator row;
  - the best evaluated hypothesis carrying the truth-closest pair (oracle tier, evaluation only);
  - label source walks.
- C4 NA19650 LOO sv:7 at 10× and 30×, and CYP2D6 NA19650 LOO sv:9 at 5× and 30× (M1-B fresh-caller runs):
  - block posterior, support and filters;
  - the M1-B structure verdict;
  - the dosage audit row;
  - label source walks.

Each case is bound to its own graph, reference, curation, profile/preset, `min_variant_bp` and model fingerprint, taken
from its recorded command line and settings.

**Also cited without new execution:** the CYP2D6 block-representation degeneracy inventory
(`cyp2d6_shortread_pilot_v1/degeneracy/`; sealed, model-only).

## Budget

- **Total: 300 s.**
- **Pre-freeze: 2 s,** disclosed. It covers small read-only extractions while gathering evidence; the figure is an
  estimate, rounded up.
- **Admission rule:** pre-freeze + spent + timeout ≤ 300.

## Disclosed defect (cosmetic, not corrected in sealed M9-A)

The committed M9-A runner's docstring still describes the M8-D runner: the substitution I made when deriving it did not
match. Its code (manifest, budget, output directory) is M9-A's, and the M9-A ledger confirms that behaviour. The M9-B
runner has its own docstring.

## Deliverables (RESULT.md)

1. **An M9-A qualification,** appended to its RESULT. It replaces "understates SV certainty 4–5×" with the conditional
   masses, defines C as joint agreement and D as agreement at one block, and states that neither is a calibrated
   probability of correctness.
2. **The output mapping from code:** prediction units, the block-to-record relation, the label-to-sequence and
   label-to-structure relations, and projection-ambiguity handling.
3. **Three targets, described:**
   - exact block-allele agreement;
   - emitted SV-record agreement, if a sound mapping exists;
   - future structural-representative agreement, with the information that must be preserved.
4. **The failure check:** which target each failure violates, and which targets remain unassessed.
5. **One proposed user-facing contract,** plus the validation still needed and the smallest next experiment, with
   development samples and held-out samples kept apart.

## Amendment 1 (declared before any further execution)

**What failed.** `fixtures` exited 1 (J1 and J2 FAIL), and the runner's gate then refused `evidence`.
- The failing gate raised after the unit directory was created, so `runs/evidence/` exists empty, with no ledger row.
  It is kept, not deleted.

**Cause: a fixture definition defect.** Label pairs are (hap1 over blocks, hap2 over blocks).
- **J1:** the competitor meant to disagree at block 0 only, `((1, 9), (3, 4))`, has block-0 genotype (1, 3) and so
  agrees there.
- **J2:** the "both blocks" competitor, `((1, 9), (3, 5))`, likewise agrees at block 0.

**Correction.** J1's competitor becomes `((9, 2), (3, 4))` (disagrees at block 0 only), and J2's becomes
`((9, 9), (3, 5))` (disagrees at both). The expected values are unchanged: marginals 0.9 and 0.9, with joint 0.8 for J1
and 0.9 for J2. The aggregation code (M9-A core) is not touched.

**Rerun.** As `fixtures.a1` and `evidence.a1` (manifest A1, same timeouts); `evidence.a1` is gated on `fixtures.a1`.
The failed `fixtures` unit stays recorded.
