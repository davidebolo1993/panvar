# Pre-registration: does fragment-level evidence fix the ranking, or only the confidence?

Written before any arm below was run. Follows `CYP2D6_ORACLE_RESULT.md`, which established that
whole-locus multiplicity restores identifiability under ideal counts and that the residual failure is
consistent with about four independent fragment observations behind 23 correlated markers.

## The distinction being tested

An effective-sample-size correction scales every candidate's score by the same factor. In this
emission `emis = baseline + rho * ll + ...` with `baseline` common to all pairs, so any pairwise
difference is `rho * (ll_A - ll_B)`: **ordering is invariant by construction.** The negative-results
ledger already records this -- "the ESS discount as the cause -- refuted, it cannot reorder within a
block".

So:

| | fixes overconfidence | can fix ranking |
|---|---|---|
| scalar ESS discount | yes | **no, provably** |
| per-clump aggregation | yes | possibly |
| fragment-level observation | yes | possibly |

Only changing **what is observed** can change which candidate wins.

## Four arms. Frozen.

| arm | construction |
|---|---|
| **F1** | whole-locus shared-marker factor, 23 marker terms, no discount |
| **F2** | F1 multiplied by ESS/markers -- the confidence-only control |
| **F3** | per-clump mean: each independent window contributes one aggregate term |
| **F4** | read-level fragment signatures: each physical fragment contributes one compatibility observation |

**F2's local ranking must be unchanged.** That is a prediction from the algebra above, not a
hypothesis; it is run as a correctness check that the implementation matches the arithmetic. A
difference there is a bug, not a finding.

## Prerequisite for F3, recorded because it currently blocks it

F3 needs clumps defined from **panel or graph coordinates and library fragment length, never truth
coordinates**. Production's existing definition cannot do this at cyp2d6 block 5:

| | |
|---|---|
| windows visible from truth coordinates (diagnostic only) | 4 |
| **clumps from `node_first_pos / fragment_len`** | **1** |
| **markers with no recorded position** | **33 of 36** |

`node_first_pos` is populated for 3 of the 36 restored markers, and stores only a marker's FIRST
occurrence, so a marker appearing four times is assigned one window. **F3 is not implementable until
that is fixed**: it would aggregate into a single window rather than four.

**F4 needs no marker positions at all.** It observes fragments directly, so it sidesteps the
prerequisite entirely. That makes it the more buildable arm, not the more speculative one, and it is
where this starts.

## F4 construction

Each read pair is one observation. For a fragment, take the ordered signature of panel markers it
carries, and score it **once** against each candidate's compatibility with that signature. Do not
also count its individual markers at full node weight -- that is the double counting `--edge-weight`
was measured to suffer from.

Design note specific to cyp2d6: a fragment carrying a shared marker frequently also carries a
block-5-only marker. That co-occurrence is what localises a shared count to a block, without any
truth input. The same mechanism is the primary route for gstm1, where a 150 bp fragment spans a
5-42 bp allele and both its boundaries; here it serves a different purpose -- preventing correlated
shared-marker votes from being counted repeatedly.

## Reported for every arm

Certified class in the top set; certified excess; top-class size; wins and losses per donor; GQ
calibration; leave-zero-out regressions. Over **multiple read seeds**, since the residual under test
is a sampling phenomenon and a single seed cannot separate it from a systematic effect.

## Gate

The improvement must hold **across seeds**, move the remaining failures toward their certified class,
and not damage the donors already exact. Specifically:

1. F2's ranking is identical to F1's -- correctness check;
2. certified excess falls against F1 by at least 10%, on the median across seeds;
3. no donor exact under F1 becomes non-exact;
4. the near-representable donors (E\* < 1% of aligned sequence) improve, since they are where the
   caller currently misses by 58-176 edits against a floor of 1-2;
5. leave-zero-out does not regress;
6. it holds on donors not used for design.

## Scope

All 25 scored donors plus the three near-representable ones misclassified earlier as
non-representable -- 28 in total. Representability is defined from certified sequence error, never
from whether a truth label exists.

## Declared outcome

- F4 improves ranking -> the physical fragment is the right observation unit; implement it.
- F4 fixes only confidence, like F2 -> correlated sampling is real but not what decides the call, and
  the remaining CYP2D6 failures need a different explanation.
- F3 becomes implementable and matches F4 -> prefer F3, it is cheaper.

No production change before F4 passes the gate.
