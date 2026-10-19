# M1-D: can active-panel units recover better structural approximations?

Branch `m1-d-unit-reachability` from `a70f71d`. This is a diagnostic experiment only:
- no production, search, prior, confidence or caller change, and no truth input to a caller;
- only `tests/genotype_probe.cpp` changes, adding the diagnostic `edit-distances` subcommand;
- the `panvar` binary is untouched.

Results: `/Users/davide.bolognini/panvar/tests/results/m1d_unit_reachability_20261019`.

## 1. M1-C scope

Four qualifications are appended to the M1-C README:
- examined alternatives only;
- 2–5.5 nats is not zero and not a calibrated equivalence criterion;
- adequacy holds for this comparison only;
- truth stays excluded by representation and prior support.

## 2. Panel-only unit inventory (`inventory.py`), frozen before any truth comparison

**Source.**
- `genotype_probe block-sequences` without `--truth-gfa` on the C4 LOO model, which contains no NA19650 path. The
  inventory asserts that no NA19650 template is among the carriers.
- Only sv:7 alleles carried by at least one active-panel template are used. No truth, no truth-derived anchor, no other
  graph sequence.

**Boundaries.**
- S is the shortest carried allele of at least 10 kb.
- An anchor is a 31-mer that:
  - occurs exactly once in S;
  - occurs in every carried allele;
  - contains no N;
  - has consecutive occurrences at least 20 kb apart.
- Phase A is the first such anchor in S. Phase B is the first one at least half of S later.
- Boundaries are reliable only if both phases give every allele the same unit count and every phase-A spacing lies in
  [20, 40] kb. Otherwise the limitation is reported and nothing is forced.

**Components** of each allele, from its phase-A anchor starts:
- prefix = `X[0:p1]`;
- period_i = `X[p_i:p_(i+1)]`;
- tail = `X[p_n:]`.

Lossless reconstruction (prefix + periods + tail = allele) is asserted for every allele.

Each component records its source allele, coordinates, orientation (the model frame), length, SHA-256 and carriers.
Identical bytes are listed once per type.

**Freeze.** The inventory files' SHA-256 are committed (`INVENTORY_FROZEN.tsv`) before any truth unit runs. The
reachability unit refuses an inventory that differs from that record.

## 3. Reachability (`reach.py`; isolated evaluator, truth allowed)

Truth's two sv:7 blocks (the M1-B projection) are partitioned with the inventory's own phase-A anchor, which was chosen
from panel bytes. Then:
- **exhaustive distances:** every truth component against every unique inventory component of the same type, and every
  truth block against every carried whole allele (`genotype_probe edit-distances`, unbanded global edlib);
- **exact unit availability:** distance 0;
- **nearest existing whole-block allele:** W(truthH);
- **compositions from existing components only** (no new bases, no consensus, no within-unit patching):
  - `compose`: nearest prefix + nearest period(s) + nearest tail;
  - `periodonly`: W with only its period(s) replaced by the nearest panel period(s).

Nearest-unit distances do not establish the best whole-diplotype likelihood.

## 4. Candidate set (`construct.py`; at most 12 per coverage; declared before scoring)

All compositions are **ORACLE** representation-feasibility tests. They were selected as nearest to truth. They are not a
truth-free generator and not evidence that a caller would find them.

Reported flanks are fixed for every candidate except `saved_2+2`, which is a whole panel diplotype. Side 1 is the
reported 3-unit homologue; side 2 is the 1-unit homologue.

| label | sv:7 side 1 / side 2 |
|---|---|
| reported | reported 3+1 (saved total as expectation) |
| saved_2+2 | M1-B panel 2+2 comparator (saved total as expectation) |
| truthblock_a | truth1 / truth2 sv:7 bytes: the truth-block reference |
| whole_nearest_a | W(truth1) / W(truth2) |
| periodonly_a | W(truthH) with the nearest period(s) |
| oracle_compose_a | compose(truth1) / compose(truth2), 2+2 |
| oracle_compose_b | the same blocks exchanged between homologues, 2+2 |
| oracle_compose_c | the composed periods exchanged, 2+2 |
| oracle_compose_d31 | P1+U1+U2+T1 / P2+T2, 3+1 |
| oracle_compose_e13 | P1+T1 / P2+U1+U2+T2, 1+3 |

The c, d31 and e13 variants are asserted to use exactly the composed a/b unit inventory. Every side is reported left
flank + block + reported right flank (asserted). Every component equals its frozen inventory record (SHA-256 asserted),
and every join's 62-bp context is recorded.

## 5. Scoring and reporting

- `score-set --explicit`, validated in M1-C, with the unchanged library, fragment model, band, exposure and options.
  There is no origin recording (likelihoods only) and one batch per coverage.
- Saved expectations for `reported` and `saved_2+2` must reproduce within |Δ| ≤ 1e-6.
- Reported per candidate (`summary.py`):
  - fragment term, exposure term and likelihood;
  - prior support, with no neutral prior and totals only where defined;
  - unordered sv:7 edit distance to truth;
  - pooled block k-mer discrepancy against truth and against the reported blocks (k31 and k151);
  - the likelihood change against `reported`, and the share of `truthblock_a`'s gain recovered;
  - the likelihood differences among the matched composed inventory (a against b, c, d31 and e13).

## 6. Decision

The report answers three questions separately:
1. Do recombined panel units give a materially better approximation?
2. Does better content resolve allocation, or do allocations stay close?
3. What remains attributable to missing unit sequence rather than missing combinations?

A better construction is positive feasibility evidence. Not finding one in this bounded set does not prove
impossibility. One next experiment is recommended; nothing is implemented automatically.

## 7. Budget

1,200 s cumulative, one process at a time, at most 4 threads, at most 300 s per unit, 6 GiB envelope. Manifests are frozen
and committed before execution, amendments are declared first, and failed units stay recorded.
