# M1-C: isolate C4 allocation from sequence-content differences

Branch `m1-c-allocation` from `5eeab46`. Diagnostic evaluator only. There are no production model, search, prior,
confidence or caller changes, and no truth input to any caller. Only `tests/genotype_probe.cpp` changes, and the
`panvar` binary is unchanged. Results: `/Users/davide.bolognini/panvar/tests/results/m1c_allocation_20261019`.

## 1. M1-B scope

Qualifications are appended to the M1-B results README:
- pooled-content superiority holds among the declared alternatives only;
- many evaluated 2+2 candidates is not exhaustive search;
- allocation was not isolated from content;
- the rebuilt-mass mismatch compared two different reductions.

## 2. Mass reconciliation (before any placement-level reading)

The scorer keys accepted origins by `StateKey` (insert, orientation, summed latent mass), counts places per class, and
reduces in key order with one `log_add_exp` per class: `log(count) + log(1/2) + log P(insert) + latent`. The M1-B check
used `log_placement_mass`, which sums origins one at a time in recorded order. That is a different reduction, so bit
identity was never expected.

`reconcile_origins` (probe) rebuilds each cell both ways and reports:
- bit agreement for each reduction;
- the maximum absolute log-mass difference;
- non-finite disagreements;
- origins, state classes and maximum multiplicity;
- duplicate coordinates, i.e. the same (start, insert, orientation) recorded twice;
- the fragment terms rebuilt from each reconstruction, minus the scorer's.

**Gates:**
- **Synthetic fixtures** (`genotype_probe mass-reconcile`, two seeds): a three-copy tandem haplotype (one copy with
  substitutions, one with a 2-bp deletion) beside a one-copy haplotype, with reads from both and both library
  orientations. Canonical reconstruction must be bit-equal in every cell, the rebuilt fragment terms bit-equal, and
  multiplicity > 1 must occur.
- **One existing C4 decomposition:** the M1-B C4 10x score-set (`reported` vs `saved_2+2`), rerun with the
  reconciliation. Its saved expectations must reproduce.

If canonical reconstruction disagrees, `placement_evidence_interpretable` is false and no placement-class split is
interpreted. Scorer-derived score comparisons are kept separately in that case.

## 3. Explicit-sequence scoring (`score-set --explicit`)

A candidate may be two complete sequences, given as FASTA records with the target-block interval of each. It is scored
with the unchanged fragment model, library, band and exposure. The prior is the panel's sequence prior. Bytes that no
allele path spells get prior −inf, reported as `prior_support none` with no total: never a neutral prior, and never a
total compared with production rankings. Each run also writes exact bits for comparison across runs:
- `fragment_term_hex` and `exposure_term_hex`;
- per unique haplotype: sequence FNV-64, exposure as hexfloat, FNV-64 of the mass column, and fragments with mass.

**Gate.** All 8 M1-B C4 10x path candidates are re-declared as explicit sequences spelled from their own paths
(`genotype_probe spell`) in a separate process, plus four controls: homologue exchange, duplicate label, a duplicate
sequence (homozygote) and the off-panel truth pair. All of the following must hold; otherwise no explicit alternative is
scored:
- candidate terms bit-equal to the path-based run, with prior and total equal;
- every on-panel haplotype's exposure and mass column bit-equal;
- the swap and duplicate-label controls bit-equal to `reported`;
- the homozygote has a defined prior;
- truth has `prior_support none` and no total.

## 4. The frozen allocation experiment (`construct.py`; C4 NA19650 at 10x and 30x only)

At most 12 explicit alternatives per coverage, fixed by rule before scoring:

| label | what |
|---|---|
| `reported` | the reported 3+1 pair: the rank-1 representative of the M1-B fresh run, full-precision saved total as expectation |
| `saved_2+2` | the M1-B selected 2+2 comparator, saved total as expectation |
| `realloc_A_p1`, `realloc_A_p2`, `realloc_B_p1`, `realloc_B_p2` | 2+2 reallocations of the reported pair's own units (below) |
| `truth` | the whole-locus truth pair (diagnostic reference; off-panel) |
| `truthblock_in_reported_a/b` | the reported flanks with the truth sv:7 bytes, in both assignments |
| `reportedblock_in_truth_a/b` | the truth flanks with the reported sv:7 bytes, in both assignments |

**Unit boundaries.** An anchor 31-mer must satisfy all of:
- 3 copies in the 3-unit block;
- 1 copy in the 1-unit block;
- 2 copies in each truth block;
- no other copy on either reported haplotype;
- consecutive spacings on the 3-unit side within 5 % of each other.

Phase A is the first anchor in the 1-unit block. Phase B is the first anchor at least half a unit later, so the moved
period is framed differently.

**Reallocation.** A full anchor-to-anchor period P_k of the 3-unit haplotype H is moved to the 1-unit haplotype G at its
anchor: H' = H[:a_k] + H[a_(k+1):] and G' = G[:q] + P_k + G[q:].
- Reconstruction of the original pair is verified before any change.
- Afterwards, base composition, total length and anchor counts of 2 + 2 are asserted.
- Cuts, moved bytes and hashes, lengths and the 62-bp context of every new junction are recorded.
- Nothing is duplicated, discarded or reversed.
- If no anchor qualifies, the limitation is recorded and nothing is forced.

**Content measured, not assumed** (`content.py`, whole haplotypes, canonical k = 31 and 151):
- pooled discrepancy from `reported`, which counts what a reallocation's new junctions changed;
- pooled and per-haplotype discrepancy from truth.

## 5. Questions and reporting

For each alternative the report gives the fragment term, exposure term, their sum (the likelihood), and the prior where
it is defined. Decompositions (scorer masses, compact per-fragment tables, top-200 origin detail):
- `reported` against each reallocation, `truth`, `saved_2+2` and `truthblock_in_reported_a`;
- `truth` against `reportedblock_in_truth_a`.

1. **Allocation evidence:** the likelihood differences between `reported` and the reallocations (same unit inventory),
   read beside their content differences from `reported`.
2. **Truth versus reported:** the likelihood of `truth` against `reported`. The block swaps separate the target-block
   content from flanking differences.

A flat result concerns these alternatives and these reads. It is not a general proof of non-identifiability. The next
experiment is recommended from the result: representation limitations, allocation evidence, or likelihood adequacy.
There is no clustering and no switch-probability tuning.

## 6. Budget

1,200 s cumulative (development runs included), one process at a time, at most 4 threads, at most 300 s per unit, and a
6 GiB envelope. The manifest is frozen and committed before execution. Amendments are declared before further execution,
and failed units stay recorded.
