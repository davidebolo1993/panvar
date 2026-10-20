# M9-E result: the three improving neighbours were never generated, because sv:5 carries no informative proposal marker, so all 6,786 allele pairs tie and the 32-pair repair shortlist is filled purely by allele-index order

**Run.** Branch `m9-e-reachability`, freeze `fafab52`, with no amendments. Results are in
`tests/results/m9e_reachability_20261020/`.
- **Execution:** 40.0 of 300 s, including 2 s of disclosed pre-freeze inspection. No failures.
- **Inputs:** saved traces and emissions-only probe runs. There was no caller rerun or new search.

## 1. Terminal contexts (saved traces)

| case | reported winner equals | final neighbourhood scan at that state |
|---|---|---|
| E1 (d1) | the trunk terminal (state after round 5 of 6) | round 6: 973 generated, 956 scored, 16 self, 1 visited, 0 duplicate, 0 prior-rejected; no move |
| E3 (d3) | the trunk terminal (after round 6 of 7) | round 7: 984 generated, 968 scored, 15 self, 1 visited, 0 duplicate, 0 prior-rejected; no move |
| E6 (d6) | **both converged tie-branch terminals** (branch 0 after 9 rounds; branch 1 after 8; both `no_improvement`). The trunk terminal is not the winner: the trunk stopped at round 5 on a tie | branch 0's terminal round: 888 candidates are recorded with `scored_in branch0_round9` |

**Branch 1's terminal scan is not resolved by the records.** `mosaic_candidates` keeps one `scored_in` per label pair, and
no candidate carries `branch1_round8`. Because branch 1 converged to the same state as branch 0, its scan may have met
only already-pooled pairs, but the records do not show it. Sol's identities are confirmed.

## 2. Shortlist reproduction (`score-set --rank-blocks … --emissions-only 1`)

**Settings match.** The probe calls the caller's own marker functions with k = 31 and s = 15, which equal the runs'
`kmer_size` and `syncmer_s`. No marker background is declared in either path; the background is estimated. The model,
curation, reads and library are the run's, and `repair_local_pairs` = 32 comes from the run's status.

**Gate.** All **15/15** recorded trunk moves (block, pair, `local_pair_rank`) are reproduced exactly by the probe's
ranking, across d1, d3 and d6 and nine different blocks.

| case | frozen improving alternative (Δtotal vs the winner) | sv:5 pair | rank of 6,786 | emission | cutoff emission (rank 32) | gap | in the 32 list | pairs tied at its emission / at the cutoff |
|---|---|---|---|---|---|---|---|---|
| E1 | A_truth_ba (+52.05) | 6/77 | **753** | 0.0 | 0.0 | 0 | **no** | 6,786 / 6,786 |
| E3 | A_truth_ab, A_truth_ba (+8.46) | 34/94 | **3,444** | 0.0 | 0.0 | 0 | **no** | 6,786 / 6,786 |
| E6 | A_truth_ab (+16.45), A_truth_ba (+17.55) | 13/21 | **1,439** | 0.0 | 0.0 | 0 | **no** | 6,786 / 6,786 |

**The deterministic tie rule:** emission descending, then `(allele_x, allele_y)` ascending (`CoordinateRepairSearch`).
- With every emission equal, the list is **(0,0), (0,1), …, (0,31)**: allele 0 paired with alleles 0–31.
- The list is built once and does not depend on the refinement context. Only the full score does.

**Other blocks are informative.** Pair emissions are almost all distinct: ctx:4 has 38,793 distinct values over 57,291
pairs, sv:13 has 57,055 over 59,685, ctx:6 has 10,781 over 13,530. **sv:5 is the only audited block with a single
emission value.**

## 3. Short alleles and markers

**How proposal markers work** (`build_marker_panel`): they are k = 31 syncmers **inside each block allele's own sequence**,
and markers occurring in several blocks are dropped. **Junction-spanning k-mers are therefore never proposal markers.**

| case | called alleles (bp; internal 31-mers) | improving alleles (bp; internal 31-mers) |
|---|---|---|
| E1 | 5 (27; 0), 71 (30; 0) | 6 (29; 0), 77 (28; 0) |
| E3 | 45 (33; 3), 107 (18; 0) | 34 (37; 7), 94 (23; 0) |
| E6 | 0 (7; 0), 21 (44; 14) | 13 (50; 20), 21 (44; 14) |

**Across the block:** sv:5 alleles span 7–52 bp, and **64 of 116 are shorter than 31 bp**, so no marker is possible for
them. For the 52 longer alleles, the internal 31-mers are few (at most 22), and **all 6,786 pair emissions are exactly
0.0.**
- **What that means:** no retained, informative proposal marker distinguishes any pair at sv:5 in this model.
- **What is not enumerated:** whether the longer alleles' syncmers were absent, fixed across alleles, or dropped as
  multi-block. There is no existing per-marker dump.

**Not length alone.** E3 and E6 involve alleles of 33–50 bp, which could carry internal k-mers. **The exclusion is caused
by the tie,** which allele-index order then resolves. It does not reflect informative marker evidence against the
improving pairs.

**This contrasts with the read evidence.** M9-D found that the full score prefers the improving pairs by +16 to +53 nats
of likelihood. That evidence sits in a few fragments **spanning the block junctions**, which is exactly what the
allele-internal marker model cannot see.

## 4. Generation and scoring

**Every frozen improving alternative, in both homologue assignments, was excluded by the shortlist.** It is never
generated at any trunk round or tie branch, because a neighbour sets sv:5 only to a listed pair.
- **E6:** this holds at the trunk and at both branches, since the list is context-independent.

**Not in the scored pool, by label or by sequence.** M9-D's inventory has no within-block duplicates and cross-block
equivalences only within sv:5..ctx:6, so only pool candidates equal outside that interval can spell the same pair:

| case | pool candidates checked (equal outside sv:5..ctx:6) | sequence-identical |
|---|---|---|
| E1 | 128 | 0 |
| E3 | 127 | 0 |
| E6 | 189 | 0 |

**Not evaluated:** the visited, duplicate and prior-rejection stages. They apply only to generated neighbours, and none
of these was generated.

**How the called sv:5 pairs entered the search:**
- **E1 (5/71) and E3 (45/107):** these pairs are outside the 32 list, so they cannot have come from a repair move at sv:5.
  They come from the union or panel start.
- **E6 (0/21):** this pair is inside the list (rank 22), and the tie that started branching was a swap at sv:5. The saved
  records do not show its origin.
- **In every case,** repair could change sv:5 only to (0, y ≤ 31).

## 5. M9-D semantics corrected (context-aware classification; fixture 6/6 PASS)

The `fixture` unit tests synthetic rows:
- context A against the reported winner: in pool and above → INCONSISTENCY; absent and above → SEARCH_GAP; below →
  no improvement;
- context B carrier versus its called insertion: not an inconsistency unless a pool member scores above the reported
  winner.

**All 52 M9-D contrasts, reclassified** (`m9d_reclassified.tsv`):

| reading | count |
|---|---|
| **INCONSISTENCY** | **0** |
| SEARCH_GAP | 5 (the five frozen improving alternatives) |
| no improvement over the winner | 25 (context A) |
| consistent | 22 (context B) |

The correction is appended to M9-D's `RESULT.md`. Its tables are unchanged.

## 6. Cause table

| case | terminal context | improving pair: rank / cutoff | proposal emission at sv:5 | stage reached | cause |
|---|---|---|---|---|---|
| **E1** | trunk terminal | 6/77: 753 / 32 | all 6,786 pairs = 0.0 | excluded by shortlist; not in pool (label or sequence) | uninformative proposal emissions at a short-allele block; the tie is resolved by allele index |
| **E3** | trunk terminal | 34/94: 3,444 / 32 | all = 0.0 | excluded; not in pool | same |
| **E6** | converged tie-branch terminal (both branches) | 13/21: 1,439 / 32 | all = 0.0 | excluded at trunk and branches; not in pool | same |

## 7. Recommended truth-free proposal improvement (one)

**For a block whose proposal emissions are tied across the local-list cutoff, derive its local-pair ranking from read
evidence the allele-internal markers miss.** Score candidate pairs with markers taken over the block **together with its
flanking junctions**: k-mers spanning block plus context, computed for the incumbent's flanks. Then the tie is broken by
reads rather than by allele index.

**Supported here by:**
1. sv:5's emissions are identically 0.0, and the list is index-ordered;
2. the missed improvements are favoured by the full likelihood through junction-spanning fragments (M9-D);
3. markers are allele-internal by construction, and 64 of 116 sv:5 alleles are shorter than k.

**Truth-free:** it uses only the reads, the model and the incumbent.

**Not proposed:**
- widening the list until the truth enters (ranks 753–3,444 of 6,786 would need a 100-fold wider list);
- reviving joint search.

**Not implemented:** no caller change.
