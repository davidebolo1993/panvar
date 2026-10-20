# M9-D result: three of the six ACOT high-confidence errors are search gaps (a one-block replacement of the reported winner scores higher and was never evaluated), two are copying-prior effects in the LOO panel, and one is a weak likelihood preference

**Run.** Branch `m9-d-attribution`.
- **Commits:** freeze `eb9a28f`, which includes the M9-C qualification; Amendment 1 `fe4cfe4`, the interval check over the
  degenerate sv:5..ctx:6; Amendment 2 `d9d6034`, a record-name fix.
- **Results:** `tests/results/m9d_attribution_20261020/`.

**Execution: 185.4 of 900 s.**

**Recorded failure:** `intervals` hit a record-name defect and was rerun as `intervals.a2`.

**Disclosed labelling limitation.** The automatic `reading` column in `contrasts.tsv` applies its rules correctly only in
context A.
- In context B, the "replacement" column holds the context itself (the truth carrier, already in the pool). The label
  "already scored higher … INCONSISTENCY" there is an artefact; the numbers are unaffected.
- The B readings are restated below.
- No constructed pair from the scored pool scores above the reported winner, so there is no real inconsistency.

**Scope:** no recalibration, tuning, caller search, residual editing, clustering, filter change or promotion. Truth-derived
alternatives are diagnostic only and appear in no caller output.

## Gates

| gate | result |
|---|---|
| model fingerprints (dump vs runs) | LOO `12cbe6be…` = all five LOO runs; LZO d1 `4606e1bd…`; LZO d2 `545ff712…` |
| spelling | 116/116: replaced block = replacement allele, flanks byte-identical, offsets recorded |
| LZO control allele mapping | by exact sequence, one hit each (LOO 5/71 → LZO 5/71; LOO 76/89 → LZO 76/89) |
| context reproduction | 14/14 contexts: fragment, exposure, prior and total within 1e-6 of the saved scores; score-set's expected totals reproduced in every case |
| decomposition integrity | `placement_evidence_interpretable` true everywhere; decomposed fragment term = scorer exactly; reconstructed bounds ≤ 2.7e-7, so every contrast below is certified |
| degeneracy inventory (model-only) | all three models: 0 within-block duplicates, **2,045 cross-block equivalences, all in sv:5..ctx:6** |

**The degeneracy inventory changes one earlier statement.** M9-C's "no representation ambiguity" held only among the
explored labellings. Globally, sv:5..ctx:6 is a degenerate representation region.

**Interval check** (Amendment 1): reported sv:5 + ctx:6 bytes (sv:9 for E4) against the truth projection, per haplotype.

| case | reported hap1 vs truth | reported hap2 vs truth |
|---|---|---|
| E1 | differs | differs |
| E2 | identical | differs |
| E3 | differs | differs |
| E4 | differs | identical |
| E5 | differs | differs |
| E6 | identical | differs |

So **none of the six errors is a label-only artefact** of the degeneracy: each involves a real sequence difference on at
least one haplotype.
- Residual differing cores are 0–4 bp, except E5 hap1 and E6 hap2. There the differences lie at two separated positions,
  because ctx:6 is also wrong and unavailable (floor 2).

## Matched contrasts (truth pair − called pair at the affected block, all other blocks fixed)

**Each cell** reports the best homologue assignment: Δlikelihood (fragment + exposure) / Δprior / **Δtotal** in nats.
- **Context A** is the reported winner: the best truth replacement minus the winner. **No truth replacement in A was in
  the scored pool.**
- **Context B** is the best evaluated truth carrier: the carrier minus the best called-pair insertion.

| case | block, cell | context A (reported winner) | context B (best truth carrier) | informative molecules | reading |
|---|---|---|---|---|---|
| **E1** | sv:5, d1 | +53.04 / −0.98 / **+52.05** | +40.00 / +37.44 / +77.44 | 9 fragments (7 > 1 nat), all spanning the block junction | **search/proposal gap**: a one-block replacement beats the reported winner by 52 nats and was never evaluated |
| **E2** | sv:5, d2 | +7.23 / −9.91 / **−2.68** | +7.23 / +18.37 / +25.60 | 3 junction fragments (1 > 1 nat) | **prior effect in context A**: likelihood favours truth by 7.2, the LOO-context copying prior favours the call by 9.9; the sign depends on context |
| **E3** | sv:5, d3 | +29.44 / −20.98 / **+8.46** | +311.31 / +39.93 / +351.24 | A: 11 junction fragments (3 > 1 nat). B: dominated by one molecule placed inside the block only on the truth haplotype (+249) | **search/proposal gap** (+8.5, not evaluated) |
| **E4** | sv:9, d3 | +3.88 / −9.20 / **−5.31** | +3.88 / +16.94 / +20.83 | 12 junction fragments (1 > 1 nat) | **prior effect in context A** (likelihood +3.9 for truth; prior −9.2) |
| **E5** | sv:5, d4 | −1.28 / 0.00 / **−1.28** | −1.02 / +21.21 / +20.19 | 13–14 junction fragments, none > 1 nat (a broad, small shift) | **weak likelihood preference for the call** (about 1 nat); context A also carries a wrong ctx:6 |
| **E6** | sv:5, d6 | +16.45 / +1.10 / **+17.55** | +16.45 / +21.21 / +37.66 | 11 junction fragments (3 > 1 nat) | **search/proposal gap** (+17.5, not evaluated) |

**LZO controls** (correct call minus LOO-called insertion, same donor and reads):

| control | Δlikelihood | Δprior | Δtotal |
|---|---|---|---|
| C1 (d1) | +53.04 | +37.87 | +90.91 |
| C2 (d2) | +7.23 | +18.38 | +25.61 |

- **The likelihood differences equal E1's and E2's to the digit**, so the block-local read preference for the truth pair
  is reproduced in an independent context.
- **The prior's sign is context-dependent.** It favours the truth pair when the donor's own haplotypes are in the panel
  (LZO) or when the context is a truth carrier (B). In the LOO reported context it favours the call (E2, E4).

## What this shows (ranked, evidence-backed)

1. **Search/proposal gaps** in three of six cases (E1, E3, E6; all sv:5).
   - Changing only sv:5 in the reported winner yields a hypothesis that scores **+8.5 to +52 nats higher** (certified).
   - It is absent from the scored pool, and the likelihood also favours it (+16 to +53).
   - The M9-C error masses (1e-38, 1.5e-24, 9e-7) were computed without it, which is why they were confidently wrong.
     They are conditional on a pool that missed a better neighbour.
   - **Why the bounded repair did not evaluate this move is unresolved.**
2. **Copying-prior effects** in two cases (E2 sv:5, E4 sv:9).
   - The reads favour the truth pair modestly (+7.2, +3.9).
   - In the LOO reported context the sequence-level copying prior favours the call (−9.9, −9.2).
   - The same prior favours truth in other contexts.
   - This is quantified only: no defect is claimed and nothing is retuned.
3. **A weak likelihood preference for the call** in one case (E5).
   - About 1 nat, spread over 13–14 junction fragments with no informative single molecule.
   - The context also contains a wrong, unavailable ctx:6 allele, so the contrast compares two imperfect contexts.
   - Its cause is unresolved, and this alone does not establish a bug.

**Common to all.** The evidence comes from **a few molecules spanning the affected block's junctions**, not from broad
read-level shifts, except E5. No contributor is attributed to sequencing error or placement: placement evidence is
interpretable, but no error model was examined.

**What remains unresolved, in order:**
1. why the repair neighbourhood missed the one-block sv:5 improvements (E1, E3, E6);
2. whether the LOO-context prior preference is expected panel-composition behaviour or a modelling limitation (E2, E4);
3. E5's weak, diffuse likelihood preference;
4. what the global sv:5..ctx:6 degeneracy implies for a block-level confidence unit.

## One next experiment (recommended, not run)

**A proposal-reachability audit at sv:5 for E1, E3 and E6,** using the existing `score-set --rank-blocks sv:5` local-rank
table (`in_local_list`) in each reported winner's context. Questions:
- Was the truth pair (both assignments) inside the caller's `repair_local_pairs` = 32 shortlist?
- If it was, why was the move not generated or kept?
- If it was not, where did it rank?

This would show whether the gap is in local pair ranking or in move generation, with no new search, tuning or caller
change.

## Correction (appended by M9-E, 2026-10-20; the text and tables above are unchanged)

**What was wrong.** The automatic `reading` column of `analyse/contrasts.tsv` was not context-aware.

**What replaces it.** M9-E reclassifies all 52 contrasts (`m9e_reachability_20261020/runs/analyse/m9d_reclassified.tsv`;
fixture 6/6). An **inconsistency is claimed only for a scored-pool member that scores above the reported winner.** A
context-B truth carrier beating its called-pair insertion is not, by itself, an inconsistency.

| reading | count |
|---|---|
| INCONSISTENCY | **0** |
| SEARCH_GAP | 5 (E1 A_truth_ba; E3 A_truth_ab and A_truth_ba; E6 A_truth_ab and A_truth_ba) |
| no improvement over the winner | 25 |
| consistent | 22 |

**What M9-E found about the search gaps.** They are shortlist exclusions. sv:5's proposal emissions are identically 0.0
across all 6,786 pairs, so the 32-pair list is allele-index order, and the improving pairs rank 753, 3,444 and 1,439.
