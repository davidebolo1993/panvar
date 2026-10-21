# M12-C result: with the six C4 blocks completed, the dosage-first rule stays INCONCLUSIVE; K1 joins K5 as a catastrophic prediction made worse; the obstacle is choosing within the correct copy state

**Run.** Branch `m12-c-dosage-first-completion`; freeze `4a99f65`; result after.
- **Results:** `tests/results/m12c_dosage_first_completion_20261021/`.
- **Budget:** 289.2 of 900 s: 60 s of disclosed development plus 229.2 s of ledger wall.
- **Units:** 8 of 8 passed first time, with no amendment and no watchdog event; peak 0.82 GiB.

**Not done:**
- no CYP2D6 rerun and no ACOT scoring;
- no reads, simulation, caller rerun, validation use, feature change, new sequence, fit, prior change, merge, push or cleanup.

**M12-B's rule, selections, results and decision are preserved.** Interpretation corrections are appended to its RESULT.

## 1. Completion of B10, B12, B13, B14, B16 and B17 (K1)

**`bind`: 20 of 20 checks pass.**
- **Executables:** `genotype_probe` `ff0932fd…` and `m10c_align_eval` `21aabc85…` equal M12-B's manifest.
- **Candidate lists:** each frozen list's count equals M12-B's selection. Hashes are in `MANIFEST.tsv`.
- **M12-B plan:** each block was admitted by the bound and NOT RUN on budget.
- **Reads:** SHA-256 equal to M12-B's.
- **Count-state selection:** not recomputed. M12-B's outputs are linked read-only.

**Scoring:** M12-B's exact score-set command, exhaustive within each selected state with both assignments. Each block took
19–43 s.

**Gates per new block** (M12-B's unchanged `evaluate`) — all pass:
- incumbent reproduced, and its distance equals the evaluator's;
- byte lengths and FNV identity;
- shared candidates (45/45 to 65/65); no alias inconsistency; coverage 100%.

**M12-B's 12 scored rows reproduce identically (12/12).**

| block | state: original → selected (truth per haplotype) | excess (bp): original → selected | outcome |
|---|---|---|---|
| B10 C4 HG02018 | (4,3,0) → (4,3,0) ((2,1,0)\|(2,2,0)) | 2 → 2 | unchanged |
| B12 C4 NA12329 | (4,3,0) → (4,3,0) | 6 → 6 | unchanged |
| B13 C4 HG00358 | (5,3,0) → (5,3,0) | 0 → 0 | unchanged |
| B14 C4 NA19331 | (5,2,0) → (5,2,0) | 2 → 2 | unchanged |
| B16 C4 HG03732 | (5,3,0) → (5,3,0) | 0 → 0 | unchanged |
| **B17 K1 C4 NA24385** | (5,4,0) → **(4,4,0) = the truth's total**, but allocation (1,1,0)\|(3,3,0) against the truth's (2,2,0)\|(2,2,0) | **39,155 → 65,462 (worse)** | RETAINED_LARGE |

**Why K1 got worse** (existing records only):
- **Within the correct total state, the objective prefers the 1+3 pair 115/48 over the evaluator-best 2+2 pair 44/34 by
  68.4 nats:** fragment +58.5 and prior +9.9 in 115/48's favour; exposure is equal.
- **The count loss also ranks a 1+3 pair (47/48) best within the state,** ahead of 44/34 (299,423 against 299,875).
- **The count state is right and the allocation is wrong.** Neither the objective nor the count loss chooses the truth's
  2+2.

## 2. K5 from the existing M12-B records (no new scoring)

| candidate | per-haplotype pooled-feature states | total state (feature-0 total) | count loss | fragment / exposure / prior | total (Δ vs selected) | excess (bp) | block bp |
|---|---|---|---|---|---|---|---|
| **selected 104/97** | (2,1,1,1,1,0,0) \| (2,1,1,1,1,0,0) | (4,2,2,2,2,0,0) (4) | **249,164** | −69,705.6 / −5,337.2 / −44.0 | −75,086.76 (0) | **24,240** | 27,469 + 27,464 |
| evaluator-best 95/10 | (3,2,2,1,1,0,0) \| (1,0,0,1,1,0,0) | **(4,2,2,2,2,0,0) — the same state** (4) | 251,675 | −69,819.7 / −5,337.6 / −41.7 | −75,199.01 (**−112.26**) | **0** (distance = floor 1,625) | 39,639 + 15,327 |
| evaluator-best 10/95 | swapped | same | 251,675 | −69,819.7 / −5,337.6 / −60.5 | −75,217.84 (−131.09) | 0 | |
| incumbent 95/97 | (3,2,2,1,1,0,0) \| (2,1,1,1,1,0,0) | (5,3,3,2,2,0,0) (5) | 258,693 | −69,535.3 / −5,489.3 / −25.1 | −75,049.73 (+37.03) | 12,075 | |
| truth (diagnostic) | (3,2,2,1,1,0,0) \| (1,0,0,**0**,1,0,0) | (4,2,2,**1**,2,0,0) (4) | — | — | — | — | 39,635 + 13,769 |

**The four quantities, kept distinct:**
- **Biological total unit count.** It is **not established by this audit:** pooled features are groups of k-mers, and none
  is identified with a biological unit.
  - Feature 0 (3,056 k-mers; 1–4 copies per allele) totals **4** in the selected pair, in the evaluator-best pair and in
    the truth.
  - Its total in the incumbent is 5.
- **Complete pooled-feature vector.** No panel pair reaches the truth's (4,2,2,1,2,0,0).
  - Every one of the 112 panel alleles carries at least one copy of feature 3 (788 k-mers), and truth haplotype 2 carries
    none.
  - **The selected state differs from the truth only in feature 3, and so does the evaluator-best pair's state.**
  - "No available pair matches the full truth vector" is therefore **not** "the correct total copy number is unavailable".
- **Per-haplotype allocation.** The evaluator-best pair's hap1 equals the truth's hap1, (3,2,2,1,1,0,0). The selected pair
  splits the copies 2|2 against the truth's 3|1.
- **Sequence excess.** Selected 24,240; evaluator-best 0.

**K5's worsening was not caused by excluding the best pair.** It was in the selected state and was scored.
- **The unchanged objective prefers 104/97 by 112.3 nats,** almost all from the fragment term (+114.1).
- **The count loss also prefers 104/97 by 2,510.**
- The rule's within-state choice made the prediction worse.

## 3. Combined coverage and decision

**All 29 original blocks** (`runs/evaluate/coverage_ledger.tsv`):

| group | blocks |
|---|---|
| scored: M12-B | 12 |
| scored: M12-C | 6 |
| not scored: ACOT, LIMITATION_STATE_TOO_LARGE (unchanged by this batch) | 11, including K2 |

**Copy state** (truth-free selection against the truth's total state; unchanged from M12-B): 27 of 29 match.
- All 24 originally correct states are kept.
- 3 of 5 catastrophic blocks are fixed: K1, K2 and K4.

**Sequence outcome, 18 scored blocks:**

| result | blocks |
|---|---|
| fixed large error | 1 (K4: 12,159 → 17) |
| large error reduced | 1 (K3: 26,210 → 6,263) |
| **already-catastrophic prediction made worse** | **2 — K5: 12,075 → 24,240 (+12,165); K1: 39,155 → 65,462 (+26,307)** |
| unchanged | 14 (all 13 near-floor controls; 0 damaged) |
| newly introduced large error | 0 |

**Zero newly introduced catastrophes is not evidence of safety.** Two of the four scored catastrophic predictions got worse.

**Decision:**
- **M12-B's historical decision:** INCONCLUSIVE (unchanged).
- **The frozen rule applied to the augmented results:** catastrophic blocks scored 4, fixed 1, introduced 0. That is
  neither MERITS (needs ≥ 3 fixed) nor FAILS (needs 0 fixed, or introduced ≥ fixed), so the result is **INCONCLUSIVE**.
- **Not claimed:** completing C4 says nothing about ACOT coverage (11 blocks unscored, including K2), and nothing here
  supports promotion.

## 4. Recommendation: **park the rule**

**The specific remaining obstacle is the choice within the correct copy state.**
- **The count-based state step is accurate:** 27 of 29 states correct, and none changed wrongly.
- **Inside the selected state, choose between allocation and content:**
  - **K1:** the unchanged objective picks a 1+3 allocation over the truth's 2+2.
  - **K5:** it picks wrong-content alleles over the floor pair.
  - **Both scores agree in both cases:** the objective by 68 and 112 nats, the count loss by 452 and 2,510.
- **Forcing the state therefore converts a dosage error into an allocation or content error.** Here that was as large or
  larger: 2 of 4 catastrophic cases got worse, and only K4 was fixed.

**No variant of this rule is supported** until a within-state selector exists that, on K1 and K5, prefers the
truth-allocated or floor content over both scores' current choice.

**Consistent with earlier findings** that C4 sv:7 allocation is weakly identifiable from reads (M1-C, M6-B, sv7
allocation).

## Closure (2026-10-21)

**M12-C is closed, and the hard dosage-first rule is parked.**

**Kept as diagnostics only:**
- the count-state outputs: pooled-feature copy states (`count-attrib --features 1`), count-state selection and the
  coverage ledger.

**Not done:**
- the selector is not integrated;
- the scoring campaign is not extended;
- no further dosage weight is tuned.

**Not claimed:** these results do not establish short-read unidentifiability of allocation or content. They show only
that, on these cases, the unchanged objective and the count loss preferred the wrong within-state candidates.

**Unchanged:**
- the caller checkpoint (M10-A + M11-C);
- the 18 ACOT validation families.

No further batch is open.
