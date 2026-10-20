# M9-F result: junction-aware proposals (JP1) break the sv:5 tie, and full scoring then finds real improvements in all three affected cells (sv:5 distance 4→2, 9→2, 43→2 bp) with both LZO controls unchanged; recommend a default-off experiment

**Run.** Branch `m9-f-junction-proposals`.
- **Commits:** freeze `2e6bf2f`; Amendment 1 `28920c3`, analysis only: distances to truth, and overlap NA where no original
  ranking exists.
- **Results:** `tests/results/m9f_junction_proposals_20261020/`.

**Execution: 173.0 of 1,200 s,** including 3 s of disclosed pre-freeze gate development.

**Failures:** none.

**Status:** a diagnostic prototype in Python, with no caller change. The proposal score is a search heuristic, not a
likelihood, confidence value or final-score term.

## Rule and gates

**JP1** (`DESIGN.md`, `jprop.py`):
- **Markers:** all canonical 31-mers in incumbent-context windows L + allele + R, with K − 1 bases each side. This is a
  declared deviation from the caller's syncmers.
- **Expected count:** rate · copies + background, using the run's own depth fit, with copies elsewhere in the incumbent
  locus included.
- **Score:** the Poisson kernel over a fixed inventory.
- **Pair ranking:** the better homologue assignment, top 32 kept, ties explicit.

**Synthetic gates: 11/11 pass** (F1–F9):
- junction-only short allele;
- empty allele;
- K-mers spanning both boundaries;
- repeats elsewhere;
- reverse-complement, homologue-swap and permutation invariance;
- same sequence under different labels;
- uninformative ties;
- the independent oracle agrees to ≤ 6e-14.

**Real-data gates per cell: all pass.**
- The context SHA equals v2's representative.
- The empirical all-K-mer rate is **0.99–1.13×** the fitted rate.
- The fast score equals the whole-haplotype oracle within ≤ 2.3e-13.
- The original sv:5 emissions are **all tied** in every model checked: LOO (M9-E) and both LZO models (1 distinct value
  each).

## Results (sv:5; shortlists frozen before truth or full scoring)

| cell | incumbent → best replacement | truth | Δtotal (Δlik / Δprior) | sv:5 distance to truth: incumbent → best | truth pair rank: index list → JP1 | M9-D improving pair in JP1 top 32 |
|---|---|---|---|---|---|---|
| **LOO d1** | 5/71 → **6/90** | 6/77 | **+16.68** (+38.86 / −22.18) | **4 → 2** | 753 → 222 | no (rank 222) |
| **LOO d3** | 45/107 → **34/93** | 34/94 | **+19.41** (+42.00 / −22.59) | **9 → 2** | 3,444 → **10** | **yes** |
| **LOO d6** | 0/21 → **18/21** | 13/21 | **+17.61** (+16.51 / +1.10) | **43 → 2** | 1,439 → **16** | **yes** |
| LOO d2 | unchanged (best −25.35) | 76/81 | — | 1 (unchanged) | — → 1,541 | — |
| LOO d4 | unchanged (best −26.59) | 59/78 | — | 6 (unchanged) | — → 3,614 | — |
| LOO d5 (sequence-tied context) | unchanged (best −6.72) | 18/102 (not attainable) | — | — | — → 13 | — |
| **LZO d1 (correct)** | **unchanged** (best −35.39) | 6/77 | — | 0 | — | — |
| **LZO d2 (correct)** | **unchanged** (best −41.40) | 76/81 | — | 0 | — | — |

**Incumbents.** Every incumbent's v2 representative score reproduced exactly; score-set's expectation held in 8/8 cells.

**Changes.** Every proposed change is a **real sequence change**: different unique haplotypes, not a label-only
relabelling.

**Number of candidates improving on the incumbent:** 3, 16 and 13 in d1, d3 and d6; 0 in the other five cells.

**Shortlist diversity:**
- 14–17 distinct alleles per shortlist;
- 3,962–5,038 distinct pair scores, against 1 originally;
- the top tie class has size 1–2, and the cutoff tie class 1–5. In d1, d2 and LZO d2 the cutoff tie crosses rank 32 and
  is reported, not hidden;
- overlap with the original index-ordered list is **0** wherever it was computed (d1, d3, d6, LZO d1, LZO d2).

**Cost per cell:**
- JP1 proposal about 1 s (Python prototype);
- full scoring of about 64 candidates (28–36 unique haplotypes) 12–17 s;
- peak memory 1.12–1.42 GiB (score-set).

## Reading

1. **Proposal.** JP1 replaces an index-ordered list with read-evidence rankings at a block where the existing emissions
   carry no information.
   - It lifts 2 of 3 known improving pairs into the top 32 (ranks 10 and 16).
   - The third rises from 753 to 222.
2. **Search usefulness under the unchanged final scorer.** In all three affected cells the new shortlist contains a
   replacement that the original objective prefers to the reported winner by 16.7–19.4 nats.
   - Each moves sv:5 **closer** to truth, ending 2 bp away in every case.
   - None of these replacements was reachable before.
3. **Controls and the other cells.**
   - Neither correct LZO control would change: the best alternative is 35–41 nats worse.
   - LOO d2, d4 and d5 are unchanged.
4. **What remains unresolved:**
   - **The objective, not the proposal, in d3 and d6.** In d3 the exact truth pair (rank 10, in the list) scores below
     34/93, which is 2 bp away: M9-D's truth replacement gained +8.46, against +19.41 here. In d6 the truth pair (+17.55
     in M9-D) and 18/21 (+17.61) are within 0.07 nats.
   - **Proposal reach in d1.** The exact truth pair (+52.05 in M9-D) is still outside the list, at rank 222.

**Not supported:**
- that this improves genotyping accuracy in general;
- that it calibrates anything.

**Limitations:**
- One block, eight inspected development cells, one locus, one coverage.
- The heuristic was designed after seeing these failures.
- A Python prototype is not the caller's incremental machinery.

## Recommendation

**Integrate as a default-off experiment.**

| criterion | evidence |
|---|---|
| recovery | 3/3 affected cells gain a real, objective-preferred improvement toward truth; 2/3 known improving pairs enter the list |
| control behaviour | 0/2 LZO controls change, and 0/3 other LOO cells change |
| measured cost | small: about 1 s per proposal in Python, and the shortlist size is unchanged at 32, so no extra neighbours |

**Scope of the experiment:**
- only for blocks whose existing proposal emissions are entirely tied;
- replacing their index-ordered local list;
- sample-scoped, with caches keyed to the flanking sequences;
- the persisted locus index untouched.

**What it needs next:**
- an identity gate showing that "off" is unchanged;
- a held-out evaluation on the remaining ACOT families, with a predeclared criterion.

**Not proposed:** a wider list, prior changes or joint search.
