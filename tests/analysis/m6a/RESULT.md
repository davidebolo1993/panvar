# M6-A result: a simple 31-mer count-fit comparator, four development cells

**Run.** Branch `m6-a-count-model`: freeze `bcfa921`, Amendment 1 (summary parsing). Results in
`tests/results/m6a_count_model_20261019/`.

**Budget.** 45.9 of 600 s. 20 units; one recorded failure (`summary`, a script parsing defect, re-run as `summary.a1`).

**No caller change.** The comparator is a diagnostic probe (`genotype_probe count-fit`) with a shared core
(`tests/count_fit.hpp`).

## Model, as frozen

**Loss.** Σ_{k∈D} (observed/depth − predicted)² over the canonical 31-mers of the active panel's complete spellings. It
is a count-fit objective: not a likelihood, never combined with the fragment model or the prior, with no posterior or
GQ.

**Counts.**
- **Predicted:** whole-haplotype multiplicities, junctions included, decomposed exactly (gate G1).
- **Observed:** both mates, canonical; non-ACGT windows skipped.
- **Depth:** mean over stable-copy k-mers / 2.

**Candidates.** Every ordered allele pair at the target block inside the caller's own reported winner as a fixed
context; identical profiles are grouped.

| cell | dictionary k-mers | stable-copy k-mers | depth (mean / median) | candidates | profile groups | wall | peak footprint |
|---|---|---|---|---|---|---|---|
| C4 LZO HG00171 5× | 242,402 | 126,643 | 1.65 / 1.5 | 14,161 | 7,140 | 3.6 s | 0.23 GB |
| CYP2D6 LZO HG02769 5× | 255,263 | 141,115 | 1.66 / 1.5 | 12,996 | 6,554 | 3.1 s | 0.23 GB |
| C4 NA19650 LOO 10× | 242,037 | 126,868 | 3.29 / 3 | 13,689 | 6,903 | 3.4 s | 0.23 GB |
| CYP2D6 NA19650 LOO 5× | 254,909 | 141,325 | 1.65 / 1.5 | 12,544 | 10,531 | 3.1 s | 0.23 GB |

- **No catalogue was truncated;** the largest was 14,161 ordered candidates, against a preflight bound of 40,000.
- **Phase is never resolved:** both assignments of an unordered pair always share one profile group.
- **The absolute loss is dominated by low-depth sampling noise** (about 150–280k over about 240k k-mers). Candidate
  differences of tens to thousands are the signal.

## Gates (synthetic, exact counts): all pass

- G1: junction-inclusive multiplicities and repeated copies.
- G2: a change in total copy content is detected.
- **G3: 2+2 and 3+1 with identical pooled content give identical scores and one group, so the allocation is
  UNRESOLVED.**
- G4: an arrangement with distinguishing junctions is separated.
- G5: an alias does not multiply evidence.
- G6: the loss equals an independent scalar calculation (relative error 1.8e-15), and candidate order has no effect.

## Results (count selections frozen before truth)

Structure is the unit pair from M1-B's rule, evaluation only. Excess is sequence excess at the target block over the
block's availability floor.

| cell | truth units | caller | count-best group (members) | count-best structure | count-best block excess | caller block excess | caller's count rank (Δloss) | first group with the truth unit pair | verdict |
|---|---|---|---|---|---|---|---|---|---|
| C4 HG00171 5× (control) | 2+2 | 2+2 | the caller (2: both assignments) | 2+2 correct | 0 | 0 | 1 (0) | 1 | **PRESERVES** |
| CYP2D6 HG02769 5× (control) | 2+2 | 2+2 | another 2+2 pair (2) | 2+2 correct | **14** | 0 | 5 (+42.9) | 1 | structure preserved; **sequence worse** |
| C4 NA19650 LOO 10× | 2+2 | 3+1 (allocation error) | the caller's 3+1 (2) | 3+1, allocation error | 65,444 | 65,444 | 1 (0) | **4** (+301.8) | **SAME_ERROR** |
| CYP2D6 NA19650 LOO 5× | 3+2 | 2+2 (dosage error) | a 3+2 pair (1) | **3+2 correct** | **0** | 12,159 | 164 (+4,056.8) | 1 | **IMPROVES** |

### Fragment-model comparison on the same candidates (score-set; the caller's total reproduced)

**CYP2D6 NA19650 LOO 5×.** The count-best 3+2 candidate scores 31.47 **below** the caller under the fragment model:
| term | 3+2 − caller |
|---|---|
| fragment | **+127.1** |
| exposure | **−151.9** |
| prior | −6.7 |

This is M1-B's picture: the fragment gain is cancelled by the exposure cost, and the prior decides. The count model's
separate depth fit sees the dosage decisively (Δ 4,056.8).

**C4 NA19650 LOO 10×.**
- **Both models reject a candidate that matches the truth at sv:7.** A 2+2 candidate in count group 4 reaches the sv:7
  block availability floor (distance 36 = floor, excess 0). The fragment model also scores it 425.9 below the caller.
- **The error is not missing block content.** The best available sv:7 content is in the catalogue, and both models
  prefer the wrong 3+1.
- **The cause stays unresolved.** It may lie outside the block, in the fixed caller context or in sequence the leave-one-out
  panel cannot spell elsewhere; the locus-wide availability floor is 113.

**CYP2D6 HG02769 5× (control).**
- The count-best pair has the correct units but block excess 14. The caller's own pair (excess 0) ranks 5th (+42.9),
  and the fragment model prefers the caller by 292.3.
- Within one structure class, the count objective does not separate the correct sequence from near neighbours.

## Decision

| question | answer |
|---|---|
| Improves structural preference? | **In one of two error cells.** It fixes CYP2D6 NA19650 dosage (3+2, sv:9 block excess 0). It does not fix C4 NA19650 allocation. |
| Preserves controls? | **Structure, yes** (2/2 unit pairs correct). **Sequence, no** for CYP2D6 HG02769: the best group is a different 2+2 pair with excess 14, where the caller is exact. |
| Cannot distinguish relevant allocations? | **In principle, yes** (G3: identical pooled content gives identical scores). On the real C4 cell the candidates are not count-tied: the model separates 3+1 from 2+2 by about 302, in the wrong direction. Phase (assignment) is never separated. |
| Biased by missing content? | **Not demonstrated at block level.** In both LOO cells the truth unit pair and floor-attaining sv content are in the catalogue. The C4 preference for 3+1 is shared with the fragment model, and its cause (context, flanks, content elsewhere) is unresolved. |

**Overall.** The simple count model is a **complementary dosage signal**, not a replacement.
- It recovers the dosage error where the fragment-likelihood combination fails (exposure vs fragment vs prior).
- It does not fix the C4 allocation error.
- It is too coarse to pick the exact sequence within a structure class.

**Limitations.**
- Four cells, one coverage each.
- Local substitution at one block in the caller's context, not whole-locus inference.
- Depth assumes the sample's stable-copy k-mers are diploid.
- The loss is uncalibrated and noise-dominated at low depth.

**Not done, per the brief:** no further loss, k change or Ctyper/Locityper replacement was added.

## Qualifications appended at M6-B (2026-10-19; Sol's review; the record above is not rewritten)

1. **The floor is the closest available sequence, not exact content.** "Floor-attaining 2+2" means the candidate reaches
   the block's availability floor: the closest available sequence (sv:7 distance 36 = floor, not 0). It does not mean the
   exact truth content is in the catalogue. "Not missing block content" above should read "the closest available block
   content is in the catalogue".
2. **M1-C/E's mechanism is not excluded.** M1-C and M1-E showed that small residual sequence changes can substantially
   change the allocation comparison: in M1-E, the 23-bp T1 correction alone moved 2+2 vs 1+3 by 7.6 at 10×. M6-A does
   not exclude that mechanism.
3. **Phase is not always tied.** Swapping the target alleles between fixed, *different* flanks changes the junction
   counts. CYP2D6 NA19650's best count group has a single member: (22,77) and (77,22) score differently. The statement
   "phase is never resolved: both assignments always share one profile group" is wrong in general. It held where the two
   context flanks are identical around the target.
4. **Swapping the two complete haplotypes is invariant** (M6-B gate G8). It is a different operation from swapping
   target alleles inside fixed contexts.
5. **A unique local count profile does not establish whole-locus phase.**

**Unchanged:** the positive CYP2D6 NA19650 dosage result (count-best 3+2, sv:9 block excess 0), and the sequence-level
control regression (CYP2D6 HG02769: count-best block excess 14 where the caller is exact).
