# M6-B result: exact grouping changes nothing; the C4 count preference is residual content, and the allocation is count-invisible

**Run.** Branch `m6-b-count-attribution`, freeze `7bec866`. Results in `tests/results/m6b_count_attribution_20261019/`.

**Budget.** 27.0 of 300 s. 9 units, all exit 0 and first-time, with no amendments.

**No caller change and no scoring rule.** The loss, dictionary, normalization and catalogue are M6-A's.

## 1. M6-A qualifications (appended to its record)

- **The floor is the closest available sequence, not exact content.** The floor-attaining 2+2 has sv:7 distance 36,
  which equals the floor; it is not 0.
- **M1-C/E's residual-sequence mechanism was not excluded by M6-A.** It is confirmed below.
- **Target swaps can change junction counts.** Swapping target alleles between fixed, *different* flanks changes junction
  counts (CYP2D6 (22,77) vs (77,22)), so "phase is never resolved" was wrong in general.
- **Swapping the two complete haplotypes is invariant** (gate G8).
- **A unique local count profile does not establish whole-locus phase.**
- **Unchanged:** the CYP2D6 NA19650 dosage result and the CYP2D6 HG02769 sequence-level control regression.

## 2. Exact grouping

**Change.** The hash is now a lookup key only; complete-profile equality decides groups.

**Gates.**
- **G7, forced collisions** (hash mask 0, every candidate in one bucket): groups = the 14 distinct complete profiles,
  the same partition as the default hash, the collision is reported, and equal 2+2/3+1/1+3 profiles stay grouped.
- **G8:** swapping complete haplotypes is invariant; a target swap under different flanks changes the profile.

**Re-verification of M6-A's four catalogues.**

| cell | candidates | loss / group-rank / group-size differences | groups (M6-A = exact) | hash-collision buckets | caller rank | evaluation selection |
|---|---|---|---|---|---|---|
| C4 HG00171 5× | 14,161 | 0 / 0 / 0 | 7,140 | 0 | 1 | identical |
| CYP2D6 HG02769 5× | 12,996 | 0 / 0 / 0 | 6,554 | 0 | 5 | identical |
| C4 NA19650 10× | 13,689 | 0 / 0 / 0 | 6,903 | 0 | 1 | identical |
| CYP2D6 NA19650 5× | 12,544 | 0 / 0 / 0 | 10,531 | 0 | 164 | identical |

**Verdict: UNCHANGED.** No collision occurred, and every M6-A grouping and selection stands.

## 3–5. C4 NA19650 LOO 10× attribution (oracle diagnostics; nothing enters a prediction catalogue)

**Setup.**
- The dictionary (242,037 k-mers), the 126,868 stable k-mers and depth (3.289) reproduce M6-A exactly.
- All 14 M1-E constructions pass the binding check: their flanks equal the fixed context byte for byte.
- `m1e_reported` reproduces `reported` exactly.
- `truth_a` is byte-identical to M1-E's truthblock.
- The truth's target block has **124** 31-mers outside the panel dictionary. They are reported, not added.

### Count loss under the three profiles

real = reads / depth; oracle 1 = noiseless truth block in the caller context; oracle 2 = noiseless complete true
diplotype.

| candidate | real | oracle 1 | oracle 2 | fragment-model likelihood (F + E, reused) |
|---|---|---|---|---|
| reported 3+1 | 149,971.3 | 1,175 | 1,971 | −153,266.7 |
| floor2x2 (catalogue 2+2) | 150,273.0 (+301.8) | 1,604 | 2,400 | −153,673.8 |
| M1-E a22_none / e13_none (content-matched from panel units) | 149,925.5 (−45.7) | 1,364 | 2,160 | −153,250.6 (both) |
| M1-E a22_T1 / e13_T1 | 149,057.6 | 279 | 1,075 | −153,099.7 / −153,107.3 |
| truth_a = truth_b = a22_all = e13_all = truthblock | 148,786.1 | **0** | 796 | −152,986.3 (a22_all) / −152,993.9 (e13_all) |

### Contrast decompositions

Contributions reconstruct every difference exactly.

| contrast | real Δ | oracle 1 Δ | differing k-mers | categories |
|---|---|---|---|---|
| reported − floor2x2 | **−301.8** | **−429** | 1,719 | **all target-interior-only**; junction 0, shared 0, outside 0 |
| reported − truth_a | +1,185.1 | +1,175 | 710 | all interior |
| floor2x2 − truth_a | +1,486.9 | +1,604 | 1,511 | all interior |
| e13_none − a22_none | **0** | **0** | **0** | identical profiles |
| e13_T1 − a22_T1 | **0** | **0** | **0** | identical profiles |
| e13_all − a22_all | **0** | **0** | **0** | identical profiles |
| reported − a22_none | +45.7 | **−189** | 1,113 | all interior |
| a22_none − a22_all | +1,139.4 | +1,364 | 1,271 | all interior |

- **Concentration is broad.** The top 10 k-mers carry about 2–6% of the summed abs Δ, and the top 1,000 carry 65–100%.
- **No single molecule or k-mer decides.**
- **Interior k-mers are not independent allocation evidence.** All candidates share the junction content because the
  blocks are divergence-trimmed, so every difference is interior, and interior k-mers in repeated modules are not
  allocation evidence.

### The unchanged catalogue under each profile

| profile | best group | reported 3+1 rank | floor2x2 rank | first group holding a 2+2 |
|---|---|---|---|---|
| real | 3+1 and 1+3 (2 members) | 1 | 4 | 4 |
| oracle 1 (truth block in context) | 3+1 and 1+3 | 1 | 4 | 4 |
| oracle 2 (complete true diplotype) | 3+1 and 1+3 | 1 | 4 | 4 |

**Normalization.** In the complete true diplotype the stable k-mers have mean multiplicity **1.998**, so the implied
depth bias is 0.999. The diploid assumption holds to 0.1% here; depth was not refitted.

## 6. Causal conclusion

**Residual sequence-content mismatch drives the count preference.**
- Even with **no sampling noise** and with the **true block** in the caller's context (oracle 1), the active catalogue
  ranks the reported 3+1 first.
- The floor-attaining 2+2 differs from the truth's counts at 1,511 interior k-mers, against 710 for the 3+1, and loses
  by 429.
- The panel's available 2+2 alleles carry module content further from NA19650's than the 3+1 alleles do. This is M1-C/E's
  mechanism, now shown for the count objective.

**Within matched content, the allocation is count-indistinguishable.**
- Content-matched 1+3 and 2+2 constructions have **identical** 31-mer profiles: zero differing k-mers, with or without
  the T1 correction or all four corrections.
- The fragment model separates them only through long-range adjacency once T1 is in: 7.6 at 10× (M1-E), and 7.62 in the
  reused scores above.
- The count model cannot see the allocation at all here.

**Not sampling or normalization.**
- The preference survives the noiseless profiles, and the normalization holds to 0.1%.
- One nuance: the content-matched a22_none (not in the catalogue) beats the reported pair on real counts (−45.7) but
  loses on noiseless truth (−189). That one contrast is sampling-sensitive, but it does not change the catalogue's
  conclusion.

**Not a fixed-context or shared-k-mer effect.**
- Every differing k-mer is target-interior-only.
- Replacing the caller context with the complete true diplotype (oracle 2) shifts every loss by the same 796 and changes
  no ranking.

**Unresolved:** nothing about the direction. The 124 truth 31-mers absent from the panel dictionary cannot inform any
candidate.

**Implication, stated, not acted on.** For C4 NA19650 sv:7 the count model faces two problems:
- **content:** the catalogue lacks NA19650's module variants;
- **observability:** the allocation is invisible to 31-mer counts given matched content.

Neither a better loss nor a different search over the existing catalogue addresses the first. The second needs
long-range adjacency evidence (M1-E). No k change, weighting, production patch or combination with the fragment model
was made.
