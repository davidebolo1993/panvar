# M6-B: count-profile correctness and C4 content-versus-allocation attribution (diagnostic; frozen before any unit)

**Ancestry.** Branch `m6-b-count-attribution` from `5ff7353` (M6-A).

**No caller change and no scoring rule.** The loss, dictionary, normalization and candidate catalogue are M6-A's,
unchanged.

## 1. Exact grouping

- **Grouping (`tests/count_fit.hpp` `group_profiles`).** The 64-bit hash is a lookup key only. Inside a bucket,
  candidates are split by equality of the complete predicted profile M1[a] + M2[b]; the base is shared. The tie order is
  M6-A's: loss, then profile hash.
- **New gates.**
  - G7: with the hash mask set to 0, every candidate shares one bucket. Groups must then equal the distinct complete
    profiles and the default partition, the forced bucket must report a collision, and equal profiles (2+2/3+1/1+3) must
    stay grouped.
  - G8: swapping the two complete haplotypes, flanks included, leaves loss and profile unchanged. A target swap under
    different flanks changes the profile.
- **`countfit.<cell>`** reruns M6-A's four cells with exact grouping. Its summary adds `hash_collision_buckets`.
- **`verify.groups`** compares every candidate's loss, group rank and group size with M6-A, and re-derives the frozen
  evaluation selection by M6-A's rule. Every difference is reported; no collision is presumed.

## 2. C4 NA19650 LOO 10× attribution (oracle diagnostics; none enters a prediction catalogue)

**Fixed context.** M1-B `call.c4_10`, rank-1 representative.

**Candidates (`candidates.tsv`, frozen).**
| label | content | oracle |
|---|---|---|
| reported | 3+1, alleles 41/47 | no |
| floor2x2 | M6-A group 4 first member, alleles 37/46 | no |
| truth_a, truth_b | exact truth sv:7 blocks substituted into the context, both assignments | **yes** |
| `m1e_*` (14) | every M1-E explicit construction (content-matched a22 2+2 / e13 1+3, single and all truth corrections, truthblock, reported) | **yes** |

- **Binding check for M1-E constructions:** bytes outside the recorded sv:7 interval must equal the context, or the
  candidate is marked unbound.
- **Reused fragment scores:** M6-A score-set for the two catalogue members; M1-E score-set for its constructions (same
  run, model, reads and library). truth_a/truth_b are matched to M1-E's truthblock by bytes. Nothing is recomputed.

**`count-attrib`.**
- The M6-A dictionary and stable-copy depth, recomputed by the same construction and checked against M6-A. Oracle
  candidates' out-of-dictionary target 31-mers are reported, never added to the dictionary.
- **Profiles:**
  - real: reads / depth;
  - **oracle 1:** noiseless counts of truth_a (truth block in the fixed caller context), unit depth;
  - **oracle 2:** noiseless counts of the complete true diplotype (`spell --truth-gfa`), unit depth.

  These are not predictions, simulated-read benchmarks or calibration data.
- **Contrasts:** reported−floor2x2, reported−truth_a, reported−truth_b, floor2x2−truth_a, e13_none−a22_none,
  e13_all−a22_all, e13_T1−a22_T1, reported−a22_none, a22_none−a22_all, reported−m1e_reported (binding identity).
  - Each is decomposed over the k-mers whose predicted counts differ, and the contributions must reconstruct the
    difference.
  - **Disjoint categories,** from the union of occurrence locations in both candidates: target_interior_only,
    target_junction_only, shared (more than one location class, i.e. target plus elsewhere, or interior plus junction),
    outside_only.
  - **Concentration:** shares of the top 10, 100 and 1,000 k-mers.
- **Catalogue under each profile:** the unchanged active catalogue (all ordered allele pairs, same context, exact
  groups) is scored against each profile.
- **Normalization check:** the multiplicity of the stable k-mers in the complete true diplotype, without refitting depth.

**Causal reading (frozen in `summary.py`):**
- SAMPLING_OR_NORMALIZATION;
- RESIDUAL_CONTENT_MISMATCH;
- CONTEXT_EFFECT;
- COUNT_INDISTINGUISHABLE;
- UNRESOLVED.

A k-mer inside the block is never read as independent allocation evidence, and the same reads are never counted twice.

## Bounds

- 300 s counting every attempt, existing envelope.
- **Declared timeouts:** gate 60 s; count-fit 60 s each; verify 30 s; spell 60 s; attrib 120 s; summary 30 s.
- Not done: no new cells, k change, robust weights, production patches or combination with the fragment model.
