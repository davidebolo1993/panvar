# M9-A: read-only uncertainty-target audit

Separates uncertainty about the exact sequence pair, the unphased block genotypes and phase, ahead of a simpler production
confidence recommendation.

**Inputs:** the saved outputs of the eight completed M8-C/M8-D runs.
- M8-C: HG03209 LZO/LOO and HG01975 LZO, each at r0.01 and r0.1.
- M8-D: HG01975 LOO at r0.01 and r0.1.

**Constraints:**
- No new placement, genotyping, model, clustering, threshold, filter, PASS rule or GQ calibration.
- Historical packages are not modified.
- Results: `tests/results/m9a_uncertainty_audit_20261020/`.

## Hypotheses and aggregation

**A hypothesis** is the caller's sequence-level group: one unordered whole-haplotype sequence pair, identified by
`sequence_candidate_id` in `mosaic_scores.tsv`.
- **Members** are the label pairs (rows of `mosaic_candidates.tsv`) that spell it.
- **Each hypothesis is counted once:**
  - a whole-haplotype swap of a label pair is the same member;
  - a label pair listed in two hypotheses is rejected as inconsistent input.
- **Mass** is the caller's own normalised score, recomputed from the saved scores and checked against
  `candidate_set_probability`. It is conditional on that run's evaluated pool.
  - Pools differ between runs (hypothesis counts are reported).
  - Differences in mass between runs are not calibration evidence.

**Classification for each target and hypothesis**, with the existing mixed-projection semantics of
`genotype_block_support.cpp`:

| class | rule |
|---|---|
| SAME | every member carries the reported state |
| DISAGREEING | no member carries it |
| MIXED | some members carry it and some do not |

- **Reported mass** is the interval [SAME, SAME + MIXED].
- **A MIXED hypothesis** is never assigned to a representative and never discarded.
- **The best genuine competitor** is the highest-scoring DISAGREEING hypothesis. For it the audit reports:
  - the total score gap (log posterior odds within the pool);
  - separately, the likelihood gap (fragment + exposure);
  - the fragment, exposure and prior components.
- **Conditional posterior mass and likelihood support are kept apart.**
- **No disagreeing hypothesis:** NO_EVALUATED_ALTERNATIVE.

## Targets

| target | reported state | predicate per member |
|---|---|---|
| A | the top sequence pair | is this hypothesis (no MIXED possible) |
| B | unordered allele pair at every block | equal at all 19 blocks |
| C | **reportable-block genotype agreement**: unordered allele pairs at the blocks with `reportable = true` (the 9 SV blocks) | equal at those blocks |
| D | each block's unordered allele pair | equal at that block |
| E | phase **given B**: label-level phase vector over the het blocks of B | carries B with the reported phase |

**C is label agreement at reportable blocks.** It is not biological structural equivalence: similar SV alleles are not
clustered, and equal gene structure is not established.

**E, phase given B:**
- **Orientation:** the haplotype carrying B's smaller allele at the first het block is 'A', and the phase vector is A's
  allele at each het block.
- **Conditional mass:** the mass of hypotheses carrying B with the reported phase vector, divided by the mass carrying B,
  as an interval when MIXED hypotheses exist.
- **Competitor:** the best hypothesis whose members all carry B with another phase.
- **Adjacent phase:** for each pair of consecutive het blocks, the conditional mass of the reported relative orientation.
- **Applicability:** requires two or more het blocks and an unambiguous winner. E is reported apart from genotype
  confidence.

**Existing phase support:** the caller's `phase_support` (READS_WITHIN_CANDIDATES vs PRIOR_ONLY, from fragment contenders)
is reproduced and reported next to E, as likelihood support rather than mass.

## Reproductions (definitions as in the code; discrepancies explained, never forced)

Normalisation:
- sum of the recomputed probabilities;
- largest absolute difference from `candidate_set_probability`.

Existing quantities:

| quantity | reproduced as |
|---|---|
| locus genotype separation (`caller_genotype_separation_gap`) | gap to the best non-winner hypothesis not carrying B only |
| mass with `carries_called_genotype = true` | the B SAME mass |
| `numeric_separation_gap` | the A gap |
| `winner_candidate_set_probability` | the A mass |
| `fragment_contenders` and `phase_support` | the rival rule in `genotype_repair.cpp` |
| per-block `conditional_genotype_posterior` | D mass (the caller uses each hypothesis's representative genotype; identical when no hypothesis is MIXED) |
| per-block `genotype_support` (READS / PRIOR_ONLY) | the fragment-contender genotype rule |

**Known sources of difference:**
- `mosaic_scores.tsv` prints 12 significant digits, about 1e-7 absolute at these magnitudes, while the locus table prints
  17.
- The separation definition counts MIXED hypotheses as "other genotype" (`!carries_only`), whereas the audit's genuine
  competitor excludes them. These are identical when no hypothesis is MIXED.

## M8-D qualification

**Values to reproduce,** from Sol's independent check:
- locus genotype separation 0.7021371507 / 0.9018471817;
- reported-unphased-genotype mass ≈ 0.405311 / 0.313111.

**Also checked:** whether the rank-2 sequence pair (the 0.04-nat gap quoted in M8-D) carries the reported unphased genotype
at every block.

**Where it is recorded:** a qualification is **appended** to `tests/analysis/m8d/RESULT.md`. The M8-D package README is
left unchanged.

## Gates (`audit.py gates`; all must pass before `audit`)

| fixture | what it checks |
|---|---|
| F1 | whole-haplotype swapping is counted once |
| F2 | same block genotypes with different phase: B = 1, and A/E split with the phase competitor |
| F3 | genuinely different block genotypes: competitor, total and likelihood gaps, per-block and C competitors |
| F4 | identical sequence under different labels: one hypothesis, mass once; agreeing blocks are not MIXED, disagreeing blocks are |
| F5a | a label pair in two hypotheses is rejected |
| F5 | mixed projection: interval [SAME, SAME + MIXED], and the competitor excludes the MIXED hypothesis |
| F6 | absent disagreeing alternatives give NO_EVALUATED_ALTERNATIVE |
| F7 | normalisation, and D equals the caller's block posterior formula |

## Execution

**Budget: 300 s of execution in total.**

**Before freezing (disclosed and counted): 98.5 s**
- a full audit dry run: 91.0 s, before the genotype precomputation made it faster;
- a one-run smoke test of HG01975 LOO r0.01: 7.1 s;
- gate dry runs: 0.4 s.

**Manifest P:**

| unit | timeout |
|---|---|
| `gates` | 10 s |
| `audit` | 120 s |

Admission rule: pre-freeze time + spent + timeout ≤ 300.
