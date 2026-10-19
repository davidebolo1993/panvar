# M6-A: bounded 31-mer count-model comparator for structural genotyping (diagnostic prototype; frozen before any unit)

**Ancestry.** Branch `m6-a-count-model` from `3efba08` (M5-A), plus the M5-A closure.

**No caller integration and no production flag.** New code is diagnostic only:
- `tests/count_fit.hpp`, the comparator core;
- `genotype_probe count-fit`;
- `tests/count_fit_gates.cpp`;
- scripts in this directory.

## The frozen observation model

**Loss.** loss(candidate) = Σ_{k∈D} (observed(k)/depth − predicted(k))². It is a count-fit objective, not a calibrated
likelihood:
- never combined with the fragment likelihood or the copying prior;
- no posterior, no GQ, and no independent-k-mer probability claims;
- no learned weights, robust loss, residual editing, sequence reconstruction or rounding.

**Dictionary D.** The canonical 31-mers of the complete spellings of every **active** (curated) panel template. There
are no truth-derived or candidate-derived additions. The marker inventory is not reused.

**predicted(k).** k's multiplicity summed over the candidate's two complete haplotypes:
- every occurrence anywhere in the locus is counted, junction-spanning 31-mers included;
- each haplotype is left + allele + right, counted as base(left, right) + middle(left[−30:] + allele + right[:30]), which
  is exact (gate G1);
- 31-mers a candidate carries outside D are reported (`middle_kmers_outside_dictionary`) and excluded from the loss.

**observed(k).** Counted over every 31-mer window of **both mates** of every read pair:
- each mate independently, with no de-duplication of mate overlap (inserts of about 350 bp exceed 2 × 150 − 31, so
  overlap is rare);
- both strands, through the canonical form;
- any window containing a base other than A/C/G/T is skipped and counted (`read_windows_skipped_non_acgt`).

**Depth: one haploid 31-mer depth per read pool,** depth = Σ observed / (2 · n) over the **stable-copy** k-mers.
- A stable-copy k-mer has exactly one occurrence in every active template's complete spelling, and no occurrence
  overlapping the target block in any template.
- Every genotype the panel can spell therefore carries it twice.
- **Assumption:** the sample also carries two copies. It cannot be verified for a sample outside the panel, such as the
  LOO donor.
- Fewer than 1,000 stable k-mers gives an explicit refusal: no stable scale.
- The median is reported beside the mean but not used.

## Candidates: local substitution only, not whole-locus inference

**Context, fixed and truth-free.** The run's reported winner, i.e. the rank-1 representative of its equivalence set:
| cells | source run |
|---|---|
| LZO | M4-B constant-arm runs |
| LOO | M1-B fresh runs `call.c4_10`, `call.cyp_5` |

**Catalogue.** Every ordered pair (a on homologue 1, b on homologue 2) of the target block's active alleles: n² ordered
candidates.
- Both assignments are always present, and identical profiles are grouped.
- Allele aliases (identical bytes) are already one canonical allele in the model.
- **Preflight:** at most 40,000 ordered candidates, otherwise a refusal (never truncated). Haploid middle vectors are
  computed once per allele and homologue; the loss is exact and incremental.

**Groups.** Candidates are grouped by identical predicted count profile (a linear 64-bit hash over the dictionary). A
multi-member group is **one** count-model state:
- it keeps every member, i.e. every sequence pair and allocation;
- it never establishes a unique allocation or genotype.

## Gates (`gate.count_fit`, synthetic, exact observed counts)

| gate | checks |
|---|---|
| G1 | multiplicities: base + middle equals whole-spelling counts for empty, 5-bp, tandem U³ and 90-bp alleles with differing homologue flanks; repeated copies are counted |
| G2 | a total copy-content change, 2+2 vs 3+2, is detected |
| G3 | 2+2 and 3+1 tandem allocations with identical pooled content have identical profiles and scores, and the best group holds both allocations: UNRESOLVED |
| G4 | an arrangement swap (UV vs VU) whose junctions differ is discriminated |
| G5 | a label alias leaves the dictionary and the scores unchanged |
| G6 | incremental loss equals an independent full-materialisation scalar loss (all candidates, noisy counts), and candidate order has no effect |

No gate assumes every structural difference is visible to 31-mers; G3 asserts the opposite.

## Cells (controls first) and order of work

| # | cell | model | target |
|---|---|---|---|
| 1 | C4 LZO HG00171 5× | `c4.gfa` | sv:7 |
| 2 | CYP2D6 LZO HG02769 5× | `cyp2d6.gfa`, curated | sv:9 |
| 3 | C4 NA19650 LOO 10× | `NA19650.loo.gfa` | sv:7 |
| 4 | CYP2D6 NA19650 LOO 5× | `loo.NA19650.gfa`, curated | sv:9 |

1. **`countfit.<cell>`** (truth-free). Writes all candidates with loss, group, profile, out-of-dictionary 31-mers and the
   caller flag. Its summary holds depth, dictionary, stable k-mers, the caller's group rank, residual split (target-window
   vs whole loss), timings and peak footprint.
2. **`freeze`** (truth-free). A digest over every count output, and the frozen evaluation selection: the caller; every
   best-group member up to 20 in file order (a larger group's remaining members are structurally evaluated but not
   fragment-scored); the first member of groups 2–5.
3. **Evaluation only:**
   - `blockseq.<cell>`: truth block bytes;
   - `units.<cell>`: repeat-unit copies by M1-B's unchanged rule, for every allele and truth haplotype;
   - `fragscore.<cell>`: score-set fragment-model F/E/prior/total and sequence excess for the frozen selection; the
     caller's total must reproduce.
4. **`summary`.** Structural verdicts per cell, rules frozen in `summary.py`:
   - best-group structure (ALL_CORRECT, TOTAL_ONLY, WRONG_TOTAL, MIXED_TOTAL, NOT_ASSESSED);
   - whether the allocation is resolved;
   - whether the truth unit pair exists in the catalogue (the missing-content check);
   - the first group carrying it;
   - total content and allocation, kept separate;
   - the verdict against the caller (IMPROVES, PRESERVES, WORSENS, SAME_ERROR, DIFFERENT_ERROR).

## Bounds

- 600 s counting every attempt, 4 threads where used, 6 GiB envelope.
- **Declared timeouts:**
  - gate 60 s;
  - count-fit: 120 / 90 / 150 / 90 s;
  - freeze 30 s; block-sequences 60 s; units 60 s; score-set 120 s; summary 30 s.
- Builds are not counted.
- Not done: no simulation, campaign, ACOT, merge or push.

## Amendment 1 (declared after `summary` failed and before any further unit; 45.7 s spent)

**What failed.** `summary` exited 1 on a parsing defect. Its arguments are `cell:block`, and block ids contain ":"
(`sv:7`), so `spec.split(":")` returned three parts.

**Fix.** `spec.split(":", 1)`. No verdict rule changes, and nothing else in the script changes. Every count, structure
and fragment unit had already completed under the frozen rules.

**Dry run.** The fixed script was run once on the recorded outputs (to the scratchpad, not counted as a unit) before the
re-run.

**Units.** Manifest A1, `summary.a1` (30 s). `summary` stays recorded as failed.
