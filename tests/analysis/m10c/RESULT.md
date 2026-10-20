# M10-C result: no existing alignment method in the repository supports every required structural distinction; the saved cases were therefore not inspected; not ready for a frozen broader evaluation

**Run.** Branch `m10-c-sequence-evaluation`.
- **Commits:** freeze `e2a467d`; Amendment 1 `326d664`; Amendment 2 `bfe862a`; result after.
- **Results:** `tests/results/m10c_sequence_evaluation_20261020/`, containing `ledger.tsv` and
  `runs/{fixtures,fixtures_b.a1,fixtures_c.a2}/`.
- **Budget:** 70.1 of 900 s, including 65 s of disclosed development (30 pre-freeze + 35 for the amendments).
- **Closed:** M10-B / B.1. The long-node skeleton is rejected as a production structural alphabet, and its outputs are
  preserved. M10-A remains the computational checkpoint.
- **Not done:** no reads, genotyping, prediction changes, calibration, priors or defaults.

## 1. Machinery inspected and reused

| need | existing implementation |
|---|---|
| whole-locus truth | `project_path_onto_locus`, documented as "EVALUATION SUPPORT, never inference", via `genotype_probe block-sequences --truth-gfa --donor`. It does **not** require truth to be model-spellable, and it gives matched chain endpoints by construction (the interior, blocks 1..n−2) |
| method A | `panvar::nw_edit_distance`: edlib NW with path, the `benchmark` round-trip-QV routine |
| method B | `panvar::minimap2_hits`: asm20, affine gaps, all regions |
| method C | `panvar::poa_msa`: abPOA, global, affine gaps, two rows |

`m10c_align_eval` is a thin diagnostic over those three. No aligner was written.

## 2. Frozen fixtures × methods (21 fixtures; F10c added by Amendment 1)

| fixture | required | A (edlib, unit cost) | B (minimap2) | C (abPOA, affine) |
|---|---|---|---|---|
| F01 identical DNA (node and boundary representation irrelevant) | CONCORDANT | yes | yes | yes |
| F02 neutral SNPs + small indels | CONCORDANT | yes | yes | yes |
| F03 55 bp deletion | DISCORDANT | **no** (gap fragmented, longest run 8) | yes | yes (D:55) |
| F03b 56 bp deletion in (CA)n | DISCORDANT | yes | yes | yes (placement-sensitive) |
| F04 copy number 3 vs 2 | DISCORDANT | **no** (longest run 24 of 600) | **no** (one hit over 2,000 of 3,400 bp) | yes (I:600) |
| F05 3+1 vs 2+2 | DISCORDANT | yes | **no** (partial hits) | yes, with `DOSAGE_BALANCED_ALLOCATION_DIFFERS` (I:600 / D:600) |
| F06 order A,B → B,A | DISCORDANT | **no** | **no** (one 601 bp block) | **no** (233 mismatches, longest gap 10) |
| F07 600 bp inversion | DISCORDANT | **no** | yes (reverse hit) | **no** |
| F07b 60 bp inversion | DISCORDANT | **no** | **no** (absorbed as mismatches) | **no** |
| F08 / F08b shifted tandem or microsatellite representation (same DNA) | CONCORDANT | yes / yes | yes / yes | yes / yes |
| F09 novel SNPs absent from any panel | CONCORDANT | yes | yes | yes |
| F10 class site declared on **one** of two identical truth haplotypes | DISCORDANT | AMBIGUOUS | AMBIGUOUS | AMBIGUOUS |
| F10b the same change at an undeclared position | CONCORDANT | yes | yes | yes |
| F10c class site declared for the locus | DISCORDANT | — | yes | yes |
| F11 assignment tie (S,S) vs (S+55,S) | DISCORDANT | **no** | yes | yes |
| F11b one of three near-identical copies deleted | DISCORDANT | yes | yes | yes (placement-sensitive) |
| F12 truth 200 bp past the anchor | UNSUPPORTED | yes | yes | yes |
| F12b truth lacks 20 bp at the start | UNSUPPORTED | **no** (noted only) | **no** | **no** |
| F13 equal-length 55 bp divergent replacement | definition required | CONCORDANT (29 edits) | CONCORDANT | CONCORDANT (16 mismatches, 7 + 7 indel bp) |
| F14 60 bp near-tandem duplication | DISCORDANT | **no** | yes | yes |

**Notes on the table:**
- F10's ambiguity is caused by my fixture, not by the method. With identical truth haplotypes, the assignment is arbitrary,
  and an ambiguous verdict is the correct response.
- F12b shows that unmatched endpoints under 50 bp are only noted by alignment. That is why endpoint matching for real cases
  is **by construction** (the same chain boundaries), not inferred.
- **Gate:** the aligner is deterministic, the edit distance is reversal-symmetric, and identical sequences give 0 edits.

**The methods against the pre-declared conditions:**
- **A:** the edit distance is exact (55, 600, 637, 315, 60) but events cannot be located.
- **B:** fails copy number, allocation, order and the small inversion. Its repeat handling is partial: a minimap2 hit set
  does not establish allocation, as warned.
- **C:** status **NONE**, because it fails order (F06) as well as inversions. Under the frozen rule, **no saved case was
  inspected**: `truth.*` and `cases` are recorded NOT RUN in P, A1 and A2.
- **Runner defect, recorded.** The A2 truth units carry the stale Amendment 1 "superseded" reason, because that check
  fired first. The governing reason is Method C's NONE status. The outcome is unchanged: nothing ran.

## 3. Simplest supported procedure, and its limits

**There is no procedure that the fixtures support for the full required set.** What each piece supports, under the
whole-haplotype assignment with Method A's exact edit distance:

| component | supported on the fixtures | not supported |
|---|---|---|
| A: full sequence distance | exact distance per homologue and per assignment; assignment ties exposed | locating any event |
| C: affine global alignment | **indels ≥ 50 bp, repeat copy-number changes, and 3+1 vs 2+2 allocation** (opposite-sign events with a balanced diploid net); located class sites; placement ambiguity flagged | **rearrangements**: order changes and inversions, both silently CONCORDANT; endpoint offsets under 50 bp |
| B: minimap2 hits | large inversions (reverse-strand hit) | repeats (partial alignment), order, small inversions |

**A narrower evaluator** ("indel / copy-number / allocation only; rearrangement-blind") could be built from A + C. **It
was not pre-declared as a target**, and declaring it after the result would be fitting the target to the method. It is
offered here only as an option for approval.

**Dependencies:**
- the truth projection, which requires the truth donor's path in a full graph whose bubbles span the model's chain (the
  `block-sequences` fingerprint check);
- the chain-interior endpoint convention;
- the 50 bp event convention;
- abPOA's alignment (affine but heuristic, placement-sensitive);
- Method A for the assignment.

**Also inherent:**
- placement within repeats is not unique (PLACEMENT_SENSITIVE on F03b, F04, F05, F11, F11b, F14);
- "which copy" allocation is undefined when copies differ only by copy-specific SNPs.

**Separation from runtime confidence:** every result here is a pairwise comparison between one prediction and one truth. It
is **not** a transitive equivalence class, and it must never pool posterior mass.

## 4. Biological definitions that still need Davide

1. **Rearrangements.** Are order changes and inversions part of the structural genotype? From what size? (The repository
   has no method that detects them reliably.)
2. **Event-size semantics.** Is the 50 bp convention the structural threshold? Is an equal-length divergent replacement (F13)
   structural?
3. **Endpoints.** Approve the chain interior as the evaluated locus (matched by construction), and the handling of truth
   whose projection is missing at a block.
4. **Class-defining sites.** Positions declared per locus (F10c supports located sites; F10 shows they must be declared for
   the locus, not for one haplotype).
5. **Copy identity.** When copies differ only by copy-specific SNPs, does it matter which copy is gained or lost (F11b)?
   Or is copy number per homologue sufficient?

## 5. Ready for a frozen broader evaluation?

**No.** The minimum missing component is a **rearrangement-aware whole-haplotype comparator**: one that aligns both strands,
allows non-collinear blocks, and handles repeats completely (an assembly-versus-assembly structural comparison). None exists
in the repository.

**The options need a decision, not more prototyping here:**
- **(a)** approve the narrower rearrangement-blind evaluator (A + C), with rearrangements reported as unsupported;
- **(b)** adopt an external assembly-comparison tool under a frozen protocol;
- **(c)** build the comparator.

No off-model parser was built. No merge, push, cleanup or promotion. Inputs, indexes, WIP and stash are preserved.
