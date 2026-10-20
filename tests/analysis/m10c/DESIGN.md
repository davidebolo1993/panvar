# M10-C: feasibility of sequence-based structural evaluation (frozen before any fixture or case was evaluated)

**Closes M10-B / B.1.** The long-node skeleton is rejected as a production structural alphabet. Its historical outputs are
kept. M10-A remains the computational checkpoint.

**Purpose.** An **offline evaluator**: compare predicted structural diplotypes with truth **without** requiring truth to
parse into the calling panel.
- **Not** a caller or confidence change.
- A pairwise, alignment-based match is **not** a transitive equivalence class, and is never used to pool posterior mass.
- **Scope:** no reads, genotyping or families.
- **Budget:** 900 s, including 30 s of disclosed development (two builds of the diagnostic; nothing evaluated).
- **Branch:** `m10-c-sequence-evaluation`. **Results:** `tests/results/m10c_sequence_evaluation_20261020/`.

## 1. Reused machinery

| need | existing implementation reused |
|---|---|
| whole-locus truth sequences | `project_path_onto_locus`, declared in the header as "EVALUATION SUPPORT, never inference", via `genotype_probe block-sequences --truth-gfa <full graph> --donor <id>`. It expresses the held-out truth path in the model's block frame with the same boundary handling and trims. **It does not require the truth to be spellable from model alleles** (truth block bytes may be novel). M1-B and M9-I used it. All blocks are requested. A block it cannot express (nullopt) makes that truth haplotype unavailable. The built model's fingerprint must equal the saved dump's |
| global alignment | `panvar::nw_edit_distance` (edlib NW + path; the `benchmark` round-trip-QV routine, including its contiguous-block split at a threshold). `m10c_align_eval` is a thin wrapper: one more walk over edlib's path records indel runs, terminal runs and located-site bases, and runs both directions |
| predictions | saved labellings (M1-B, M9-H, M10-A) spelled from the M10-B.1 model dumps, which are byte-identical to M10-B's |

**Matched endpoints.** The primary scope is the chain **interior**: block 1 through block n−2. The leading and trailing
flank blocks hold each path's own assembly extent, so they are excluded **identically** for truth and prediction (the truth
projection uses the same trims). The flank-inclusive scope is also reported, as an endpoint-compatibility measurement.

**Assignments.** Only the two whole-haplotype assignments, never a per-block choice.

## 2. Frozen fixtures (`m10c.py` `fixtures_def`)

| fixture | content | required verdict |
|---|---|---|
| F01 | identical DNA (node splits and boundaries are representation only) | CONCORDANT |
| F02 | 4 neutral SNPs, a 2 bp deletion and a 3 bp insertion | CONCORDANT |
| F03 | 55 bp deletion (many short nodes in a graph) | DISCORDANT |
| F03b | 56 bp deletion inside a (CA)n repeat | DISCORDANT |
| F04 | copy number 3 vs 2 | DISCORDANT |
| F05 | 3+1 vs 2+2 | DISCORDANT, with the dosage-balanced flag |
| F06 | order A,B vs B,A | DISCORDANT |
| F07 | 600 bp inversion | DISCORDANT |
| F07b | 60 bp inversion | DISCORDANT |
| F08 | shifted representation of one tandem expansion (identical copies, same DNA) | CONCORDANT |
| F08b | microsatellite expansion placed left vs right | CONCORDANT |
| F09 | 6 novel SNPs absent from any panel | CONCORDANT |
| F10 | 1 bp change at a declared (located) class site | DISCORDANT |
| F10b | the same change at an undeclared position | CONCORDANT |
| F11 | assignment tie: (S, S) vs (S+55, S) | DISCORDANT |
| F11b | one of three near-identical copies deleted (which copy is ambiguous) | DISCORDANT |
| F12 | truth extends 200 bp past the right anchor | UNSUPPORTED |
| F12b | truth lacks the first 20 bp | UNSUPPORTED |
| F13 | equal-length 55 bp divergent replacement | DEFINITION_REQUIRED (reported, not scored) |
| F14 | 60 bp near-tandem duplication in truth only | DISCORDANT |

**Copies** carry copy-specific SNPs. **50 bp** is used only as the existing event-size convention. Nothing is fitted to
development errors.

**Gate.** The aligner is deterministic on a repeat run, its edit distance is symmetric under reversal, and identical
sequences give 0 edits.

## 3. Measurements and the experimental verdict

**Measurements** (reported always, per homologue, under the chosen assignment, for both assignments' totals):
- edits, mismatches, inserted and deleted bases;
- E1: indel runs ≥ 50 bp, with sizes and positions;
- E2: the benchmark's non-match blocks ≥ 50 bp;
- terminal runs, located-site bases, net length difference;
- all of the above in both directions.

**Verdict:**
- UNSUPPORTED: truth unavailable, or a terminal run ≥ 50 bp;
- AMBIGUOUS: an assignment tie with different structure, or the two directions disagree on whether an E1 event exists;
- DISCORDANT: an E1 event, or a declared-site difference;
- CONCORDANT: otherwise.

**Flags:** E1_E2_DIFFER, PLACEMENT_SENSITIVE, E1_SIZES_DIRECTION_DEPENDENT, ASSIGNMENT_TIE, ENDPOINT_NOTE,
DOSAGE_BALANCED_ALLOCATION_DIFFERS, CLASS_SITE_DISCORDANT.

**Method failures** are results. A fixture whose required verdict is not met is reported as a failure of the method on that
distinction, not hidden.

## 4. Inspected saved cases (only if the fixtures support the method; development only; no accuracy claim)

| case | predictions | truth |
|---|---|---|
| C4 NA19650 | reported 3+1, saved 2+2 (M1-B) | `c4.gfa` donor NA19650 |
| C4 HG00096 control | M10-A `D.c4_5` | `c4.gfa` donor HG00096 |
| CYP2D6 NA19650 | 5× 2+2, 30× twin 3+2, saved 2+3 (M1-B) | `cyp2d6.gfa` donor NA19650 |
| ACOT L4 HG02083 | M9-H JP on / off | `acot.gfa` |
| ACOT L2 HG02074 | M9-H JP on / off | `acot.gfa` |

Provenance (truth GFA SHA-256, fingerprint match) and endpoint compatibility are recorded per case.

**Units (manifest P):** `build`, `gate.aligner`, `fixtures`, `truth.<5>`, `cases`. Admission: spent + declared timeout ≤
900 s.

## Amendment 1 (after `fixtures`; no truth extraction or case has run)

**What `fixtures` showed (recorded, unchanged).** Method A (unit-cost global alignment) gets the **edit distance** right,
but **cannot locate events**. Unit-cost gaps fragment freely:
- the F03 55 bp deletion has a longest indel run of 8;
- the F04 600 bp copy insertion, 24;
- inversions and order changes appear as scattered mismatches and indels.

Neither E1 nor the benchmark's E2 detects them. Required distinctions missed: F03, F04, F06, F07, F07b, F11, F14 (all read
as CONCORDANT) and F12b (a 20 bp endpoint offset, only noted).

F10 came out AMBIGUOUS. The cause is my fixture: it declared the site on one of two identical truth haplotypes, so the
assignment is arbitrary. The evaluator's ambiguity is the correct response. F10 stays recorded.

**Method B (added, existing implementation).** `panvar::minimap2_hits` (asm20, affine gaps, all regions), in both roles,
with the frozen event rule in `m10c.py` (`events_b`):
- reverse-strand hit → INVERSION;
- chain order → ORDER;
- target overlap → DUPLICATION;
- inter-hit gap → GAP;
- CIGAR I/D ≥ 50 → INDEL;
- uncovered ends → endpoint;
- hits overlapping by ≥ 50 query bases → ALTERNATIVE_PLACEMENT (flag).

The assignment and located-site bases still come from Method A. **Method A's results stay recorded as its failures.**

**Fixture added:** F10c, the class site declared for the locus on both truth haplotypes (required DISCORDANT).

**Condition for inspecting saved cases.** They run only if Method B meets **every** structural and concordance fixture:
F01, F02, F03, F03b, F04, F05, F06, F07, F07b, F08, F08b, F09, F10b, F10c, F11, F11b, F14.
- Excluded: F12/F12b and F13, because endpoint matching for real cases is **by construction** (the same chain boundaries
  for truth and prediction), and F13 needs a definition.
- If the condition fails, no case is inspected.

**Runs (manifest A1):** `build.a1`, `fixtures_b.a1`, `truth.<5>.a1`, `cases.a1`. The P units `truth.*` and `cases` are
superseded (recorded NOT RUN). Disclosed amendment development: 30 s.

## Amendment 2 (after `fixtures_b.a1`; no case inspected)

**What Method B showed (recorded).**
- **Met:** F03, F03b, F07 (600 bp inversion as a reverse hit), F11, F11b, F14, F10c and the concordance fixtures.
- **Failed:**
  - F04 and F05 (copy number; 3+1 vs 2+2): minimap2 returned **one** hit over the first 2,000 bp and dropped the rest of a
    sequence with repeated units. The verdict read that as an unmatched endpoint.
  - F06 (order): only one 601 bp block was aligned.
  - F07b: the 60 bp inversion was absorbed into the hit as mismatches.

  So a minimap2 primary or secondary hit set does not establish copy number or allocation (as the brief warned). Under the
  frozen condition, no case was inspected (`truth.*.a1` and `cases.a1` recorded NOT RUN).

**Method C (added, existing implementation).** `panvar::poa_msa` (abPOA: global, affine gaps) as a two-row alignment,
forward and on plain-reversed sequences.
- **Events:** gap runs ≥ 50 bp.
- **Verdict rule:** Method A's.
- **Assignment:** Method A's exact edit distance.
- **Known in advance:** a single-strand global alignment cannot represent an inversion.

**Condition for inspecting saved cases (declared now, before running C).**
- **FULL:** C meets every fixture in the Amendment 1 list.
- **INVERSION_BLIND:** C fails only F07/F07b. Cases are then inspected under the explicitly **narrower, named
  "inversion-blind"** evaluation, with inversions listed as unsupported.
- **NONE:** otherwise; no case is inspected.

**Runs (manifest A2):** `build.a2`, `fixtures_c.a2`, `truth.<5>.a2`, `cases.a2`. Disclosed development for Amendments 1–2:
35 s.
