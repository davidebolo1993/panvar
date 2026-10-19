# M3-C: likelihood-guided shortlisting with incremental scoring (diagnostic)

**Question.** Does observation-likelihood ranking recover useful two-block proposals at a practical computational cost?

**Ancestry.** Branch `m3-c-likelihood-shortlist` from `226b4fb` (M3-B). There is no library change since baseline
`282e64e`, and M2 stays closed.

**Scope.** No caller integration, prior or default change, dosage repair, residual editing, novel bases or FASTA.

## 1. Cells and frozen rule

**Cells.** The same three as M3-B, with the same input hashes (`INPUTS.tsv`, copied unchanged) and incumbents (the
rank-1 representative row of `equivalence_set.tsv`):

| cell | role |
|---|---|
| `lzo_5x` | correct LZO control, **run first** |
| `mosaic_5x` | correct twin; its budget is reserved before the positive case |
| `mosaic_2x` | the positive case |

**Rule (truth-free).**
1. **Singles.** For every block b and homologue h, score **every** single-block alternative with the other homologue
   fixed.
2. **Ranking.** Rank by **fragment term + exposure term** (the context's `log_likelihood`, prior excluded), descending,
   with ties broken by allele index.
3. **Kept alleles.** Keep the top **4** non-incumbent alternatives.
4. **Catalogue.** The kept singles plus every same-homologue two-block combination across all block pairs, including
   combinations whose intermediates lose. Deduplication is by unordered sequence pair, through the context's
   byte-confirmed slot identity, with every equivalent move kept as provenance.
5. **Truth exclusion.** No truth, error coordinates, allele identities or M3-B observations enter generation.

**Bounds** (separate, checked before any delta):
- at most **2,048** single-block alternatives per cell;
- at most **2,048** catalogue candidates per cell.

On a breach the cell stops with nothing scored. The extra single-block work is reported apart from the joint work.

## 2. Incremental exact scoring (existing machinery only)

All scoring uses `PreparedScoringContext` (stage-5 reference path, the default):
- **Incumbents:** `add_haplotypes` with censuses, one full placement each.
- **Every single:** `add_repaired_haplotype` from the incumbent, a block replacement reusing every origin outside the
  affected window. No census is retained.
- **Kept singles (≤ 4 per block and homologue):** promoted by repeating the same call with `retain_census = true`. This
  fills the existing slot's census by delta and materialises it, as the repair search promotes a winner. Censuses are
  therefore retained selectively: at most 88, never one per trial haplotype.
- **Joints:** `add_repaired_haplotype` from the promoted single at the **earlier** block. The later block's interval is
  shifted by the first replacement's length change.
- **Totals:** `score_candidates` (combination from retained columns) plus the sequence-level pair prior at r = 0.01.

**Reported:**
- full placements beyond the two incumbents, as the fallback count (expected 0);
- origins reused and recomputed, DP executions and window bases, per phase;
- census count and bytes;
- context capacity;
- timings for context and incumbent, singles, promotion, joints, combination and prior, plus the unit's wall time and
  peak RSS.

No library code changes.

## 3. Gates (frozen tolerances)

1. **`gate`, before any cell.** Tiny fixtures: an 8 kb random haplotype with ~1,800 simulated fragments.
   - **Cases:** distant, adjacent and overlapping-window replacements, each at equal length and with length changes
     (+37/−23, −41/+29, +55/−30).
   - **Variation:** both edit orders, and the edited haplotype as either homologue.
   - **Requirements:**
     - |Δfragment| and |Δexposure| ≤ **1e-9** against cold placement of the complete joint sequence in a fresh context;
     - the joint slot equals the slot the context assigns to the cold joint bytes (byte identity);
     - both orders give one slot;
     - 0 full placements on the delta path.
2. **Incumbent reproduction.** Each cell's incumbent total must reproduce the saved total within **1e-6** before
   anything else, or the cell refuses.
3. **`compare.m3b`.** Every M3-C candidate whose unordered allele-path pair also appears in M3-B's cold reference scores
   (2× complete; 5× merged parts) must agree on fragment, exposure, prior and total within **1e-6**. LZO has no M3-B
   scores and is reported as such.

## 4. Selection, evaluation and interpretation

- **Selection.** `rank` (truth-free; M3-B's rank.py with its Amendment-2 parsing) freezes the following, with a digest,
  before any evaluation:
  - the winner, ties, the improvement and the incumbent's rank;
  - the best joint move with its intermediates, decomposed.
- **Evaluation.** `eval.<cell>` (truth) reports:
  - exact sequence match;
  - block projection;
  - unordered distance;
  - the phase diagnostic;
  - whether the truth was generated;
  - truth alleles' **likelihood ranks** (from the frozen singles table), inspected only afterwards.
- **No retuning.** 4 is never changed if the correct allele stays excluded.
- **A higher-scoring but less accurate prediction** is a model finding, not a correction.
- **Incomplete controls.** If the controls cannot complete within budget, that is reported as an engineering limitation,
  and no partial positive result substitutes for the gate.

## 5. Units and budget

**Units, in order:**

| unit | declared timeout |
|---|---|
| `gate` | 60 s |
| `lik.lzo_5x` | 180 s |
| `lik.mosaic_5x` | 180 s |
| `lik.mosaic_2x` | 180 s |
| `rank` | 30 s |
| `compare.m3b` | 30 s |
| `eval.<cell>`, ×3 | 60 s each |

**Limits:** 600 s cumulative, counted as spent time plus the declared timeout. At most 4 threads, 6 GiB, one process at a
time, at most 180 s per unit.
