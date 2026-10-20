# M9-I result: compound interval targets separate block-boundary disagreements from sequence disagreements; on saved evidence they remove up to 236 bp of apparent error that was pure labelling, with no case of compound projection ambiguity

**Run.** Branch `m9-i-compound-targets`, freeze `a4b2098`, with no amendments. Results are in
`tests/results/m9i_compound_targets_20261020/`.
- **Execution:** 234.7 of 600 s, including 3 s of disclosed pre-freeze inspection. No failures.
- **Read-only:** no caller rerun, read, fit or output change.
- **M9-H:** closed as SUPPORTIVE bounded validation. JP1.1 stays default off.

## 1. Contract

**Diagnostic units:**
- **ACOT sv:5..ctx:6;**
- **CYP2D6 sv:3..sv:5,** including ctx:4;
- **controls:** single blocks, ACOT sv:13 and CYP2D6 sv:9.

These are declared diagnostic units, not claimed minimal or exhaustive.

**Model evidence (`inventory`), each model having exactly one degenerate region:**

| locus | models | within-block duplicates | cross-block equivalences | region |
|---|---|---|---|---|
| ACOT | five (LOO, both M9-C LZO, both M9-H LZO; M9-H's re-run here) | 0 | 2,045 each | sv:5..ctx:6 |
| CYP2D6 | seven, sealed pilot | 0 | 41–42 each | sv:3..sv:5 |

**The compound target** is the unordered pair of concatenated interval sequences.

| property | behaviour |
|---|---|
| homologue exchange | invariant |
| synonymous block labellings | invariant |
| any byte change | sensitive |
| a different within-interval arrangement | sensitive |

**It is an exact-sequence target:**
- **no** tolerance of SNPs or indels;
- **not** structural equivalence;
- its probability is not "structural confidence".

**Pool-conditional bounds.** For every explored sequence group, every recorded labelling is projected, and each group's
mass is counted once.

| bound | definition |
|---|---|
| support_lower | groups whose projections all equal T |
| support_upper | groups with any projection equal to T |
| error bounds | the direct sums of the complements, never 1 − rounded support |

- **Ambiguity:** a group with several projections is ambiguous, and its mass is never split.
- **Scope:** the bounds cover explored candidates and recorded labellings only. They certify nothing about unsearched
  alternatives.
- **The prediction is kept** whatever the annotation says.

## 2. Gates

**Synthetic gates: 8/8 pass,** against an oracle that spells whole haplotypes and slices by offsets.

| gate | checks |
|---|---|
| G1 | synonymous labellings form one projection |
| G2 | a real difference is error |
| G3 | homologue exchange and permutation invariance |
| G4 | mixed projections give an interval, and the mass is not split |
| G5 | duplicate labellings create no mass |
| G6 | an absent target gives support 0; a wholly supporting pool gives NO_EVALUATED_ALTERNATIVE |
| G7 | a single block reproduces the existing per-block bounds |
| — | normalisation and lower ≤ upper |

**Real-data control.** The single-block computation reproduces the saved `sequence_blocks` `error_mass_lower/upper` in
**32/32 v2 runs** (maximum relative difference 4.4e-6, within the 6-digit printing). The 2 legacy pilot runs have no saved
bounds.

## 3. Saved cases (34 runs; `compound.tsv`, `controls.tsv`, `l4_top_set.tsv`)

**Where degeneracy changes the picture** (block view = existing per-block evaluator distance; compound = exact edit
distance of the interval, allowing homologue exchange):

| case | block labels: exact? (distance) | compound sequence: exact? (distance) | compound error bounds | labellings in the reported group → compound projections |
|---|---|---|---|---|
| **CYP M1-B cyp_5** (NA19650 LOO 5×) | sv:3 no (118), ctx:4 yes, sv:5 no (118) | **yes (0)** | 1.4e-3 / 1.4e-3 | 1 → 1 |
| **CYP M1-B cyp_30** (30×) | sv:3 no (118), sv:5 no (118) | **yes (0)** | 1.0e-7 | 1 → 1 |
| **CYP pilot lzo.HG02769.30** (legacy distance_tie no-call) | the two top-set labellings disagree at sv:3 and sv:5 | **yes (0)** | 7.5e-12 | **2 → 1** (2 block-label variants) |
| **CYP pilot loo.HG02953.30** (legacy no-call) | same | **yes (0)** | 9.3e-6 | **2 → 1** |
| **ACOT M9-H on L4** | sv:5 no (3), ctx:6 no (3) | **yes (0)** | 0.018 | 1 → 1 (both decision-set groups) |
| ACOT M9-H off L2 | sv:5 7, ctx:6 3 (block sum 10) | no (5) | 0.12 | 1 → 1 |
| ACOT M9-C LOO d3 | sv:5 9, ctx:6 4 (13) | no (6) | 2.1e-6 | 1 → 1 |

**Wherever the block view and the compound view disagree, the block view overstates.** Across the ACOT cases, the
compound distance is never larger than the per-block sum. Examples: 0 against 6 (L4 on); 5 against 10 (L2 off); 6 against
13 (d3).

**Controls.**
- **LZO:** all exact in both views.
- **Non-degenerate single blocks:** the two views coincide by construction.

**The M9-H L4 check, against truth bytes:**
- **on arm:** both recorded top-set groups (33557 and 84904, both labelled 102/103 | 92/33) spell the truth interval
  **exactly**;
- **off arm:** both (108/45 | 92/92) do not.
- So the claimed exact interval holds for every recorded top-set labelling.

**Projection ambiguity.**
- **ACOT:** none, because every group has one recorded labelling.
- **CYP pilot:** the two label-ambiguity no-calls each have one sequence group with 2 labellings and 2 block-label
  variants, which project to **one** compound sequence.
- **Overall:** no saved group had more than one compound projection. The mixed-projection path is therefore validated by
  the fixtures only.

**Confidence (pool-conditional, not calibrated).** The compound error mass is not evidence of correctness: the target is
still exact-sequence, and the pool is still conditional.
- **Inexact calls can still look confident:** M9-C d2 has compound error 5e-28 at 1 bp from truth.
- **Others are honest:** M9-C and M9-G off d6 sit at 0.50, M9-G on d1 at 0.20, M9-H on L1 at 0.33.

**Locus excess.** The existing measure sums per-block distances. It is block-aggregated evaluation, not a whole-locus
sequence edit distance, and it inherits the labelling overstatement shown above. It should not be relabelled.

## 4. Recommendation

**Can compound annotations be implemented safely?** **Yes, as an additive annotation layer.**
- **What:** for a declared interval, one row with the compound status, the decision-set compound projections, and the
  support and error bounds defined above.
- **Where:** beside, never replacing, `sequence_blocks`. Predictions are unchanged.
- **How:** the computation only needs data the v2 report already holds (groups, every explored labelling, probabilities,
  allele bytes).
- **Interval source:** the existing representation-interval machinery or a model-level degeneracy inventory. The
  inventory is cheap (under 4 s per model here) and exhaustive within its stated search.

**Still needed for a simple user-facing confidence decision:**
1. **An interval-declaration rule** used by the caller: from the model inventory, merged with the per-run
   representation intervals. It must be conservative and must say that it is not minimal.
2. **Evaluation and reporting in the compound unit wherever an interval exists,** so labelling disagreements are never
   counted as genotype errors. Block-level exactness and "locus excess" should be reported only as diagnostics there.
3. **Independent calibration,** on held-out data and in the compound unit, before any threshold. The inspected cases
   here are not a basis for a PASS rule, and the remaining overconfidence (M9-C d2) shows that pool-conditional mass
   alone is not a correctness probability.
4. **A separate, explicitly named target for structure, if wanted.** Exact compound sequence remains strict on SNPs and
   indels.

**How to avoid exposing representation differences as errors:**
- report the prediction's sequence (FASTA or extraction references, already v2 policy);
- give block labels only as provenance;
- describe confidence and errors over compound intervals wherever the model has degeneracy.

**Not done:** no implementation, threshold, filter or default change.
