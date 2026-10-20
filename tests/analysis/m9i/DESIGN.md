# M9-I: representation-aware compound confidence targets (read-only, saved evidence)

**Where M9-H stands:** closed as SUPPORTIVE bounded validation. JP1.1 stays default off, and no families are consumed
here.

**Not done:** no caller rerun, new reads, fitting, tuning, production-output change, filter, model merging, clustering,
default change, merge, push or cleanup.

**Results:** `tests/results/m9i_compound_targets_20261020/`.

## 1. Target contract (frozen)

**Diagnostic intervals:**
- **ACOT sv:5..ctx:6;**
- **CYP2D6 sv:3..sv:5**, including the intervening ctx:4.

These are declared diagnostic units, not claimed minimal, and not a claim that every degeneracy has been found.

**Controls:** ordinary single blocks, ACOT sv:13 and CYP2D6 sv:9.

**Model evidence (`inventory`):**
- **ACOT:** M9-D's inventory of the LOO and M9-C LZO models (2,045 cross-block equivalences, all within sv:5..ctx:6), and
  the same exhaustive search re-run here on M9-H's LOO, Z1 and Z2 models.
- **CYP2D6:** the sealed pilot inventory (sv:3..sv:5; 41–42 equivalences per model).

**The compound target T** is the **unordered pair of concatenated interval sequences**: the allele bytes of every block in
the interval, per homologue.

| property | behaviour |
|---|---|
| exchanging homologues | invariant |
| synonymous block labellings (different labels spelling the same bytes) | invariant |
| any genuine byte change | sensitive |
| a different within-interval arrangement | sensitive |

**It is an exact-sequence target:**
- it tolerates **no** SNP or indel;
- it is **not** structural equivalence;
- its probability is not "structural confidence".

## 2. Bounds (`m9i.py`, read-only)

**Groups:** the explored sequence groups of `mosaic_scores`.

**Members:** **every** recorded labelling of each group in `mosaic_candidates`, not only the representative. Each is
projected to T.

**Mass:** each group's `candidate_set_probability` is counted **once**.

| bound | definition |
|---|---|
| support_lower | mass of groups whose recorded projections all equal T |
| support_upper | mass of groups with at least one projection equal to T |
| error_lower | mass of groups with no projection equal to T, summed directly (never 1 − rounded support) |
| error_upper | mass of groups with some projection unequal to T, summed directly |

**Ambiguity.** A group with several distinct projections is reported as projection-ambiguous; its mass is never divided
or assigned.

**Scope.** The bounds cover explored candidates and recorded labellings only. They certify nothing about unsearched
alternatives.

**Predictions are kept** even when the annotation is ambiguous or unassessed.

**Reported labelling:**
- **v2 runs:** the representative labelling (`sequence_locus`), with the decision-set groups taken from
  `sequence_alternatives`;
- **legacy pilot runs:** the rank-1 equivalence-set group's first recorded labelling. Every member projects identically
  in the cases used, which is checked.

## 3. Gates (`gates`; synthetic, with an independent oracle that spells whole haplotypes and slices by offsets)

| gate | checks |
|---|---|
| G1 | two labellings spelling the same compound sequence are one projection |
| G2 | a real sequence difference is error mass |
| G3 | homologue exchange and allele-index permutation are invariant |
| G4 | mixed projections give an interval, and the mass is not split |
| G5 | duplicate synonymous labellings create no mass |
| G6 | an absent target gives support 0; a pool that wholly supports the target gives NO_EVALUATED_ALTERNATIVE |
| G7 | a single block reproduces the existing per-block error-mass calculation |
| — | normalisation (support + error = total, both pairings) and lower ≤ upper |

**Real-data control.** On v2 runs the single-block control must reproduce the saved `sequence_blocks`
`error_mass_lower/upper` within their 6-significant-digit printing.

## 4. Saved cases (explicit list)

**ACOT** (30 runs):
- M9-C: 8 cells;
- M9-G: 10 runs (off and on × LOO d1, d3, d6, LZO d1, d2);
- M9-H: 12 runs (off and on × L1–L4, Z1, Z2).

**CYP2D6** (4 runs):

| run | model | case |
|---|---|---|
| pilot `lzo.HG02769.30` | original | a sequence-identical label-ambiguity no-call |
| pilot `loo.HG02953.30` | loo.HG02953 | a sequence-identical label-ambiguity no-call |
| M1-B `cyp_5` | loo.NA19650, v2 | NA19650 LOO, 5× |
| M1-B `cyp_30` | loo.NA19650, v2 | NA19650 LOO, 30× |

**Truth bytes** (`truth.<…>`): the existing `block-sequences` truth projection for the interval blocks, per (model
curation, donor), 17 units.

**Compound distance.** The reported interval against the truth interval, as an exact edit distance (`genotype_probe
edit-distances`) minimised over homologue exchange.

**Reported side by side:**
- block-label exactness and distances (the existing evaluator);
- compound exactness and distance;
- projection ambiguity;
- compound support and error bounds;
- the existing block-level bounds.

**Specific check.** Every recorded labelling of every top-set group in M9-H L4 (both arms) must spell the claimed interval
sequence.

**Locus excess.** The existing "locus excess" is a block-aggregated evaluation, the sum of per-block distances, and is not
whole-locus sequence edit distance. It is not relabelled.

**No PASS threshold** is proposed from these inspected cases.

## 5. Budget

**Total: 600 s,** including 3 s of disclosed pre-freeze inspection.

**Admission:** spent + timeout + 120 s (the analysis). One process at a time.
