# Phase 0: where the emission's preference actually comes from

lpa block 13, leave-one-out, 15 donors with a contrast to measure (the other 5 have the certified pair
AS the emission optimum, and are exactly the 5 with zero excess). Noiseless counts, full block cap.
Table: `results/genotype_kiv2_phase0.tsv`. Script: `cohort/ph0_analyse.py`.

## Why this experiment and not the factorial

The present/absent factorial decomposed `C(noiseless) - C(observed)` -- the CHANGE the injection
causes. At this block truth-absent markers carry ~500 counts across ~6,700 markers, almost all already
zero, so rewriting them was a near no-op (~0.1% of the block's 435,000 counts). Its "absent = 1.5%"
therefore said nothing about whether absent markers decide the ranking. **This decomposes the contrast
itself**, at fixed counts.

## Method, and the check that makes it usable

For every marker: `t` = truth multiplicity, `c` = certified pair's, `r` = the emission optimum's.
Per-marker contribution is `nb_h(o, lambda*c+mu) - nb_h(o, lambda*r+mu)` at the injected count `o`,
mirroring `genotype.cpp` including its Poisson branch when the fitted phi degenerates to 0 (which it
does on 4 of 15 donors -- modelling those as negative binomial would score a different likelihood
than the one under test).

**Harness check.** The per-marker sum times the block's clump discount must reproduce the binary's own
probe delta. Implied rho is 0.05456-0.05488 across all 15 donors -- one block constant. The offline
model is measuring the production emission, not an approximation of it.

> **Two different clump counts, and they are not interchangeable.** The rho above implies ~479 units.
> That is production's `marker_clumps`, which bins (allele, position) pairs over every allele. The
> grouping used by the Phase 1 scores is `node_first_pos / fragment_len`, one position per marker
> slot, and it yields **381-383** clumps on these donors. Both are panel-derived and independent of
> the candidate and the truth, so either satisfies the contract a production score must meet -- but
> they count different things and an earlier draft of this file quoted 479 for both.

## Result: support mismatch dominates

Counts below are **donor-marker observations pooled over 15 donors**, not unique marker sequences. The
same marker appears once per donor and a single sequence difference can supply many markers.

Gross negative -15,964.4; gross positive +1,735.5; **net -14,228.9**. Both are reported because a
share of the gross negative alone would hide that two classes actively favour the certified pair.

| class | obs | sum dLL | % gross neg | % net |
|---|---:|---:|---:|---:|
| `t>0 c=0 r>0` truth supported, competitor covers it, certified does not | **263** | **-11,872** | **74.4%** | **83.4%** |
| `t=0 c>0 r=0` certified predicts content truth lacks | **295** | **-3,589** | 22.5% | 25.2% |
| `t>0 both>0` pure multiplicity | 10,871 | -291 | 1.8% | 2.0% |
| `t=0 both>0` both predict absent content, differing | 19 | -213 | 1.3% | 1.5% |
| `t>0 c>0 r=0` certified covers truth, competitor does not | 6 | **+213** | - | -1.5% |
| `t=0 c=0 r>0` competitor predicts content truth lacks | 133 | **+1,523** | - | -10.7% |
| `c=r` uninformative | 119,678 | 0 | - | - |

**558 observations of presence mismatch outvote 10,871 of dosage evidence.** The multiplicity term
contributes 1.8% of the gross loss and is close to neutral. The certified pair is penalised in BOTH
presence directions -- for lacking sequence the truth has, and for carrying sequence it does not --
and is helped in both mirror directions, which is why any bound must be symmetric.

## And the loss is concentrated in a few clumps

Adjacent syncmers share reads, so N disagreeing markers are not N independent observations. Per donor,
over the informative markers:

| | range across 15 donors |
|---|---|
| clumps carrying the negative gap | 1 - 45 (median ~15) |
| top-1 clump share | 15.7% - 62.7% |
| **top-5 clump share** | **53.6% - 99.2%** |

HG00321's entire dominant class is **5 markers in 1 clump** worth -393. HG00097: 12 markers, 2 clumps,
top-5 = 99.2%. `rho` scales the total afterwards, but only after each correlated marker has already
voted.

> **What concentration does and does not show.** Markers sharing a clump sit within one fragment, so
> they are certainly correlated. It is CONSISTENT with a few dozen underlying sequence differences
> counted many times each, and that is the natural reading -- but it has not been demonstrated. A
> clump could equally hold several genuinely distinct differences. The clump-to-sequence mapping is
> what settles it, and until then this file claims concentration in the evidence, not a count of
> sequence differences.

## Heterogeneity: one mechanism, one secondary direction

13 of 15 donors have `t>0 c=0 r>0` dominant. HG00350 and HG00128 are dominated by the opposite
direction (`t=0 c>0 r=0`). HG00320 has pure multiplicity dominant but its total loss is -31.9, the
smallest in the cohort. This is not adaptive-scoring territory: one bound, applied symmetrically.

## What this does and does not say

**Says:** the interface between off-panel marker support and the current likelihood geometry is
misaligned with sequence-nearest projection. A small, highly concentrated set of presence mismatches
outweighs a large body of dosage evidence. How much sequence those mismatches actually represent is
not established here -- calling them "small differences" would assume the mapping's answer.

**Does not say:** anything about whether the retained markers are well chosen. An earlier draft said
"the markers are fine". That is withdrawn -- it was never demonstrated. This experiment localises the
defect to the interface between marker support and likelihood geometry; some markers may still be
poorly selected, or may redundantly encode one sequence difference many times over. Mapping the
dominant clumps back to sequence spans is what would settle it, and has not been done.

**Clump membership must stay usable at genotyping time.** `node_first_pos / fragment_len` is
panel-derived and so satisfies that contract. A grouping taken from held-out truth coordinates would
be diagnostic only, and must not leak into a production score.

**Does not license `--all-kmers --marker-outlier`.** That combination's held-out experiment failed
(Gate B: one pair fixed, one broken 4x, six untouched, 9-point call-rate cost) and that result stands.
Bounding returns only as a new, pre-declared, clump-level model.

## Branch taken

Support mismatch dominates consistently AND is concentrated in few clumps, which selects:
**a per-clump bounded likelihood, or a block-wide off-panel residual.** That is Phase 1.
