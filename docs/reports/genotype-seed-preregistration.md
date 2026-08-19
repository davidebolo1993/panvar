# Pre-registration: does the GC association reproduce across read draws?

Written and committed before any seed beyond the first is analysed. Four claims in this thread were
corrected because the conclusion was drawn from a check chosen after seeing the data. Fixing the
statistic, the model and the decision rule first removes that freedom.

## What is being tested

On one lpa leave-one-out sample, markers were found to show a weak positive association between the
GC of their k-mer and their observed count per copy. It survived dosage normalization, a block-specific
linear position adjustment, and sign agreement across most blocks, at r-squared about 0.001. That was
one read draw. This asks whether the association is a property of the panel and the simulator, or of a
draw.

## The model

One fit per seed, over markers the sample actually carries:

    log E[y_k] = alpha_b + beta * GC_k + f_b(u_k) + log(m_k)

| term | definition |
|---|---|
| `y_k` | the observed count of marker k, the raw count and not a ratio |
| `m_k` | truth multiplicity, this sample's own copy number from the spelled truth sequences; enters as a fixed offset so dosage is divided out inside the model rather than in the response |
| `alpha_b` | a free intercept per block |
| `GC_k` | GC fraction of the marker's own k-mer, the term of interest |
| `f_b(u_k)` | a natural cubic spline in normalized truth position, BLOCK-SPECIFIC |
| `u_k` | `truth_pos_mean`, the mean occurrence position within the truth sequence, normalized to [0,1] by that sequence's own length |

Pinned choices, so none of them can be selected afterwards:

- **Family**: Poisson with a log link, fitted by IRLS. Dispersion is not estimated and no marker-level
  standard error is used or reported, because markers covary through shared fragments -- measured at
  roughly 9-fold for anchors and 20-fold for markers -- and any marker-level interval would be invalid.
- **Spline**: natural cubic, 4 degrees of freedom, interior knots at the 25th, 50th and 75th percentiles
  of `u` within each block, boundary knots at 0 and 1. Block-specific, not pooled.
- **Inclusion**: `truth_mult > 0`, `truth_pos_mean` present, and the block has at least 200 such
  markers. Blocks below that are excluded entirely rather than fitted without a spline.
- **Weighting**: unweighted; each marker is one observation and multiplicity enters only through the
  offset. A multiplicity-weighted refit is reported as secondary.
- **Pooling**: one `beta` shared across blocks per seed. Per-block `beta_b` from separate fits is
  reported as secondary heterogeneity, not as the primary quantity.
- **Seeds**: 50, `SEED = 42 + 466k`, which holds the genotype fixed and changes only the reads.

## Evidence: the 50 seed-level values, and nothing below them

From the 50 estimates `beta_1 .. beta_50`, report exactly:

- the mean slope;
- the across-seed standard deviation;
- the 95 percent t interval for the mean, on 49 degrees of freedom;
- the number of seeds with a positive slope.

Excess spread is NOT to be attributed to any particular source. An earlier draft of this file proposed
comparing the across-seed spread against the mean within-seed OLS standard error, which is incoherent:
that standard error is invalid for the same clustering reason the marker-level statistics are, so it
cannot diagnose whether spread exceeds counting noise. Separating counting noise from other variation
requires a fragment-level bootstrap or cluster-robust uncertainty, and neither exists yet.

## Decision rule, and its limits

**Primary, statistical.** The association is called reproducible when the 95 percent t interval for the
mean slope excludes zero AND at least 45 of 50 seeds are positive.

Fewer than 45 positive seeds means **no reproducible positive effect was demonstrated in this wgsim
experiment**. It does not close the real-library hypothesis, and nothing measured here could. An
earlier draft said the efficiency line "closes" below 90 percent, which is both too absolute and out of
proportion to the difference between 44 and 45 of 50.

**Co-primary, practical.** A tiny slope reproducing in 50 of 50 seeds can be real and irrelevant. So
alongside the slope, report the quantity that decides whether it matters:

    predicted lambda bias at block 13 = beta_mean * (weighted GC delta) / 2

against the all-correct window half-width of 0.025, where the weighted GC delta is the
multiplicity-weighted GC of block 13's truth-carried markers minus the mean GC of the anchors. On the
first seed those were 0.4366 and 0.3991. An effect is practically relevant only if this predicted bias
is an appreciable fraction of 0.025.

**Secondary, reported without a threshold.** The same slope on anchors alone, where multiplicity is 1
by construction; the fraction of (block, seed) fits positive; and at block 13 the per-seed mean AND
median of `count_per_copy`, because on the first seed the mean sat 1.39 percent above the anchors while
the median sat 0.03 percent above. That gap says the excess is a tail rather than a shift, and only a
mean-based lambda would follow it.

## What this run cannot establish

These reads come from wgsim, which samples fragments approximately uniformly and models no library GC
selection. A slope that reproduces here is a reproducible property of the marker set and the simulator.
It is not evidence of GC bias in a real library, and no number of seeds makes it so. The real-library
check with an external single-copy depth control remains the one that decides whether any of this
generalizes.

At 50 seeds, a proportion near 0.9 carries a standard error of about 4 percentage points. That is
enough to separate the bands above and not enough to calibrate anything finer.
