# Pre-registration: does the GC association reproduce across read draws?

Four claims in this thread were corrected because the conclusion was drawn from a check chosen after
seeing the data. Fixing the statistic, the model and the decision rule first removes that freedom.

**Amendment, 2026-08-19, made before any confirmatory slope was fitted or inspected.** The first draft
claimed this was written "before any seed beyond the first is analysed". That was false. Seeds 42, 508,
974, 1440 and 1906 had already been analysed at length -- they are the five-seed tables in section 4.7
and 8.3 of the verification report -- and they are exactly `k = 0..4` of this schedule. Those five are
therefore **exploratory** and are excluded. The confirmatory set is `k = 5..54`, fifty seeds none of
which has been looked at. The exploratory five may be reported separately and carry no weight in the
decision rule.

A second amendment in the same revision corrects the practical co-primary, whose formula did not match
the model it was supposed to summarise. Details in the co-primary section.

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
- **Seeds**: `SEED = 42 + 466k` for `k = 5..54`, which holds the genotype fixed and changes only the
  reads. `k = 0..4` are exploratory and excluded, see the amendment above.

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
alongside the slope, report the quantity that decides whether it matters.

**Corrected formula.** The model has a log link with GC in the exponent, so the multiplicative effect
of a GC difference is `exp(beta * delta_GC)` and the predicted change in the per-copy rate is

    predicted lambda bias = lambda_anchor * (exp(beta_mean * delta_GC) - 1)

where `lambda_anchor` is the mean `count_per_copy` over anchors in the same seed, and `delta_GC` is the
multiplicity-weighted GC of block 13's truth-carried markers minus the mean GC of the anchors.

The first draft wrote `beta_mean * delta_GC / 2`. That is wrong twice over: it treats a log-link
coefficient as if it were additive on the count scale, and the factor of two was carried across from
the anchor relation `count = 2*lambda`, which has nothing to do with a regression coefficient. On the
exploratory seed 42, with beta 0.0301025, delta_GC 0.0373547 and lambda_anchor 11.63235, the two give
**0.000562 against 0.013088** -- a factor of 23, and the difference between 0.02 and 0.52 of the
half-window. The wrong formula would have reported the effect as negligible.

Reported as the **ratio** to the 0.025 half-width, with no binary practical verdict attached. Naming a
threshold for "appreciable" after seeing a 23-fold formula error would be choosing the bar to suit the
number.

**Note from the exploratory seeds, recorded rather than acted on.** Fitting `k = 0..4` showed the
pooled `beta` scattering around zero (+0.030, +0.012, -0.037, -0.024, +0.019) while the anchor-only
slope was consistently positive and an order of magnitude larger (+0.41, +0.38, +0.29, +0.30, +0.35).
The likely reason is mechanical: Poisson IRLS weights each observation by its mean, and an array marker
at multiplicity 20 carries roughly ten times the weight of a single-copy anchor, so a single shared
`beta` is dominated by high-multiplicity markers and is closer to an array-marker slope than to a
global one.

The primary is deliberately NOT changed on the strength of this. The anchor-only slope is already a
pre-registered secondary and will be reported whatever the pooled figure does, and rewriting the
primary after looking at data -- even exploratory data -- is the freedom this document exists to
remove. It is recorded here so the reading of the confirmatory result is fixed in advance: if the two
disagree in the same way, the anchor-only slope is the one that speaks to the depth denominator,
because anchors are what set lambda.

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
