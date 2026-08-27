# CYP2D6 background oracle: six arms, 24 donors

Scored against `CYP2D6_BACKGROUND_PREREGISTRATION.md`, unchanged since it was committed before any arm
ran. Simulated reads, leave-one-out, cyp2d6 blocks 3 and 5.

## Result

| arm | counts | truth in top class | median \|top\| | median rank |
|---|---|---:|---:|---:|
| A1 production, 10 markers | observed | **24/24** | **28** | 1 / 105 |
| A2 unfiltered, 36, no background | observed | 0/24 | 1 | 53 / 105 |
| A3 conditional, truth background | observed | 0/24 | 1 | 25 / 105 |
| A2n unfiltered | noiseless | 1/24 | 1 | 38 / 105 |
| **A3n conditional** | **noiseless** | **24/24** | **1** | **1 / 105** |
| A4 joint 3+5 | observed | 0/24 | 4 | 773 / 6,930 |
| **A4n joint 3+5** | **noiseless** | **24/24** | **2** | **1 / 6,930** |

A4n gets block 5 right in 23/24 and block 3 in 21/24, out of 6,930 combinations, without being handed
either block's answer.

Marker restoration verified per donor, not assumed: production 10 markers, unfiltered 36, block 3 23 --
identical across all 28 donors. All 16 markers that fail both filters reach scoring.

## What this establishes

**Production is locally UNDETERMINED at block 5, not actively wrong.** The truth pair is in the
emission's top class in 24 of 24 donors, but that class has a median of **28 members**. The emission
does not prefer a wrong pair; it cannot separate. Something downstream -- the chain, candidate
ordering, tie-breaking -- then selects one of the 28. The final call is still wrong, and that remains
a caller problem; but the locus of the failure is identifiability, not a mis-ranking.

**Relaxing both filters restores identifiability under ideal expected counts.** A3n and A4n recover
the truth 24/24 with a median top class of 1 and 2. The filter diagnosis is confirmed as far as it
goes.

## What this does NOT establish

**That these markers can safely be used on observed reads.** The filter diagnosis is necessary and
not sufficient. On observed counts every restored-marker arm fails: A2 at median rank 53, A3 at 25,
A4 at 773 of 6,930. A2 is *worse than production*, so restoring markers without accounting for their
counts is actively harmful.

## The discrepancy, stated narrowly

On the 16 shared markers, pooled across 24 donors:

| | counts |
|---|---:|
| observed | 25,954 |
| predicted from blocks 3 and 5 truth | 17,567 |
| **unexplained by the current block model** | **+8,387 (+48%)** |

Median +47% per donor, range +1% to +89%.

> **The source of this discrepancy is not known.** An earlier draft called it "read acquisition /
> count contamination". That is premature and is withdrawn. These are simulated reads generated
> directly from the truth haplotypes, so mapping is not implicated. Candidates still open: marker
> occurrences in other bubbles, backbone or flanks; canonical k-mer collisions; sequencing-error
> matches; depth or dispersion misspecification; sampling variance; or an error translating paths into
> block multiplicities.
>
> Two candidates ARE excluded by construction: slot-to-code is strictly 1:1 across all 36 markers, so
> no observed count is shared between slots; and per-allele multiplicity reaches 10, so repeated
> occurrences inside one allele are already represented.

**Pooled excess does not establish causation.** The extra mass could sit on markers that do not
separate the candidates. Which markers actually drive the wrong decision is not yet measured, and the
pooled figure cannot answer it.

## The pre-registered reading

"Pass noiseless, fail observed." The honest statement of it:

**CYP2D6's discarded markers contain sufficient information in the ideal model, but the current block
decomposition cannot explain their observed counts.**

## Consequence for the plan

The route is not closed, but it has an intermediate step that was not in the plan: before shared
markers can be restored, the caller must account for every place in the locus that can generate their
counts. Building a joint block-level shared-marker factor now would be building on a model that
explains half the signal.

## Exclusions

24 of 28 donors analysed. Three -- HG00171, NA19036, NA19317 -- have **no representable truth at block
5** under leave-one-out and cannot be scored on recovering it; they are separated here rather than
counted as failures, and they are owed a certified-floor analysis of their own. HG00096's first run
was made with `-q` and carried no model line; it is recovered from its non-quiet run.


---

# Check 1: whole-locus multiplicity. The block model was missing real marker sources.

`--dump-truth-marker-counts` scans both truth haplotypes' complete walks and counts every panel
marker with the panel's own syncmer code, so the two cannot drift. Background is then measured from
sequence rather than reconstructed from block annotations.

| background for block 5 | median count error | truth at rank 1 | median rank |
|---|---:|---:|---:|
| blocks 3 + 5 only | **+47%** | 0 / 25 | 22 |
| **whole-locus minus block 5** | **+3.8%** | **16 / 25** | **1** |

The surplus is structural: on the shared markers each truth haplotype carries the marker once more
outside blocks 3 and 5, a consistent **+2 per marker**.

Counting the pre-registered way -- top **equivalence class**, not exact pair:

| | |
|---|---|
| rank 1 | 16 / 25 |
| **rank <= 2** | **21 / 25** |
| rank <= 5 | 21 / 25 |

Five of the nine misses sit within **1.5 log units**, two of them essentially tied at 0.07 and 0.08.
Only four miss by a real margin.

## Residual MAGNITUDE predicts failure. Sign does not.

| stratum | n | rank 1 | median \|err\| |
|---|---:|---:|---:|
| negative residual | 5 | 4/5 | 8.7% |
| positive, < 10% | 11 | 9/11 | 2.6% |
| positive, 10-20% | 4 | 3/4 | 16.4% |
| **positive, >= 20%** | **5** | **0/5** | 23.5% |

> **CORRECTION.** An earlier draft said "sign matters -- the negative cases don't break the genotype
> the way the positive ones do", from a single donor at -19.2%. Withdrawn. Across five negative
> donors the pattern does not hold: NA19347 succeeds at -19.2% while HG02554 fails at -12.4%. What is
> supported is **magnitude**: every donor above +20% fails, and below it 16 of 20 succeed.

> **A positive residual does not prove an additional sequence source.** Every exact occurrence in the
> complete truth haplotypes is already counted, so the remainder could equally be sequencing-error
> matches, depth underestimation, positional coverage (lambda x multiplicity assumes every occurrence
> gets identical coverage, and occurrences near a boundary have fewer possible fragment starts),
> canonicalisation collisions, or duplicated counting. An earlier draft called it "an unmodelled
> source still contributing counts". That is one hypothesis of several and is not established.

# Check 2: the deficit is concentrated, and entirely on shared markers

For the truth pair against the pair the score chose, per marker delta log-likelihood, ranked by
absolute value, over the nine failing donors:

| | |
|---|---|
| markers carrying 50% of the deficit | median **4** |
| markers carrying 90% | median **9** |
| single worst marker's share | median 17% |
| **deficit on markers shared with block 3** | **100%** |

Not diffuse miscalibration: a handful of shared markers decide it. That is the tractable case, and it
means a locus-wide shared-marker factor is aimed at the right objects.

> **Known defect in this table.** The filter-fate breakdown came out 56% `both` / 44% unattributed,
> because ledger slot indices come from a different panel build than the scored run. The fate column
> is unreliable; the shared-vs-block-5-only split is computed from the data itself and is sound.
> Keying on marker code rather than slot would fix it.

# Precise conclusion

**Whole-locus multiplicity explains most of the CYP2D6 failure and validates the architecture of a
locus-wide shared-marker factor. The remaining observed-count residual must be attributed before that
factor can be implemented safely.**

# Still owed, before any production change

1. Exact simulated-read accounting: trace every counted marker occurrence to its source read and
   coordinate, and partition into expected-position occurrence, sequencing-error-created match,
   canonical reverse-complement collision, count reused across slots, otherwise unexplained. Include
   marker-specific read opportunity, since `lambda x multiplicity` assumes uniform coverage per
   occurrence.
2. Multi-seed repeat with haplotypes fixed, to separate persistent marker-specific residuals
   (systematic model error) from residuals that change sign across seeds (sampling).
3. Certified-floor analysis for the three donors with no representable truth: measure
   `error(selected) - error(certified nearest)` and whether the top class contains a certified-optimal
   pair, rather than asking whether an unavailable pair ranks first.
4. The fate-attribution fix above.

# The production construction this points at

One unique marker -> one likelihood factor `y_j ~ NB(lambda * SUM multiplicity over every block state
carrying j, + mu)` -> connected to every block in which that marker occurs. Each observed marker
scored **once**, never once per block. The marker-to-block incidence graph decomposes into connected
components: exact enumeration for small ones, beam or message passing for larger. At cyp2d6 block 5
the component is just blocks 3 and 5.


---

# Closing the analyses: two of my explanations were wrong, and the third is supported

## Retraction 1: the read-opportunity computation was circular

`opportunity()` counted the **realised** wgsim fragments covering each occurrence, taken from the read
names. That predicts observed counts using the observed reads, so the reported collapse of the
residual to about -3% was near-tautological. **Withdrawn.** Positional read opportunity is not
established as the explanation, and the flanking-context experiment designed to test it is not the
priority it looked like.

A production-compatible opportunity term would have to be geometry-derived -- expected fragment
starts from position, sequence length, read length and insert distribution -- never from realised
reads.

## Retraction 2: per-donor lambda bias, refuted by multi-seed

The per-marker ratio `observed / (lambda * multiplicity)` looked like a tight per-donor constant
(median IQR 0.11) varying between donors (0.87 to 1.26), which suggested a systematic depth
misestimate. Re-simulating the same haplotypes under three further read seeds refutes it:

| donor | implied lambda across four seeds | sd |
|---|---|---:|
| HG03452 | 14.00, 9.00, 14.00, 14.00 | 2.17 |
| NA20509 | 13.94, 12.77, 12.38, 9.00 | 1.84 |
| HG00358 | 11.50, 12.00, 11.38, 13.33 | 0.78 |

The implied lambda swings by more within one donor across seeds than it varies between donors. It is
not a stable property of the donor, so it is not a depth bias.

## What the residual actually is: four independent observations, not 23

At block 5 the markers occupy almost no independent evidence:

| | |
|---|---|
| markers with occurrences | 23 |
| total occurrences over both haplotypes | 68 |
| **distinct 350 bp fragment windows** | **4** |
| occurrences per window | median 18, max 36 |
| per-marker ratio sd | 0.243 |
| Poisson sd for a single count of that size | 0.209 |

The per-marker spread matches the noise of a **single** count. The markers are not averaging down --
they sit inside four fragments and move together. Effective sample size is about 4, not 23, so a
+-25% swing in the pooled ratio between read realisations is exactly what should happen.

That also explains why residual magnitude predicted failure (0/5 above +20%) without any of it being
a fixable count-model defect: an unlucky realisation of four fragment windows inflates the residual
AND breaks the call, because the same four windows carry the decision.

**This is the clump/effective-sample-size problem, at a locus where it was not previously located.**
The emission already discounts by `rho`; what this shows is how severe the discount should be here --
the evidence behind block 5 is four observations, and treating it as 23 overstates certainty by
roughly sqrt(23/4), about 2.4x.

# Certified floor: the three excluded donors were misclassified

They were dropped because the truth ALLELE INDEX is absent from the reduced panel. Measured against
the certified nearest panel pair instead:

| donor | E* certified | E called | excess | best identity | stratum |
|---|---:|---:|---:|---:|---|
| HG00171 | 2 | 60 | 58 | 0.9977 | near-representable |
| NA19036 | 1 | 177 | 176 | 0.9987 | near-representable |
| NA19317 | 2 | 119 | 117 | 0.9969 | near-representable |

**All three are near-representable, not panel-limited.** The panel reconstructs them to within 1-2
edits; the caller misses by 58-176. Excluding them understated the failure rate, and the exclusion
criterion itself was wrong -- exactly what the release plan's three strata were written to prevent.
Representability must be defined from certified sequence error, never from whether a label exists.

# Denominator

25 donors, not 24. The earlier figure predates the `params.log` fallback that recovered HG00096's
model line. With the three near-representable donors added the scored set should be 28.

# Where CYP2D6 stands

**Established:** the block decomposition misses real marker occurrences, and whole-locus multiplicity
fixes it -- median count error +47% to +3.8%, recovery 0/25 to 16/25 exact and 21/25 within a
two-member class. The deficit that remains is concentrated (median 4 markers carry half) and sits
entirely on markers shared with block 3.

**Not established:** any count-model defect beyond that. Two candidate mechanisms were tested and both
failed. The residual is consistent with sampling over four independent fragment windows.

**Consequence for the design:** a locus-wide shared-marker factor remains the right architecture, and
it does NOT need a positional coverage term. What it does need is an honest effective sample size --
scoring 23 markers that occupy four fragments as four observations, not 23.
