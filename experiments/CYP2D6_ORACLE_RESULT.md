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
