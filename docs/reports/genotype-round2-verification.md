# `genotype` round 2 — verification pass and first experiment

Date: 2026-08-18. Status: item 1 of the agreed ordering is implemented (commit 9ac2318, a reporting
change with no behaviour change); nothing in the emission, scorer or band logic has been touched.
Sections 1 to 7 are the record of the verification pass and the first experiment. Section 8 is the
live status and the open questions.

Inputs to this report: an external defect review of `panvar genotype`, an accompanying redesign
specification, and a proof-of-concept script (`experiments/multiplicity_genotyper_poc.py`) proposing a
dosage-first model. Every checkable claim was read in the tree before being accepted or rejected.

---

## 1. Verification of the defect review

**17 claims confirmed, 2 needing correction.** Confirmed, with locations:

| claim | location |
|---|---|
| index omits `all_kmers` and `marker_clumps` | `genotype_index.cpp:84-96` ends after `anchor_slots`; both fields live on `ReadPanel` |
| indexed route is a reduced pipeline | `genotype_command.cpp:319-358` returns before coverage evidence, PanGenie and joint refinement |
| `--evidence coverage` and `--model-pangenie` accepted and ignored under `--index` | same |
| oracle ranks by syncmer-set Jaccard | `genotype_command.cpp:1478-1492`; `best_identity` aligns only the top 16 at `:1518` |
| haplotype pair and allele pair decoded separately, can disagree | `genotype.cpp:781-805` |
| PanGenie mode overwrites 6 fields of a panvar row | `genotype_command.cpp:1244-1252` |
| joint refinement then discards the PanGenie result | `genotype_command.cpp:1292` |
| `counts.edge` computed, never read | no reference in `genotype.cpp` |
| adjacency keys sorted, conflating A to B with B to A | `genotype_reads.cpp:31`, `genotype_markers.cpp:30` |
| clump window hard-coded at 350 while `--fragment-len` is an option | `genotype_markers.cpp:998` |
| transition is `recomb_rate / n_blocks` | `genotype.cpp:643` |
| pruned allele receives the all-background likelihood | `genotype.cpp:635` |
| no numeric range validation | only enums and `--kmer-size` are checked |
| help says "Only --audit is implemented so far" and "No genotypes are called yet" | `genotype_command.cpp:41,60` |
| no genotype test registered in CTest | 9 other modules have one |
| `synthetic_bench.sh` shares one wgsim seed across both homologues | `tests/synthetic_bench.sh:80` |
| both genotype documentation pages missing | `docs/modules/`, `docs/algorithms/` |

### Two corrections

**Pair-aware assignment already exists.** `genotype_command.cpp:1419-1426` aligns all four
called-to-truth combinations and takes the assignment with the smaller total edit distance. A length
rank is also already reported. The genuine gap is a certified nearest-pair oracle, not pair-awareness
itself.

**Bubble id order is asserted, not assumed.** `genotype_blocks.cpp:96-103` throws when ids are not
ascending. The residual risk is narrower than stated but real: `node_pos` silently returns 0 for a
non-numeric node name, so a graph from another constructor would compare every boundary as equal
rather than raising an error.

---

## 2. Recommendations already answered by prior measurements

These were measured previously in this project and are recorded so they are not re-run.

**Boundary fragments at tandem arrays are absent information, not a missing feature.** Sweeping read
length at the LPA KIV-2 block, at 350 bp the order evidence favours the wrong allele (86.2 percent
against 87.6 percent). The crossover is one repeat unit, between 2 and 5 kb against a 5,547 bp unit.
Fragment evidence remains worth building for blocks a fragment can span.

**Distance-aware transitions are not binding.** `--recomb-rate` swept over 100,000-fold gives
byte-identical output, and at KIV-2 the emission alone ranks the target pair 4270th of 104,653. The
reproducibility argument for physical distance stands; the power argument does not.

**Fragment-window marker clumping is already implemented** as `marker_clumps`, and is what took the
paralogous, segdup and near-twin fixtures to 36/36.

**A non-negative mixture over panel alleles is already falsified** because the alleles are collinear.
The recommendation to make the basis unit variants rather than whole alleles matches what this project
reached independently from the graph side.

---

## 3. Literature check

The redesign cites five methods. Their relevant properties were read rather than recalled.

**Ctyper.** Uniqueness is enforced at the family boundary: k-mers occurring outside the homologous
group are discarded, while k-mers shared between paralogs inside it are retained and allocated by the
matrix. That is the implementable form of not discarding multi-mapping markers. Ctyper also states it
does not provide confidence values for genotypes, and its depth normalization is an acknowledged
limitation. It reports a 6.7 percent F1 degradation under leave-one-out.

**danbing-tk.** Per-marker bias correction is worth batch-r-squared 0.531 to 0.820 for VNTR length,
the largest single effect reported. It estimates composition, not unit order.

**PanGenie.** The 16 and 32 unique-k-mer per-allele caps are a memory optimization of the
`UniqueKmers` structure, not a statistical choice. A fair comparison should report capped and uncapped
results so the marker rule is isolated from the memory budget.

**Paragraph.** Genotypes breakpoints independently and infers the variant genotype from them. Prior
work in this project measured that composing independent records describing one allele fails badly,
with records-per-carrier as the predictor. The local-realignment idea transfers; the
breakpoint-independent decomposition does not.

---

## 4. Experiment: auditing the proof-of-concept

Run first, deliberately, because the redesign's motivation leans on its result.

### 4.1 It reproduces exactly

Through the existing harness unmodified, `LOO=1 SEED=42 tests/genotype_sim.sh lpa 1 30 0.001`. SEED=42
is the default, and its pair at p=0 is the one the review used (path indices 42 and 55).

    15/20 blocks exact, bubbles 7/10, 5 unrepresentable          matches the review
    block 13 = bubble 7, block_class=array, 457 alleles, 8754 markers
    production called 53,91 at 257,814 bp                         matches the PoC README
    PoC identity 0.915611 / 0.917052, production 0.825155 / 0.824934   to six decimals

### 4.2 It is not circular

Data flow traced end to end. `--truth-bp` reaches only a printed column and the final error line. It
never enters the latent dosage, the pair scores, the graph length calibration, the band construction
or the band selection. The claim in its README holds.

### 4.3 A prediction of ours that was wrong

We predicted the rescue would break at lambda 11.70, our own implied value against the 11.50 supplied.
It does not. The correct band spans 11.50 through 11.70.

| lambda | selected band | error |
|---|---|---|
| 11.00 | 279081..280232 | +16,719 |
| 11.07 to 11.21 | 273513..274689 | +11,114 |
| 11.30 to 11.40 | 267975..269137 | +5,588 |
| **11.50 to 11.70** | **262442..263594** | **+43** |
| 11.80 to 11.94 | 256896..258048 | -5,571 |
| 12.00 | 251333..252492 | -11,103 |

Errors are exactly 0, 1, 2 and 3 repeat units of 5,547 bp. The value 11.50 sits at the bottom edge of
the good window, which is why the script's own minus-1-percent audit fails while plus-1-percent does
not.

### 4.4 The band is the mechanism, not the loss

| lambda | unconstrained rank 1 | constrained rank 1 |
|---|---|---|
| 11.50 | 178,312, error +5,588 | 180,393, error **+43** |
| 11.70 | **53,91, error +5,571** (production's exact call) | 180,393, error **+43** |

The robust count loss alone never selects the right pair at either lambda, and at 11.70 it independently
reproduces the production NB call. The novelty is the placement of the window, not the loss. Separation
is small throughout: unconstrained ranks 1 and 2 differ by 44.9 on a loss of 5,513, or 0.8 percent.

### 4.5 Within a fixed band, length and identity are decoupled

    rank 1  180,393  263,428  +43 bp  identity 0.9156
    rank 2   53,294  263,387   +2 bp  identity 0.8419

Two bp from truth and 7.4 identity points worse than the pair 43 bp away. An earlier claim in this
project that length and sequence conflict was withdrawn as inferred rather than checked. This is the
checked version, measured inside a band where total length is no longer the free variable.

### 4.6 The depth it rests on is a fallback, not a measurement

`estimate_depth` skips a block with fewer than `--min-anchors` anchors, defaulting to 20. Block 13 has
13 anchors. The array block therefore never estimates its own depth and inherits the region-wide value;
its reported `anchor_median` and `anchor_mad` are that fallback rather than a measurement. The
high-anchor blocks around it read 11.0743 (1145 anchors), 11.2097, 11.7253, 11.9130 and 11.9373 (1395
anchors), a spread of 3.8 percent across a locus where the band window is about 2 percent.

### 4.7 It does not survive independent read draws

Five read simulations of the same pair. SEED and SEED+466k select the same path indices at p=0 while
producing different wgsim seeds, so only the reads change and the harness needs no modification.

Production is perfectly stable and stably wrong: 53,91 at 257,814 bp on all five, one repeat unit
short, with lambda 11.5 from 13 anchors every time.

The PoC at that same lambda is right 3 times in 5.

| seed | mass_bp | selected band | error |
|---|---|---|---|
| 42 | 266,198 | 262442..263594 | **+43** |
| 1440 | 265,098 | 262442..263594 | **+43** |
| 1906 | 264,864 | 262442..263594 | **+2** |
| 508 | 268,328 | 267975..269137 | **+5,588** |
| 974 | 268,356 | 267975..269137 | **+5,588** |

**Superseded, see section 7.1.** The paragraph below concluded that depth was not the problem. That is
wrong: lambda was held constant, which shows mass noise is sufficient to move the band, but lambda 11.5
is itself a quantized fallback and is biased low. Both effects matter. The paragraph is kept as written
so the correction is legible.

The failure mode is not depth. Lambda was identical in every run. What moves is the latent marker mass:
the two seeds with the highest mass are exactly the two selecting one repeat unit too long, and the
correlation across all five is perfect. The observed cross-seed spread is 3,492 bp, or 0.63 repeat
units, against a band spacing of 5,547 bp. The production `mass_bp_sd` of 2,892 is larger than the
observed spread, so that uncertainty estimate is honest and conservative.

The scorer alone is right 1 time in 5; with the band, 3 times in 5.

**The proof-of-concept trades a stably wrong answer for an unstable one that is right 60 percent of the
time.** It refutes its own headline and vindicates its design recommendation, which is a good outcome
for an experiment and the reason it was run first. The script is being kept.

### 4.8 A plug-in band diagnostic centred on a biased estimator, and how it fails

**Partly superseded, see section 7.3.** As first written this section was titled "an honest posterior
does not rescue it" and concluded that a posterior cannot help. That conclusion was an artifact of the
width chosen: it used the cross-seed empirical standard deviation from five draws rather than the
model's own, larger `mass_bp_sd`. At the model's own width the posterior is never confidently wrong.
What survives is the narrower claim that a posterior centred on a biased estimator can be confidently
wrong, which is what the numbers below show. The word "honest" was also an overclaim; this is a plug-in
Gaussian band diagnostic.

The natural next step, and the remedy both the redesign and we proposed, is to carry an uncertainty on
the span instead of a point estimate, and decline when it straddles a band boundary. That was built and
measured on the same five dumps.

    cross-seed point estimates: 265472, 267625, 267608, 264357, 264110
      mean 265835, sd 1706 bp (0.308 repeat units)
      independent-Poisson sd(span) about 387 bp
      variance inflation 19.4-fold

| seed | span | P(correct band) | P(top band) | top correct | decision at 0.80 |
|---|---|---|---|---|---|
| 42 | 265,472 | 0.638 | 0.638 | yes | **DECLINE** |
| 508 | 267,625 | **0.033** | 0.966 | **no** | **CALL** |
| 974 | 267,608 | **0.034** | 0.965 | **no** | **CALL** |
| 1440 | 264,357 | 0.931 | 0.931 | yes | CALL |
| 1906 | 264,110 | 0.955 | 0.955 | yes | CALL |

The posterior declines a correct call and confidently accepts both wrong ones at 96.6 percent. Any
threshold that catches seeds 508 and 974 must exceed 0.966, which declines everything.

The reason is that the error is bias, not variance. Errors against truth are +2,087, +4,240, +4,223,
+972 and +725. All five are positive, with mean +2,449 bp, or 0.44 repeat units, against a noise sd of
1,706. The bias is 1.4 times the noise, and a posterior centred on a biased estimator cannot correct
for that however well its width is calibrated.

The 19.4-fold variance inflation is itself a result: it is correlation between markers sharing
fragments, not per-marker overdispersion, and a negative binomial with independent per-marker
dispersion cannot represent it.

### 4.9 Incidental

Genotype quality for the identical call on the same pair, with only read noise differing: 70.8, 49.9,
99.0, 55.7, 18.1. A 5.5-fold range on an unchanged answer.

---

## 5. What this changes

**Debiasing the span estimator is now the blocking item**, ahead of the solver it was meant to feed.
Candidate sources, all already in scope:

1. Per-marker mappability weighting, which danbing-tk measures as the largest single effect in their
   method. This moves from a refinement to a prerequisite.
2. The depth value is a region-wide fallback at the one block whose answer is a copy number. Dividing
   by a low lambda inflates mass, which is the correct sign for the observed bias.
3. The error background is hand-set to 1. Too low a background also inflates mass, again the correct
   sign.
4. The length calibration is fitted on panel alleles and applied to a sample whose true alleles were
   removed from the panel. That is a bias risk inherent to the leave-one-out condition being tested.

**A block that inherits its depth should say so.** Nothing in the current output distinguishes a
measured lambda from a fallback, at precisely the block where copy number is the answer.

**Acceptance criteria need rethinking** if the honest output is sometimes a set of adjacent copy-number
states rather than one.

---

## 6. Questions for review (round one, all answered)

These were the questions posed after section 4.7. They were answered in the round-two review; the
answers are folded into sections 7 and 8. Kept for the record. Live questions are in section 8.4.

### The questions as posed

1. The remedy proposed for the depth sensitivity was a depth posterior. Depth was not the problem here:
   lambda was stable and correct across all five draws, and the failure is a positive bias in the mass
   estimate that exceeds the noise. Does that change the solver design, and should the effort go into
   debiasing before any posterior is built?

2. A Poisson-independent variance under-predicts the observed spread by 19.4-fold. That is correlation
   between markers sharing fragments, which an independent per-marker dispersion parameter cannot
   represent. Does this make the fragment-level likelihood load-bearing for the tandem-array case as
   well, rather than only for boundary evidence? The existing effective-sample-size discount is an
   approximate version of the same correction, which suggests the right idea is being applied in the
   wrong place.

3. What is the acceptance criterion when the honest answer is a set? If the correct output on 2 of 5
   draws is two adjacent copy-number bands with no decision, accuracy is the wrong gate. Calibration of
   the band posterior, expected cost, or something else?

4. Is there prior art on calibrated posteriors over integer copy-number states for a locus of this
   kind? Ctyper explicitly does not report genotype confidence, and danbing-tk reports r-squared rather
   than per-sample uncertainty.

5. The planned order is measurement integrity first (certified oracle, index and direct parity,
   haplotype and allele consistency, harness redesign, input contracts), then an ideal-multiplicity
   oracle with no reads to establish whether the truth is representable in the node basis at all, and
   only then the solver. Given sections 4.7 and 4.8, would you reorder anything?

---

## 7. Round-two corrections

An external review of sections 4.7 and 4.8 raised four corrections. All four are accepted. Three were
re-derived here before acceptance; two of the reviewer's own measurements are recorded as theirs.

### 7.1 "The failure mode is not depth" was wrong

Holding lambda fixed shows that marker-mass sampling noise alone can move the band. It does not show
that lambda is correct, and lambda is not correct.

The region fallback is a plain lower median of raw integer anchor counts (`median_of` returns
`v[v.size()/2]`), so the fallback value is quantized to half-integers: it can report 11.5 or 12.0 while
the band decision responds to shifts of 0.1 to 0.2.

An independent check, using only the simulation parameters and no reference to truth:

    haploid depth                          15
    31-mer coverage factor (150-31+1)/150   0.8
    error-free probability 0.999^31         0.9695
    predicted marker lambda                11.6335

Sweeping lambda and counting how many of the five draws select the correct band:

| lambda | correct |
|---|---|
| 11.500 (the fallback) | 3 of 5 |
| 11.600 | **5 of 5** |
| 11.607 (fitted to truth, diagnostic only) | **5 of 5** |
| **11.633 (simulation theory, no truth used)** | **5 of 5** |
| 11.650 | **5 of 5** |
| 11.700 | 3 of 5 |

The theory value lands inside the all-correct window without using truth. The fallback lambda is
biased low by about 1.2 percent, and that is sufficient to move the array by one repeat unit.

The corrected statement is: fixed-lambda experiments prove mass sampling noise exists, but lambda 11.5
is itself a quantized fallback and is biased low. Both effects matter, and because the target is a
ratio their covariance matters too. A shared coverage fluctuation raises anchors and array markers
together and cancels in the ratio, which independent perturbation of numerator and denominator would
miss.

One qualification to add: the all-correct window is 11.60 to 11.65, about 0.43 percent wide. Debiasing
alone therefore does not make this a solved point estimate. It removes the systematic component; the
covariance-driven spread still has to be carried.

### 7.2 A reporting defect, not documentation debt

`BlockDepth::anchor_median` is documented in the header as "the raw anchor-count median, never
rewritten by a depth model", with `median` holding the model's fitted value. `write_read_audit` writes
`d.median` under the column header `anchor_median`, and `d.mad` under `anchor_mad`.

For a block below `--min-anchors` the first loop is skipped entirely, so `d.anchor_median` and `d.mad`
keep their defaults of 0, while `d.median` is later overwritten with the region fallback. Block 13
therefore prints `n_anchor=13 anchor_median=23 anchor_mad=0 lambda_hap=11.5`, which reads as a precise
local measurement and is an inherited quantized constant. The correct value is already computed and
stored; only the writer conflates them.

Required output instead: raw anchor median and MAD as distinct columns from the fitted value, plus an
explicit `depth_source` of LOCAL, SHRUNK, REGION_FALLBACK, QUANTILE, BASES or JOINT.

### 7.3 The band posterior in 4.8 used the wrong width

Section 4.8 used the cross-seed empirical standard deviation, 1,706 bp estimated from five draws, as
the posterior width. Using the model's own `mass_bp_sd` instead, carried to the calibrated-span scale:

| seed | sd used | P(top band) | top correct | P(correct band) | plug-in 95% set size | truth in set |
|---|---|---|---|---|---|---|
| 42 | 2,884 | 0.535 | yes | 0.535 | 2 | yes |
| 508 | 2,908 | 0.716 | no | 0.216 | 3 | yes |
| 974 | 2,908 | 0.715 | no | 0.218 | 3 | yes |
| 1440 | 2,872 | 0.685 | yes | 0.685 | 2 | yes |
| 1906 | 2,870 | 0.709 | yes | 0.709 | 3 | yes |

At the model's own conservative width the posterior never becomes confidently wrong, and every 95
percent credible set contains the true band. So the conclusion in 4.8 that a posterior cannot rescue
this was an artifact of the width chosen, and section 4.7 had already noted that `mass_bp_sd` was the
larger and more conservative estimate before 4.8 declined to use it.

The finding that survives is narrower and still useful: a posterior centred on a biased estimator can
be confidently wrong, which is what the 1,706 bp version demonstrates.

One number to record alongside: the call rate at P(best) at least 0.90 is **0 of 5**. Abstention
converts a point caller that is right 60 percent of the time into one that is honest and never calls.
That is the correct behaviour if the evidence does not identify the band, and it is why the deliverable
has to be a calibrated set rather than a point.

### 7.4 Two smaller corrections

The table in 4.7 shows production `mass_bp`, but band selection uses the proof-of-concept's
regression-derived target span. For seed 42 those differ across the decision boundary: production
`mass_bp` 266,198 against a target span of 265,472, with the boundary at 265,785. The rank ordering
argument survives, since the two highest values are the two failures under either quantity, but the
table should show both and name the one that drives selection.

"Honest posterior" in 4.8 is an overclaim. The width comes from five draws of one genotype, it is
evaluated on those same draws, it assumes a Gaussian estimator with a flat prior over supported
intervals, it conditions away the gaps between bands, and it models neither depth, background, marker
bias, calibration uncertainty nor an off-panel state. It is a plug-in Gaussian band diagnostic.

### 7.5 Reviewer measurements recorded but not re-derived here

Raising the proof-of-concept background from 1.0 to 1.5 reduces the mean bias only from 2,450 to
1,872 bp, and roughly 3.1 would be needed to remove it at lambda 11.5, which is less plausible than the
measured depth correction. Leave-one-allele-out cross-validation of the marker-mass to length
regression across all 457 alleles gives a residual standard deviation of 342 bp with no length trend,
far below the observed 2,450 bp bias. Together these remove background and panel-calibration
extrapolation as leading explanations, leaving the depth fallback as the primary one.

### 7.6 Revised ordering

1. Fix measurement provenance: correct the mislabelled depth columns, report `depth_source`, separate
   raw from shrunk from fallback.
2. Replace the quantized median fallback with a continuous estimator, for example
   `anchor_count_k ~ NB(2 * lambda * q_k + background_k, phi)`, initially a robust or trimmed mean.
3. Re-run the five seeds without truth tuning, comparing median fallback, pooled-anchor mean, robust
   mean, NB maximum likelihood, with the simulation-theory value as a diagnostic control only.
4. Expand to 50 to 100 independent read draws before any calibration claim. Five found the defect and
   are far too few to calibrate against.
5. Add fragment to marker incidence to the diagnostic dump, so depth and array mass can be bootstrapped
   from the same fragments and retain their covariance.
6. Only then build the copy-number posterior, integrating depth and mass jointly.
7. Continue with the certified nearest-pair oracle and the ideal-multiplicity oracle, which remain
   correctly ordered before the full solver.

### 7.7 An adopted implementation shortcut

For the certified nearest-pair oracle: with A candidate alleles and two truth haplotypes, only 2A exact
alignments are needed. Cache each allele's alignment against each truth haplotype, then score all
A(A+1)/2 pairs by the better of the two assignments, read from the cache. At this block that is 914
alignments rather than 104,653 pairs, so the oracle can be exact with no Jaccard preselection and no
shortlist at all. This removes the approximation the original review objected to rather than bounding
it.

---

## 8. Implementation status

Live status of the ordering agreed in section 7.6. Sections 1 to 7 are the record of how we got here;
this section is what is true now. Resolved questions are removed rather than accumulated.

### 8.1 Done

**Item 1, depth provenance.** Commit `9ac2318`, a reporting-only change verified byte-identical on lpa.
`DepthSource` and a shrinkage coefficient on `BlockDepth`; the audit separates the block's own raw
observations from the model's fitted value; a log line names the bubble blocks taking the region's
depth whole. `write_read_audit` also moved after joint refinement, which replaces every fitted depth,
so the audit no longer describes a state the emission never used.

It made visible that 9 of 12 lpa bubble blocks have no accepted local estimate, and that 61 anchors
across six blocks exist but are discarded by the `min_anchors` cliff.

**Item 2, continuous depth estimator.** `DepthEstimator` of Median, Mean or TrimmedMean behind
`--depth-estimator`, with all three reported side by side on every run whichever is selected. The
default remains `median`, so nothing has changed behaviour; the flag exists so the comparison below is
measured rather than assumed.

### 8.2 The quantization defect: diagnosis confirmed, without using truth

The previous revision predicted that a pooled mean would land near an anchor count of 23.27 and give
lambda about 11.63. Measured on lpa, 19,330 region anchors:

| estimator | anchor centre | lambda | correct band, 5 draws |
|---|---|---|---|
| median (current default) | 23.0000 | 11.5000 | 3 of 5 |
| trimmed mean, 10 percent | 23.1396 | 11.5698 | 3 of 5 |
| **arithmetic mean** | **23.2505** | **11.6252** | **5 of 5** |
| simulation theory | 23.2670 | 11.6335 | 5 of 5 |

The mean agrees with the independently derived theory value to 0.07 percent and uses no truth at any
point. **The trimmed mean does not work**, exactly as the review predicted before the measurement
existed: a trimmed mean estimates the trimmed centre, not the expected count, and for a skewed count
distribution the difference is enough to miss the window. The requirement is Fisher consistency for the
expected count, not robustness as such. Had this been adopted on intuition it would have looked like a
fix and been none.

### 8.3 What the estimator change buys, and a correction

`mass_bp` is the quantity to judge, not `called_bp`: at `block_class=array` the allele pair is the
closest content match and its length is not the copy number, while `mass_bp` is the documented answer
and depth is the denominator that produces it. Five read draws of one pair, truth 263,385 bp:

| seed | median lambda | median error | mean lambda | mean error | allele call moves |
|---|---|---|---|---|---|
| 42 | 11.5 | +2,813 | 11.6252 | -67 | no |
| 1440 | 11.5 | +1,713 | 11.6809 | -2,499 | yes |
| 1906 | 11.5 | +1,479 | 11.6311 | -1,609 | yes |
| 508 | 11.5 | +4,943 | 11.5672 | +3,382 | no |
| 974 | 11.5 | +4,971 | 11.5485 | +3,835 | no |

|  | bias | sd | RMSE | mean absolute |
|---|---|---|---|---|
| median | +3,184 | 1,695 | 3,526 | 3,184 |
| mean | +608 | 2,878 | 2,645 | 2,278 |

**Correction to a first reading of this experiment.** Seed 42 alone gives +2,813 to -67, and that was
quoted as a 42-fold improvement before the other four draws had run. Seed 42 is the favourable draw.
Across all five the honest figure is a 25 percent RMSE reduction: the mean removes 81 percent of the
bias and adds 70 percent to the variance.

Three consequences:

Lambda is no longer constant across draws, ranging 11.5485 to 11.6809 against a flat 11.5 under the
median. That spread was always present and the integer lattice concealed it. Removing the bias converts
it into variance that now has to be carried rather than suppressed, which is the argument for treating
depth and mass jointly rather than fixing depth and moving on.

The allele call changes on 2 of 5 draws, and blocks exact improves from 15/20 to 16/20 on those two,
with bubbles unchanged at 7/10. So this is not a reporting-only change and cannot be committed as one.

A plain mean is the right direction and not the destination. It is Fisher consistent, which is why the
bias goes, and it is also the least robust choice available, which is the likely source of the added
variance. A negative binomial or Huber estimator targeted at the expected count should keep the
correction without inheriting the tail sensitivity.

**The default therefore stays `median`** until the synthetic ladder and a second real locus agree. Five
draws on one locus is not a basis for changing every locus.

### 8.4 Defects found in our own item-2 change

Four, all raised in review and all fixed before commit. Recorded because three of them are the same
class of error this project keeps making.

1. **The flag did not reach the default model.** Joint refinement read `anchor_median`, which is
   deliberately pinned to the median so the audit column stays comparable across runs. A run with
   `--depth-estimator mean` would have changed the first pass and left the final joint value
   median-based. A separate `local_center` now follows the estimator and joint uses it. This is the
   project's own standing rule, broken in a new form.
2. **The indexed route did not receive the estimator**, recreating exactly the indexed-versus-direct
   drift that is already on the defect list.
3. **`region_weight` was a shrinkage coefficient, not a fraction of the fitted value.** Renamed.
4. **Zero was doing double duty for "not computed".** Zero is a legitimate anchor count; a
   `local_available` flag now drives an explicit NA in the audit.

### 8.5 In flight

Synthetic ladder A/B, median against mean, across clean, nested-del, abutting, paralogous, segdup,
near-twin, folded-consensus and per-design VNTR. The paralogous case is the one that matters most: with
one allele stripped of markers it can separate heterozygous from homozygous only by the absolute count
level, which is precisely what the depth estimator sets, so a regression should appear there first.

Generator flags for this were verified against the tool rather than guessed. An earlier draft used
`--nested`, `--segdup` and `--vntr`, none of which exist; it would have run the same default locus
eight times under eight labels and reported no difference anywhere.

### 8.6 Not started

The `min_anchors` cliff, which should become a continuous shrinkage over every non-empty block with
only n=0 as a pure fallback, and with the arbitrary 200 pseudo-anchors replaced by a precision
estimated from within-block and between-block variability. A negative binomial or Huber estimator
consistent for the expected count. Tests: commit `9ac2318` added none and there is still no registered
`genotype_stats.sh`; the needed assertions are a below-threshold block reporting NA with
REGION_FALLBACK, a shrunk block whose raw and fitted values differ by the expected coefficient, the
audit carrying the final joint depth rather than the first-pass value, direct and indexed runs honouring
the same estimator, and continuous behaviour across 19, 20 and 21 anchors. Items 3 through 7 of section
7.6 are untouched, as are both oracles.

The `--explain-pair` and `--deconvolve` diagnostic paths request the historical median explicitly; that
needs to be either deliberate and labelled, or changed.

### 8.7 Open questions

**Estimator.** Does a negative binomial or Huber estimator targeted at the expected count keep the bias
correction while recovering the variance the plain mean gives up? The five-draw variance increase is
the specific thing to beat.

**The cliff and the pseudo-count.** Is 200 defensible at all, or should the shrinkage precision be
estimated? A block with 13 anchors currently contributes nothing where it should contribute about 6
percent.

**Scale of the default decision.** The ladder plus one more real locus is what we plan to require
before flipping. Is that enough, given the call moves on 2 of 5 draws at a single block?

**Certified oracle, a correctness point not yet acted on.** Minimum summed edit distance over the two
assignments is not equivalent to maximum mean identity when the two alignment lengths differ. The
caching shortcut is still valid, but the pair assignment must optimise whichever identity definition is
actually reported.
