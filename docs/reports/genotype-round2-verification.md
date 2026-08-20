# `genotype` round 2 — verification pass and first experiment

Date: 2026-08-19. Status: seven commits landed, in depth estimation, reporting, the marker dump and the
module's first registered test. Nothing in the emission, the scorer or the band logic has been touched
and no default has changed. Verified defects V1 (index omissions) and V8 (hard-coded clump window) are
closed. Sections 1 to 7 are the record of the verification pass and the first experiment; section 8 is
the live detail; section 9 is the summary and the recommended order. The seed run is pre-registered in
`genotype-seed-preregistration.md` and has not been started.

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

### 8.5 Synthetic ladder: no regression, and a weak test

Median against mean, 4 pairs, seed 11, 30x. Every case identical in both arms.

| case | leave-zero-out | leave-one-out |
|---|---|---|
| clean, nested-del, paralogous, segdup, folded, per-design-vntr | 36/36 blocks, 16/16 bubbles | 36/36, 16/16 |
| abutting | 32/32, 16/16 | 32/32, 16/16 |
| near-twin | 36/36, 16/16 | 14/16, 14/16 (panel ceiling) |

The paralogous case is the one with the most at stake, because with one allele stripped of markers it
can separate heterozygous from homozygous only by the absolute count level, which is exactly what the
depth estimator sets. It does not move.

**Stated as a limitation rather than a result:** nearly every case sits at 36/36 under both arms, so
the ladder is saturated. It can detect a regression and cannot detect a gain, and passing it is
necessary but weak evidence. A version at lower depth, where depth estimation is noisier and the
lattice is relatively coarser, would have discriminating power that this one does not.

The generator flags were verified against the tool rather than guessed. An earlier draft of this
experiment used `--nested`, `--segdup` and `--vntr`, none of which exist; it would have run the same
default locus eight times under eight labels and reported no difference anywhere.

### 8.6 Robustness: one hypothesis removed, the question not settled

An earlier revision of this section claimed robustness was not needed. **That was overstated and is
withdrawn.** Matching three summary statistics does not establish that 19,330 counts follow a
homogeneous Poisson distribution; marker-specific efficiency, GC, repeated markers and mixtures can all
preserve mean, median and trimmed mean simultaneously.

The check also used the wrong trimming definition. It compared against a Poisson expectation computed
by retaining complete integer states between the 10th and 90th percentiles, which holds 85.4 percent of
the mass rather than 80 percent. The implemented estimator drops exactly 10 percent of observations
from each end. Simulating that estimator directly:

| statistic | Poisson(23.2505) target | observed |
|---|---|---|
| median | 23 | 23.0000 |
| 10 percent trimmed mean, as implemented | 23.1503 | 23.1396 |
| mean | 23.2505 | 23.250491 |

The corrected gap on the trimmed mean is 0.011 rather than the 0.062 previously reported, so the data
is closer to Poisson than the flawed check suggested. The error made the data look less Poisson than it
is, which happens to be the conservative direction, but the method was wrong either way.

What this does support, narrowly: the median's downward bias is fully accounted for by the discreteness
and skew of the count distribution, so **a gross contaminating tail is not what drives the median-mean
gap**. That removes one hypothesis. It does not establish the distribution.

What would settle it, none of it yet done: variance-to-mean ratio and fitted negative binomial
dispersion, observed against Poisson tail frequencies, zero frequency, quantile residuals, the same
summaries stratified by block, GC, uniqueness and marker clump, and held-out predictive deviance for
Poisson against negative binomial. All of this needs per-anchor counts exported, which they currently
are not.

**A separate correction that changes what "build a robust estimator" would even mean.** A negative
binomial with one common mean has essentially the same point estimate as the arithmetic mean; the
negative binomial changes the variance model, not the location estimate. So negative binomial is not
the robustness lever. Robustness would require marker-specific efficiency q_k, calibrated weighting, or
a Fisher-corrected robust score. The earlier conclusion "do not build the negative binomial estimator"
was right, but for the wrong reason: it is not that robustness is unnecessary, it is that a negative
binomial would not have provided it.

### 8.7 The variance objection, reframed with two corrections

Section 8.3 recorded that the mean increases the spread of `mass_bp`, which read as a cost. Two
corrections to how that was stated, then what it means.

**It is a 70 percent increase in standard deviation, not in variance.** Standard deviation goes 1,695
to 2,878; variance goes up about 188 percent, a factor of 2.88. The earlier wording conflated the two.

**Calling it "genuine between-draw depth variation" was wrong.** All five simulations share the same
nominal depth parameter, so nothing about the underlying depth changed. The variability is
cluster-correlated sampling of overlapping markers in the realized coverage. It must still be
propagated, but it is sampling variability, not a moving parameter.

    region lambda across 5 draws:   sd 0.0532
    independent-Poisson prediction: sd 0.0173
    ratio 3.07, variance inflation 9.4-fold, 95 percent CI 3.4 to 77.6-fold

**The 9.4-fold figure is not a calibrated factor.** It comes from five draws, and its interval spans an
order of magnitude in each direction. It is evidence that the anchors covary, which matters, and it is
not a number a posterior can carry as though it were known. An earlier revision said a depth posterior
must carry the 9.4-fold inflation "exactly"; that is withdrawn.

The direction of the conclusion survives, with the wording corrected once more: the added spread is
cluster-correlated estimator sampling error under a fixed lambda. It is not independent Poisson noise,
but it is still uncertainty in the estimator rather than something outside it. What holds is that
neither the depth nor the mass uncertainty can be derived from independent per-observation counting.
The magnitude needs far more draws.

### 8.8 What this does not do

The estimator change improves the depth denominator and moves two allele calls on lpa. Bubble accuracy
is unchanged at 7/10, and cyp2d6 is identical in both arms. This is a calibration repair. It is not
evidence that the model chooses the nearest available haplotype, and it does nothing for mosaics. The
certified nearest-pair oracle and the ideal-multiplicity oracle remain the decisive next experiments.

A reporting defect of exactly the class this pass exists to remove was found in the pass itself: the
region summary logged its selected anchor centre as though it were the lambda in force, which is false
under the Bases, Quantile and Joint models where the fitted depth comes from elsewhere. It now prints
`selected_anchor_center` and does not name it lambda.

### 8.9 The min_anchors cliff, removed

`estimate_depth` skipped any block below `--min-anchors` before computing a local estimate, so the
shrinkage pass had nothing to shrink and those blocks took the region value whole. That contradicted
the code's own stated rationale, which is that a hard threshold is the wrong shape and shrinkage
replaces it. At lpa it discarded 61 anchors across six bubble blocks.

Every block with any anchors now produces a local estimate and is shrunk by its own weight; only a
block with none is a pure fallback. `--min-anchors` now flags thin local evidence rather than gating
it.

    depth provenance: 4 REGION_FALLBACK, 21 SHRUNK        (was 10 and 15)

    block 11  bubble 6   n=17  local=27  fitted=23.3134  shrink=0.922
    block 13  bubble 7   n=13  local=23  fitted=23.0000  shrink=0.939
    block 15  bubble 8   n=1   local=15  fitted=22.9602  shrink=0.995
    block 17  bubble 9   n=13  local=19  fitted=22.7559  shrink=0.939
    block 19  bubble 10  n=12  local=28  fitted=23.2830  shrink=0.943
    block 21  bubble 11  n=5   local=21  fitted=22.9512  shrink=0.976

The local values are genuinely varied, 15 to 28, and the shrinkage does its job: block 15's single
anchor reading 15 is pulled to 22.96 rather than trusted. `mass_bp` moves by 1 bp, the truth check is
unchanged at 15/20 and 7/10, and graded accuracy is identical to six decimals. A consistency fix with
essentially no behaviour change.

### 8.10 The distributional checks, run

`--dump-anchors` writes one row per (block, marker) with the count, the marker k-mer's GC and the
confinement audit. On lpa that is 19,330 rows, and it makes section 8.6's checks answerable.

    variance-to-mean ratio   1.0075        (1.0 under Poisson)
    NB dispersion alpha      0.000322      (0 means Poisson adequate)
    zero frequency           0.000155 observed against ~0 Poisson, i.e. 3 true zeros
    beyond +/-2 sd           1.05-fold
    beyond +/-3 sd           1.43-fold
    randomized quantile residuals vs uniform: KS D 0.0087, p 0.104

The pooled distribution looks Poisson-like. **The Kolmogorov-Smirnov p-value is not a valid acceptance
test here and should not be read as one**: it assumes independent observations while these anchors
covary through shared fragments, which this same report measures at roughly 9-fold; the Poisson mean
was estimated from the data being tested; and the randomized residuals add a further source of
randomness. Read the whole block as descriptive evidence that outlier robustness is not the immediate
problem in this clean simulation, not as a formal acceptance of homogeneous Poisson sampling.

The three zeros are all in block 0, a flank, which is itself evidence of spatial heterogeneity that a
homogeneous model does not predict.

The stratified views show the structure the summaries hid, which is what the review predicted:

| GC quintile of the marker k-mer | n | mean count |
|---|---|---|
| 0.000 to 0.290 | 2,666 | 22.714 |
| 0.290 to 0.355 | 3,351 | 23.267 |
| 0.355 to 0.419 | 4,379 | 23.216 |
| 0.419 to 0.484 | 3,590 | 23.298 |
| 0.484 to 0.839 | 5,344 | 23.504 |

There is a positive GC-associated trend, and an earlier revision of this section called it monotonic
and treated it as an established effect. **Both were overstated.** The quintiles are not monotonic --
23.267 falls to 23.216 in the middle -- and re-binning on equal-count quintiles gives a different shape
again, peaking in the fourth rather than the fifth. The association is real but weak:

    pooled Pearson r(count, gc)   0.0399     r-squared 0.0016
    within-block centred r        0.0410
    excluding both flanks r       0.0387

So it is not an artifact of differing coverage between blocks or of the flanks, and it explains about
0.16 percent of the variance. Per-block slopes are mixed: of eleven blocks with at least 200 anchors,
two are negative and one substantially so at -6.759.

The arithmetic linking it to lambda is correct as far as it goes -- the extreme-bin difference is about
3.4 percent of the mean, which is 0.40 in lambda against an all-correct window 0.43 percent wide, so
roughly eight-fold wider. But an extreme-bin difference from one seed with binning-sensitive quintiles
is not a stable estimate of an effect size, and that comparison should not be quoted as though it were.

**And a mechanism qualification that matters more than the statistics.** These reads come from wgsim,
which samples fragments approximately uniformly and does not model library GC selection. There is
therefore no GC-bias mechanism in this data for the trend to be. GC is more likely acting as a proxy
for position, sequence context, syncmer ascertainment or a realized local coverage fluctuation. The
right phrase is GC-associated marker efficiency, not a GC effect, and certainly not a GC-dependent one.

The comparison with danbing-tk is directionally consistent and should be stated no more strongly than
that: this independently motivates the same class of correction. It does not reproduce their result,
because panvar has not shown that correcting the association improves copy-number or genotype accuracy.

**The honest conclusion is that marker-specific efficiency is now the highest-priority hypothesis to
test, not that it is the dominant remaining bias source.**

By block, the variance-to-mean ratio runs 0.826 to 1.418, with both flanks near 1.4 and interior
backbones below 1.16. Flanks are the locus edges where coverage falls away, so the region estimate
probably should not weight them equally; they currently contribute 1,841 of 19,330 anchors.

**Two defects in the dump itself.** It covers anchors only, so it cannot compare anchor GC against the
array block's marker GC, which is the decisive test for whether the trend biases the array at all. And
its confinement columns are useless as written: `expected` is 0 for all 19,330 rows because
`dbg_expected` is populated from the informative-marker expectation map while anchors are filtered
against a separate anchor expectation, and `actual` is uniformly 464, which is just the panel size
after exclusions and is the definition of an anchor rather than a fact about any particular one.

Extending the dump to informative markers is not sufficient on its own either. Without per-allele
multiplicity it would compare array markers carrying twenty copies against single-copy anchors and
confound efficiency with dosage. The extended dump needs, at minimum: anchor against informative,
block and within-block position, clump membership, per-allele multiplicity and carrier count, the
simulated truth multiplicity where truth is supplied, and the observed count normalized by that
multiplicity.

### 8.11 The dump extended, a bug in it found, and the GC question answered

`--dump-markers` replaces `--dump-anchors` (which remains as an alias) and covers every marker, not
only the invariant ones, with `marker_class`, per-allele multiplicity and carrier count, the marker's
offset and clump, and -- when `--truth-haplotypes` is supplied -- this sample's own copy number and the
count divided by it. That last pair is what separates efficiency from dosage. The anchor `expected`
field is fixed: it read 0 on every row because anchors are absent from `by_block`, and now carries the
anchor's own bound, which equals `actual` when the contract holds.

**A bug in the first version of this dump, found by diagnosing rather than reporting.** On lpa,
informative markers appeared to read 26 percent higher per copy than anchors, with 2.5 times the
spread. Stratifying by multiplicity showed it at once:

    truth_mult == 1   n=2364   mean count_per_copy 20.829   (about 2*lambda)
    truth_mult == 2   n=4922   mean count_per_copy 11.739   (about lambda)

Where one truth haplotype's allele was unrepresentable the loop skipped it silently, so `truth_mult`
was half the sample's actual dosage and `count_per_copy` doubled. It now requires both haplotypes or
reports NA. After the fix, anchors read 11.684 and informative markers 11.670, agreeing to 0.1 percent.
On the synthetic fixture the same quantity is exactly lambda for every marker of both classes, which is
what the regression test now asserts.

**The GC association, with dosage and position controlled.** The falsifiers were written down before
the numbers: dosage confounding, position, and sign disagreement across blocks.

    r(gc, raw count)                             +0.0349
    r(gc, count_per_copy)                        +0.0318
    r(gc, truth_mult)                            +0.0189
    r(gc, count_per_copy) less position in block +0.0469
    per-block sign agreement                     9 of 11 positive

It survives all three, so it is not dosage confounding, not position alone, and not a sign accident.
It is also very small: r-squared about 0.001.

**Whether it bites at the array. Superseded, see 8.12** -- the figures below were computed over every
candidate panel marker in block 13, because the held-out truth alleles were unrepresentable and the
block therefore carried no truth dosage at all. The target of the whole analysis was excluded from it.

    all anchors        mean GC 0.3991
    block 13 (KIV-2)   mean GC 0.4301     delta +0.0310
    implied lambda bias at array  +0.0261 against a window half-width of 0.025

That was an extrapolation from the anchor slope and panel-wide marker composition, and the claim that
it accounted for roughly a quarter of the span bias was not supported by it.

### 8.12 Dosage from the truth SEQUENCES, and the array measured directly

Truth multiplicity was computed from panel allele indices, which do not exist for a held-out array
allele. On lpa that left all 8,754 informative markers at block 13 -- and every marker at blocks 6, 11,
18 and 20 -- with `truth_mult` NA. The conservative both-haplotypes rule was right and its consequence
was that the block under study contributed nothing to the analysis about it.

Dosage now comes from `ts1` and `ts2`, the spelled truth sequences, which are retained whether or not
the allele is representable. Markers are counted straight out of the sequence against
`panel.node_codes`. Block 13 now reports dosage for all 8,754 rows, of which 2,007 are carried, at a
mean multiplicity of 20.34.

Two related contract fixes: an anchor's dosage was hard-coded to 2, which is wrong wherever a truth
haplotype takes a bypass -- on lpa 7 of 19,330 anchors have only one traversing haplotype -- and a
marker with dosage 0 reported `count_per_copy` 0, which put background observations into the efficiency
population. It is NA now.

Note that `--depth-calibration` skips blocks where either truth allele index is negative, so that
diagnostic has the same blind spot and has never covered the array either.

**The array, measured rather than extrapolated:**

| | value |
|---|---|
| anchors, mean count_per_copy | 11.628 |
| block 13 carried markers, mean count_per_copy | 11.790 |
| block 13 carried markers, MEDIAN count_per_copy | 11.632 |
| multiplicity-weighted GC delta against anchors | +0.0375 |
| GC-predicted lambda bias | +0.0316 |
| lambda bias implied by the mean per-copy excess | +0.0808 |

The mean per-copy rate at the array is 1.39 percent above the anchors, which implies a lambda bias of
+0.081 against a window half-width of 0.025, and GC predicts about 39 percent of it.

**But the median excess is 0.03 percent, not 1.39.** The array's typical marker reads exactly what an
anchor does; the mean is pulled by a right tail. So there is no uniform per-copy inflation to correct,
and whether the tail biases lambda depends on which estimator consumes it -- a mean follows the tail, a
median does not. That interacts directly with section 8.2, which recommended the mean on the grounds
that anchors are clean Poisson with nothing to be robust against. Anchors are; array markers, on this
evidence, are not. A single estimator may not be right for both populations.

This is one seed and the tail is not characterized. It is the first direct measurement of the quantity
that matters and it is not yet a conclusion.

### 8.13 Contract work cleared before the seed run

**The clump window is no longer hard-coded** (verified defect V8). `--fragment-len` reached only the
emission while marker clumping used a literal 350, and the clump count IS the effective sample size the
emission divides by, so a run with a different library silently kept the default. `MarkerOptions` now
carries it and `ReadPanel` records what it was computed at. Verified reaching the model rather than
assumed: quadrupling the fragment length takes the distinct clump count from 3 to 1 on the fixture.

Worth recording how that fix nearly failed. The assignment was first placed at the declaration of
`MarkerOptions`, which is BEFORE the argument parse loop, so it copied the default and the flag would
have gone nowhere -- the identical failure this module has now hit four times. It is after parsing, and
`--fragment-len` is range-checked.

**The index serializes what inference depends on** (verified defect V1, the highest-severity item in
the original review). `marker_clumps`, `node_first_pos`, `all_kmers` and `fragment_len` join the format
at version 5. Absent `marker_clumps`, an indexed run fell through to `span / fragment_len` for its
effective sample size and produced a different GQ from the same reads and the same calls, which is
exactly the symptom that was measured. Direct and indexed genotype rows now agree byte for byte
including GQ, and the test asserts that rather than only lambda.

**Truth-derived positions.** `node_first_pos` is one occurrence from one allele, which for a repeat
marker carried twenty times is one of twenty places it sits. `truth_pos_mean` and `truth_pos_span` are
computed from the truth sequence itself, so every occurrence is seen and the covariate is the sample's
own.

**The position control in section 8.10 was overstated by a small amount.** Rows in blocks too small to
residualize were left unadjusted and still included, mixing adjusted and unadjusted data. Dropping them
gives about +0.044 rather than +0.0469. The adjustment also removes only a block-specific linear trend,
so it does not exclude nonlinear spatial coverage.

The seed-level statistic is pre-registered in `SEED_STATISTIC.md` before any seed beyond the first is
analysed: a block-adjusted, spline-position-adjusted slope of `count_per_copy` on GC, with the decision
rule and the failure reading fixed in advance, and reproducibility judged from the distribution of the
per-seed slope rather than from marker-level p-values.

### 8.14 Contract corrections found in the previous round's own work

Four defects in the fixes themselves, all raised in review.

**The GQ parity assertion did not compare GQ.** It cut columns 1 to 12 while GQ is column 13, so it
compared `hap_posterior` under a message claiming GQ was included. The whole file is compared now, and
it is byte-identical between routes.

**An indexed run could still mix fragment lengths.** The index records the length its clumps were built
at, but the emission was assembled from the runtime value, so an index built at 1400 and run at the
default took its clumps from one and every other term from the other. `--fragment-len` is now inherited
from the index when unspecified and refused when it contradicts, both asserted.

**The 350-versus-1400 clump check was run by hand**, not registered. It is an assertion now.

**A double bypass was treated as unknown dosage.** When neither truth haplotype traverses a block, the
dosage is a known zero, and calling it unknown discards the one case where a marker's absence is the
informative observation.

**Truth positions pooled raw coordinates from two sequences of different lengths.** At an array those
differ by tens of kilobases, so a pooled raw offset correlates with copy number rather than with
location. Occurrences are normalized within their own sequence to [0,1] before pooling.

### 8.15 The seed experiment, held twice before starting

Three further defects, all found before any confirmatory data was fitted.

**The practical co-primary formula did not match the model.** It read `beta * delta_GC / 2`, which
treats a log-link coefficient as additive on the count scale and carries a factor of two across from
the anchor relation `count = 2*lambda`, which has nothing to do with a regression coefficient. The
correct form is `lambda_anchor * (exp(beta * delta_GC) - 1)`. On the exploratory seed those give
**0.000562 against 0.013088**, a factor of 23, and the difference between 0.02 and 0.52 of the
half-window. The wrong formula would have reported the effect as negligible.

**The first five seeds were not confirmatory.** Seeds 42, 508, 974, 1440 and 1906 had already been
analysed at length in sections 4.7 and 8.3, and they are exactly `k = 0..4` of the schedule, so the
claim that the plan predated any analysed seed was false. They are exploratory now and excluded; the
confirmatory set is `k = 5..54`.

**`--fragment-len 0` reopened the index route.** Zero means "disable clumping", and inheritance was
guarded on `> 0`, so an index built at 0 and run without the flag reverted to 350: GQ 33.03 against the
direct route's 99, files differing. Inheritance is unconditional now and both non-default cases, 1400
and 0, are asserted as whole-file parity. Until this, V1 and V8 were not actually closed.

Also: the truth contract is tri-state rather than overloading an empty string, since a haplotype that
spells nothing is a known dosage of zero only when an allele was found for it, and unknown otherwise.

**Two implementation defects in the analysis itself**, found by running it on the exploratory seeds
before touching the confirmatory ones. The IRLS convergence tolerance sat at 1e-10, below the numerical
noise floor of the least-squares solve, so every fit was flagged non-convergent while being stable to
eight decimals from iteration 10. And the pooled slope scatters around zero across the exploratory five
while the anchor-only slope is consistently positive and an order of magnitude larger, which is
mechanical: Poisson weights each observation by its mean, so a marker at multiplicity 20 carries about
ten times the weight of a single-copy anchor and a single shared slope is dominated by array markers.
That is recorded in the pre-registration and the primary is deliberately not changed on the strength of
it.

### 8.16 CONFIRMATORY RESULT: the GC association does not reproduce

50 seeds, `k = 5..54`, analysed exactly as pre-registered. The exploratory five are excluded.

    PRIMARY   mean beta -0.00641   95% t [-0.01450, +0.00167]   positive in 16 of 50
              -> no reproducible positive effect demonstrated in this wgsim experiment

    CO-PRIMARY  predicted lambda bias -0.00278, ratio 0.11 of the 0.025 half-width

    SECONDARY, with intervals
      weighted refit         -0.01680   95% CI [-0.02880, -0.00479]   19 of 50
      anchor-only            +0.00282   95% CI [-0.00850, +0.01413]   27 of 50
      per-block fits         356 of 750 positive, 47.5 percent, a coin flip
      block 13 mean excess   +0.203 percent  95% CI [-0.359, +0.765]  25 of 50
      block 13 median excess +0.158 percent  95% CI [-0.244, +0.560]  25 of 50
      mean minus median      +0.045 percent  95% CI [-0.261, +0.351]  21 of 50

The pre-registered rule required the interval to exclude zero and at least 45 of 50 seeds positive.
Neither holds, and the slope is positive in under a third of draws.

**Correction to a first statement of this result.** It said every estimator agrees at zero. That is
false: the multiplicity-weighted refit is **significantly negative**, with an interval excluding zero.
It does not rescue the positive-GC hypothesis -- the sign is opposite -- and the primary is untouched,
but the summary was wrong and the weighted figure it quoted, -0.00442, came from a defective
implementation described below. The accurate statement is that the pooled, anchor-only and per-block
estimators are consistent with zero while the weighted one is negative.

All 50 seeds converge at full rank.

**Two earlier findings of this report are withdrawn as single-draw artifacts.**

Section 8.10's GC-associated marker efficiency does not reproduce. It was one read draw, and across 50
the slope is centred slightly below zero.

Section 8.12's array per-copy excess does not reproduce. On seed 42 the array read 1.39 percent above
the anchors per copy; across 50 seeds it reads **0.203 percent with a 95 percent interval of -0.359 to
+0.767**, positive in 25 draws of 50, which is exactly a coin flip. The median excess behaves the same
way. The mean-minus-median gap that seemed to indicate a right tail is +0.045 percent with an interval
spanning zero, so the tail does not reproduce either, and the consequence drawn from it -- that a
mean-based lambda would follow an array tail a median-based one would not -- has no support and is
withdrawn.

**Two bugs in the analysis.** The first was found before the result was reported; the second after,
in review, which is worth recording as the less good outcome.

The multiplicity-weighted refit scaled the design and the response by `sqrt(m)`. That is not a weighted
Poisson GLM -- it fits a different model. Observation weights belong in the IRLS weight as `mu * w`
with X, y and the offset untouched. Corrected, the estimate moves from -0.00442 to **-0.01680** with an
interval excluding zero, which is what makes the "every estimator agrees at zero" sentence above wrong.

Separately, the IRLS formed the normal equations `X'WX`, squaring the condition number. That spurious
ill-conditioning was enough to stall one seed of fifty at the iteration cap. Solving the weighted least
squares directly on `sqrt(W) X` converges it in five iterations and moves its coefficient by 3e-11, so
no reported number depended on it -- but the earlier "k=46 failed to converge" note was an artifact of
the solver, not of the data.

The anchor-only secondary first read +0.327, positive in 50 of 50 seeds, which contradicted every per-block fit
including the array's own -0.017. Pooling with block intercepts cannot exceed its constituents, so the
estimate was wrong rather than interesting: 336 anchor rows of 19,330 belonged to blocks below the
200-anchor threshold and therefore had no intercept and no spline, so their entire count level was
absorbed by the shared GC coefficient. Restricting to modelled blocks, as the primary already did,
moves it from +0.4605 to +0.0367 on that seed and to +0.00282 across the 50. **1.7 percent of the rows
were driving the coefficient by a factor of twelve.**

That also disposes of the note added to the pre-registration from the exploratory seeds, which
attributed the anchor-versus-pooled gap to Poisson weighting by the mean. The mechanism was wrong; the
gap was the orphaned rows.

**What survives.** The depth-estimator work rests on evidence independent of any of this: agreement
with simulation theory to 0.07 percent, the ladder null across eight cases, the cyp2d6 null, and the
lattice argument that a median of integers cannot express a value between half-integers. None of it
depended on GC.

**What is closed, stated at its actual scope.** Marker-specific efficiency is **not supported as an
explanation for the array residual in this one fixture: lpa, one genotype pair, 30x, one error rate,
and a simulator with no GC-selection mechanism to detect**. That is narrower than "closed for
simulated data" and much narrower than closed: one locus, one genotype pair, one depth, one error rate,
and a simulator with no GC-selection mechanism to detect. A real library, a different depth, or a
different locus could all reopen it, and the real-library check with an external single-copy depth
control remains untouched by this experiment.

The fifty per-seed values are committed as `experiments/genotype_seed_slopes.tsv` rather than existing
only as console output.

### 8.17 The depth ladder: the fixture was trivial, and the estimator is neutral

Run against a criterion fixed in the script before any result existed -- flip, hold, or revert -- with
noise defined in advance as a one-block difference on one case at one depth.

**First finding, and it is about the harness rather than the estimator.** At 30x, 10x and 5x, with
exact twins retained, every case scored 36/36 blocks and 16/16 bubbles under both estimators. Read
counts confirmed depth reached the simulator (5,035 / 1,677 / 838 reads), so the ladder genuinely
scores perfectly at 2.5x per haplotype. The reason is defect V15: the generator emits each design twice
as an EXACT twin, so holding out a sample leaves an identical haplotype in the panel and leave-one-out
has a perfectly representable answer that barely needs coverage to find.

**The ladder was never saturated by generous depth. It was saturated because the task is trivial**, and
no reduction in depth can make a trivial task discriminate. Every ladder figure in this project's
history was measured in that regime, including the 36/36 results quoted as validation for the
paralogous, segdup and folded-consensus fixes. Those remain valid as evidence that the fixes break
nothing. They were never evidence that the model performs well off-panel, and the earlier observation
that routing "never loses and never gains" reads differently once it is clear the ladder could express
neither.

**With `--twin-divergence 40`, leave-one-out becomes an off-panel test and behaves like one.**
Leave-zero-out stays 36/36 everywhere, correctly, since the answer is still in the panel.

| case | 30x | 10x | 5x |
|---|---|---|---|
| clean | 14/16 | 14/16 | 11/16 |
| paralogous | 13/16 | 11/16 | 11/16 |
| segdup | 15/16 | 14/16 | 12 to 13/16 |
| folded | 12/16 | 11 to 12/16 | 11/16 |
| per-design-vntr | 16/16 | 14/16 | 13/16 |

**Second finding: the estimator is neutral even where it should matter most.** 28 of 30 cells are
identical. The two that differ are one block each and in opposite directions -- folded at 10x has
median 12 against mean 11, segdup at 5x has median 12 against mean 13 -- so both fall under the noise
definition fixed in advance, and neither reproduces across depths or cases.

The lattice argument predicted this was where the mean should help: one median step is about 4 percent
of lambda at 30x and roughly 25 percent at 5x. It does not.

**Verdict: HOLD. The default stays `median`.** The case for the mean now rests entirely on theory --
Fisher consistency for the expected count, the impossibility of a median of integers expressing a value
between half-integers, and agreement with an independently derived simulation value to 0.07 percent.
That is a real argument and it is not what this sweep was run to provide. Three empirical tests have
now been neutral: the 30x ladder, cyp2d6, and this.

**Done in `8d15ef3`:** `--twin-divergence` is the harness default at 40 SNPs, with `TWIN_DIVERGENCE=0`
retained as the exact-representation regression mode so figures measured under the old default remain
reproducible by asking for them.

**But the divergence is narrower than it sounds, and the ladder numbers above must be read
accordingly.** `--twin-divergence` perturbs SNP sites only. Every SV site pushes both `#a` and `#b`
onto the same node, so **twins still share every structural bubble allele**. Holding out a haplotype
therefore removes an identical whole haplotype while leaving the bubble alleles in the panel, and the
exact-allele check skips any block whose truth allele is absent. The 11 to 16 of 16 figures measure
genotyping of KNOWN bubble alleles on an unseen haplotype background. That is a real and useful test.
It is not nearest-haplotype accuracy, and the earlier description of it as "a genuine off-panel test"
was too strong.

Measured on the same run, the two endpoints diverge: bubbles score 14/16 exact, 88 percent, while the
call is the most similar AVAILABLE allele in only 55 of 72 haplotype-blocks, 76 percent, with a mean
identity shortfall of 0.004663. The exact score overstates the behaviour the real-data problem is
about.

Two harness defects were fixed alongside. Both homologues shared one wgsim seed (defect V14), making
their coverage fluctuations identical and removing exactly the per-haplotype depth noise that decides
heterozygous from homozygous; `genotype_sim.sh` was fixed for this long ago and this harness was not,
so every absolute ladder figure in this report predates the fix. And the oracle summary read
`gt_lo.accuracy.tsv` after the pair loop, where each pair had overwritten the last, so it reported ONE
pair under a header that said otherwise; it accumulates now.

### 8.18 Not started

A negative binomial or robust weighted estimator, now understood to require marker-specific efficiency
rather than a change of likelihood family. Tests: there is still no registered `genotype_stats.sh`, and
neither depth commit added one. The needed assertions are a zero-anchor block reporting NA with
REGION_FALLBACK, a shrunk block whose raw and fitted values differ by the expected coefficient, the
audit carrying the final joint depth rather than the first-pass value, direct and indexed runs
honouring the same estimator, and continuous behaviour across 19, 20 and 21 anchors. Items 3 through 7
of section 7.6 are untouched, as are both oracles.

The `--explain-pair` and `--deconvolve` diagnostic paths request the historical median explicitly; that
needs to be either deliberate and labelled, or changed.

### 8.19 Agreed plan before the default moves

Retain `--depth-estimator mean` as the leading experimental mode, default unchanged. Then: 50 to 100
seeds on the fixed lpa pair; lower depths and higher error rates so the saturated ladder becomes
discriminative; per-anchor counts and fragment or clump membership exported; uncertainty estimated by
fragment-level bootstrap rather than from five draws; several real loci, then at least one real
sequencing library with an external single-copy depth control. The default does not move until a
registered genotype regression test exists and the real-library distribution checks have been run.

### 8.20 Open questions

**The cliff's replacement constant.** 200 pseudo-anchors is arbitrary. Should the shrinkage precision
be estimated from within-block and between-block variability instead, and is there enough data per
locus to do that?

**Certified oracle, still unacted.** Minimum summed edit distance over the two assignments is not
equivalent to maximum mean identity when the alignment lengths differ. The 2A caching shortcut remains
valid; the assignment must optimise whichever identity definition is reported.

**Depth on real data.** The simulation-theory control that validated the mean exists only because depth
and error rate are known by construction. What plays that role on a real library, before per-marker
efficiency has been learned from a cohort? The suggestion of unique single-copy flanking markers
normalized by cohort behaviour, with mean q_k constrained to 1, is the current best answer and is
untested here.

---

## 9. Summary and recommended order

### 9.1 What landed

| commit | change |
|---|---|
| `9ac2318` | depth provenance: the audit stops reporting a region fallback as a local measurement, and is written after joint refinement rather than before |
| `5812e48` | `--depth-estimator median\|mean\|trimmed`, with the estimator reaching the joint model and the indexed route rather than only the first pass |
| `7d9bd65` | the `min_anchors` cliff removed, so every block with anchors contributes and is shrunk by its own weight; the region summary stops naming its anchor centre lambda |
| `6c01d4c` | `--dump-anchors`, and the distributional checks it makes possible |
| `fdc8b3d` | `genotype_stats.sh`, the module's first registered test: 30 assertions, verified to fail against the pre-fix build |
| `d0e012c` | `--dump-markers` over every marker with dosage, position and clump; anchor `expected` fixed |
| `3b7500f` | dosage from the truth sequences, so a held-out array is measurable at all; anchor and zero-dosage contracts fixed |
| `ae0da53` | `--fragment-len` reaches marker clumping (V8); index v5 carries `marker_clumps`, `node_first_pos`, `all_kmers` and the construction fragment length (V1) |
| `0c29557` | GQ actually compared; index fragment-length refused when contradicted; double-bypass and normalized truth positions |
| `<this>` | co-primary formula corrected; exploratory seeds excluded; `--fragment-len 0` inheritance; truth tri-state; analysis committed |

**No default option changed. One intended default-path consistency fix can change output slightly:**
commit `7d9bd65` alters the default median path for every block with 1 to 19 anchors, whose anchors now
contribute through shrinkage instead of being discarded. On lpa that moves `mass_bp` by 1 bp. An
earlier revision of this section said default behaviour was unchanged throughout, which contradicted a
measurement recorded two sections earlier in the same document.

### 9.2 What the depth thread found, stated at its actual strength

The median-to-mean correction moves lambda by about 1.1 percent, which took the proof-of-concept's band
selection from 3 of 5 draws to 5 of 5, and cut the bias in `mass_bp` by 81 percent across five draws
for a 25 percent RMSE improvement. That part is measured and holds.

The anchor dump then showed a positive GC-associated trend in the counts. An earlier revision of this
section promoted that into the claim that marker efficiency is the dominant remaining bias source.
**That is withdrawn.** What the data supports:

- a weak association, pooled r 0.0399, r-squared 0.0016;
- surviving within-block centring and flank exclusion, so not a between-block artifact;
- not monotonic across quintiles, and sensitive to how the quintiles are cut;
- mixed per-block slopes, two of eleven negative, one at -6.759;
- from a single wgsim draw, in which no library GC-selection mechanism exists at all, so GC is likely
  standing in for position, context, syncmer ascertainment or a local coverage fluctuation.

The extreme-bin difference is about eight times the width of the lambda window, and that arithmetic is
right, but an extreme-bin difference from one seed is not a stable effect size and should not be
quoted as one.

**Tested and not reproduced.** The pre-registered 50-seed experiment in section 8.16 returns a mean
slope of -0.00641, positive in 16 of 50 draws, with anchors, weighted and per-block estimators all at
zero. The association was a property of one read draw. Marker-specific efficiency is not the
explanation for the array residual in this fixture, and the array per-copy excess reported in 8.12 does
not reproduce either.

Recording plainly that section 9.4 below is a lesson this section had already broken when it was
written: an anchor-only, single-seed association was used to support a claim about the dominant bias in
genotyping. It is enough to reorder the next experiment. It is not enough to establish the conclusion.

### 9.3 What is blocking

**The decisive test is not runnable, and the dump written to enable it is why.** `--dump-anchors`
covers invariant depth anchors only, so it cannot compare them against the informative array markers.
It also has a real defect: its `expected` column is 0 on every row, because it is populated from the
informative-marker expectation map while anchors are filtered against a separate anchor expectation,
and its `actual` column is uniformly 464, which restates the panel size rather than saying anything
about a given marker.

Extending it to informative markers is necessary and not sufficient. Without per-allele multiplicity,
array markers carrying twenty copies would be compared against single-copy anchors and efficiency would
be confounded with dosage.

**Resolved.** `genotype_stats.sh` is registered with 37 assertions, verified to fail against the
pre-fix build, and the dump covers every marker with dosage, position and clump. V1 and V8 are closed.
What blocks the default now is replication, not machinery: one seed, one locus, and the tail at the
array described in 8.12 is uncharacterized. The seed run is pre-registered and not yet started.

### 9.4 A pattern in this pass worth recording

Four claims made in this thread were corrected: that depth was not the problem, that a posterior could
not rescue band selection, the Poisson trimming comparison, and the strength of the GC finding. The
common shape is that each conclusion was drawn from a check narrower than the claim it supported -- one
lambda held fixed, one choice of posterior width, one wrong definition of trimming, one seed and one
marker class.

The stratified descriptive analysis in section 8.10 is what exposed structure that four agreeing
summary statistics had missed. An earlier revision credited that to the distributional tests having
power to discriminate, which overstates it: the Kolmogorov-Smirnov calibration is invalid under marker
dependence, so it was the stratification and not the formal test that did the work.

### 9.5 Recommended order

1. ~~Write and register `genotype_stats.sh`~~ **done**, `fdc8b3d`, now 37 assertions.
2. ~~Fix the anchor `expected` field and extend the dump~~ **done**, `d0e012c` and `3b7500f`.
3. ~~Run the pre-registered seed experiment~~ **done**; negative, see section 8.16. The GC line is
   not supported in this fixture; see the scope statement in section 8.16.
4. ~~Lower-depth ladders~~ **done**, section 8.17: the ladder was trivial rather than saturated, and
   with diverged twins it discriminates while the estimator remains neutral. Higher error rates
   untested.
5. **Fit a cross-validated efficiency model** and test whether it improves held-out mass and
   copy-number accuracy and calibration.
6. Only then reconsider the default.

The certified nearest-pair oracle and the ideal-multiplicity oracle remain ahead of any solver work and
neither has been started. Nothing in this pass constitutes progress on choosing the nearest available
haplotype: bubble accuracy on lpa is unchanged at 7 of 10, and cyp2d6 is identical under both
estimators.
