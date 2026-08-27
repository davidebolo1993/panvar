# Noiseless probe and the present/absent factorial at KIV-2

lpa block 13, leave-one-out, simulated reads (wgsim 30x, e=0.001), 20 held-out donors, full block cap
(`--max-alleles-block 13:512`) so no candidate is pruned. Tables:
`results/genotype_kiv2_noiseless.tsv`, `results/genotype_kiv2_factorial.tsv`.

## What "noiseless" means here, precisely

`--noiseless-counts 13` rewrites block 13's marker counts to `lambda*(m1+m2)+mu`, the emission's own
mean for a stated pair. **Only that block's counts change.** Depth, lambda and dispersion are still
estimated from the reads, and the injected counts are built from that same lambda -- which is what
makes generation and scoring self-consistent, and also means the run is not noiseless end to end.
What it removes is acquisition and sampling noise **at the target block, conditional on the estimated
nuisance parameters**.

Counts are integers and mu is ~0.23, so a marker the source pair lacks is written to **0**, not to mu.
That is the sharp zero-background limit, sharper than any real read pile-up.

## Experiment 1: does removing that noise fix the ranking?

Three arms per donor: `ctl` (counts from a panel pair, probing that same pair -- the control),
`nl` (counts from the truth haplotypes, off-panel under LOO), `obs` (no injection).

| | result |
|---|---|
| on-panel controls at rank 1 | **20/20** |
| failing donors with the certified pair NOT at rank 1 under noiseless counts | **12/12** |
| their ranks | 2, 2, 6, 10, 21, 25, 40, 46, 170, 460, 462, 506 |

The control passing 20/20 shows the injection and the likelihood are internally self-consistent:
counts generated from a panel pair recover that pair, uniquely. So the off-panel result is a
statement about the model, not about the instrument.

**Concordance, 20/20 with no exceptions:** the certified pair is the noiseless emission's unique top
choice **if and only if** the call has zero excess edits. Five donors have rank 1 and all five have
excess 0; fifteen have rank > 1 and all fifteen have nonzero excess, including three near-misses at
3, 4 and 18 edits.

> **Scope.** This predicts EXACT equality, not error magnitude. Spearman against excess is +0.79 over
> all 20 donors but only **+0.51** (rank) and **+0.28** (delta) among the fifteen that fail. The
> strong number is almost entirely the zero/nonzero split.

### What this establishes

The sequence-nearest panel pair is generally **not** the likelihood-nearest pair once the truth is
projected into the retained marker space -- even with target-block noise removed and all 457 alleles
available. Combined with the attribution result (pruning -11,045 edits, i.e. expanding the candidate
set makes it worse), the accurate statement is:

**The current retained-marker representation and emission objective are not aligned with
sequence-optimal off-panel approximation at KIV-2.**

It is NOT principally read sampling, candidate pruning or the HMM. It does not yet distinguish marker
selection, multiplicity representation, the absence background, diploid summing, or the count-derived
nuisance parameters.

## Experiment 2: the present/absent factorial, and its limit

2x2 over which half of the block's markers the injection rewrites, with 3-4 FIXED reference pairs
probed in every arm (certified, the observed optimum, the all-injected optimum, production's call).
Fixed-pair raw-score contrasts are the statistic; ranks and per-arm deltas move with the arm's own
optimum and are secondary.

The partition identity `C(all) - C(present) - C(absent) + C(observed) == 0` holds to **4.7e-10** over
all 25 contrasts, so the halves partition the marker set exactly.

| stratum | present | absent |
|---|---:|---:|
| all contrasts (n=25) | 98.5% | 1.5% |
| worst half of donors | 97.7% | 2.3% |
| milder half | 100.0% | 0.0% |

> **CORRECTION -- this does NOT show that absent markers are unimportant.** The factorial decomposes
> `C(noiseless) - C(observed)`: the CHANGE the injection causes. Truth-absent markers carry ~500
> counts across ~6,700 markers (0.08 each) and were already almost all zero, so rewriting them to zero
> was very nearly a no-op -- about 0.1% of the block's 435,000 marker counts. The absent half could
> dominate the ranking through the veto term and this experiment would not see it.
>
> The claim it supports is only: **removing block-level read noise acts almost entirely through
> markers the truth carries.** An earlier version of this document concluded "the work is the
> multiplicity model for carried markers". That did not follow and is withdrawn; the absolute
> decomposition in `PHASE0_*` is what answers it.

Five donors produce no contrast at all -- certified pair, emission optimum and production call
coincide -- and they are exactly the five with zero excess. The instrument has nothing to measure
where nothing is wrong.

## Two hypotheses of mine that these runs killed

- **"Noiseless is systematically worse, so clean zeros sharpen the veto."** Reported as 8/8. The
  donors run worst-first by excess, so that was the first eight rows of a sorted list; across all 20
  it is **9/20**. A sorted prefix is not a cohort result.
- **"The work is multiplicity scaling of carried markers."** See the correction above.

## Standing caveats

- one block, one locus, simulated reads, leave-one-out;
- cyp2d6 and gstm1 fail through a different bucket (the emission cannot separate), so nothing here
  should be read across to them;
- all 20 donors are dev under the frozen split in `results/genotype_donor_split.tsv`.
