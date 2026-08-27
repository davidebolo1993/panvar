# Pre-registration: does the linkage constraint survive end to end?

Written 2026-08-27, **before any arm was scored**. Frozen here so the criteria cannot be chosen
after the numbers are seen.

## What is already settled, and what is not

`S5_LINKAGE_RESULT.md` measured that linkage moves calls off a unique block-local emission optimum
93 times over 6 loci x 20 donors: 20 to the right answer, 73 to a wrong one. A threshold rule was
proposed, swept, and implemented as `--max-linkage-emission-loss`, which excludes losing states
**before** forward-backward so the posterior and GQ are recomputed under the constraint.

That experiment's own stated limit is the reason for this one:

> The counterfactual treats blocks independently, and the chain does not. Forcing block *i* to its
> local optimum changes the forward-backward messages and therefore the calls at neighbouring
> blocks. The -62 is what happens if you revert each block in isolation; the real end-to-end effect
> could be larger or smaller. **This is a reason to implement the rule and measure it, not to trust
> the number.**

So the question here is not "is linkage sometimes wrong" -- that is measured. It is whether a real
constrained run, with messages propagating, still wins.

## Frozen before the run

**tau = 0.25.** Taken unchanged from the prior sweep. It is not re-tuned here, and the sweep is
reported only as context, never as the basis of the verdict.

**Arms.** Two, differing in one flag: baseline (flag absent) and `--max-linkage-emission-loss 0.25`.
Same reads, same seeds, same index-free direct route in both.

**Cohort.** 6 loci (cyp2d6, gstm1, lpa, acot, ankrd36c, c4) x 20 donors x {LZO, LOO}. Donors are the
first 20 individuals carrying both haplotypes, in sorted name order -- a rule fixed by the panel, not
chosen. Reads: wgsim, 150 bp paired, insert 350, e=0.001, 30x diploid, per-haplotype seeds derived
deterministically from locus and donor name.

**PRIMARY METRIC: sequence error mass**, summed over graded bubble blocks:

```
error_bp = sum over blocks of (1 - identity) * true_bp
```

`identity` is edit-distance derived. `dbp` is a LENGTH difference and is NOT the metric -- the
project has already published one wrong conclusion by reading it as sequence distance.

**Why this metric is the honest test.** tau was chosen on the prior cohort by counting *allele-index
mismatches*. Sequence error mass was never used to select it, so it cannot have been tuned to. This
is a held-out **metric**, not a held-out donor set: the same donors and graphs are reused, and that
limit is stated here rather than discovered later.

## Declared criteria

The rule is accepted only if **all three** hold:

1. **Cohort sequence error mass decreases.** Total `error_bp` at tau=0.25 is below baseline.
2. **Leave-zero-out is not broken.** LZO total `error_bp` at tau=0.25 must not exceed baseline by
   more than 1%. A setting that damages the regime where the answer is provably in the panel is
   disqualified whatever it does elsewhere -- the rule that retired `--no-region-unique`.
3. **cyp2d6 and gstm1 both improve.** These are the two loci blocking release, and the prior
   experiment found linkage never rescues a unique-optimum block at either (0 rescues against 21 and
   6 overrides). If the rule does not help where it was predicted to help most, the mechanism is not
   the one claimed.

## Declared in advance so it is not read as a surprise

**lpa is expected to be the worst arm and may regress.** It holds every long-range rescue in the
prior data (7 rescues, max move 1.96), which is what a tandem array should look like, since order
there comes from context rather than content. An lpa regression is consistent with the mechanism and
does not by itself reject the rule; criterion 1 is the aggregate test.

**Negative result is a real outcome.** If criteria fail, the rule stays default-off and this file
records the failure. The prior arm S2/S3 at KIV-2 was rejected on exactly this basis.

## What this does NOT establish

- Simulated reads only; real-read confirmation is a separate run.
- Same donors as the prior experiment, so this holds out the metric, not the sample.
- One global tau. The prior per-locus table argues for a per-block-class tau; that is a further
  experiment and is explicitly out of scope here.
