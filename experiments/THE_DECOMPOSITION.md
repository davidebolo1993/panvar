# The one experiment: decompose a failure into its stages

Written 2026-08-25.

## Why nothing has worked

Count the things tried: denser markers, all k-mers, no region filtering, smaller k, adjacency
weighting, mapped coverage, inverse-mean node weighting, bounded outliers, mass windows, marker-rule
variants, per-block rule selection, candidate expansion. **Twelve interventions, zero successes.**

They have one thing in common: every one was an **intervention**, not a **diagnosis**. Each was a
guess about which part of the pipeline was broken, applied end to end, and scored on whether accuracy
moved. When a guess fails it teaches almost nothing, because the failure is consistent with the guess
being wrong *and* with the target stage being fine.

**We still cannot say, for a typical failing call, which stage failed.** Three anecdotes exist —
cyp2d6 block 5 (retained markers cannot separate the answer), KIV-2 (the emission ranks the right
pair 144th on noiseless counts), gstm1 block 11 (no markers exist at all) — and they point at three
*different* stages. There is no distribution, so there is no basis for choosing where to work.

## The pipeline has four places to fail, and they are separable

Between the reads and the reported call:

| stage | question | fails if |
|---|---|---|
| **S1** | do any markers survive selection for this block? | selection strips it bare |
| **S2** | do the surviving markers *distinguish* the best available pair from every other pair? | **retention/selection** discarded the discriminating ones |
| **S3** | fed the *noiseless* counts that pair would produce, does the emission rank it first? | **the likelihood** is wrong |
| **S4** | fed the *observed* counts, does the emission still rank it first? | **noise, depth, or read acquisition** |
| **S5** | does the final call equal the block-local argmax from S4? | **the HMM chain** overrode a correct block answer |

These are mutually exclusive and exhaustive: every failing call lands in exactly one bucket, and each
bucket maps to exactly one place in the code. S2 involves no emission and no reads, so it cannot be
wrong about the model. S3 removes read noise by construction. S5 has **never been tested anywhere** —
every diagnostic to date has been block-local, and the module is a chain.

## The experiment

For every failing (sample, block) in the 15 blocks carrying 80% of the error:

1. Get the **certified** best available pair with `--certified-oracle`. It does exact 2A alignments
   and reports `best_a`, `best_b`, `called_rank`. **Do not use `best_identity`** — that is shortlisted
   by top-16 syncmer Jaccard and is only a lower bound, so it would give the wrong target.
2. Get the retained panel and the observed counts with `--dump-block` (`conf.tsv` carries
   allele, slot, mult, obs).
3. Compute, in order:
   - **S2** — is `vec(best_a, best_b)` unique among all diploid pairs? (pair-vector classes, the
     instrument already used at cyp2d6 block 5)
   - **S3** — score every pair under the production emission on counts `lambda * vec(best_a, best_b)`.
     Rank of the best pair?
   - **S4** — same, on the observed counts. Rank?
   - **S5** — does production's `allele1/allele2` equal S4's argmax?
4. Assign the bucket at the **first** stage that fails.

Run on **simulated reads under leave-one-out first** — that is where the error lives, and clean reads
make S4 a control rather than a confound. Repeat on real reads to size read acquisition separately.

Roughly 15 blocks x 20 donors = 300 decompositions. No new modelling, no new reads, and every
instrument already exists.

## Predictions, recorded before running

So the result cannot be rationalised afterwards:

- **S2 dominates** at cyp2d6 3/5 and lpa 9/15 — blocks where the filter takes retention to near zero
  while all-informative markers separate 103/105.
- **S3 dominates** at lpa 13 (KIV-2), where a noiseless rank of 144 is already measured.
- **S5 is non-zero** at zero-marker blocks. 92% of their PASS calls are exact, so the chain is
  carrying them; if the chain can carry a block to the right answer it can carry one to the wrong
  one, and nobody has checked.
- **S4 is small on simulated reads** and materially larger on real ones.

If S2 dominates overall, the work is marker retention — and specifically a background-aware count
model, since globally keeping contaminated markers is already falsified. If S3 dominates, the work is
the likelihood, and the presence-veto geometry is the standing suspect. If S5 is large, the work is
the chain, and none of the marker work matters until it is fixed.

**If the buckets come out roughly even, that is also an answer**: it means there is no single defect
and the module needs staged rebuilding rather than a fix.

## What makes this different from the twelve

It cannot fail to produce a result. An intervention either moves accuracy or does not, and "does not"
is nearly uninformative. A decomposition returns a distribution over four mutually exclusive causes
whatever the outcome, and that distribution is exactly the missing input to every decision currently
blocked.

It is also the generalisation of the only two diagnostics that ever produced a real finding — Gate 2
at KIV-2 and 9A at cyp2d6 — both of which were single-block, single-stage versions of it.

## Do this before

- the background-aware filter (needs S2 to be the dominant bucket to be worth building);
- the boundary-transition factor (blocked on a second locus regardless);
- GQ calibration (worth doing anyway, but its target depends on which stage produces the confident
  errors).
