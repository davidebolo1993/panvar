# Pre-registration: does constraining linkage reduce genotyping error?

Written 2026-08-25, **before any constrained run was scored**. The donor split and the threshold grid
below are frozen; the split file is `results/genotype_donor_split.tsv`.

## What is being tested

`--max-linkage-emission-loss <tau>` excludes any diploid state losing more than `tau` to its block's
emission optimum, **before** forward-backward, so posterior and GQ are computed under the constraint.
Default is infinity and reproduces current output byte for byte.

## Why

Exploratory measurement over 6 loci x 20 donors (simulated, leave-one-out): the chain moved off a
**unique** block-local optimum 93 times — 20 rescued the call, 73 broke it — and the two separate by
how far it moved (rescues median 0.15, max 1.96; overrides median 1.57, max 6.99). A per-block
counterfactual suggested tau ≈ 0.25 would remove 62 of 419 errors.

**That number is exploratory and is not the claim under test.** It reverts blocks independently while
the chain is coupled: constraining block *i* changes the forward-backward messages and therefore the
calls at its neighbours. The end-to-end effect is what this experiment measures, and it may be larger
or smaller.

## Frozen split

42 donors, assigned by `sha256(donor_id) % 2` so a donor lands in the **same arm at every locus** —
splitting by row or by locus would let one donor train at one locus and validate at another.

- **dev** (20): tuning and inspection.
- **val** (22): scored once, at the end, at the tau chosen on dev.

## Frozen grid

`tau in {0, 0.1, 0.25, 0.5, 1, 2, inf}`, plus a **local-emission-only** arm (no chain) as the opposite
extreme. Seven constrained points and two references.

## Primary outcome

**Certified excess sequence-error mass**: summed edit distance of the called pair to truth, minus the
same for the best *reachable* pair, over the certified 2A oracle rather than the top-16 Jaccard
shortlist. Allele-index mismatch is **secondary** — sequence-equivalent alleles carry different
indices, so index accuracy can move without any sequence changing.

Also reported, all pre-declared:

- PASS errors and call rate (a rule that improves accuracy by declining more is not an improvement);
- rescues and overrides counted separately, so a net figure cannot hide a bad trade;
- LZO and LOO separately;
- simulated and real reads separately;
- representable blocks separately from the certified off-panel set;
- by depth if the cohort supports it, since emission-loss units are not depth-invariant — a fixed tau
  may mean different things at 15x and 45x, which would make a single global value the wrong shape.

## Gate

Accept a finite tau only if, **on val**:

1. certified excess sequence-error mass falls;
2. call rate does not fall by more than 1 point;
3. LZO does not regress (it is 99.7% correct and has almost nothing to gain);
4. no locus regresses on the primary outcome by more than it gains elsewhere;
5. the direction agrees between real and simulated reads.

## Declared risks

- **Coupling.** The exploratory estimate is per-block; the real one is not. If the end-to-end gain is
  much smaller than 62 errors, the coupling absorbed it, and that is a result, not a failure to
  explain away.
- **tau is not scale-free.** Emission loss is in log-likelihood units, which grow with marker count
  and depth. A block with 8,754 markers and one with 7 do not have comparable margins. If the best
  tau differs sharply by marker count, a *relative* rule (a fraction of the block's own spread) is
  the right shape and this absolute one should be retired rather than tuned.
- **Internal, not external, validation.** This cohort has been explored extensively across many
  experiments. A held-out donor split controls for donor-specific overfitting; it does not make the
  result independent evidence.

## Explicitly out of scope

**Per-locus thresholds are not fitted.** The exploratory data show cyp2d6 and gstm1 with zero rescues
and lpa holding every long-range rescue, which makes a per-class tau tempting. Fitting it on six loci
would overfit them. A biologically defined class rule (array vs simple) may be pre-registered as a
**separate** experiment after a global tau is validated, and only with enough held-out rescue and
override events per class to support it.
