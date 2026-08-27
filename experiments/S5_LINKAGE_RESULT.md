# S5: linkage, measured in both directions

6 loci x 20 donors, leave-one-out, simulated reads. 2,300 calls, 1,623 representable, 419 wrong.

## Why both directions

Linkage demonstrably earns its place somewhere: blocks with **zero markers are called correctly 92%
of the time**, which can only be the chain carrying them. A measurement that counted only the cases
where linkage lost would conclude "linkage is harmful" from a one-sided ledger.

`called_emission_rank` / `called_emission_delta` (added for this) give the other side: how far the
reported pair sits below the block-local emission optimum. Delta 0 means the chain kept the local
optimum; a large negative delta means it overrode strong local evidence.

## Result: at blocks with a UNIQUE local optimum

Restricted to unique optima because that is the only place "revert to the local optimum" names a
single pair. Where the optimum is tied it names a set and reverting picks arbitrarily — simulating
that as a win would flatter the rule. This restriction matters: on the full moved set the rescue tail
ran further and the trade-off looked better than it is.

**93 calls moved off a unique local optimum: 20 rescued, 73 overrode.**

| | n | median move | p90 | max |
|---|---:|---:|---:|---:|
| linkage **rescued** (moved to the right answer) | 20 | 0.15 | 0.80 | 1.96 |
| linkage **overrode** (moved to a wrong answer) | 73 | **1.57** | 4.40 | 6.99 |

Rescues cluster at near-ties; overrides at real margins. That separation is what makes a threshold
possible at all.

Of the 73 overrides, **71 are genuinely fixable** — the unique local optimum IS the truth pair, so
reverting corrects the call. (The other 2 have `truth_rank > 1`: the local optimum is also wrong, and
reverting would not help. Checked rather than assumed; without the check the numbers below would have
been overstated by 2.)

## Candidate rule and its trade-off

*Linkage may move off the block-local optimum only when the local margin is under tau; above it, the
emission stands.*

| tau | rescues lost | overrides genuinely fixed | net errors |
|---:|---:|---:|---:|
| 0.10 | 14 | 69 | -55 |
| **0.25** | **6** | **68** | **-62** |
| 0.50 | 4 | 64 | -60 |
| 1.00 | 1 | 52 | -51 |
| 2.00 | 0 | 29 | -29 |

At tau = 0.25 this removes **62 of 419 leave-one-out errors (15%)**, and 62 of the 77 errors that
occur at unique-optimum blocks (80%).

## Per locus, and the part that decides the design

Genuinely-fixable overrides, with rescues at the same blocks:

| locus | rescued | overrode | max rescue move | median override move |
|---|---:|---:|---:|---:|
| cyp2d6 | **0** | 21 | — | 1.52 |
| gstm1 | **0** | 6 | — | 2.65 |
| lpa | 7 | 23 | **1.96** | 1.68 |
| ankrd36c | 8 | 10 | 0.45 | 1.41 |
| c4 | 2 | 8 | 0.52 | 2.52 |
| acot | 3 | 5 | 0.14 | 0.89 |

**At cyp2d6 and gstm1 linkage never rescues a unique-optimum block — it only ever loses.** Those are
the two loci with the largest model gaps. lpa is the opposite: it holds every long-range rescue
(max 1.96), which is what a tandem array should look like, since order there genuinely comes from
context rather than content.

That argues for a **per-block-class** rule rather than one global tau: an array block earns a larger
tau, a simple bubble a small one. The data support that, but it is a second experiment, not this one.

## What this does NOT establish

**The counterfactual treats blocks independently, and the chain does not.** Forcing block *i* to its
local optimum changes the forward-backward messages and therefore the calls at neighbouring blocks.
The -62 is what happens if you revert each block in isolation; the real end-to-end effect could be
larger or smaller. **This is a reason to implement the rule and measure it, not to trust the number.**

Other limits:
- simulated reads only; the real-read cohort is not re-run;
- tau is chosen on the same data it is scored on, so it needs a held-out split before it is a value
  rather than an illustration;
- errors are counted as allele-index mismatches, not sequence-error mass, so a fixed call and a
  trivial one weigh the same;
- covers only the 70% of leave-one-out blocks whose truth pair is reachable.

## Next

1. Implement behind a flag (`--linkage-margin`), default off so it cannot change current output.
2. Registered tests: at 0 the rule is inert and output is byte-identical; direct and indexed agree.
3. Score end-to-end on a held-out donor split, by sequence-error mass as well as index mismatch, on
   both real and simulated reads.
4. Only then consider the per-block-class tau the locus table above suggests.
