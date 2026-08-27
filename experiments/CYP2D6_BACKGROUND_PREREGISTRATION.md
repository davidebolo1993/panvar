# Pre-registration: does shared-marker background explain CYP2D6?

Written before any arm below was run. cyp2d6 blocks 3 and 5, simulated reads, leave-one-out.

## Why this locus and this block

Largest model-only gap of any locus: simulated off-panel mean 0.0223, against lpa's 0.0100.

Measured at block 5, 14 alleles, 105 diploid pairs, identical across 8 samples and both regimes:

| marker set | distinct node markers | distinct pair-vectors | ambiguous pairs |
|---|---:|---:|---:|
| production (retained) | 10 | 10 / 105 | 102 |
| relax confinement alone | 10 | 10 / 105 | **102 -- no change** |
| relax over-expected alone | 20 | 28 / 105 | 90 |
| **all informative** | **36** | **103 / 105** | **3** |

Marker fates, freshly confirmed: **10 retained, 10 fail over-expected only, 16 fail both, 0 fail
confinement only.** Relaxing either rule alone is useless. Every dropped marker names **block 3** and
no other, so the coupling is a single pair of blocks.

## The four arms. Frozen.

| arm | markers | target block state | other blocks |
|---|---|---|---|
| **A1 production** | 10 retained | inferred | as production does |
| **A2 unfiltered** | all 36 | inferred | ignored -- raw observed counts |
| **A3 conditional** | all 36 | block 5 inferred | contribution supplied from truth, block 3 included |
| **A4 joint** | all 36 | blocks 3 and 5 inferred **together** | truth supplied for the rest only |

A3 and A4 are each run on **noiseless** and on **observed** counts. A2 is the control that says
whether restoring markers without a background model is enough on its own -- it should not be, since
`--no-region-unique` already broke leave-zero-out.

**The joint arm is not optional.** A conditional oracle can pass merely by being handed block 5's true
state while genotyping block 3, or the reverse. Production knows neither.

## Background enters the MEAN, never the counts

```
mean_j(candidate) = mu_j
                  + lambda * multiplicity_j(candidate at the target block)
                  + lambda * multiplicity_j(truth at every other block carrying j)
```

Subtracting the external contribution from the observed count can drive it negative and destroys the
count likelihood. An earlier draft of the plan said "subtract"; that is wrong and is corrected here.

## Interpretation, declared in advance

| A3 | A4 | reading |
|---|---|---|
| pass | pass | strong evidence for a joint shared-marker factor -- proceed to prototype |
| pass | fail | the information exists; **multi-block inference** is the binding problem |
| fail (noiseless) | - | marker filtering is **not** the complete diagnosis; the route is closed |
| pass noiseless, fail observed | - | read acquisition / count contamination is binding, not filtering |

## Success requires all of

1. all 36 markers, **including the 16 that fail both filters**, actually reach scoring -- asserted,
   not assumed;
2. on-panel noiseless controls remain correct;
3. the **top equivalence class contains the certified pair**;
4. certified excess improves by the release plan's threshold (>= 10%);
5. the improvement is not carried by a single donor;
6. leave-zero-out controls do not regress;
7. it repeats on donors not used for design.

### On the 3/105 ceiling

Three pairs remain in one indistinguishable equivalence class even with all 36 markers. That is **not
three unavoidable wrong calls.** The oracle succeeds if it identifies the correct equivalence class
and reports all three members where the evidence cannot separate them. Gating on a single exact pair
would fail a correct answer.

## Out of scope here

- Any production change before A4 passes on untouched donors.
- Estimating the background from the caller's own calls and feeding it back as truth -- constraint
  C-1, learned from two separate failures, and a route to a stable wrong fixed point.
- Adjacency evidence at this block: closed. Retained edges give the same 10 / 105 and are filtered at
  the same rate as nodes, 11 of 41 against 10 of 36. That closes "add the currently retained adjacent
  syncmer edges". It does **not** falsify whole-fragment evidence built from unfiltered markers, which
  stays a lower-priority option here and the primary direction for gstm1.
