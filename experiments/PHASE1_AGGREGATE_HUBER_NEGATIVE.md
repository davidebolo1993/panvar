# Negative result: the aggregate-Huber construction

Kept so it is not re-attempted, and named distinctly so it is not confused with the clump-normalized
score in `PHASE1_ORACLE_PREREGISTRATION.md`.

## What it was

Multiplicities summed within a clump, then a Huber loss on the aggregate residual:

```
score(a,b) = - sum over clumps k of Huber_2.0( (sum_k mhat - sum_k m(a,b)) / sqrt(max(1, sum_k m(a,b))) )
```

plus, in the second variant, a total-dosage term at `w = 1.0`.

## Why it does not do the job

**It never removed marker-count replication inside a clump.** Summing multiplicities means a clump
holding 20 markers still contributes 20 times the mass of a clump holding one, which is the precise
thing Phase 0 says has to stop. It bounds the aggregate, not the influence.

## Measured, 15 development donors, noiseless counts

Certified pair ranked first: **S1 0/15, aggregate-Huber 2/15, +dosage 2/15**.

Ranks improved substantially -- 462 to 137, 506 to 127, 46 to 17 -- but **rank improvement is not
sequence improvement**, which is why excess and not rank is the pre-registered outcome. Of the eight
donors whose selected pair could be priced from existing measurements:

- it reproduced production's exact call in **five**;
- on HG00344 it selected the pair worth 22,063 excess against production's 16,534;
- it recovered the certified pair on two (HG00350, HG00320).

The dosage term was inert: about 5 log units against Huber sums in the hundreds, so the two variants
agreed almost everywhere. That is consistent with the standing advice not to add a dosage term before
the support geometry is fixed.
