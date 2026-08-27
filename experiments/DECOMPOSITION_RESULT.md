# Decomposition result (COMPLETE: 6 loci x 20 donors x 2 regimes)

Method frozen in `THE_DECOMPOSITION.md`. Simulated reads (wgsim, 30x, e=0.001, per-homologue seeds),
6 loci x 20 donors x {LZO, LOO}, direct route (no index). Every failing call is assigned to the
**first** stage that fails, so the buckets are mutually exclusive.

## Stages

| | question | fails ⇒ blame |
|---|---|---|
| S1 | was the truth pair pruned before scoring? | candidate pruning |
| S2 | is the truth pair at a **non-unique maximum** of the emission on observed counts? | see the correction below |
| S3/S4 | does the emission rank truth below the top? | **likelihood** (S3) or **read noise** (S4) |
| S5 | did the chain override a *uniquely* correct emission? | **linkage / HMM** |

> **CORRECTION (2026-08-25, after review). S2 as measured is NOT a pure evidence failure, and the
> headline "evidence, not likelihood" is premature.**
>
> The frozen method said S2 would compare raw retained-marker pair vectors, independently of the
> likelihood. **That stage was not run.** What is measured instead is that the production emission has
> several pairs at its maximum on observed counts, which lumps together:
> (a) pairs with *identical* retained-marker vectors — a genuine identifiability failure;
> (b) pairs with *different* vectors the likelihood happens to score equally;
> (c) ties induced by mass constraints, numerical structure or count sampling.
>
> Only (a) is a marker problem. Since S3/S4 counts only cases where truth scores strictly *below* the
> maximum, every likelihood-induced tie is currently charged to S2. **"The likelihood accounts for
> only 4%" is therefore a lower bound, not a finding.** The established claim is narrower and still
> useful: *in 77% of representable leave-one-out errors the truth pair sits at a non-unique maximum of
> the current emission.*
>
> Completing the promised separation needs: raw retained-vector comparison, then the production
> likelihood on noiseless expected counts, then on observed counts, then the chain.

**S3 and S4 are NOT separated here.** wgsim at e=0.001 and 30x carries sequencing and sampling noise,
so this bucket mixes "the likelihood is wrong" with "noise overturned it". Separating them needs a
re-run on synthetic noiseless counts and has not been done.

## The instrument that made this possible, and the bug it caught

`truth_rank` counts only STRICTLY better pairs, so a block whose markers separate nothing reports
rank 1 for **every** pair. A first pass reading rank 1 as "the emission got it right" concluded that
**94% of leave-one-out errors were the chain overriding a correct emission**. That was an artifact.

`truth_emission_ties` (added for this) counts pairs within 1e-9 of the best, over the same triangle
the rank loop uses. At gstm1 block 11 it reads 231 for 21 alleles — every pair in the block — because
the block has no markers at all. With it, "uniquely preferred" and "said nothing" separate.

Registered in `tests/genotype_stats.sh`: zero markers ⇒ all A(A+1)/2 pairs tie; separating markers ⇒
exactly one pair at the top; every count within those bounds.

## Result

### LZO — 2,175 calls, 99.7% correct

7 wrong: 6 S2, 1 S3/S4. The control behaves.

### LOO — 419 wrong calls at representable blocks

| bucket | n | share of wrong |
|---|---:|---:|
| **S2 — emission cannot separate** | 323 | **77.1%** |
| **S5 — chain overrode a unique optimum** | 71 | **16.9%** |
| S3/S4 — emission ranked truth below top | 18 | 4.3% |
| S1 — pruned before scoring | 7 | 1.7% |

Only **37% of wrong calls are reported** (`PASS`); the rest are declined. But the reported share
differs sharply by bucket, and it inverts the priority a raw count would give:

| bucket | n | reported as PASS |
|---|---:|---:|
| S1 | 7 | 0% |
| S2 | 323 | 34.7% |
| S3/S4 | 18 | 38.9% |
| **S5** | 71 | **53.5%** |

**S5 failures are the most likely to reach a user.**

Per locus:

| locus | S1 | S2 | S3/S4 | S5 | total |
|---|---:|---:|---:|---:|---:|
| lpa | 0 | 81 | **12** | 21 | 114 |
| gstm1 | 0 | **103** | 0 | 6 | 109 |
| cyp2d6 | 0 | **61** | 0 | 21 | 82 |
| ankrd36c | 1 | 44 | 2 | 10 | 57 |
| acot | 6 | 30 | 4 | 5 | 45 |
| c4 | 0 | 4 | 0 | **8** | 12 |

lpa is the only locus where the likelihood contributes materially (12 of 18 S3/S4 cases across the
whole cohort), which matches the standing KIV-2 result. c4 is mostly S5 — its few errors are linkage,
not evidence.

## What this changes

**The likelihood is the smallest bucket that is cleanly attributable (4.3%)** — but see the
correction above: likelihood-induced ties fall into S2, so this is a floor. What can be said is that
no intervention aimed at the emission's *ranking* can address the 77% where the emission has no
unique maximum to rank.

**S2 is not one disease.** Splitting by how tied the block is:

The tie count is over the candidates actually SCORED (pruned to `--max-alleles`, default 64), not
over every allele in the block. A first pass divided by the full allele count, which understates how
tied a large block is; 110 of the 323 S2 observations come from blocks above the cap. Corrected:

| S2 stratum | n | median markers | median tied pairs |
|---|---:|---:|---:|
| ALL scored pairs tied (no evidence) | 59 | **0** | 231 |
| **10-100% tied** | **110** | **11** | 210 |
| 1-10% tied | 85 | 398 | 45 |
| under 1% tied (near-miss) | 69 | 1,219 | 3 |

> **CORRECTION.** With the wrong denominator this read 49 / 114 / 55 / 105 and I reported that the
> largest stratum was a near-miss group at ~1,000 markers. **It is not.** The largest stratum is
> blocks that are 10-100% tied with a median of **11 markers**, and 169 of 323 S2 cases (52%) are at
> least 10% tied. S2 is dominated by *low-evidence* blocks, not by last-mile resolution failures.
> `n_scored_alleles` is now emitted so the denominator never has to be inferred again.

**S5 is real and unexamined.** 71 calls where the emission had a *unique* optimum and the chain chose
otherwise. Median 1,236 markers, median GQ 12.9 — well-supplied blocks, and the module was
appropriately unconfident. This is a **floor**: ties are charged to S2 by the first-failing-stage
rule, so any case where the chain also had a hand is counted earlier.

## Caveats

- Covers only the ~70% of LOO blocks where the truth pair remains in the reduced panel. The other 30%
  need the certified 2A oracle to establish the best reachable pair; not yet run.
- `best_identity` elsewhere in the cohort is a top-16 Jaccard shortlist and is a lower bound — it
  falls below the called identity in 4 observations. Not used in this decomposition, but it is used
  in the per-block plots' "panel floor".
- S2 conflates identifiability with likelihood-induced ties (see the correction above); the raw
  pair-vector stage is still owed.
- Sequence-equivalent alleles are scored as wrong if the index differs. The certified-optimal *set*,
  not a single best pair, is the correct target; not yet applied.
