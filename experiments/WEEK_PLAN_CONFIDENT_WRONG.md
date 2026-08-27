# One week: make the module's confidence mean something

Written 2026-08-25. Data: `results/genotype_confident_wrong.tsv`, `results/genotype_per_block_summary.tsv`.

## The target

A call is **confidently wrong** when it is reported (`PASS`), materially wrong (`1 - identity > 0.001`),
and wrong beyond what removing the donor forces (`err - panel_floor > 0.001`).

**192 of 6,549 PASS calls (2.93%).** 16 under leave-zero-out, 176 under leave-one-out. They sit in 36
blocks; the top 10 blocks hold 60% of them.

This is the right target because the module's other errors are already declined. Among simulated
leave-zero-out calls that PASS, the wrong-pair rate is **0.00%** — every failure there is flagged
`LOWGQ`. The module mostly knows what it does not know. Where it does not, users cannot tell.

## Two things measured this week that change what to build

### 1. GQ is monotone but its scale is meaningless

GQ is a Phred: GQ = q claims an error probability of 10^(-q/10).

| GQ bin | claimed error | observed confidently-wrong | observed wrong-pair |
|---|---:|---:|---:|
| 10-20 | 3.16% | 8.16% | 41.2% |
| 20-30 | 0.32% | 3.59% | 19.1% |
| 30-40 | 0.03% | **3.21%** | 12.6% |
| 40-60 | ~0.001% | **1.00%** | 6.6% |
| 60-80 | ~0% | 0.49% | 5.4% |
| 80-100 | ~0% | 0.57% | 5.6% |

At GQ 40-60 the observed rate is about **1,000x** what the number claims, and the curve **saturates**
near 0.5% — beyond GQ 60 a higher GQ buys nothing. GQ ranks calls correctly (AUC 0.747 for separating
wrong from right) but its absolute value is unusable, so no threshold a user picks means what they
think it means.

### 2. There is no threshold or filter fix. Both obvious ones were tested and both are bad trades

| intervention | confidently-wrong removed | correct calls lost | ratio |
|---|---:|---:|---:|
| decline blocks with zero markers | 39 | 435 | **1 : 11** |
| decline under 2 markers/allele | 124 | ~1,260 | 1 : 10 |
| require GQ >= 40 | 166 | 3,186 | **1 : 19** |

**92% of the PASS calls at zero-marker blocks are exact.** A block with no markers is usually decided
correctly by the chain's linkage from its neighbours, so "no evidence in this block" does not mean
"cannot be genotyped". Declining them throws away far more than it saves. This kills the intuitive
fix — refuse to call what has no local evidence — and it is worth stating because it is the first
thing anyone proposes.

## The week

### Day 1-2 — Calibrate GQ (the whole point)

Build the calibration curve above on held-out data and fit a monotone map from raw GQ to observed
error rate (isotonic regression, or binned with a floor). Report the calibrated value.

This converts an unusable number into the one users need, and it is what makes the statement
*"this block cannot be genotyped in this sample"* possible — today the module cannot say it, because
GQ 40 and GQ 90 have nearly the same real error rate.

- **Gate:** in a held-out split, observed error in each calibrated-GQ decile is within 2x of claimed.
  Ranking must not change (calibration is monotone, so AUC is invariant by construction — assert it).
- **Deliverable:** `GQ_CAL` alongside `GQ`, plus the curve in the docs. Do not replace `GQ` in place;
  the raw value stays for comparison and for anyone with a tuned threshold.

### Day 2-3 — Find out WHY the model is confident where it should not be

The calibration curve saturating near 0.5% says a floor of irreducible confident error exists. Two
candidate causes, both testable with instruments that already exist:

1. **The posterior is over the wrong space.** GQ comes from the marginal over allele pairs given the
   markers. At a block whose markers cannot separate the truth pair from 34 others (measured at
   cyp2d6 block 5: 10 retained markers, 10 of 105 pair-vector classes), the posterior can be sharp
   and wrong because it is sharp over a class, not over a pair. **Use `--ledger-block` plus the
   pair-vector analysis to compute, per block, the size of the truth pair's equivalence class, and
   test whether GQ falls when that class is large.** If it does not, GQ is ignoring the one thing
   that should cap it.
2. **Linkage confidence is being counted twice.** At zero-marker blocks the call comes from the
   chain, and 92% are right — but the 8% that are wrong carry GQ up to 60.7. Check whether the chain's
   contribution to GQ is discounted by the same clump/ESS factor as the marker contribution.

- **Gate:** a stated mechanism for the saturating floor, with the measurement behind it. Not a fix.

### Day 3-4 — The evidence-rich confidently-wrong

Splitting the 192 by evidence: **124 sit below 2 markers/allele** (GQ should have been lower — day 1-2
handles them) and **68 sit at or above it**, where the evidence is adequate and the model still picks
wrong. Those 68 are the real modelling defect and they concentrate:

| locus | block | n | median err | median GQ | markers | alleles |
|---|---:|---:|---:|---:|---:|---:|
| gstm1 | 5 (backbone) | 34 | 0.0088 | 33.1 | 42 | 9 |
| lpa | 13 (KIV-2 array) | 14 | 0.0426 | **51.5** | 8,754 | 457 |
| ankrd36c | 19 (bubble) | 9 | 0.0083 | 27.6 | 17 | 9 |
| ankrd36c | 17 (bubble) | 7 | 0.0081 | 19.7 | 76 | 56 |
| cyp2d6 | 11 (bubble) | 7 | 0.0021 | 19.2 | 25 | 34 |

**gstm1 block 5 is the best new lead in the cohort**: 9 alleles, 42 markers (4.7 per allele, well
above the threshold where blocks are normally exact), 34 confidently-wrong calls, and the error is
identical at 0.0088 every time. A deterministic, reproducible, evidence-rich failure — the cleanest
possible case for debugging the emission, and it has never been looked at.

lpa block 13 is the known KIV-2 problem and is the highest-GQ failure in the cohort (GQ 99.0). It is
the one place where the module is maximally confident and wrong, and it has resisted every
intervention in the ledger. Do not start there.

- **Gate:** an explanation of gstm1 block 5 that predicts the identical 0.0088, verified against a
  second block.

### Day 4-5 — Surface it

- Per-block reliability in the output: the truth-class size where computable, the calibrated GQ, and
  an explicit `NOEVIDENCE`-style reason string so a user reads *why* a call is declined.
- Update the plots: the filter split is in as of today (filled = reported, hollow = declined), so a
  reader can see which error the module already refuses.
- `docs/modules/genotype.md` and `docs/algorithms/genotype.md` — still owed, and now they have
  something worth saying.

## What this week does NOT do

The three structural workstreams in `BAD_BLOCKS_PLAN.md` — background-aware filtering, the
boundary-transition factor, the KIV-2 emission — are all multi-week and two of them are blocked
(the transition factor needs a second short-allele locus; the filter needs an oracle prototype).
Calibration is worth doing first regardless of how those land, because every one of them will be
evaluated against a confidence number that currently means nothing.
