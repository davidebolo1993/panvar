# KIV-2 attribution: complete, 20 donors

lpa block 13, leave-one-out, simulated reads (wgsim 30x, e=0.001, per-homologue seeds), 20 held-out
donors x 2 arms = 40 runs, all complete. Per-donor table in
`results/genotype_kiv2_attribution.tsv`.

Metric is **certified excess edit distance**: `called_total_edits - best_total_edits`, where `best`
comes from the certified 2A oracle over all 457 block-13 alleles. It is not the top-16 Jaccard
shortlist and it is not a length difference.

## The three arms

| arm | what it removes |
|---|---|
| **base** | nothing (production; `--max-alleles 64` at every block) |
| **b13** | **pruning**, at block 13 only (`--max-alleles-block 13:512`) |
| **b13t0** | pruning **and** linkage (`--max-linkage-emission-loss 0`: the chain may not move off the block-local emission optimum) |

Only block 13 is expanded. Expanding the whole locus would also change the neighbouring emissions the
chain sees, and a difference at block 13 could then have come from a neighbour. An earlier version of
this experiment made exactly that mistake.

## Result

| | edits |
|---|---:|
| certified best (the floor) | 2,542 |
| **baseline excess** | **126,664** |
| pruning (base → b13) | **−11,045** |
| linkage (b13 → b13t0) | +16,307 |
| **residual (b13t0 excess)** | **+121,402 = 96%** |

**Pruning is negative: expanding the candidate set makes it worse.** It changed the call in 3 of 20
donors and only one of those improved. HG00344 goes 16,691 → 27,747 — the larger set contains a pair
the emission likes more and truth likes less.

Linkage helps, but almost all of the +16,307 is two donors (HG00146 −10,762, HG00344 −5,527), and
HG00344's share is repairing damage the expanded cap did in the same arm. Its constrained result
(22,220) is still worse than production's (16,691).

**Even taking the best of all three arms per donor with full hindsight** — not a rule anyone could
apply — the total falls only 126,664 → 115,869, recovering **9%**.

This retracts an earlier claim of mine. I had observed that the certified-optimal pair is often
outside the 64-candidate set and read that as the cause. The observation is true and the causal
reading was wrong: restoring those candidates does not fix the call.

## What the residual is made of

Eight of 20 donors are already essentially exact (excess ≤ 18). The other twelve fail, and **every one
of them fails by an integer number of KIV-2 repeat units**, u ≈ 5,539.5 bp by least squares over the
twelve nonzero excesses.

> **CORRECTION.** A first version of this section argued from that alone that the module "never picks
> a wrongly-composed array of the right size", using **diploid** called-vs-truth length. That argument
> does not work and GPT was right to reject it: at three donors the diploid total is correct to within
> 33 bp while the excess runs to 33,457 edits, which on the diploid instrument looks exactly like a
> composition error. The conclusion survives, but only on the **per-homologue** instrument, and the
> reason it survives turns out to be the more useful finding.

Per-homologue length error under the better haplotype assignment (`called_h1_lenerr`,
`called_h2_lenerr` from the full certified oracle; per-donor table in
`results/genotype_kiv2_mechanism.tsv`):

| donor | excess | h1 units | h2 units | diploid units | mechanism |
|---|---:|---:|---:|---:|---|
| HG00146 | 33,457 | **+3.02** | **−3.01** | 0.01 | dosage-split |
| HG00344 | 16,534 | −3.01 | −0.00 | −3.01 | total-dosage |
| HG00140 | 16,502 | −1.00 | +2.00 | 1.00 | mixed |
| HG00350 | 11,156 | **+1.00** | **−1.00** | −0.00 | dosage-split |
| HG00329 | 11,079 | **−1.00** | **+1.00** | 0.00 | dosage-split |
| 7 donors | 5,112–5,537 | ±1 or 0 | ±1 or 0 | ±1.00 | total-dosage |
| 8 donors | 0–18 | 0 | 0 | 0 | (correct) |

Every per-homologue error is an integer number of units — the largest deviation from a whole multiple
is 0.02.

### Splitting the excess into length and non-length, exactly

> **CORRECTION.** A first version computed "edits remaining once repeat count is accounted for" as
> `excess − (|h1_lenerr| + |h2_lenerr|)`. That subtracts an **absolute** quantity from a **difference**
> and is not a decomposition of anything: `excess` is measured against the certified pair, so the
> certified pair's own length error has to come off too. It produced negative residuals — −415 edits,
> −375% — which I printed and did not flag. Impossible values in an output are the cheapest error
> signal there is and I walked past them.

Excess is `D_called − D_certified`. Writing `L = |Δlen_h1| + |Δlen_h2|` for each pair:

```
length excess      =  L_called − L_certified
non-length excess  = (D_called − L_called) − (D_certified − L_certified)
```

| component | edits |
|---|---:|
| length contribution | **127,583** |
| non-length contribution | **−919** |
| **sum (= excess)** | **126,664** ✓ |

The non-length term is **negative**. Beyond the unavoidable length difference the called pairs carry
1,353 edits and the certified pairs carry 2,272 — so relative to the certified solution, composition
contributes *no* positive excess at all. (The certified pair minimises **total** edits, which can
trade a little non-length divergence for a better length match; it is not the composition-optimal
pair.)

So the honest statement is: **virtually all excess relative to the certified pair is per-homologue
length. Non-length sequence differences exist, but are no greater — and collectively smaller — than
the certified pairs' own.** The module does not misassemble arrays; it miscounts them.

The mechanism buckets below classify *failing calls* and are a separate cut from this decomposition:
the "composition at correct dosage" bucket holds 25 edits of excess because those donors are
essentially correct, not because composition was measured to cost 25 edits.

### The split that the diploid total hides — 44% of the loss

| mechanism | edits | share |
|---|---:|---:|
| **dosage-split** (diploid count right, allocation between homologues wrong) | 55,692 | **44.0%** |
| total-dosage (diploid count itself wrong) | 54,445 | 43.0% |
| mixed | 16,502 | 13.0% |
| composition at correct dosage | 25 | 0.0% |

At HG00146 the two homologues are +3 and −3 units: the diploid array is the right total size, made of
the wrong two pieces. Same at HG00350 and HG00329 at ±1.

**This changes what §2c and C2 of the plan are worth.** Those say the continuous dosage estimator
(`mass_bp`, accurate to ~0.42 repeat units) should be the copy-number answer rather than the selected
pair's length. That is a statement about the **diploid total** — and the diploid total is *already
correct* in 44% of the loss. A better diploid estimator cannot touch those cases, because nothing
about the diploid number is wrong. **Allocating a correct total between two homologues is a separate
problem and no part of the module currently addresses it.**

## What this does NOT establish

The residual is measured on **observed simulated reads**. It removes pruning and linkage; it does not
remove read acquisition, sequencing error, count sampling or depth estimation. So "the likelihood
prefers the wrong pair" is established *for the emission evaluated on observed counts* — not yet for
the likelihood in the abstract. The noiseless-count probe is what narrows that, and it narrows it
rather than closing it: only the **target block's** counts are made noiseless. Depth, lambda and
dispersion still come from the reads, and the injected counts are built from that same lambda, so what
it removes is acquisition and sampling noise *at that block, conditional on the estimated nuisance
parameters*.

Other limits:

- one block at one locus; the decomposition's S3/S4 bucket is 12 of 18 cases at lpa, so this block is
  where the likelihood contributes at all, and the finding should not be read across to cyp2d6 or
  gstm1, whose errors are S2 (the emission cannot separate).
- `u` is fitted from these 12 observations, not measured from the panel's allele lengths. It is an
  empirical step consistent with the array's repeat, not an independent measurement of it.
- "composition is not the problem" is measured **at this block, at these dosages, and relative to the
  certified pair**. It says the called pairs are no worse than the certified ones in non-length terms;
  it does not say composition would stay right if the counting were fixed and harder cases became
  reachable.
- **repeat count is inferred from length ÷ fitted unit, not counted.** The KIV-2 unit structure makes
  that reading compelling, but two alleles of equal length could differ in unit composition. Counting
  units explicitly from the walk would settle it and has not been done.
- simulated reads only.

## Next

The interpretation split stands, but the target is now specific: the emission's preference is a
function of repeat count, and the failure is in **how many copies each homologue gets**, not which
sequence they are made of.

1. **Cohort-wide noiseless probe.** Rank every certified-optimal pair under the full block cap on
   noiseless expected counts from the true haplotypes, and again on observed counts. *Correction
   accepted: this is not a new experiment for the motivating pair* — `GENOTYPE_NEXT_STEP.md` already
   records the certified pair ranking **144** on held-out truth multiplicities while the production
   pair ranked 1. A genuine noiseless emission failure is therefore already proven for at least one
   donor. What is unrun is the 20-donor version, i.e. whether that generalises.
   (`--noiseless-counts <BLK|BLK:A,B>` is written for this; the `BLK:A,B` form is the on-panel control
   a proper likelihood cannot fail.)
2. Record, alongside the rank, the **total multiplicity difference and per-homologue length
   difference** between called and certified pairs — the quantities the mechanism table says the
   preference actually tracks.

Given the mechanism split, the architecture this points at is **three separated outputs**, not one
argmax over pairs:

| output | current state |
|---|---|
| diploid total dosage | already good — `mass_bp` reaches ~0.42 units; 43% of the loss is here |
| **allocation between homologues** | **not modelled at all; 44% of the loss** |
| composition given dosage | essentially solved — 0.0% of the loss |

The middle row is the gap, and it is the one nobody has looked at. It is also what would make mosaics
identifiable, since a mosaic is precisely a statement about allocation rather than total.
