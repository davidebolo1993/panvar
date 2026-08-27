# panvar `genotype`: six-locus held-out cohort, real and simulated reads

Date: 2026-08-23. Self-contained report for external review. No production source was changed;
`tests/genotype_sim.sh` gained a self-pairing default and sequence-distance reporting (commit
`cc00cd2`), and `src/genotype_command.cpp` had its help text corrected (`5fe5f36`).

---

## 1. Design

Six real loci. Twenty held-out individuals each, run in **two regimes from identical reads**:

- **LZO** (leave-zero-out): the individual's two assembly haplotypes remain in the panel. An exact
  call is reachable, so this is an implementation ceiling — anything lost here is a defect, not a
  limitation.
- **LOO** (leave-one-out): both haplotypes excluded via `--exclude-haplotypes`. Generalisation.

Both regimes score against those same two haplotypes (`--truth-haplotypes`). Default options
throughout; no sweeps, no per-sample tuning, no HMM changes.

Run twice over the same individuals with two read sources:

- **real**: 1000G 30x CRAMs, region pulled remotely with `samtools view <url> <region>` over the
  locus +/- 20 kb derived from each graph's own PanSN GRCh38 path name. Region-mapped reads only;
  no unmapped-read rescue.
- **simulated**: `wgsim` at 30x directly from the two haplotype sequences spelled out of the graph.
  No mapping step, therefore no unmapped or mismapped fraction.

234 sample-locus runs, 4,764 scored blocks.

### Metric definitions, because two of them are easy to confuse

- `identity` — edit-distance identity of the called allele against the truth haplotype (NW, edlib).
- `best_identity` — the best any panel allele reaches, computed over a top-16 syncmer-Jaccard
  shortlist. It is therefore a LOWER bound on the true ceiling: real headroom is at least what is
  shown.
- **`gap` = `best_identity - identity`** — model headroom. This is the headline.
- `dbp` — `|len(called) - len(truth)|` summed per homologue. **A LENGTH difference, not a sequence
  distance.** `best_dbp` is the minimum of that over all panel alleles, and the source comment warns
  it "only says some allele matched the truth's LENGTH, which among hundreds of alleles happens by
  coincidence".

An earlier draft of this analysis read `dbp` as sequence distance and drew two wrong conclusions from
it (section 5).

---

## 2. Result: real reads

| locus | regime | n | blocks | identity | panel best | gap | exact |
|---|---|---:|---:|---:|---:|---:|---:|
| cyp2d6 | LZO | 18 | 342 | 0.99846 | 1.00000 | 0.00153 | 93% |
| cyp2d6 | LOO | 18 | 342 | 0.97700 | 0.99809 | **0.02109** | 44% |
| c4 | LZO | 20 | 220 | 1.00000 | 1.00000 | 0.00000 | 94% |
| c4 | LOO | 20 | 220 | 0.99801 | 0.99935 | **0.00134** | 46% |
| gstm1 | LZO | 19 | 376 | 0.99739 | 1.00000 | 0.00261 | 91% |
| gstm1 | LOO | 19 | 376 | 0.98132 | 0.99964 | **0.01832** | 43% |
| acot | LZO | 20 | 380 | 0.99972 | 1.00000 | 0.00028 | 93% |
| acot | LOO | 20 | 380 | 0.99721 | 0.99929 | **0.00208** | 55% |
| ankrd36c | LZO | 18 | 378 | 0.99698 | 1.00000 | 0.00302 | 80% |
| ankrd36c | LOO | 18 | 378 | 0.99508 | 0.99996 | **0.00487** | 50% |
| lpa | LZO | 19 | 475 | 0.98210 | 1.00000 | **0.01790** | 92% |
| lpa | LOO | 19 | 475 | 0.97077 | 0.99941 | **0.02864** | 50% |

## 3. Result: simulated reads, same individuals

| locus | regime | n | blocks | identity | panel best | gap | exact |
|---|---|---:|---:|---:|---:|---:|---:|
| cyp2d6 | LZO | 20 | 380 | 0.99978 | 1.00000 | 0.00022 | 99% |
| cyp2d6 | LOO | 20 | 380 | 0.97586 | 0.99815 | 0.02230 | 44% |
| c4 | LZO | 20 | 220 | 1.00000 | 1.00000 | 0.00000 | 100% |
| c4 | LOO | 20 | 220 | 0.99876 | 0.99935 | 0.00059 | 47% |
| gstm1 | LZO | 20 | 396 | 1.00000 | 1.00000 | 0.00000 | 100% |
| gstm1 | LOO | 20 | 396 | 0.98108 | 0.99962 | 0.01854 | 45% |
| acot | LZO | 20 | 380 | 0.99994 | 1.00000 | 0.00006 | 99% |
| acot | LOO | 20 | 380 | 0.99699 | 0.99929 | 0.00230 | 57% |
| ankrd36c | LZO | 20 | 420 | 1.00000 | 1.00000 | 0.00000 | 100% |
| ankrd36c | LOO | 20 | 420 | 0.99790 | 0.99996 | 0.00206 | 60% |
| lpa | LZO | 20 | 500 | 1.00000 | 1.00000 | **0.00000** | 100% |
| lpa | LOO | 20 | 500 | 0.98947 | 0.99944 | **0.00996** | 56% |

## 4. Paired difference: what real reads cost

| locus | LZO real | LZO sim | read cost | LOO real | LOO sim | read cost |
|---|---:|---:|---:|---:|---:|---:|
| cyp2d6 | 0.00153 | 0.00022 | +0.0013 | 0.02109 | 0.02230 | −0.0012 |
| c4 | 0.00000 | 0.00000 | 0.0000 | 0.00134 | 0.00059 | +0.0008 |
| gstm1 | 0.00261 | 0.00000 | +0.0026 | 0.01832 | 0.01854 | −0.0002 |
| acot | 0.00028 | 0.00006 | +0.0002 | 0.00208 | 0.00230 | −0.0002 |
| ankrd36c | 0.00302 | 0.00000 | +0.0030 | 0.00487 | 0.00206 | +0.0028 |
| **lpa** | **0.01790** | **0.00000** | **+0.0179** | **0.02864** | **0.00996** | **+0.0187** |

### Three conclusions

**(a) The implementation is sound.** Simulated LZO is 0.00000 at four of six loci and <= 0.00022 at
the other two, with 99–100% of blocks exact. With perfect reads and the answer in the panel, the
model finds it. This excludes a broad class of suspected defects — indexing, chaining, pruning,
emission scale — as explanations for anything below.

**(b) LPA's LZO anomaly is entirely read acquisition.** Real 0.01790, simulated **0.00000**. This was
previously flagged as possibly the highest-value bug in the module; it is not a model defect at all.
LPA's real LOO gap of 0.02864 decomposes as **0.0187 read acquisition (65%) + 0.0100 model (35%)**.
The plausible mechanism is that reads from the most divergent KIV-2 repeat copies are exactly those
that fail to map to GRCh38, so region-restricted extraction loses them preferentially. Unmapped-read
rescue (`kfilt`-style) therefore becomes the largest single lever for LPA.

**(c) Read acquisition is negligible everywhere else.** |read cost| <= 0.003 at the other five loci,
and negative (simulated marginally worse, i.e. noise) at three of them. This is an LPA-specific
effect driven by that locus's repeat structure, not a general property of real data.

### Ranked true model gap (simulated LOO, free of read effects)

| | gap |
|---|---:|
| cyp2d6 | **0.0223** |
| gstm1 | **0.0185** |
| lpa | **0.0100** |
| acot | 0.0023 |
| ankrd36c | 0.0021 |
| c4 | 0.0006 |

Three loci are effectively solved. cyp2d6 and gstm1 are now the largest genuine targets, and neither
is affected by read quality.

---

## 5. Pre-registered question, closed

`experiments/ARRAY_ESTIMATOR_PREREGISTRATION.md` froze this before the cohort was run: replace
`called_bp` with the continuous `mass_bp` as the reported array copy number only if, per block class,
`mass_bp` is strictly closer in **B > 0.60** of observations AND its mean error is lower. The declared
prediction was that it would fail everywhere.

**Outcome: failed, decisively.** Over 476 array-block observations across 25 array blocks and six
loci, B = 0.00 at 22 blocks; the maximum is 0.37 at LPA block 13 (KIV-2). Mean absolute length error:
`called_bp` **412 bp** against `mass_bp` **7,981 bp** — pair selection is 19x better.

`called_bp` remains the copy-number answer; `mass_bp` stays a diagnostic. This retires the
estimator-selection line of work.

### Why an earlier analysis pointed the other way

On GPT's simulated chimera pairs 16–23 at LPA KIV-2 only, `mass_bp` reached 0.42 repeat units while
pair selection exceeded 1.5, which suggested the swap. That comparison was doubly selective: one
block, and chimeric diploids (two haplotypes from unrelated donors) where pair selection at KIV-2
fails much harder than it does on real individuals. Across all array blocks and real individuals the
ordering reverses.

---

## 6. Corrections to earlier claims, recorded

1. **"`dbp` is a sequence distance."** It is a length difference. Two conclusions built on it —
   "the LPA panel ceiling is 80 bp so the panel represents everyone" and "CYP2D6 is essentially
   solved" — were wrong. On sequence identity both loci have real headroom.
2. **"Chimeras are harder than real individuals."** From n=1. At n=8 per arm the self-paired
   synthetic run lost MORE (94,982 bp) than the chimera control (78,300 bp).
3. **"Simulation is much harder than reality."** Also from n=1, also false: real mean 7,735 bp/sample
   against synthetic 11,873, overlapping ranges.
4. **"Report `mass_bp` instead of `called_bp` at arrays."** Refuted by section 5.
5. **"LPA's LZO gap is a model bug and the highest-value target."** Refuted by section 4b — it is
   read acquisition.

Two earlier GPT results were also re-scored against their own pre-registrations: the mapped-coverage
Poisson method passes its stated criteria only because of one sample (removing pair 13 turns a 24%
improvement into 18% worse), and the inverse-mean weighting fails all three criteria on its
validation set.

---

## 7. Known limits of this cohort

- **Path-holdout, not graph rebuild.** The graph, bubbles and block chain were built using the
  held-out individual. Measured leakage for one individual at LPA: 6 private nodes (6 bp of 198.8 kb,
  0.003%) and 14 private oriented adjacencies (0.17%), with median 462 of 464 retained haplotypes
  supporting each adjacency. Believed small, verified once.
- `best_identity` is shortlisted by top-16 Jaccard, so all gaps are lower bounds.
- One ancestry-skewed cohort (1000G individuals that are also HPRC assemblies), 30x short reads,
  one simulator (`wgsim`, uniform error, no GC or fragment bias).
- Simulated reads come from the same graph the panel is built from, so they cannot expose reference
  or assembly error — only mapping and coverage effects.
- 18–20 individuals per locus; per-block estimates on 18–20 observations.

---

## 8. Questions for review

1. Does the LZO/LOO decomposition support the claim in 4(a) that implementation defects are excluded,
   or is simulated LZO too easy a test to carry that weight?
2. Is the read-acquisition interpretation in 4(b) the best explanation of an LPA-only +0.0179, or
   should mismapping within the region (rather than unmapped reads) be separated first — and what
   measurement would separate them?
3. cyp2d6 (0.0223) and gstm1 (0.0185) now carry the largest true gaps, and neither is a tandem array
   of KIV-2's kind. Is there a reason to expect a common cause, or should they be attacked separately?
4. Given the exact-block fraction halves under LOO (80–94% -> 43–57%) while identity moves 0.2–2
   points, is `gap` the right headline, or is there a better single statistic?
5. The estimator-selection line is retired (section 5) and the dense-marker channel was falsified
   earlier. What remains as the strongest candidate for the cyp2d6/gstm1 gap?
