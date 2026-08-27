# panvar `genotype`: experiment summary and current status

Date: 2026-08-25. Self-contained for external review. Supersedes the results sections of
`COHORT_REPORT_2026-08-23.md`, which remains valid for its methods description.

---

## 0. One-paragraph status

Internal self-consistency is strong: with the answer in the panel and reads simulated from it, the
model retrieves it (simulated LZO gap 0.00000 at four of six loci, 99–100% exact). That does NOT
validate off-panel discrimination, where a third of blocks fail even when the truth remains
representable. Three of six real loci are effectively solved. The residual splits into **read
acquisition at LPA** (~64% of its loss, paired, positive in 19/19 samples) and **evidence-poor simple
blocks at cyp2d6/gstm1**. Six interventions were pre-declared and gated; one passed, worth 0.5–4% of
the gap it targets. The evidence-poor blocks are best described as approaching an identifiability
limit OF THE CURRENT globally region-unique, marginal syncmer-count representation — not as a
short-read information limit.

---

## 1. Experimental design

Six real loci, twenty held-out individuals each, two regimes from identical reads:

- **LZO** (leave-zero-out): the individual's two assembly haplotypes stay in the panel. An exact call
  is reachable, so this is an implementation ceiling.
- **LOO**: both excluded. Generalisation.

Run twice, over the same individuals:

- **real**: 1000G 30x CRAMs, region pulled remotely (`samtools view <url> <region>`), locus +/- 20 kb
  from each graph's own PanSN GRCh38 path name. Region-mapped reads only, no unmapped rescue.
- **simulated**: `wgsim` 30x directly from the two haplotype sequences. No mapping step, so no
  unmapped or mismapped fraction.

234 sample-locus runs, 8,934 scored blocks. Per-block table:
`experiments/cohort_2026-08-23/all_blocks.tsv` (36 columns, haplotype names not indices).

### Metric

**`gap` = `best_identity` − `identity`**, per block. `best_identity` is the best any panel allele
reaches, over a top-16 syncmer-Jaccard shortlist, so every gap quoted is a LOWER bound on true
headroom.

`dbp` is a LENGTH difference (`|len(called) − len(truth)|`), not a sequence distance. An earlier
analysis conflated the two and drew two wrong conclusions; see section 5.

---

## 2. Where the module stands

### Real reads

| locus | LZO gap | LOO gap |
|---|---:|---:|
| c4 | 0.00000 | **0.00134** |
| acot | 0.00028 | **0.00208** |
| ankrd36c | 0.00302 | **0.00487** |
| gstm1 | 0.00261 | **0.01832** |
| cyp2d6 | 0.00153 | **0.02109** |
| lpa | **0.01790** | **0.02864** |

### Simulated reads, same individuals

| locus | LZO gap | LOO gap |
|---|---:|---:|
| c4 | 0.00000 | 0.00059 |
| acot | 0.00006 | 0.00230 |
| ankrd36c | 0.00000 | 0.00206 |
| gstm1 | 0.00000 | 0.01854 |
| cyp2d6 | 0.00022 | 0.02230 |
| lpa | **0.00000** | **0.00996** |

**Simulated LZO establishes strong internal SELF-CONSISTENCY, not overall soundness.** It is 0.00000
at four of six loci, <= 0.00022 at the other two, 99–100% of blocks exact. But it only shows that when
the exact allele is in the panel and the reads were simulated from it, the system retrieves it. It does
NOT validate off-panel discrimination, which is where the module actually operates. Under simulated
LOO, among observations whose truth genotype IS still exactly representable:

| locus | representable observations | exact |
|---|---:|---:|
| c4 | 115 | 89.6% |
| acot | 263 | 82.9% |
| ankrd36c | 307 | 81.4% |
| lpa | 395 | 71.1% |
| cyp2d6 | 251 | 67.3% |
| gstm1 | 288 | 62.2% |

So selection fails on a third of blocks whose answer is present. Candidate pruning, emission
specification and linkage are NOT excluded by LZO.

**LPA's LZO anomaly is read acquisition, not a model defect.** Real 0.01790 against simulated 0.00000.
Paired over the 19 samples with both arms: LZO penalty +0.01790 and LOO penalty +0.01820, **positive in
19/19 samples**, leaving a model component of ~0.01043 — about **64% acquisition, 36% model**.

**Calling it 'unmapped reads' is NOT yet supported.** Reads are counted from sequence without using
their mapping positions, so a read misplaced onto the wrong KIV-2 copy but still inside the extraction
interval is still available to the genotyper. The loss could be reads that are unmapped, placed
outside the interval, dropped by mate handling, or affected by real-library coverage bias. These imply
different fixes and the experiment in section 9B is designed to separate them.

Read acquisition is negligible at the other five loci (|cost| <= 0.003, negative at three).

**Independent confirmation** from a different statistic: the fraction of blocks where the model copies
from the individual's own haplotypes under LZO is 97.3% simulated against 79.8% real. Under LOO it is
0.0% in both, which verifies `--exclude-haplotypes` across all 234 runs.

---

## 3. Where the remaining gap lives (M1)

Simulated LOO, so free of read effects.

**Two distinct failure modes, and the loci sort into them:**

| locus | loss concentrated in | array-block gap | simple-block gap |
|---|---|---:|---:|
| cyp2d6 | **simple blocks (98% of bp)** | 0.00003 | **0.02647** |
| gstm1 | **simple blocks (81%)** | 0.00007 | **0.02062** |
| lpa, acot, ankrd36c, c4 | array blocks (95–99%) | — | 0.0001–0.0024 |

cyp2d6's and gstm1's array blocks are essentially perfect — roughly 1000x smaller gaps than their
simple blocks. **The tandem-array line of work is irrelevant to the two largest remaining gaps.**

**The predictor is a marker-density ratio with a sharp threshold near 2.** Naming matters here: the
`n_markers` column is the count of DISTINCT MARKER SLOTS over the whole block (the union across its
alleles), so `n_markers / n_alleles` is block-wide marker slots per candidate allele — NOT the mean
number of markers that discriminate each allele from its neighbours. It mixes evidence density with
block complexity. It is a good empirical warning signal and a poor mechanistic quantity. The rows below
are simple (non-array) blocks only, and are repeated sample-block observations over roughly 8–51
UNIQUE blocks per bin, so they are correlated:

| markers/allele | n | mean gap | blocks losing > 0.01 |
|---|---:|---:|---:|
| 0 – 0.5 | 236 | **0.04483** | **47%** |
| 0.5 – 2 | 240 | **0.03882** | **35%** |
| 2 – 5 | 160 | 0.00119 | 1% |
| 5 – 20 | 982 | 0.00107 | 0% |
| 20 – 100 | 178 | 0.00002 | 0% |

A 30–40x drop across the threshold. The share of marker-poor blocks predicts the locus ranking almost
monotonically: gstm1 44%, cyp2d6 38%, lpa 29%, acot 23%, ankrd36c 7%, **c4 0% with the smallest gap**.

Worst individual blocks: cyp2d6 block 5 chooses among **14 alleles on 10 markers** (gap 0.138);
gstm1 blocks 9 and 11 among **16 and 21 alleles on ZERO markers** (gaps 0.058, 0.085).

At marker-poor simple blocks under simulated LOO: **94.3% still have an exactly representable truth
genotype**, mean `best_identity` is **0.99896**, and exact recovery is only **55.5%**. The panel
usually contains the answer and the answer is clearly distinct from what is chosen.

**This does not establish a short-read information limit.** It shows these blocks approach an
identifiability limit of the CURRENT representation: globally region-unique markers, pooled as
marginal counts. The dense-marker failure (4a) falsifies adding more markers under that same global
filtering and pooling; it does not falsify recruiting fragments locally first and then using ambiguous
internal markers conditionally. See section 9A.

---

## 4. Interventions, all pre-declared and gated

**Gate, fixed in advance for every one:** blocks with markers/allele < 2 must improve AND blocks at or
above the threshold must not regress. Reported per stratum, paired per block, never as an aggregate.

### 4a. Evidence SUPPLY — dense markers (`--all-kmers`). FALSIFIED.

| locus | stratum | n | gap before → after | better / worse |
|---|---|---:|---|---|
| cyp2d6 | < 2 | 120 | 0.06217 → 0.05853 | 11 / 9 |
| cyp2d6 | >= 2 | 200 | 0.00505 → 0.00506 | 7 / **14** |
| gstm1 | < 2 | 39 | 0.04249 → **0.04731** | 0 / 3 |
| gstm1 | >= 2 | 50 | 0.00016 → **0.00024** | 2 / 7 |

Poor blocks move by a coin flip at one locus and worsen at the other; well-supplied blocks regress.
Dense k-mers pass through the same filters and land in the same place. **A two-pass per-block refactor
was planned and never written, because the premise failed first.**

### 4b. Evidence RETENTION — drop region-uniqueness (`--no-region-unique`). FALSIFIED.

| locus | regime | n | gap before → after | samples worse |
|---|---|---:|---|---|
| cyp2d6 | LZO | 13 | **0.00000 → 0.10227** | **13/13** |
| cyp2d6 | LOO | 13 | 0.01428 → **0.10585** | **13/13** |

A 7x worse LOO gap, and it **breaks LZO**, which was exactly 0. A setting that fails when the answer
is provably in the panel is disqualified on logic, so gstm1 was not run. Mechanism (already in the
ledger): a syncmer occurring in several places accumulates counts from all of them, so its count stops
reflecting the block being genotyped. **Region uniqueness is doing necessary work; the markers it
removes are genuinely ambiguous.**

### 4c. Adjacency evidence (`--edge-weight`, NEW). PASSES, small.

Adjacency multiplicities are computed for every allele and were **never read by the emission** — a
whole evidence class carried and discarded. Now scored with the same negative binomial against the
same depth, discounted by the same `rho` as nodes, behind its own weight (default 0).

All four arms complete, 20 samples each:

| weight | locus | markers/allele < 2 | >= 2 |
|---|---|---|---|
| **0.25** | cyp2d6 | 0.06217 → 0.06185 (5b/2w) | 0.00505 → 0.00504 (**12b/0w**) |
| **0.25** | gstm1 | 0.04684 → **0.04497** (6b/2w) | 0.00016 → 0.00016 (6b/2w) |
| **1.0** | cyp2d6 | 0.06217 → **0.06961** (6b/**11w**) | 0.00505 → 0.00504 (21b/3w) |
| **1.0** | gstm1 | 0.04684 → **0.04174** (**13b/2w**) | 0.00016 → 0.00016 (17b/8w) |

**Weight 0.25 passes the gate at both loci**, gaining 0.5% (cyp2d6) and 4% (gstm1) of the targeted gap.

**Weight 1.0 splits the loci, and this is the important result.** cyp2d6 gets 12% WORSE (6 better
against 11 worse) while gstm1 gets 11% BETTER (13 better against 2 worse) — its best result of any
arm. An earlier reading of the cyp2d6 arm alone concluded that full weight fails and that this
confirmed the ledger's double-counting mechanism. **The full data does not support that**: the ledger's
prior holds at cyp2d6 and is reversed at gstm1.

**There is therefore no single safe weight.** The optimum is locus-dependent, and picking it per locus
from the data is exactly the self-selecting tuning that design constraint C-1 forbids — the same shape
as the per-block rule selection that picked the permissive set in 23/23 blocks because a noisy marker
set is confidently wrong.

Conservative reading: **0.25 is the only weight that is safe everywhere tested**, and it buys 0.5–4% of
a 0.019–0.022 gap. Recommendation unchanged: leave the default at 0, document 0.25 as tested-safe on
two loci, and do not adopt a higher weight without a mechanism that explains why gstm1 tolerates it and
cyp2d6 does not.

### 4d. Earlier, same gate

- **Dense markers + bounded outlier** (8 held-out pairs): fixed the motivating pair (49,901 → 5,588 bp)
  and broke another 4x (5,550 → 22,185), six unchanged. Excluding the motivating pair, 60% worse.
- **Array copy number from `mass_bp` instead of `called_bp`** (pre-registered, 476 array observations,
  25 blocks, six loci): B = 0.00 at 22 blocks, max 0.37 against a 0.60 threshold. Mean absolute length
  error `called_bp` **412 bp** vs `mass_bp` **7,981 bp**. Declared prediction held.
- **GPT's mapped coverage** (BWA → gfainject → gafpack), re-scored against its own pre-registrations:
  the Poisson arm passes only because of one sample (removing pair 13 turns a 24% improvement into 18%
  worse); the inverse-mean arm fails all three criteria on its validation set.

---

## 5. Claims retracted during this work

1. **"`dbp` is a sequence distance."** It is a length difference. Two conclusions built on it — "the
   LPA panel ceiling is 80 bp so the panel represents everyone" and "CYP2D6 is essentially solved" —
   were wrong.
2. **"Chimeras are harder than real individuals."** From n=1. At n=8 per arm, self-paired lost MORE
   (94,982 bp) than the chimera control (78,300 bp).
3. **"Simulation is much harder than reality."** Also n=1, also false: real 7,735 bp/sample against
   simulated 11,873, overlapping ranges.
4. **"Report `mass_bp` instead of `called_bp` at arrays."** Refuted by the pre-registered test.
5. **"LPA's LZO gap is a model bug and the highest-value target."** It is read acquisition.
6. **"The region-uniqueness filter is too aggressive."** It is necessary; removing it breaks LZO.

The recurring error is a comparison that is selective in a way not noticed at the time — one block,
one sample, one metric. Pre-registration with a declared prediction and per-stratum reporting caught
the later ones.

---

## 6. Current priorities

1. **R2: unmapped-read rescue at LPA** (`kfilt`-style). Quantified at 65% of LPA's loss. Input
   plumbing, not modelling. The only read-side item worth doing.
2. **cyp2d6 / gstm1 marker-poor blocks.** No intervention has moved them materially. Present belief:
   an information limit under the current marker paradigm.
3. **Surface the indicator.** `markers/allele < 2` predicts a 30–40x higher gap, and `evidence=linked`
   marks blocks called from neighbours with no local evidence (152 of 2,175 blocks). Neither is
   documented; `docs/modules/genotype.md` and `docs/algorithms/genotype.md` have never been written.

---

## 7. Known limits

- **Path-holdout, not graph rebuild.** Measured leakage for one individual at LPA: 6 private nodes
  (0.003% of 198.8 kb) and 14 private oriented adjacencies (0.17%), median 462 of 464 retained
  haplotypes supporting each adjacency. Believed small, verified once.
- `best_identity` is shortlisted by top-16 Jaccard, so all gaps are lower bounds.
- One ancestry-skewed cohort (1000G individuals that are also HPRC assemblies), 30x short reads, one
  simulator (`wgsim`: uniform error, no GC or fragment bias).
- Simulated reads come from the same graph the panel is built from, so they cannot expose reference or
  assembly error — only mapping and coverage effects.
- 18–20 individuals per locus; per-block estimates rest on 18–20 observations.
- Six real samples failed download, so real is 114/120 sample-loci against simulated 120/120; paired
  real-vs-simulated figures use each arm's own mean rather than strict per-individual pairing.

---

## 8. Questions for review

1. Is simulated LZO too easy a test to carry the "implementation is sound" claim? It excludes a broad
   class of defects in one stroke, and simulated reads cannot expose reference or assembly error.
2. LPA's +0.0179 is attributed to unmapped reads, but reads that mapped to the WRONG repeat copy would
   look identical in this measurement. What would separate mismapping from non-mapping?
3. Marker-poor blocks have `best_identity` 0.99896 and a 0.0418 gap: the right allele is present and
   clearly different, yet neither more markers nor keeping more markers helps. Is "information limit"
   the right conclusion, or is there a class of evidence not yet considered?
4. Adjacency at weight 1.0 helps gstm1 by 11% and hurts cyp2d6 by 12%, while 0.25 is safe at both and
   worth only 0.5–4%. Is there a principled way to set that weight, or does a locus-dependent optimum
   mean the term should stay off by default?
5. Given three of six loci are solved and the rest split into a plumbing problem and an information
   limit, what would constitute "done" for this module?

---

## 9. Next experiments, as proposed in review

### 9A. Marker-poor oracle — tests the most questionable conclusion here

For every failing marker-poor block, score all candidate pairs on the NOISELESS expected
multiplicities of the current markers:

- truth wins -> sampling, calibration or linkage is the problem;
- truth loses -> the emission formulation is the problem;
- truth ties -> the current representation is genuinely non-identifying.

Then repeat with exact paired-FRAGMENT signatures generated from the alleles at the observed insert
distribution. If fragments identify the truth where marginal markers do not, the indicated build is a
two-stage alignment-free fallback: recruit fragments to a bubble using unique flank or junction
syncmers, then inside that recruited set use the ambiguous internal markers plus marker-pair/order
evidence, score complete candidate haplotypes, and return an equivalence set or no-call when the
evidence remains non-identifying.

This is the same shape as the Gate 2 oracle that localised the array failure, and it decides whether
"identifiability limit" is a property of the data or of the representation.

### 9B. LPA acquisition funnel — separates loss from misplacement

Using truth-labelled SIMULATED LPA reads run through the real mapping/extraction pipeline, genotype at
five stages: original reads; after mapping; region-extracted only; all mapped reads recovered by name
plus mates; unmapped added back. Stratify reads by correct-copy, wrong-copy-inside-region,
outside-region, and unmapped.

This distinguishes reads lost from reads misplaced, and says whether the fix is extraction, mate and
unmapped rescue, or the emission.

## 10. Status of `--edge-weight`

Experimental, not finished. Present state: builds, validated to reject negative/NaN/infinite values,
and refuses to combine with `--compositional` (the compositional emission is a multinomial shape score
and the adjacency term a negative binomial over counts -- and the edge term reached only the
non-compositional branch, so the combination would previously have been silently ignored).

Still missing before the option should be retained: a registered test proving the flag changes the
scorer, that weight 0 reproduces current output byte-for-byte, that indexed and direct routes agree,
and documentation in `--help`. The favourable gate covers only simulated cyp2d6/gstm1 and rests on
correlated sample-block observations. Default remains 0.

## 11. Reproducibility gap

`experiments/` is gitignored, so the intervention tables, commands and aggregation scripts behind this
report are NOT durably versioned. Before this is treated as the authoritative record it needs: commit
and binary hashes, exact commands and options, the sample list and the six failed downloads, paired
baseline/intervention tables, the aggregation scripts, and the completed edge-weight results.

---

## 12. 9A RESULT (2026-08-25): the information is present and MARKER SELECTION discards it

**Revised 2026-08-25 after review.** The first version of this section attributed the whole loss to
confinement, on the strength of an audit column that has since been found to be vacuous. The central
claim survives and is now measured directly; the attribution does not. Corrections are marked below.

Run on the worst marker-poor blocks at cyp2d6 and gstm1, LOO, using `--dump-block` so the panel
matrix analysed is the one the model actually used, and `--ledger-block` (added for this) so the
filter decisions are read out rather than inferred.

### Two distinct causes, and neither is a short-read limit

**(a) Alleles shorter than k carry no marker at all.** Across all six loci, every allele under 31 bp
has zero markers; everything longer has 0-1%:

| allele length | n | zero markers |
|---|---:|---:|
| **< 31 bp (shorter than k)** | 62 | **100%** |
| 31-99 | 390 | 1% |
| 100-999 | 315 | 1% |
| 1000+ | 2,194 | 0% |

A 21 bp allele cannot contain a 31-mer, so no syncmer exists to keep. `nearest_sibling_jaccard` reads
1.0000 because `jaccard()` returns 1.0 for two empty sets. This is a parameterisation mismatch, not a
filter and not an information limit.

> **CORRECTION 1.** The first version added "`n_nodes_lost_region` is 0 for all of them: nothing was
> filtered". That number was meaningless. The bubble-audit path wrote
> `if (mult > options.max_multiplicity) continue`, and the option's contract is that **0 means no
> cap**, so at the default every candidate was discarded before region filtering ran and the loss
> column read 0 whatever the panel contained. Fixed, with a regression test that fails against the
> old code (0 markers retained, against 433 after the fix) and a new `n_retained_nodes` column --
> the audit previously reported only raw-inventory counts, so an emptied retained set left no trace
> in the table at all. Production genotyping was never affected: `build_block_marker_panel` guards
> the cap correctly. **The conclusion for case (a) is unchanged** -- it rests on allele length
> against k, which the bug does not touch.

**(b) Where alleles ARE long enough, selection is what loses the signal.** cyp2d6 block 5, 14 alleles
of 261-851 bp, 105 possible diploid pairs. Marker sets computed from the dumped allele sequences with
an independent syncmer implementation; equivalence classes computed from the production ledger.

| marker set | markers | distinct pair-vectors | pairs in non-singleton classes | largest class |
|---|---:|---:|---:|---:|
| all canonical 31-mers | 215 | 103 / 105 | 3 | 3 |
| closed syncmers (s=11, the default) | 36 | **103 / 105** | 3 | 3 |
| **what the model retains** | **10** | **10 / 105** | **102** | **35** |

**Syncmer sampling costs no noiseless pair-vector identifiability on this block** -- 36 syncmers
separate these alleles into exactly the same classes as all 215 k-mers. The model keeps 10 of those
36, and 102 of the 105 pairs then collapse into classes no scorer can resolve.

That claim is deliberately narrow. 36 observations are not 215 observations at finite depth, so
syncmers may still cost robustness to read noise even where they cost no identifiability. Nothing
here measures that.

> **CORRECTION 2.** The syncmer row was labelled s=20. The production default is
> `min(11, (k+2)/3)` = **11** at k=31 ([syncmer.cpp:27](../src/syncmer.cpp#L27)); the 36 informative
> markers is the s=11 count. Review reports s=20 gives 55 informative markers and the same 103/105
> classes, so the conclusion is insensitive to s -- but the row was mislabelled.
>
> **CORRECTION 3.** "Colliding pairs: 35" conflated two quantities. 35 is the LARGEST equivalence
> class, not the number of colliding pairs; there are seven ambiguous classes holding 102 of the 105
> pairs. All three numbers are now reported.

### The failure is not a scoring failure

HG00096's truth pair at this block is (2,4); production called (4,7).

| marker set | pairs sharing the truth's vector | is the called pair among them? |
|---|---:|---|
| the 10 retained | **35** | **yes** |
| retained + over-expected (20) | 5 | no |
| all 36 informative | 1 (unique) | no |

Under the markers it was given, the model could not have distinguished the truth from the pair it
chose. This is an evidence failure, not an emission failure, and it is why every reweighting tried in
section 4 moved nothing.

### WHICH filter: measured, not inferred

Two independent rules drop a marker
([genotype_markers.cpp:992](../src/genotype_markers.cpp#L992)): it varies in more than one block
(confinement), or the panel shows it more often than the blocks account for (over-expected). A marker
can fail both. The ledger records each candidate's fate before the erase. cyp2d6 block 5, identical
across 8 samples x both regimes (one sample differs in LOO, 33 informative / 7 retained):

| fate | distinct markers |
|---|---:|
| retained | **10** |
| fails confinement ONLY | **0** |
| fails over-expected ONLY | 10 |
| fails BOTH | 16 |

> **CORRECTION 4.** The first version concluded "the remaining filter is CONFINEMENT". It is not the
> only one, and on its own it is not the operative one. Relaxing confinement alone would rescue
> **zero** markers here, because every multi-block marker at this block also fails the over-expected
> test. This was attributed by elimination from the vacuous audit column above; it is now instrumented.

What each relaxation would actually buy, on the same 105 pairs:

| markers restored | markers | distinct vectors | pairs ambiguous | largest class |
|---|---:|---:|---:|---:|
| none (production) | 10 | 10 / 105 | 102 | 35 |
| confinement relaxed alone | 10 | 10 / 105 | 102 | 35 |
| over-expected relaxed alone | 20 | 28 / 105 | 90 | 20 |
| both relaxed (all informative) | 36 | 103 / 105 | 3 | 3 |

**Neither single relaxation is sufficient.** The 16 markers that fail both rules carry most of the
separating power. Locus-wide the picture agrees: the existing confinement-by-context audit reports
98-99% of varying markers already confined to one block at cyp2d6 and gstm1, so confinement was never
where the volume of loss was.

> **CORRECTION 5.** The first version said `--no-region-unique` "tests a different filter". False.
> That flag disables the entire guarded block
> ([genotype_markers.cpp:826](../src/genotype_markers.cpp#L826)), confinement included, so M2b did
> test it. What M2b falsified is narrower than "retention does not help": restoring every ambiguous
> marker globally, with no model of where the extra counts came from, breaks LZO. It did not falsify
> selective restoration, longer context, or joint scoring of the blocks involved. A further defect:
> with the flag off the audit reports 0/0 confinement because the instrumentation is skipped, which
> should print as disabled/NA.

### Across blocks: the two loci fail for DIFFERENT reasons

Ledger run on the worst blocks of both loci plus a marker-rich control, and on the block-marker audit
for allele lengths. Identical across every sample and both regimes unless noted.

| block | kind | alleles | allele bp | informative | retained | multi only | over only | both |
|---|---|---:|---|---:|---:|---:|---:|---:|
| cyp2d6 5 | bubble | 14 | 261-851 | 36 | 10 | 0 | 10 | 16 |
| cyp2d6 3 | bubble | 14 | — | 23 | 7 | 0 | 0 | 16 |
| cyp2d6 7 | — | — | — | 1,329 | **1,253** | 25 | 30 | 21 |
| gstm1 9 | backbone | 16 | **14-42** | 7 | **0** | 0 | 0 | 7 |
| gstm1 11 | backbone | 21 | **5-21** | **0** | **0** | — | — | — |

Two separate readings:

**Filtering is cheap where markers are plentiful and total where they are scarce.** The marker-rich
control keeps 94% and its losses split evenly between the two rules. At the poor blocks the dropped
markers almost all fail BOTH rules, so no single relaxation reaches them — 0 markers at cyp2d6 5 and
gstm1 9 would be rescued by relaxing confinement alone.

**gstm1's worst blocks are dominated by case (a) — but only block 11 is purely case (a).** Block 11's
alleles are 5-21 bp, every one shorter than k, including an empty bypass allele, so it has zero
informative markers of EITHER kind before any filter runs (`n_informative_nodes` and
`n_informative_edges` both 0, all 21 alleles unseparable). That one is purely a parameterisation
failure. **Block 9 is BOTH**: 11 of its 16 alleles are under k, but seven informative markers do exist
and all seven are then removed for failing both filters. Calling gstm1's worst blocks "not a selection
problem" was too strong — at block 9 the two mechanisms are stacked, and fixing either alone leaves
the other binding.

**cyp2d6 block 5, by contrast, has 261-851 bp alleles and is purely a selection failure.** The two
largest remaining locus gaps still have different dominant causes, and section 4's single "marker
scarcity" framing conflated them — but the clean dichotomy is only clean at the extremes.

Note also that both gstm1 blocks are BACKBONE blocks carrying 16 and 21 alleles. A backbone is meant
to be the invariant chain between bubbles; one with 21 distinct alleles, all under k, suggests the
block partition is putting genuine variation where the marker model cannot see it. Worth a look
independently of markers.

### Correction to section 3

The conclusion recorded there -- that these blocks approach an identifiability limit of the current
representation -- is **withdrawn**. The representation is adequate: syncmers retain full separating
power at 36 markers. What fails is marker SELECTION on top of it.

This also explains 4a: `--all-kmers` enlarges the candidate pool, but both filters still apply to the
enlarged set, so density alone cannot help.

> **CORRECTION 6.** The first version explained `--edge-weight`'s consistent sign by saying adjacency
> markers "span junctions and are therefore reachable where internal k-mers of a short allele are
> not". They are not. Adjacencies are built between consecutive syncmers **inside** `allele_seq`
> ([genotype_markers.cpp:58](../src/genotype_markers.cpp#L58)), so an allele with no syncmers has no
> adjacency either. Whatever the edge term is worth, it is not rescuing sub-k alleles. The
> boundary-spanning construction below is still the right idea -- it just does not exist yet.

### Adjacency evidence at cyp2d6 block 5: FALSIFIED before implementation

The proposed next architectural step was longer internal context — score adjacencies, which are a
different evidence class and might escape the filters that remove single syncmers. The ledger now
records edge candidates alongside nodes, so this is answerable without writing an emission.

| marker set | markers | distinct pair-vectors | pairs ambiguous | largest class |
|---|---:|---:|---:|---:|
| retained nodes (production) | 10 | 10 / 105 | 102 | 35 |
| retained edges | 11 | **10 / 105** | 102 | 35 |
| retained nodes + edges | 21 | **10 / 105** | 102 | 35 |
| ALL informative nodes | 36 | 103 / 105 | 3 | 3 |
| ALL informative edges | 41 | 103 / 105 | 3 | 3 |
| ALL informative nodes + edges | 77 | 103 / 105 | 3 | 3 |

**Adding every surviving adjacency to every surviving node changes nothing** — same 10 classes, same
35-pair largest class, and the truth pair still shares its vector with the pair production actually
called. Longer context does not escape the filters either: edges are dropped at the same rate as
nodes (11 of 41 retained against 10 of 36), with the same fate profile.

So the adjacency line is closed for this block, and it explains why `--edge-weight` measured at only
0.5-4%: the adjacencies that survive carry no information the surviving nodes do not already carry,
and the ones that would help are removed by the same two rules. **The binding constraint is the
filter, not the context length.**

### What this makes testable

1. **`--kmer-size` sweep: DONE, and it answers the question. See
   `KMER_SIZE_PREREGISTRATION.md`.** Pre-registered with the sub-k stratum declared from panel
   geometry (gstm1 blocks 3, 9, 10, 11), then run at k=21 and k=15.

   The gate's letter passes — sub-k mean gap 0.04533 -> 0.04149, control +0.00073 — but **1 of 30
   sub-k observations changed at all**, and blocks 9 and 11 (69% and 100% sub-k, the two the mechanism
   was identified in) are byte-identical. The declared prediction is falsified. The ledger says why:

   | block | k=31 informative / retained | k=21 | k=15 |
   |---:|---|---|---|
   | 3 | 11 / **0** | 16 / **0** | 30 / **0** |
   | 9 | 7 / **0** | 11 / **0** | 11 / **0** |
   | 10 | 19 / 4 | 19 / 3 | 15 / **0** |
   | 11 | 0 / 0 | 1 / **0** | 0 / 0 |

   Lowering k does create the candidates it was meant to — block 3 goes 11 -> 30 — and **every one is
   then filtered away.** At these blocks both mechanisms operate at once, and because filtering is
   downstream it is the one that binds. Case (a) and case (b) are not alternatives here; they are
   stacked.

2. **Boundary-spanning markers: the above is a design constraint on them, found before building.**
   A k-mer anchored in a flank and crossing into a short allele contains flank sequence by
   construction, and that sequence is shared with the neighbouring blocks' alleles. It therefore
   varies in more than one block and is rejected by confinement — already the dominant fate at these
   blocks (gstm1 block 3: 7 of 11 dropped for multi-block alone, at every k). **As specified, these
   markers would be filtered out on arrival.** They need to be attributed to a JUNCTION rather than to
   a block, or paired with a rule that expects a boundary-spanning marker in both blocks it spans.
   A second complication: at gstm1 block 11 the neighbouring blocks carry 11 and 48 alleles, so there
   is no invariant flank to anchor in and the junction k-mer depends on which neighbouring allele the
   haplotype carries.
2. **Instrument before relaxing.** The ledger now exists; extend it past one block before any filter
   is loosened. Any experiment on the two rules needs **separate** switches -- `--no-region-unique`
   moves both at once and cannot attribute a result.
3. **Over-expected is the rule to understand first**, since it is the only one that rescues anything
   alone. But its markers are exactly the contaminated ones M2b showed cannot simply be restored, so
   the useful form is one that keeps the discrimination without importing foreign counts: longer
   context that is confined where the single syncmer is not, or fragment-local recruitment.
4. Same stratified gate throughout: evidence-poor blocks improve AND the rest do not regress, LZO
   does not break, direct and indexed agree, on real and simulated cohorts.

### Caveat

Case (b) rests on **one block at one locus**, measured across 8 samples and both regimes. The 36 /
10 split and the equivalence classes are now read from the production ledger rather than inferred, and
the marker-set comparison is computed from allele sequences with an independent syncmer implementation
-- but the generality across blocks is not established. The gstm1 zero-marker blocks are running.
This section has now revised its own conclusion twice (information limit -> allele length -> marker
selection), each time on looking one layer deeper into the same data; treat the attribution as the
current best-instrumented reading, not as settled.
