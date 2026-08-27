# Finishing `genotype`: point-by-point to a clean release

Rewritten 2026-08-27, replacing the earlier draft, which was too centred on KIV-2 and on "honest
uncertainty". Status of each item: DONE, RUNNING, READY (unblocked, not started), BLOCKED (waiting on
something named), or CLOSED (tested and rejected).

**Headline.** The caller has good selective accuracy and poor coverage, with a residual of confident
errors. Simulated leave-zero-out has zero PASS errors above 0.01 at any locus; real leave-zero-out has
none above 0.05. But call rate is 48-76% off-panel and 57-92% on real leave-zero-out, and 1.7-7.9% of
simulated off-panel PASS calls and up to 12.2% of real ones exceed 0.01. The release target is the
**risk-coverage trade-off**, especially on near-representable off-panel blocks -- not accuracy alone,
and not KIV-2.

---

## 0. Freeze the release dashboard. Do this first.

One immutable baseline, generated before any further algorithm change, recording the exact commit,
commands, donor lists and denominators.

### Representability has three strata, not two

Defined from **certified sequence error**, not from whether the truth allele's label exists:

```
E*      = error of the certified best panel pair
Ecall   = error of the reported pair
excess  = Ecall - E*
```

| stratum | definition |
|---|---|
| exact | the true pair exists in the panel |
| **near-representable** | exact pair absent, but a panel pair reconstructs it closely: **E\* < 1% of aligned sequence** |
| panel-limited | even the certified nearest pair has substantial error |

The threshold is frozen here, before any intervention result. Continuous `E*` is reported alongside
it so no conclusion depends on the cut.

**Leave-one-out exists to test approximation when the exact allele is absent.** Panel-limited
observations are therefore NOT required to be declined -- that would reward silence and discard the
most interesting part of the experiment. What is required is that they are not confidently much worse
than the certified floor, and that panel limitation is exposed in the output.

### Per stratum, report together

Observations; `E*`; `Ecall`; excess; all-call accuracy; PASS-only accuracy; call rate; PASS calls
above 0.01 and above 0.05; exact-pair and equivalence-set accuracy; GQ reliability.

Stratify also by ordinary blocks, short-allele blocks, paralog/shared-marker blocks, and tandem
arrays. Define `leave-zero`, `leave-one` and `leave-all/off-panel` explicitly in the table.

**No accuracy figure is ever reported without call rate beside it.**

### Risk-coverage

Error as progressively lower-confidence calls are declined, and the coverage attainable while keeping
PASS error below 0.01 and below 0.05.

### Gates, frozen now

| stratum | gate |
|---|---|
| simulated LZO | >= 99% call rate, no PASS error above 0.01 |
| simulated off-panel, exact + near-representable | >= 90% call rate, mean excess < 0.005, < 1% PASS above 0.01 |
| real LZO | >= 85% call rate per locus, no PASS catastrophe above 0.05, mean PASS excess < 0.005 |
| real off-panel, exact + near-representable | >= 75% call rate per locus, < 5% PASS above 0.01 |
| KIV-2 short reads | certified answer inside the reported equivalence set, allocation flagged unresolved |

The off-panel call-rate gates apply to the exact and near-representable strata. **Overall off-panel
call rate stays visible** so conditioning cannot hide poor practical utility. The simulated LZO gate
is unconditioned: the truth is in the panel there, so widespread abstention is a real defect.

> Certified representability uses truth and is **evaluation-only**. Production must infer reliability
> from observable fit residuals, evidence counts and held-out calibration. This distinction has to
> survive into the code or someone will ship the oracle.

Where we stand against these today, from `results/real_data/*/genotype/*per_block_observations.tsv`:
simulated LZO call rate 90.5-100%; real LZO 57.4-91.8%; simulated off-panel 48.7-73.6%; real
off-panel 46.0-75.9%. Accuracy passes almost everywhere; coverage fails almost everywhere.

---

## 1. KIV-2 likelihood tuning -- CLOSED

Phase 1 ran to completion and every arm failed (`PHASE1_RESULT.md`). S1 362,913 certified excess; S2
473,054 (+30.3%); every S3 bound +102% to +146%. The 50-donor holdout was **not** opened and stays
sealed.

**Scoped precisely:** clump-mean normalization and symmetric clump clipping do not improve KIV-2
off-panel projection. No further tuning in that score family without a new mechanistic hypothesis --
not another tau, not a median, not a per-locus variant. This is not a proof that every likelihood
change is invalid.

---

## 2. CYP2D6 shared-marker retention -- READY. Start here.

The largest model-only gap of any locus: simulated off-panel mean 0.0223, against lpa's 0.0100.

Measured at block 5, identical across 8 samples and both regimes:

| marker set | markers | distinct pair-vectors | ambiguous pairs | largest class |
|---|---:|---:|---:|---:|
| production (retained) | 10 | 10 / 105 | 102 | 35 |
| relax confinement alone | 10 | 10 / 105 | **102 -- no change** | 35 |
| relax over-expected alone | 20 | 28 / 105 | 90 | 20 |
| **all informative** | **36** | **103 / 105** | **3** | 3 |

Two independent rules drop markers -- confinement (`blocks_with > 1`) and over-expected
(`actual > expected`). At this block **0 markers fail confinement only, 10 fail over-expected only,
and 16 fail both**, and those 16 carry most of the separating power. Relaxing either alone is
useless; the oracle must supply background for both.

Already closed here: adjacency evidence. Retained edges give the same 10 / 105, and edges are
filtered at the same rate as nodes (11 of 41 against 10 of 36).

Also closed: restoring every ambiguous marker globally. `--no-region-unique` broke leave-zero-out at
cyp2d6, 0.00000 -> 0.10227, worse in 13 of 13 samples. What that falsified is narrower than
"retention does not help": it falsified restoring everything **with no model of where the extra counts
came from**. Selective restoration and joint scoring were never tested.

### 2a. Truth-background oracle -- do this before any production change

On simulated data every other block's true state is known. For each restored shared marker:

```
observed count = contribution from the target block
               + contribution from other blocks and paralogs
               + background
```

Supply the known other-block contribution, then genotype blocks 3 and 5 with all informative markers.
One decisive question: **if shared-marker background were known perfectly, is the correct pair
recovered?** No means the filter diagnosis is incomplete. Yes means background handling is the
missing mechanism.

Ceiling to expect, not exceed: 3 of 105 pairs remain ambiguous even with all 36 markers.

### 2b. Replace truth with joint inference -- BLOCKED on 2a passing

A shared marker must not be assigned independently to one block. Score it once from its total
predicted multiplicity:

```
expected count = lambda * SUM over every block carrying the marker + background
```

which is a factor connecting several block states. Candidate implementations: joint or beam inference
across the involved blocks; a local factor over only the blocks where the marker occurs; or
marginalising over their possible contributions.

**The background must never be estimated from the caller's own calls and fed back as truth.** That is
constraint C-1, learned twice, and it risks a stable wrong fixed point.

### 2c. Gate

Blocks 3 and 5 materially improve; certified excess falls; representable-block call rate rises;
simulated LZO preserved; other cyp2d6 blocks not damaged; still beneficial on donors not used for
design; and no gain obtained by declining more calls.

---

## 3. GSTM1 fragment evidence -- READY after 2

Blocks 9 and 11 hold alleles of 5-42 bp. Internal 31-mers cannot represent them. `--kmer-size 21` was
tested and falsified: of 30 affected observations, exactly **one** changed.

A 150 bp read spans such an allele and both boundaries, so fragment evidence is potentially decisive
here -- and this is the locus where **half the calls are currently declined** (48.7% call rate).

This is **not** `--edge-weight 1`. That is co-occurrence of consecutive syncmers, a few to tens of bp,
and it is off by default because nodes and adjacencies come from the same reads and double-count.
Treat each read or read pair as **one observation**:

1. represent each neighbouring-state transition by an ordered syncmer signature spanning left
   context, the short allele and right context;
2. extract the corresponding signature from each fragment;
3. score the fragment once against compatible transitions;
4. add it as a transition factor between neighbouring block states;
5. do not also count its internal syncmers at full node weight;
6. apply fragment-level effective-sample-size discounting.

Test on gstm1 blocks 9 and 11; a synthetic locus with short indel alleles and variable flanks; both
regimes; and depth and error ladders. Success = materially higher call rate and lower error at gstm1
without harming cyp2d6 or ordinary blocks.

---

## 4. Real-read acquisition -- READY, independent of 2 and 3

Real LZO at lpa loses 0.0179 to read acquisition, and real call rates run 12-20 points below
simulated. Run a funnel on truth-labelled simulated reads:

```
generated -> maps to locus -> maps to the correct copy -> survives extraction
          -> yields markers/fragments -> reaches the model
```

separating unmapped divergent reads, reads mapped to the wrong repeat copy, discarded multi-mappers,
mates outside the region, and reads present but markerless. Only then implement rescue -- probably
locus or graph remapping that preserves alternative placements and retrieves mates. Re-run real LZO
afterwards.

---

## 5. KIV-2, two honest modes -- BLOCKED on 2-4

### Short-read mode

Report total array dosage only if independently calibrated; report composition and homologue
allocation separately; emit an equivalence set for unresolved allocation; suppress confident arbitrary
diplotypes. **Do not reuse `mass_bp`** -- its pre-registered replacement rule failed at B = 0.00 on
six of seven array blocks, mean length error 7,981 bp against `called_bp`'s 412.

Justification that this is the right answer and not a consolation: the emission's mean is
`lambda * (m_a + m_b) + mu`, a function of the SUM alone, so two allocations with the same summed
vector are **exactly tied by construction**, and 44% of KIV-2 loss is exactly that.

### Long-read mode

Long reads must not be reduced to the same k-mer counts -- the model never uses the fact that two
markers came from one read, so longer reads produce the same counts. Build a separate proof of
concept that preserves read order: align or chain each read against candidate allele walks; compute a
per-read compatibility score per walk; infer two haplotype paths or an equivalence class from all
reads; retain per-read linkage across repeat copies; connect blocks with the chain only after
read-to-path evidence works.

Decisive demonstration: a matched read-length ladder at 150 bp, 350 bp fragments, 2 kb, 5 kb, 10 kb
and 20 kb, where the same donors fail with short reads and approach their certified panel floor once
reads span the repeat structure. This directly tests the measured 2-5 kb crossover against a ~5.5 kb
unit.

This is the largest remaining item and is genuinely a second evidence mode.

---

## 6. Uncertainty and outputs -- BLOCKED on the evidence models

Calibrate GQ on untouched donors (currently ~1000x optimistic at GQ 40-60, saturating near 0.5% above
GQ 60). Separate confidences for total dosage, content and allocation. Emit equivalence-set size and
posterior mass. Reliability fields: marker count, independent fragment/clump count, evidence mode,
panel representability, no-call reason. The chain's confidence must not conceal a locally unresolved
emission. **No apparent accuracy obtained by raising the no-call rate.**

---

## 7. Correctness and test-contract defects -- READY, before new algorithms

- Reject non-numeric node names, or stop deriving coordinates from numeric names
  (`genotype_blocks.cpp:113` returns 0 for a malformed name).
- Validate every option range and every incompatible combination.
- Fix `n_scored_alleles` at non-representable blocks and the remaining diagnostic-contract errors.
- Require a genotype-enabled release-candidate build; genotype is `OFF` at `CMakeLists.txt:16`.
- **Fail configuration if genotype is expected but its test is not registered.** `genotype_stats` is
  registered only inside the enabling branch, so a green CTest can mean genotype was absent.
- Reassert direct/index byte identity, fragment-length provenance, 1-vs-8-thread determinism.
- Settle the depth estimator from existing held-out evidence; if the continuous one is neutral, keep
  the default and remove the public experimental choice.

---

## 8. Block-level diagnosis completion -- READY

Run the marker ledger and pair-vector identifiability calculation over **all** blocks carrying 80% of
excess error, not the three motivating blocks. Assign each to exactly one of: filtering discarded
sufficient evidence; allele shorter than the marker; retained evidence remains tied; emission ranks
truth poorly; candidate pruning; linkage override; read acquisition; panel limitation.

No algorithm work proceeds from marker count alone -- cyp2d6 block 5 already showed retained count is
a misleading proxy.

Outstanding from Phase 0: map dominant clumps to sequence spans. Partly done -- median 26 markers over
268 bp, 19.6% recurring elsewhere, and 15 clumps implicated in more than one donor carrying 44.6% of
the dominant-class loss. The sequence-level count of distinct differences per clump is still missing.

---

## 9. Code and CLI -- BLOCKED on flags being frozen

Two builds: a release genotype with stable user-facing options only, and a diagnostic build behind
`PANVAR_ENABLE_GENOTYPE_DIAGNOSTICS`. Behind the diagnostic build: truth and oracle inputs, noiseless
injection, pair probes, ledger and dump options, debug explanations.

**Remove** refuted algorithm branches from the production parser and scorer rather than hiding them --
git history and these reports preserve reproducibility. Of 44 current flags, ~9 are pure
instrumentation and ~5 are refuted experiments.

The public help should describe one coherent short-read algorithm, and if it succeeds one long-read
evidence mode -- not dozens of abandoned experiments.

---

## 10. Documentation -- LAST

Only once flags and outputs are frozen: `docs/modules/genotype.md`, `docs/algorithms/genotype.md`, and
a walkthrough covering one ordinary bubble, one short allele, one ambiguous array, and one
no-call/equivalence result. The algorithm then explains simply:

```
panel paths -> informative markers or fragment signatures -> read evidence
            -> local diplotype/equivalence likelihood -> linkage between blocks
            -> calibrated call or explicit uncertainty
```

---

## 11. Final release gate

Genotype enabled in the release build; clean build and install smoke; all genotype tests definitely
registered; sanitizers and thread determinism pass; the frozen dashboard satisfies its gates; results
regenerate from a manifest; public help free of diagnostic flags; documentation matching the binary.

---

## Discipline

Every experiment ends in **pass**, **fail-and-close**, or **production integration** -- never another
moving interpretation. No further KIV-2 likelihood experiment, and no opening of the 50-donor holdout,
until CYP2D6 filtering, GSTM1 fragment evidence and real-read acquisition have each had one gated
implementation pass.

## Closed routes -- do not retry without a new mechanism

`--all-kmers`; `marker_outlier`; `mass_window`; `mass_bp` as the dosage answer; `compositional`;
`robust_c`; read-derived block-class routing; `--no-region-unique` as a global restore; `--kmer-size`
for sub-k alleles; adjacency evidence at cyp2d6 block 5; clump-mean normalization and symmetric clump
clipping at KIV-2; the cosigt mapped-coverage route.
