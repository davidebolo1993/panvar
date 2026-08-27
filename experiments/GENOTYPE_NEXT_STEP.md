# Genotype next-step experiment: ideal multiplicity to an array-specific rescue

Date: 2026-08-21. This is an experiment, not a proposed default. Production C++ was not changed.

## Question and causal ladder

The LPA KIV-2 failure is evaluated in this order:

1. Is the sequence-certified pair present in the panel?
2. Is its marker-multiplicity vector unique?
3. Does the current likelihood recover it from its own noiseless vector?
4. Does it recover it from the held-out truth's exact marker multiplicities?
5. What changes when read sampling is introduced?
6. Does candidate pruning or the HMM change the local answer?

`ideal_multiplicity_oracle.py` reproduces the default negative-binomial emission and the 50/50
detected-marker/length candidate pruning independently of `src/`. It scores all 104,653 unordered
pairs for the 8,754-syncmer dump. The ESS factor and full-universe baseline are omitted because both
are constants within a block and cannot alter ranks. `channelized_emission_poc.py` tests a bounded
support channel separately from NB multiplicity deviance.

One terminology correction is important: 57,958 is the final called pair `(53,91)`'s rank under the
sequence oracle. It is not the certified `(42,178)` pair's emission rank.

## Results on the original 8,754-marker KIV-2 block

The target `(42,178)` has exactly one equivalent diploid marker vector: itself. From its own noiseless
expected counts it ranks first, and both alleles survive the 64-candidate pruning. Thus the model can
recognize an exactly represented sample and pruning is not the motivating failure.

From the held-out haplotypes' exact multiplicities, before read noise, the result reverses:

| evidence | full-space best | certified rank | certified rank after pruning | called pair rank |
|---|---:|---:|---:|---:|
| panel-self noiseless | 42,178 | 1 | 1 | 2,923 |
| held-out truth noiseless | 53,91 | 144 | 25 | 1 |
| observed reads | 53,294 | 212 | 30 | 2; 53,91 is best after pruning |

The default local emission therefore already selects the production pair. The HMM is not the cause at
this block, and neither sampling nor candidate loss creates the main error.

The per-marker contrast explains it. Relative to `(53,91)`, the certified pair gains about 249 log
units on 1,166 shared markers carrying 38,194 truth copies at different multiplicities. But 61 markers
carried only by `(53,91)` contain just 58 truth copies and cost the certified pair about 1,974 log
units. Eleven independent fragment-scale clumps of presence mismatch overpower the multiplicity
signal spread across 54 clumps.

Three plausible likelihood-only repairs were falsified:

- the existing bounded outlier improves the certified truth rank from 144 to approximately 10--15,
  but never makes it best;
- inverse-clump marker weighting makes the result worse because the unsupported markers occupy small
  clumps and are upweighted;
- an explicit bounded support channel plus shared-support NB deviance also leaves the target far from
  first.

This is an off-panel projection problem in the selected marker geometry, not merely an overconfident
negative-binomial tail.

## Dense-marker falsifier

The same graph, reads, exclusions, depth model, 64-candidate budget and HMM were rerun with
`--all-kmers`. The KIV marker set grows from 8,754 to 52,492; the whole locus grows from 54,453 to
317,290 retained markers.

Dense markers alone do not change the call: production still selects `(53,91)`. In the independent
candidate oracle the certified pair improves to rank 7 by raw L1 multiplicity and rank 4 by squared
error, but the default NB emission ranks it 26th on noiseless truth and 30th on reads. The same support
veto has simply become denser: shared multiplicities favour the target by about 1,291 log units, while
340 unsupported markers carrying 306 truth copies cost about 10,195.

Dense markers plus a small bounded outlier produced a material rescue on the motivating pair. At `marker_outlier=0.001`,
the observed local emission and the end-to-end HMM select `(178,312)`:

| pair | mean exact identity | length-weighted identity | diploid bp |
|---|---:|---:|---:|
| certified `(42,178)` | 0.996820 | 0.996354 | 263,428 |
| production `(53,91)` | 0.825155 | 0.824934 | 257,814 |
| dense + outlier `(178,312)` | 0.972433 | 0.976058 | 268,973 |

Across the whole locus, length-weighted identity rises from 92.41% to 98.98%, and total sequence
distance falls from 49,966 to 5,658 bp. The new KIV call has GQ 8.69 and `LOWGQ`, versus GQ 69.9 and
`PASS` for the wrong production pair. That uncertainty is desirable: short reads support a much
better content class but do not identify the unique sequence optimum.

The setting is not a global default. Applied to every block it reduces exact scored blocks from 15/20
to 12/20 and creates 13 low-GQ calls. Dense, bounded-off-panel evidence should therefore be an
array-specific hypothesis to test, not a replacement for the ordinary-bubble model. The outlier value
was examined on the motivating truth and therefore cannot establish generalisation.

## Gate B result: dense plus outlier rejected

The required held-out comparison was subsequently run on eight LPA pairs. The motivating pair alone
improved (49,901 to 5,588 bp of KIV-2 sequence distance). One independent pair worsened fourfold
(5,550 to 22,185 bp), and the remaining six were essentially unchanged at KIV-2. Across all pairs the
headline distance improved by 36%, but excluding the motivating pair it worsened by 60% (27,796 to
44,429 bp). Locus-wide exact scored blocks changed from 106/143 to 105/143, while call rate fell from
75.5% to 66.3%. The number of correct PASS calls stayed at 92; conditional PASS accuracy rose only
because the PASS denominator fell from 123 to 108.

This falsifies implementation-recipe steps 2--3 below for the tested `--all-kmers` plus
`marker_outlier=0.001` design. They are retained as historical context, not as recommendations. The
separate conclusions in steps 1 and 4--8 remain useful: keep diagnostic ranks distinct, separate
dosage from content, expose equivalence sets, retain block switching, add per-homologue truth, and
reserve path synthesis for genuinely unrepresentable truth.

## Implementation recipe

1. Correct the report language: keep sequence-oracle rank, local emission rank, candidate survival and
   final HMM call as four separate fields.
2. **Rejected by Gate B:** add an array-only marker mode behind a flag. Build dense canonical k-mers only for blocks already
   classified as arrays/high-redundancy; retain syncmers for simple bubbles and backbone blocks. This
   contains memory and avoids the measured ordinary-block regression.
3. **Rejected by Gate B:** give the array channel an explicit bounded off-panel component. Do not silently enable a value.
   Sweep a predeclared grid on training samples and select it without using the motivating pair.
4. Keep dosage independent of content. `mass_bp` (with its depth uncertainty) estimates total array
   amount; the dense marker posterior estimates local content. Do not turn the selected panel pair's
   length into copy number.
5. Emit an array equivalence set/posterior, not only its argmax: local top pairs or allele marginals,
   posterior mass, calibrated content GQ, sequence-length band, and the independent `mass_bp` interval.
   A low-GQ but 97%-identity class is more truthful than a high-GQ 82%-identity hard call.
6. Preserve block-level switching for mosaics. Flanks and ordinary blocks can phase local states; an
   array whose order exceeds fragment reach must be allowed to remain an equivalence class. Do not
   force one locus-wide panel pair.
7. Extend the truth-only diagnostic dump with `truth_mult_h1` and `truth_mult_h2`. The current diploid
   sum shows that `(42,178)` fails but cannot localize the failure to one homolog or expose compensating
   mixtures.
8. Only after this gate, test graph-flow synthesis for genuinely unrepresentable blocks. KIV-2 already
   contains a unique sequence-nearest pair; enlarging its candidate space is not its first remedy.
   Flow remains useful for novel mosaics if it is constrained to two realizable walks and may return
   multiple equivalent walks.

## Required validation before a default changes

Run four arms on many held-out LPA pairs: sparse/default, sparse/outlier, dense/default, and
dense/outlier. Include the motivating pair only as a final test. Predeclare the outlier grid and choice
rule. Report per-haplotype exact identity, length-weighted identity, total edit distance, candidate
survival, local rank, final HMM rank, GQ calibration, no-call rate, runtime and peak memory.

Repeat on every real array locus and on ordinary-bubble, deletion, paralog, folded-consensus and mosaic
fixtures at several depths. The acceptance condition is not merely a better KIV argmax: array identity
must improve without damaging ordinary calls, and low GQ/equivalence coverage must be calibrated.

The committed result tables are:

- `kiv2_ideal_multiplicity_results.tsv`
- `kiv2_outlier_sweep.tsv`
- `kiv2_clump_sweep.tsv`
- `kiv2_channelized_sweep.tsv`
- `kiv2_all_kmers_candidate_oracle.tsv`
- `kiv2_all_kmers_outlier_sweep.tsv`
