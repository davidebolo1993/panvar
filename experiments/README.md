# Genotype experiments — NOT part of the release

Everything here belongs to `genotype`, the one module excluded from the default build
(`-DPANVAR_ENABLE_EXPERIMENTAL_GENOTYPE=ON` restores it). These scripts neither change nor import the
production modules. They are kept in the tree because their results are cited in
`docs/reports/genotype-round2-verification.md` and re-deriving them is expensive; they are not
maintained to the standard the shipped modules are held to, and nothing in the pipeline runs them.

Untracked `.tsv` output alongside these scripts is run output, not input.

---

# Multiplicity-aware genotyper proof-of-concept

This experiment is deliberately outside `src/` and `tests/`: it neither changes nor imports the
production genotyper. It tests the smallest useful version of a dosage-first model before that design
is proposed for implementation.

Run the controlled assertions:

```bash
python3 experiments/multiplicity_genotyper_poc.py synthetic
```

They demonstrate three distinct claims:

1. presence/absence cannot identify tandem copy number, while multiplicity identifies the diploid
   total;
2. internal repeat evidence may be ambiguous while source/allele and allele/sink junction features
   identify the bubble genotype;
3. a block HMM can switch copying states and recover a mosaic that no single locus-wide panel pair
   represents.

## Real held-out KIV-2 reproduction

First ask the current binary for its block diagnostic (the paths named below are the held-out sample):

```bash
build/panvar genotype \
  -i results/real_data/lpa/bubble/bubble.sorted.gfa \
  -b results/real_data/lpa/bubble/bubble \
  -r 'GRCh38#0#chr6:160500000-161000000' \
  -o /tmp/lpa_loo \
  -R READ1.fq.gz -R READ2.fq.gz \
  --truth-haplotypes 'HG00408#1#CM085956.1:162480704-162774500,HG00597#2#CM085784.1:162389888-162733646' \
  --exclude-haplotypes 'HG00408#1#CM085956.1:162480704-162774500,HG00597#2#CM085784.1:162389888-162733646' \
  --dump-block 13 -q
```

Then run the independent projection:

```bash
python3 experiments/multiplicity_genotyper_poc.py real \
  --conf /tmp/lpa_loo.block13.conf.tsv \
  --alleles /tmp/lpa_loo.block13.fa \
  --haploid-depth 11.5 --background 1 \
  --truth-bp 263385 --current-pair 53,91 \
  --gfa results/real_data/lpa/bubble/bubble.sorted.gfa \
  --bubbles results/real_data/lpa/bubble/bubble.bubbles.csv --bubble-id 7 \
  --truth-path 'HG00408#1#CM085956.1:162480704-162774500' \
  --truth-path 'HG00597#2#CM085784.1:162389888-162733646'
```

`--truth-bp` is printed only after ranking; it is not used by the estimator. The estimated span comes
from read-derived latent marker dosage and a calibration learned from graph paths. Exact sequence
identity is computed by the same vendored edlib kernel used by panvar.

The command also perturbs haploid depth by ±1%. If that changes the selected length band, the output
marks the dosage result unstable. A production implementation must integrate a depth posterior learned
from local flanks, and decline or report multiple CN states when that posterior crosses band boundaries;
it must not turn the point estimate into a confident hard call.

### Observed result on the motivating held-out sample

Using the block-13 dump generated from the current tree (`515f920`) with both sample haplotypes removed:

| result | production pair | dosage-first PoC |
|---|---:|---:|
| diploid span | 257,814 bp | 263,428 bp |
| absolute error against the graph truth (263,385 bp) | 5,571 bp | 43 bp |
| mean exact NW identity to the two truth walks | 0.825155 | 0.915611 |
| length-weighted NW identity | 0.824934 | 0.917052 |

The PoC learned `bp = 6.445981 * marker_mass - 148.596` from the graph paths (`R2 = 0.99985577`)
and selected the 262,442--263,594 bp supported band. This is a strong rescue of the concrete failure,
but not yet a production-ready call: changing haploid depth from 11.5 to 11.385 (-1%) selects the next
band, 267,975--269,137 bp. That sensitivity is a positive design result from the experiment: model
depth jointly with integer dosage and expose ambiguity instead of adding another hard mass window.

The matching in-panel control was regenerated at 30x with `tests/genotype_sim.sh` and the same two
truth haplotypes retained. Production selected alleles `41,54`; the PoC also ranked `41,54` first,
with 0 bp span error and exact NW identity 1.0. Excluding those two paths from that same run reproduced
the motivating production call `53,91`; its diagnostic FASTA and count table were byte-identical to
the held-out dump used above. Thus this one real-locus A/B test shows a rescue under exclusion without
regressing the exactly representable control. It is still one sample and one array, not cohort-level
validation.

## Intended production design, if the gate passes

- Build marker-to-graph-component incidence, including source/first-step and last-step/sink junctions.
- Infer local non-negative integer component dosage from read multiplicities with a bounded count
  loss. Keep a residual and decline when no small-residual explanation exists.
- For an ordinary bubble, project that dosage onto legal source-to-sink allele pairs.
- For a tandem array only, enable a copy-number stratum when marker mass predicts path length with a
  predeclared high R2 and supported lengths form a separated ladder. Do not use this as a generic mass
  window.
- Join local pair states with a diploid Li-Stephens-style HMM. Permit switches, report posterior switch
  probability, and call a mosaic only when several adjacent blocks support the new state.

The PoC does **not** yet establish calibration, GQ, robustness across samples/loci, or read mapping in
paralogous sequence. Those require cohort-level validation before replacing any production path.
