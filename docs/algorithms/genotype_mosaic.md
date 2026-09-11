# Mosaic genotyper: development contract

This document is the implementation contract for the replacement short-read genotyper. It describes
work in progress; `panvar genotype` is not exposed until the end-to-end caller meets the acceptance
gates below.

## Objective

Call one unordered pair of alleles per block while permitting the two sample haplotypes to be mosaics
of panel templates. The caller has two evidence stages:

1. A PanGenie-like syncmer emission and Li-Stephens model proposes a small, deterministic set of
   whole-locus mosaic diplotypes.
2. A Locityper-like paired-read likelihood rescores those fixed candidates against their spelled
   sequences. Proposal scores are not added a second time; the final prior is recomputed exactly for
   each fixed mosaic.

The proposal may be approximate. Candidate scoring must be exact or return a certified interval and
refuse to choose a winner when intervals overlap.

## Components

New code is split by contract rather than accumulated in a command file:

- `genotype_model`: immutable panel, read, library and candidate types.
- `marker_model`: canonical syncmer observations and per-block proposal emissions.
- `ls_hmm`: Li-Stephens transitions and the exact prior of a fixed mosaic.
- `mosaic_beam`: deterministic top-K proposal with backpointers.
- `mosaic_spelling`: strict block-allele-to-sequence construction.
- `fragment_likelihood`: candidate-neutral paired-read placement and likelihood.
- `genotype_command`: parsing and orchestration only.

The removed factor caller is available at Git tag
`archive/genotype-factor-2026-09-11`. Code is extracted from it only one primitive at a time, with a
focused test; old source files, commands and diagnostic flags are never restored wholesale.

## Model contracts

### Library model

Read lengths are observed per fragment. For a pair with mate lengths `r1` and `r2`, the physical
insert floor is `max(r1,r2)`, so overlapping mates remain valid. Insert mean, standard deviation and
support width are estimated or explicitly supplied once, validated, and recorded with provenance.
The insert distribution is a discrete Gaussian normalized over that fragment's physical and declared
support. No call site may carry its own insert constants.

### Mosaic prior

At a block boundary the haploid template transition is

```
T = (1-r) I + (r/n) 11^T.
```

For a fixed allele mosaic, its prior marginalizes over every compatible template history in
`O(blocks * templates)`. It is not the score of the best retained beam history. The two homologues are
independent under this prior; an unordered pair of distinct mosaic sequences aggregates both orders
and receives exactly a `log(2)` multiplicity.

### Fragment score

For candidate haplotypes `a,b`, a paired fragment contributes through the existing signal/background
model only after each haplotype's placement mass has been computed. Homozygotes carry two copies. The
implementation must preserve physical-origin multiplicity: origins with identical edit and insert
statistics contribute more mass, while rediscovery of one origin through several seeds contributes
once.

Version one is paired-end, substitution-only and Hamming-certified. Non-ACGT sequence, unsupported
indel scoring, or an exhausted work budget refuses rather than silently changing model.

## Implementation sequence and gates

1. **Model foundation.** Fragment-specific insert support and exact fixed-mosaic Li-Stephens prior
   agree with independent exhaustive oracles. Implemented by `panvar_genotype_mosaic_core`.
2. **Panel/block model.** Parse the chain into immutable block alleles. Validate every panel template
   has exactly one allele per block and round-trip panel paths through the strict speller.
3. **Marker proposal.** Score all alleles at C4 (no 64-allele truncation), then construct deterministic
   top-K mosaics. Retain flat arrays and backpointers; candidate identity must not depend on threads.
4. **Strict spelling.** Spell every candidate across all blocks with a complete coordinate map. A
   panel-carried candidate must reproduce its panel sequence byte-for-byte.
5. **Fragment likelihood.** Index read pieces, scan each unique candidate haplotype once, verify full
   mates, join only valid FR placements, and compress statistically identical physical origins with
   an exact log-multiplicity. Score only the requested candidates, never their full Cartesian square.
6. **Certified refinement.** Start with bounded placement work; deepen only candidates whose score
   intervals can still win. Emit a call only when the winning lower bound exceeds every competing
   upper bound.
7. **Command integration.** Add one small `panvar genotype` surface. Experimental probes remain unit
   tests or separate developer tools, not public flags.
8. **Performance gate.** One frozen C4 donor at `K=64` must finish in at most 90 seconds and remain
   under 4 GB RSS before any cohort run. Report stage timings and counters; profile the real slow stage
   before optimizing it.
9. **Accuracy.** Run the K ladder, then an 8-donor development pilot, then held-out LZO donors at
   nested coverages. Per block, plot true edit distance (edlib) for both predicted haplotypes. LOO uses
   the closest retained alternative as its availability floor and is evaluated only after LZO.

## Deliberate exclusions from version one

- higher-order fragment factors, ownership scopes and factor planning;
- the `genotype-frag` command and hybrid/audit/probe flags;
- experimental depth/CN likelihoods;
- indel-aware fragment emission without an exact placement oracle;
- calibrated global posterior or GQ claims from a truncated candidate set.

Until omitted candidate mass is bounded, output probabilities are labelled candidate-set weights and
the runner reports the score gap to the second candidate rather than a globally calibrated posterior.
