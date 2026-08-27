# Inverse-multiplicity node validation

This note freezes the next experiment before reads or coverage are generated
for pairs 16--23.

## Development result

Mapped node coverage from pairs 8--15 (plus the original mapped pair 9) showed
that node error increases strongly with panel multiplicity.  Hard filtering is
undesirable at KIV-2 because high multiplicity is real biology.  Three
panel-only soft weights were explored with absolute diploid dosage L2:

`1 / (c + mean haploid panel multiplicity)`, for `c = 0.25, 1, 4`.

`c = 1` was best on those eight development samples: 33,574 bp aggregate
per-homologue length error and mean identity 0.982046.  A leave-pair-out weight
learned from simulated mapping residuals reached 44,461 bp, so the simpler
panel-only rule was also better than the supervised alternative.

## Frozen method

For each new sample:

1. Simulate and map reads exactly as in the mapped-node experiment: BWA primary
   alignment to all panel haplotypes, gfainject projection, gafpack
   length-scaled node coverage.
2. Estimate haploid depth as half the median coverage on >=100 bp nodes that
   every panel path traverses exactly once.
3. Within LPA bubble 7, exclude the two truth paths and enumerate all unordered
   remaining panel pairs.
4. For every node whose multiplicity varies in that remaining panel, compare
   observed diploid dosage with the pair's summed traversal multiplicity using
   squared error weighted by `1 / (1 + panel_mean_multiplicity)`.
5. Choose the minimum-loss pair deterministically.  There is no hard mask,
   fitted mapping-error model, sequence oracle, HMM, or sample-specific tuning.

## Independent validation

Pairs 16--23 are untouched validation data.  Compare against the production
caller on the identical reads using total absolute per-homologue length error
and mean sequence identity.  Advance the channel only if:

- aggregate error is lower than production;
- mean identity is not lower than production; and
- no more than one pair worsens by over one KIV-2 repeat (~5.55 kb).

Report every pair and equivalence/tie information.  Passing would justify a
second-locus and real-library test, not a production default: BWA is an
approximation to the historical bwa-mem2 route and primary placement among
equally scoring panel alignments can be implementation-dependent.
