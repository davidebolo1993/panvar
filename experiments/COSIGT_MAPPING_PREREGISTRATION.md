# Mapped-node validation rule

This note freezes the next comparison before mapped coverage is inspected for
any pair other than pair 9.

## Development observation

On deterministic LPA KIV-2 pair 9, reads were aligned to the 466 panel
haplotypes with BWA, projected to graph nodes with gfainject, and summarized by
gafpack.  Primary-alignment node dosage correlates 0.9951 with the exact truth
multiplicity.  Of the tested scores, an unweighted Poisson deviance over
variable bubble nodes selected a held-out-panel pair with 19 bp total absolute
per-homologue length error and mean identity 0.998661.  The production caller
had 22,255 bp error on this sample.

This result selected the method; pair 9 is therefore development evidence and
is excluded from the validation aggregate.

## Frozen method

For each sample:

1. Estimate haploid depth as half the median gafpack coverage on graph nodes
   of at least 100 bp traversed exactly once by every panel haplotype.
2. Divide bubble-node coverage by that depth to obtain observed diploid node
   dosage.
3. Exclude the two truth paths from the candidate panel.
4. Score every unordered panel pair with Poisson deviance between its summed
   traversal multiplicity and observed dosage, using every node whose
   multiplicity varies in the remaining panel.  Nodes receive equal weight;
   there is no length weighting, masking, clipping, tuning, or HMM.
5. Select the lexicographically first pair only if numerical scores tie.

The mapping arm uses primary alignments.  Pair 9 showed that projecting all
equal-quality alternate mappings and weighting them by occurrence both cost
orders of magnitude more storage and changed the good result into a 27,755 bp
miss.  This remains a useful negative control, not the frozen method.

## Validation set and outcome rule

Run unchanged on deterministic pairs 8 and 10--15.  Report total absolute
per-homologue length error (`dbp`) and mean sequence identity, alongside the
production caller on the same reads.  The method advances only if:

- aggregate `dbp` is lower than production;
- mean identity is not lower than production; and
- it does not worsen more than one validation pair by over one KIV-2 repeat
  unit (~5.55 kb).

All per-pair values remain visible even if the aggregate passes.  This is a
simulated BWA experiment approximating the historical bwa-mem2 pipeline; it
cannot establish performance on real libraries.
