# Cosigt cross-pair validation rule

Written before evaluating deterministic LPA pairs 8--15.

## Development evidence

Pairs 0--7 use the deterministic selection in `tests/genotype_sim.sh` with seed 42. Pair 0 is the
original motivating KIV-2 failure. On pairs 1--7, an exact node-vector loss composed of 25% absolute
traversal-multiplicity L2 and 75% per-haplotype node-presence-dosage L2 reduced summed diploid length
error from 27,796 bp for the production sparse caller to 16,754 bp. The coefficient was selected after
seeing those pairs and is therefore development evidence only.

## Frozen validation

The primary method for pairs 8--15 is `node_mult_presence_l2_0.25`, unchanged:

```
loss = 0.25 * normalized_L2(node traversal multiplicity)
     + 0.75 * normalized_L2(per-haplotype node presence dosage)
```

Both terms compare the diploid sum of two panel paths with the diploid truth vector. The two truth
paths are excluded from the candidate panel. No node is selected or weighted using pairs 8--15.

Primary outcomes, in order:

1. summed per-homologue absolute KIV-2 length error;
2. mean exact sequence identity after the edit-distance-optimal homologue assignment;
3. number of validation pairs improved, unchanged (within 100 bp), or worsened against production.

The method passes only if it improves summed length error without reducing mean identity and without a
single validation pair worsening by more than one KIV-2 repeat unit (approximately 5.5 kb). Conditional
accuracy is never reported without call rate.

Raw Cosigt, absolute multiplicity L2, presence-only scores, node-length weighting, masks, centering,
and other mixture coefficients are secondary diagnostics. They cannot replace the primary result.

## Mapping arm

The same frozen 0.25/0.75 score will be evaluated after `gfainject` plus
`gafpack --len-scale --weight-queries`. Coverage is divided by an independently measured haploid flank
depth before scoring. Exact node vectors remain the representation oracle; mapped coverage measures
whether alignment and multi-hit projection preserve that signal. Any coverage-based node filtering is
trained only on pairs 0--7 or on panel-only statistics.
