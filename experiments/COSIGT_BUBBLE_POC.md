# Bubble-local Cosigt experiment

Date: 2026-08-22. This is an independent experiment; no production source was changed.

## Question

Can Cosigt's summed-haplotype coverage vectors and cosine similarity be moved from a whole regional
graph to individual panvar bubbles, and can the vectors be populated alignment-free by counting the
same syncmers/k-mers that `panvar genotype` already indexes?

The motivating falsifier is the held-out LPA KIV-2 pair. The panel contains the unique
sequence-nearest allele pair `(42,178)` (length-weighted identity 0.996354), but the production sparse
syncmer emission selects `(53,91)` (identity 0.824934).

## What was run

`cosigt_bubble_poc.py` applies Cosigt's weighted cosine exactly to the 457 allele x 8,754 marker
multiplicity matrix exported by `panvar genotype --dump-block 13`. It scores all 104,653 unordered
pairs from four query vectors: the target pair's own vector, the held-out walks' exact marker
multiplicities, observed read counts, and background-subtracted read counts. It also reproduces the
old pipeline's panel-mean IQR coverage mask in marker coordinates and tests inverse-clump weights.

`cosigt_node_poc.py` is a stricter ideal oracle. It constructs the original Cosigt node-traversal
vectors directly from all graph paths within bubble 7, removes the two truth paths, and scores their
exact summed graph vector against all remaining path pairs. No read mapping or sampling can explain a
failure in this arm. The test includes Panplexity's actual mask and continuous weights, generated from
the released LPA graph with Cosigt's historical settings (`k=16`, window 100, linguistic complexity,
automatic threshold).

Both scripts report exact top-score equivalence counts. This matters because cosine differences do
not by themselves define a calibrated genotype quality.

## Result

Cosigt is mechanically applicable per bubble, and marker counting reproduces its node-vector answer,
but cosine does not solve KIV-2.

| evidence | mask | best pair | certified rank | best pair identity |
|---|---|---:|---:|---:|
| exact held-out marker multiplicity | none | 53,284 | 1,120 | 0.858242 |
| observed marker counts | none | 53,284 | 1,308 | 0.858242 |
| exact held-out node traversal | none | 53,284 | 2,402 | 0.858242 |
| exact held-out node traversal | Panplexity | 130,365 | 2,683 | 0.846159 |

The marker and node formulations independently choose `(53,284)`. That is a useful positive result
for an alignment-free implementation: syncmer multiplicity vectors retain the signal that drives the
node-coverage cosine. It is not an accuracy result. The cosine winner improves on production by about
3.3 identity percentage points. A dense-marker plus bounded-outlier experiment reached 0.976058 on
this motivating pair, but the held-out validation below subsequently rejected that method.

The ideal node cosine is 0.999889 for a pair whose sequence identity is only 0.858242. Cosine is
therefore badly uncalibrated as a hard-call confidence at this array. Its scale invariance is also
fundamental: if two array candidates differ only in total copy number, they have cosine 1. Cosine can
describe relative content or shape; it cannot be the total-copy-number channel.

The corresponding exact whole-locus node oracle does not rescue this held-out pair. Without masks it
selects the two HG00609 paths; their KIV alleles are `(53,54)` and reach 0.876107 identity. The global
coverage mask improves the certified KIV content class from rank 768 to rank 65, but the winner is
still `(53,54)`. Thus the usual good whole-locus behaviour is plausibly the result of distinctive
flanks and neighboring variants identifying a globally similar panel haplotype. It is not evidence
that the KIV node subvector identifies local sequence content. In this example the global evidence
actually prefers two paths from one donor while the desired local KIV alleles live on other paths.
That is precisely the trade that a switching block HMM is meant to avoid.

Several panel-only retention rules were also tested. None rescued `(42,178)`:

| local node/edge rule | certified rank | winner identity |
|---|---:|---:|
| all nodes | 2,402 | 0.858242 |
| inverse panel-mean node multiplicity | 666 | 0.876170 |
| node IDF | 519 | 0.876170 |
| all oriented edges | 1,946 | 0.876170 |
| inverse-mean oriented edges | 344 | 0.876170 |
| multiplicity residualized against bubble length | 2,377 | 0.876170 |

Node weights based on the fraction or number of graph-unique 16/31/51-mers were worse again (best
certified rank 1,833). These negatives matter: “retain unique nodes” and “retain variant nodes” are
not sufficient recipes for this array. Neither a hard subset of graph nodes nor the tested dense
sequence-marker likelihood is a validated rescue.

## Subsequent held-out falsification

The dense-marker plus bounded-outlier setting was then compared with the sparse baseline on eight
held-out LPA pairs. It fixed the motivating pair (49,901 to 5,588 bp of KIV-2 sequence distance), but
worsened another pair fourfold (5,550 to 22,185 bp) and left the other six essentially unchanged at
KIV-2. Excluding the motivating pair, total KIV-2 distance increased by 60% (27,796 to 44,429 bp).
Across the locus, exact scored blocks changed from 106/143 to 105/143 and call rate fell from 75.5%
to 66.3%. Correct PASS calls remained exactly 92; the apparent PASS accuracy increase from 74.8% to
85.2% was entirely the smaller denominator (123 to 108 calls).

This rejects the particular dense-marker plus `marker_outlier=0.001` channel as an implementation
target or default. It does not invalidate the independent diagnosis that the sparse emission chooses
the wrong pair from noiseless held-out multiplicities, nor does it rule out every possible robust
likelihood. Any future rescue must be selected and tested without the motivating pair and must report
call rate alongside conditional accuracy.

## Masks

The suspicion that Panplexity removes most KIV-2 nodes is not supported on the original graph:

| filter, computed over the locus | KIV nodes removed | KIV bases removed |
|---|---:|---:|
| Panplexity | 306 / 3,789 (8.08%) | 1,539 / 16,695 (9.22%) |
| Cosigt panel-coverage IQR | 2,106 / 3,789 (55.58%) | 6,627 / 16,695 (39.69%) |
| combined hard mask | 2,285 / 3,789 (60.31%) | 7,948 / 16,695 (47.61%) |

The dangerous part is the global coverage-spike mask. High panel traversal is the expected biology
inside KIV-2, not an artefact. Recomputing the IQR within bubble 7 removes zero nodes, but still does
not rescue the certified pair. In the alignment-free marker analogue, the IQR mask removes 1,907 of
8,754 markers and worsens the certified rank. Neither raw high multiplicity nor low sequence
complexity should therefore be a hard exclusion in an array. They should inflate uncertainty or
down-weight a channel only after conditioning on the expected panel multiplicity.

## Recommended hybrid

Do not replace panvar's emission with Cosigt. Reuse its useful geometric idea as one explicitly
limited channel:

1. Count canonical syncmers/k-mers directly from FASTQ and project each count to its pre-indexed
   bubble node, edge/junction or repeat component. No read-to-graph alignment is required.
2. Keep three local channels separate:
   - **mass**: absolute background-corrected marker count, calibrated by flanking depth, for total
     array amount/CN;
   - **content shape**: cosine or an angular residual on normalized internal marker multiplicities;
   - **junction/support**: bounded evidence for allele-specific edges and flank-to-allele junctions.
3. Use a per-bubble robust residual mask: compare observed count with the count expected from panel
   multiplicity and depth. Do not mask a coordinate merely because its mean multiplicity is high.
   Treat Panplexity complexity as variance inflation or a soft weight, not automatically as zero.
   For the content channel, prioritize multiplicity residuals and observable oriented junctions, but
   require held-out validation: the simple inverse-mean/IDF versions improved rank without solving the
   block.
4. Calibrate the relative channel weights and cosine temperature only on held-out samples. Emit an
   equivalence set/no-call when several pairs have indistinguishable shape. Never convert cosine
   directly to GQ.
5. Feed the resulting local likelihoods to the existing block HMM. The HMM may switch panel states
   between bubbles, so A-A-B-B and A-B-A-B mosaics remain representable. It should not force an order
   inside KIV-2 that fragments shorter than a repeat unit cannot observe.
6. Use cosine as a shortlist or posterior-shape feature, not as an expansion of the hypothesis space.
   KIV-2 already contains the unique sequence-nearest pair and cosine ranks it poorly.

The synthetic `multiplicity_genotyper_poc.py synthetic` assertions still establish the surrounding
architecture: absolute multiplicity identifies CN where presence cannot, junction markers identify a
bubble when internal markers cannot, and a block HMM recovers an A-A / B-B mosaic that one locus-wide
pair cannot represent.

## Engineering note on the original executable

The historical Go implementation should not be embedded verbatim. Its result collector ranges over
`resultsChan` but is also the only goroutine that attempts to close that channel, while the main
goroutine starts reading the result map immediately after workers finish. Up to a channel buffer of
scores can therefore be omitted nondeterministically. The mathematical cosine tested here is only a
few lines and should be implemented inside panvar with deterministic pair ordering. For repeated
samples, precompute the weighted panel Gram matrix; each sample then needs one panel-query product
plus O(A^2) pair assembly rather than materialising every pair-by-marker vector.

The numeric outputs are in `kiv2_cosigt_bubble_results.tsv` and `kiv2_cosigt_node_results.tsv`.
