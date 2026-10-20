# M9-F: bounded diagnostic prototype of junction-aware proposals

**Question:** can read evidence spanning block boundaries give useful repair shortlists where the allele-internal proposal
emissions are completely uninformative? This is a diagnostic, not caller integration.

**Results:** `tests/results/m9f_junction_proposals_20261020/`.

## Frozen rule JP1 (one version; `jprop.py`)

**Scope.** It applies to a block whose existing allele-pair proposal emissions are all equal. sv:5 qualifies in the LOO
model (M9-E: 6,786 pairs, one value). For the LZO models this is re-checked by `ranks.LZO.d1` and `ranks.LZO.d2`, and the
result is reported.

**Rule:**
- **Context:** the authoritative v2 representative (`sequence_locus.tsv` `representative_labelling_hap1/2`), verified by
  spelling it and matching `representative_hap1/2_sha256`. This includes the sequence-tied d5, which is used as a frozen
  context and never treated as a strict call.
- **Window:** W(h, a) = L_h + allele a + R_h, where L_h and R_h are the last and first K − 1 = 30 bases of incumbent
  homologue h around the block. They are spelled through neighbouring blocks and shortened only at locus ends.
  - Every K-mer of W overlaps the allele or spans its boundary.
  - For an empty allele the K-mers span L | R.
  - For an allele shorter than K − 1, some span both boundaries.
- **Markers: all canonical 31-mers,** the lexicographic minimum of the K-mer and its reverse complement, which equals the
  caller's 2-bit minimum. K-mers with a non-ACGT base are skipped.
  - **Declared deviation:** the caller's proposal markers are closed syncmers (s = 15). JP1 uses every K-mer, because
    syncmer thinning leaves junction windows of short alleles nearly empty.
  - Counting is verified by an independent oracle, and the per-copy rate is verified empirically (gates below).
- **Inventory:** the distinct canonical K-mers over W(h, a) for every allele a and both homologues. It is fixed for the
  comparison. Multiplicities are kept per window.
- **Copies:** copies(m | x on h1, y on h2) = n(m, W(1, x)) + n(m, W(2, y)) + e(m).
  - e(m) counts m in both incumbent haplotypes **outside** the two window spans, so occurrences elsewhere in the locus
    are kept.
  - Occurrences outside the locus are not modelled (declared).
- **Reads:** y(m) is the number of occurrences of canonical m in every read of both mates. This is the caller's
  `count_marker_observations` convention applied to all K-mers.
- **Expected count:** λ(m) = rate · copies(m) + background, using the run's own fit: `mosaic_status`
  `marker_depth_rate_per_copy` and `marker_background`.
- **Score:** S = Σ over the inventory of [y(m) · log λ(m) − λ(m)], the Poisson kernel with constants dropped. It is
  computed exactly as the delta over the K-mers the candidate's windows touch.
  - **Consequence, declared:** an allele predicting more K-mers that are unobserved scores lower even when no read
    touches the block. Equal-length alleles without read evidence tie.
- **Pair ranking:** the unordered pair {x, y} is ranked by max(S(x on h1, y on h2), S(y on h1, x on h2)). Both
  assignments are kept for full scoring.
- **Shortlist:** **32** unordered pairs.
- **Ties:** |ΔS| ≤ 1e-9 is a tie, chained within a class. Tie classes are reported, including whether the cutoff class
  crosses rank 32. The order inside a class is (x, y) ascending, which is **not** evidence that reads distinguish the
  pairs.

**What JP1 is not:** a likelihood, a confidence value or a final-score term. There is no parameter grid, no cutoff tuning
and no iteration.

## Gates (`gates`, synthetic, K = 31; the oracle is an independent whole-haplotype and read scan)

| gate | checks |
|---|---|
| F1 | a 5-bp allele distinguishable only through junctions ranks first, untied; fast delta = oracle |
| F2 | an empty allele; oracle |
| F3 | a 3-bp allele: window K-mers = 33 = len + K − 1, the 27 K-mers spanning both boundaries are present |
| F4 | the junction K-mer repeated elsewhere in the locus is counted in e(m); oracle |
| F5 | reverse-complement invariance of every pair score |
| F6 | homologue-swap invariance |
| F7 | allele-index permutation with sequences preserved |
| F8 | the same sequence under two labels scores and ties identically |
| F9 | genuinely uninformative (equal-length alleles, no read touching the block) is one tie class |

**Real-data gates per cell, in `propose`:**
- the context SHA matches v2;
- the empirical all-K-mer rate (mean read count / 2 over K-mers present exactly once in each incumbent haplotype) lies
  within [0.75, 1.25] × the fitted rate;
- the fast delta matches the oracle within 1e-6 for three candidates at ranks 1, middle and last.

A failed gate sets `propose` to exit 1, and that cell is not scored.

**Disclosed development before the freeze:** three synthetic gate runs, all under 0.2 s. The third followed two fixes:
- F3's expected count had an off-by-one;
- F9 used unequal allele lengths, which the declared model treats as informative.

The rule was not changed. Counted as 3 s.

## Units and order (budget 1,200 s; spent + timeout + 30 s for the analysis)

| step | timeout |
|---|---|
| `gates` | 60 s |
| `ranks.LZO.d1`, `ranks.LZO.d2`: the original probe sv:5 ranking for the LZO models (M9-E's covers the LOO cells) | 60 s each |
| `propose.<cell>` for LOO d1, d3, d6, LZO d1, d2, then LOO d2, d4, d5 | 120 s each |
| `freeze`: SHA-256 of every shortlist, and the score-set candidate files, before any truth or full scoring | 10 s |
| `score.<cell>`, same order | 200 s each |
| `analyse` | 30 s |

**What `score` does.** It runs the existing `genotype_probe score-set`, with the original likelihood, exposure and copying
prior and the run's model, curation, reads and library.
- **Candidates:** the incumbent, with its v2 `representative_score` as the expected total (a non-reproduction refuses
  the cell, exit 4), plus both assignments of every shortlisted pair, deduplicated by unordered label pair.
- **Never** all 6,786 pairs. **No** caller output is altered, and no replacement prediction is emitted.

**What `analyse` reports, after the freeze:**
- ranks and inclusion of M9-D's frozen improving pairs (d1 6/77, d3 34/94, d6 13/21);
- the truth pair's rank (old and new) and exact block agreement where attainable;
- diversity, tie sizes, and overlap with the original 32;
- the best improvement over the incumbent, with Δlikelihood and Δprior;
- whether the change is correct;
- whether either LZO control would change;
- a real sequence change versus a label-only one (score-set unique haplotype identity);
- proposal and scoring time, and peak memory.

**Status of the samples:** development samples, inspected. Not held-out validation.

**Boundaries:** no prior change, joint search, wider shortlist, confidence remapping, filter, default promotion or caller
integration. The persisted locus index is untouched. All JP1 state is in-process and sample-scoped.
