# M3-B: bounded two-block search audit (diagnostic)

**Question.** Can generic two-block proposals recover the demonstrated C4 mosaic search failure under the unchanged
scorer and prior?

**Ancestry.** Branch `m3-b-two-block` from `c8dd50b` (M3-A). That is `6d85cce` (M1-E) plus M3-A diagnostics; there are
no library changes since baseline `282e64e`, and M2 stays closed.

**Scope.** Genotype search only. There is no production integration, prior change, dosage repair, residual sequence
editing, novel-base reconstruction or published FASTA. Panel alleles are combined into trial mosaics for diagnostics
only.

## 1. Frozen cells and inputs (`INPUTS.tsv`, SHA-256)

| cell | run | graph and curation | reads and declared library |
|---|---|---|---|
| `mosaic_2x` | C4 `mosaic/one_homologue_once_1/call.2` (M1-A SEARCH error) | `staged.gfa` + `mosaic_exclude.tsv` | `cov2` |
| `mosaic_5x` | its correct twin `call.5` | same | `cov5` |
| `lzo_5x` | correct LZO control `lzo/HG00096/call.5` | `c4.gfa` | `reads/HG00096/cov5` |

- **Settings.** All three runs use profile `genotype_c4_v6_campaign_local` (sha `74a3723b…`), switch probability 0.01,
  mismatch band 6, k = 31, s = 15. All share model fingerprint `d2b8c6c2e0b4079d`.
- **Incumbent.** The rank-1 representative row of `<run>.equivalence_set.tsv` (full-precision total and allele paths).

## 2. Proposal rule (`two-block-catalogue`; truth-free)

1. **Ranking.** For every block b and homologue h, the alternative alleles a ≠ incumbent are ranked by the caller's own
   read-derived proposal emission with the **other homologue's allele fixed**: `emission.at(a, other)` descending, ties
   by allele index.
   - The marker panel, observations, depth fit and emissions are built as in `score-set --rank-blocks`, with k and s
     from the run.
   - The incumbent is kept separately.
2. **Kept alleles.** The top **4** non-incumbent alternatives.
3. **Candidates.**
   - Every single-block move.
   - Every combination changing **two different blocks on the same homologue**, for every block pair, whether or not
     either intermediate loses.
   - No block is selected by truth, known error coordinates or evaluator output.
4. **Deduplication.** By unordered complete sequence pair, keyed by two independent 64-bit hashes plus length per
   haplotype. Every equivalent move and its provenance is kept, and candidates equal to the incumbent are counted.
5. **Limit.** **2,048 raw candidates per cell.** If the preflight exceeds it, the unit stops before writing a catalogue,
   and the score unit is recorded as not run. Nothing is trimmed and the limit is never raised.
6. **Preflight report.** Unrestricted one-block changes (2·Σ(A_b−1)) and two-block same-homologue changes
   (2·Σ_{b<c}(A_b−1)(A_c−1)), computed but never materialised.

**Scope of "exhaustive".** The catalogue is exhaustive only within the declared top-4, same-homologue two-block
neighbourhood, not over all genotypes.

## 3. Scoring (`score-set`, unchanged; at most 180 s per unit)

- **Model.** The full fragment likelihood (paired-read placement summation), the exposure model and the sequence-level
  prior at r = 0.01, with sequence-equivalent paths aggregated by bytes.
- **Reproduction first.** The incumbent's saved total must reproduce (|Δ| ≤ 1e-6) or score-set refuses before
  reporting.
- **Recorded:** candidate counts, distinct haplotypes, wall time, peak RSS and score-set's work line.

## 4. Selection (`rank`; truth-free; frozen before evaluation)

- **Outputs:**
  - every catalogue member ranked by total, with exact ties reported;
  - the winner, the incumbent's rank and the count above the incumbent;
  - for the best two-block candidate, the incumbent, both single-block intermediates and the joint candidate, each
    decomposed into fragment, exposure and prior.
- **Freeze.** A SHA-256 digest of the selection is written. The evaluation units refuse to run without it.

## 5. Evaluation (`two-block-eval`; truth enters only here)

**Truth.**
- **Mosaic:** per block, the first carrier template listed in `truth.tsv` that is present in the model.
- **LZO:** the donor's two panel templates.

**Reported per candidate:**
- exact sequence-pair match;
- block projection (blocks whose unordered allele pair differs);
- unordered per-block edit distance;
- a phase diagnostic: blocks, heterozygous in both candidate and truth, where per-block pairing disagrees with the
  global homologue assignment.

**Also reported:**
- whether the truth pair was generated, and its rank;
- whether its intermediates would have been rejected (they score below the incumbent);
- if truth is absent, the proposal ranks of the truth alleles, read from the frozen truth-free ranking table;
- whether either correct control deteriorates.

## 6. Interpretation rules (declared now)

- **A more accurate winner** recovers the failure within this neighbourhood.
- **A higher-scoring but less accurate winner** is a scoring or model finding, not a repair.
- **Truth absent:** the proposal stage limits recovery. The 4 is never changed after seeing the answer.
- **Equal saved fragment and exposure terms** are score ties on these data, not proof of fundamental
  indistinguishability.
- **Not run here:**
  - the distance-prior representation issue stays on the roadmap; its matched-hazard scale was fixed for the
    subdivision test, and recomputing it after repartitioning needs another assessment;
  - the flank-truncation ablation.

## 7. Units and budget

**Units:**

| unit | declared timeout |
|---|---|
| `cat.<cell>` | 60 s |
| `score.<cell>` | 180 s |
| `rank` | 30 s |
| `eval.<cell>` | 60 s |

**Limits:** 600 s cumulative, counted as spent time plus the declared timeout. At most 4 threads, 6 GiB, one process at a
time.
