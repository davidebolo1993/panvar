# M9-D: isolated attribution of the M9-C ACOT high-confidence errors

**Question:** do these errors come from missing candidate combinations, from the copying prior, or from the
fragment/exposure objective?

**Not done:** recalibration, tuning, caller search, residual editing, clustering, confidence fitting, filter change,
promotion, merge, push or cleanup.

**Results:** `tests/results/m9d_attribution_20261020/`.

**First step, before any execution:** two qualifications are appended to `tests/analysis/m9c/RESULT.md`:
- excluding sv:5 after inspection does not establish calibration;
- the actual sv:5 allele lengths are 7–50 bp.

All selected cells are now development evidence.

## Cases (frozen)

| case | cell | block | kind |
|---|---|---|---|
| E1 | LOO d1 (NA19909) | sv:5 | floor-zero high-confidence error |
| E2 | LOO d2 (NA18747) | sv:5 | 〃 |
| E3 | LOO d3 (HG02178) | sv:5 | 〃 |
| E4 | LOO d3 (HG02178) | sv:9 | 〃 |
| E5 | LOO d4 (HG00639) | sv:5 | 〃 |
| E6 | LOO d6 (HG03831) | sv:5 | 〃 |
| C1 | LZO d1 (NA19909) | sv:5 | control |
| C2 | LZO d2 (NA18747) | sv:5 | control |

No other donors, reads or coverages.

## Matched comparisons (`m9d.py construct`)

**Contexts for errors E1–E6:**
- **A:** the reported winning sequence pair (`locus_genotype` block paths, `total_score`).
- **B:** the highest-scoring evaluated hypothesis that carries the exact truth pair at the block (`mosaic_candidates`
  and `mosaic_scores`).

**Within each context:**
- the block is set to the called pair and to the exact truth pair;
- each pair is scored in both homologue assignments (`ab`, `ba`);
- every other block is unchanged.

**Controls C1 and C2:**
- **context:** the LZO winner, whose sv:5 call is correct;
- **compared:** the correct call against the LOO-called sv:5 pair of the same donor;
- **mapping:** LOO alleles map into the LZO model by **exact sequence** (one hit required), never by numeric id.

**Dedupe.** Identical unordered **sequence** pairs (SHA-256 of the spelled haplotypes) are scored once. Duplicates are kept
in `membership.tsv`.

**Status of the constructed pairs.** Truth supplies diagnostic alternatives only: they are not deployable predictions and
never enter caller outputs. Every constructed pair is checked for membership in the run's scored pool (by source: `panel`,
`repair`, or `mosaic` = the 128-candidate proposal shortlist) or absence from it.

## Identity gates

| gate | how |
|---|---|
| Model fingerprints | `models`: the relinked `cyp2d6_block_alleles` dumps the LOO and both LZO models with each run's curation; each dump's fingerprint must equal every corresponding run's `model_fingerprint` |
| Replacement bytes | `construct`: for every constructed haplotype, the replaced block equals the replacement allele's sequence, and the flanks are byte-identical to the context's; offsets may shift, and the shift is recorded |
| Context reproduction | `score`: `score-set --candidates` with each context's saved total as the expectation. A total that does not reproduce within 1e-6 refuses the case (exit 4). `analyse` also compares fragment, exposure and prior with the saved `mosaic_scores` row within 1e-6 |
| Run settings | same model, curation, reads, library, `--switch-probability 0.01` and `copying_transition` (the run's transition-vector digest is checked by score-set) |

**If reproduction fails, that case is not interpreted.**

**Degeneracy inventory** (`degeneracy`, model-only, non-gating): the exhaustive search for cross-block equivalences and
within-block duplicates in the three models (`degeneracy_inventory.py`). It records whether label identity implies sequence
identity there. Pool membership is checked by labels.

## Decomposition (`analyse`)

**For each matched contrast,** truth (or inserted) minus called within the same context and every assignment pairing, it
reports:
- **Δfragment, Δexposure, Δlikelihood, Δprior** (the full sequence-level copying prior, `log_unordered_sequence_pair_prior`)
  and **Δtotal**;
- **a reconstructed numerical bound.** It applies the caller's per-candidate formula, (fragments + 16)·ε·(Σ|terms| +
  |prior|), floored at 0.5e-9, to each side and sums them; Σ|terms| comes from the per-fragment table. A contrast is
  certified if |Δtotal| > bound.

**Per-fragment differences** come from `score-set --decompose A,B --decompose-block <block>`:
- the scorer's own masses;
- placement classes inside the block, across its junction, or outside;
- whether `placement_evidence_interpretable` holds;
- concentration: the top-5 and top-20 share of |Δ|, and the number of fragments with |Δ| > 1 nat.

**Reading** (one rule per contrast, as frozen):

| finding | reading |
|---|---|
| wins total and was not evaluated | search/proposal gap in that context |
| wins total and was already scored | inconsistency; investigate first |
| wins likelihood, loses total | prior effect, quantified and not retuned |
| loses likelihood | inspect the evidence; not by itself a bug |
| within the bound | tied |

**Context dependence** is read by comparing contexts A and B. No universal block-local preference is claimed. No
attribution is made to sequencing error, placement or a prior defect without supporting evidence.

## Execution

- **Budget: 900 s,** counted as spent + the unit's timeout + 30 s (the analysis).
- **Order:** `build`, `models`, `construct`, then `score` for E1–E6, C1, C2 (90 s each), then `degeneracy`, then `analyse`.
- **One process at a time, four threads,** existing resource settings unchanged.
- **Refused units** are recorded NOT RUN.

## Amendment 1 (declared after `analyse`, before any further execution)

**Trigger.** The non-gating degeneracy inventory found **2,045 cross-block equivalences**, all within **sv:5..ctx:6**, in
all three models. There are no within-block duplicates. So label identity at sv:5 does not determine the sequence over
sv:5 + ctx:6.
- In d1 and d2, ctx:6 already agrees with truth (evaluator distance 0).
- In d3, d4 and d6 it does not (distances 4, 6 and 5).

**Question.** Whether each reported haplotype's sv:5 + ctx:6 sequence equals or differs from the truth's. This separates
representation ambiguity from a real sequence difference.

**Added units, read-only and truth used diagnostically:**
- `interval.<cell>` for LOO d1, d2, d3, d4 and d6: the existing `genotype_probe block-sequences --blocks sv:5,ctx:6,sv:9
  --truth-gfa acot.gfa --donor <donor>` (truth projection bytes; LOO model and curation), 60 s each;
- `intervals`: per haplotype, the reported sv:5 + ctx:6 bytes (sv:9 alone for E4) against the truth projection's, under
  the better pairing. It reports identity and the residual differing core after trimming common prefix and suffix (not an
  edit distance), 30 s.

**Unchanged:** the matched contrasts, gates and readings.
