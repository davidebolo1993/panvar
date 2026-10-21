# M12-B: bounded, offline dosage-first selection experiment (frozen before execution)

**What this is:** an experimental point-selection decision rule.

**What this is not:** a posterior, a production dosage repair, a confidence model or a copying-prior change.

**Status of earlier batches:** M12-A is closed, with its K1 qualification appended to its RESULT.
- **Branch:** `m12-b-dosage-first`, from `m12-a-rate-exposure-audit` `dadd5f2`.
- **Results:** `tests/results/m12b_dosage_first_20261021/`.
- **Budget:** 1,800 s, including 200 s of disclosed development: code reading, two probe builds (about 50 s of compute), and
  selftest iterations.
- **Resources:** 6 GiB; one scoring process at a time; 4 threads.

**Not done:**
- no new reads, simulation, full caller rerun, production integration, cap increase, repaired FASTA, prior change, fit,
  clustering, merge, push or cleanup;
- the original calls and every earlier output are untouched;
- the 18 ACOT validation families are untouched.

## 1. Scope (frozen manifest: `MANIFEST.tsv`, from `m12b.blocks_manifest()`)

**All 29 M11-D `repeated_unit` blocks** that are predicted and truth-defined, in M11-D cell order:
- CYP2D6 sv:9, 6 blocks (B01–B06);
- C4 sv:7, 12 blocks (B07–B18);
- ACOT sv:13, 11 blocks (B19–B29; ACOT NA21106 is a caller no-call and is not a block here).

The catastrophic blocks are B04 (K4), B06 (K5), B17 (K1), B18 (K3) and B21 (K2).

**Candidates:** only the existing active-model alleles at the target block, placed in the reported winner's fixed context
(rank-1 representative of `equivalence_set.tsv`). No other block changes.
- No truth-derived construction, discovery, recombination or post-audit reconstruction.
- The spelled sequences are diagnostic scoring inputs only.

## 2. The frozen rule

**1. Copy state.**
- An allele's copy state is its copies of each **pooled copy-count feature** of the block. These come from
  `panvar::dosage::pooled_features`: the dosage audit's own graph-only census (default `DosageSettings`: k = 31, at least
  500 k-mers per feature).
- A pair's state is the per-feature sum.
- The probe exposes it as `count-attrib --features 1` (diagnostic; it also writes the feature k-mers).
- **No new structural alphabet.** A block with no pooled feature, or with no count loss, is **UNSUPPORTED**.

**2. Count loss.** The unchanged M6-A/count-attrib objective:
- loss = Σ (observed / depth − predicted)² over the active-panel 31-mer dictionary;
- depth = the mean count of the stable-copy k-mers / 2;
- every **ordered** allele pair at the block in the fixed context is scored (count-attrib's catalogue, real-read profile);
- feature inventory, normalization and scope are all unchanged.

**3. State selection.**
- Each state's best count loss is computed, and the state(s) with the lowest are selected.
- Exact ties are retained; a tie means equal printed loss at count-attrib's 12 significant digits.
- No probability or interval is formed.

**4. Within the selected state(s).**
- **Every** ordered pair is scored: both homologue assignments, nothing shortlisted.
- Scoring uses the **unchanged** `score-set` objective (fragment + exposure + copying prior).
- The highest total is the selection. Ties within 1e-9 are reported as ties and each one is evaluated; none is broken.

**5. Not done:**
- count loss is never added to the score;
- overlapping evidence is never multiplied as though independent;
- no dosage weight is tuned.

**6. Content sensitivity** (recorded):
- state margin = (the runner-up state's best loss) − (the selected state's best loss);
- the block is flagged **CONTENT_SENSITIVE** when that margin is smaller than (median − min) of the losses inside the
  selected state. In that case a typical candidate of the selected state would lose to the runner-up state, so the count
  preference is not a pure copy-number measurement.

**7. Bounds.**
- **Preflight:** score-set on the first block of each locus, with the incumbent plus the first 24 alleles in the hap1
  context. It measures seconds per unique haplotype and the fixed overhead; no outcome is read.
- **Estimate:** overhead + 1.5 × (unique haplotypes in the selected state(s)) × seconds per haplotype.
- **Admission:**
  - a block is scored iff its estimate is ≤ 240 s, with declared timeout ceil(1.5 × estimate + 60);
  - otherwise it is **LIMITATION_STATE_TOO_LARGE_TO_SCORE_EXHAUSTIVELY** — never a shortlist;
  - the runner's global admission (spent + timeout ≤ 1,800) applies, in manifest order;
  - a block it refuses is a recorded budget limitation.

## 3. Gates (before any truth comparison)

| gate | check |
|---|---|
| `gate.selftest` | the selection logic on synthetic data: total-state grouping (3+1 and 2+2 share a state), both assignments, exact ties retained, order invariance, missing features or losses → unsupported, multi-feature sums |
| `gate.order` | score-set totals and fragment-term bits are identical with the candidate list reversed (the CYP2D6 preflight list) |
| `gate.features` | the probe's pooled features equal the caller's own `dosage_feature_copies.tsv` wherever the audit wrote one. On real C4 alleles a 3+1 and a 2+2 pair share one total state, so total dosage alone cannot resolve allocation |
| per scored block (evaluate) | incumbent total reproduced (score-set refuses otherwise) and incumbent distance = the evaluator's; every scored haplotype's length = flanks + the declared allele, with full FNV identity for the incumbent and the first and last candidates; candidates the run itself scored carry the same total (≤ 1e-5); aliased alleles (identical bytes) give identical totals |

## 4. Evaluation (truth opened only in `evaluate`, after every selection is frozen in `select` and the score units)

**Per block:**
- original and selected states;
- count losses: incumbent, selected best, state margin, content sensitivity;
- full score components against the incumbent;
- the selected pair(s): distance to truth (exact edit distance, minimum over pairings — the evaluator's convention, checked
  on the incumbent), the existing floor, excess and normalized excess;
- the truth's per-feature copy state (median count of each feature's k-mers in the truth block). It is assessable when
  ≥ 90% of a feature's k-mers sit at that integer median;
- total-copy recovery, allocation match (per-haplotype states), coverage, time and memory.

**Allocation, sequence identity and phase are kept separate.** Total-copy recovery does not establish correct haplotypes.

**Outcome per block:**
- original ≥ 1,000 and selected < 1,000: `FIXED_LARGE`;
- both ≥ 1,000: `RETAINED_LARGE`;
- original < 1,000 and selected ≥ 1,000: `INTRODUCED_LARGE`;
- otherwise `IMPROVED`, `UNCHANGED` or `WORSENED`.

A near-floor control is a block with original excess ≤ 10; it is damaged when the selected excess is > 10.

**Decision (frozen):**
- **`MERITS_DEFAULT_OFF_CALLER_EXPERIMENT`:** at least 3 catastrophic blocks scored, at least 3 fixed, 0 introduced large
  errors and 0 damaged near-floor controls;
- **`FAILS`:** at least 3 catastrophic blocks scored and (0 fixed, or introduced ≥ fixed);
- **`INCONCLUSIVE`:** otherwise.

**Status of the evidence:** inspected development evidence, not validation or calibrated confidence.

**No post-hoc changes** to features, loss, copy-state rules or tie handling.
