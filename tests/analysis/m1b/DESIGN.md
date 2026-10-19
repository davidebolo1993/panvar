# M1-B: structural-error attribution on the frozen current caller

Branch `m1-b-structural` from `233cb65`. The performance baseline and every scientific setting stay frozen; no parameter
tuning, distance-dependent prior, clustering, confidence change or production integration. Results:
`/Users/davide.bolognini/panvar/tests/results/m1b_structural_20261019` (ledger, manifests, runs).

## 1. Audit set (frozen before execution; `AUDIT_SET.tsv`, `SELECTION.tsv`, `freeze_audit_set.py`)

| case | locus | sample | cov | target | why |
|---|---|---|---|---|---|
| c4_10 | C4 | NA19650 (LOO) | 10 | sv:7 | the selected sample at 10x |
| c4_30 | C4 | NA19650 (LOO) | 30 | sv:7 | largest retained sv:7 excess at 30x |
| cyp_5 | CYP2D6 | NA19650 (LOO) | 5 | sv:9 | the known sv:9 error (M1-A) |
| cyp_30 | CYP2D6 | NA19650 (LOO) | 30 | sv:9 | its nearly correct 30x twin |

**C4 selection rule.** Frozen package `genotype_c4_production_11fe681_dosage_repair/final_26401b4/package/blocks.tsv`:
experiment `loo`, block `sv:7`, coverage 30, `recommended == True`; sort by excess descending, then subject ascending
(tie-break). 45 rows match. Top: `loo.NA19650.30`, excess 65,444 bp (distance 65,480, floor 36); next `loo.HG03807.30`,
13 bp. No tie. (`fresh_loo` has no rows.) In that package the dosage repair was not applied at this block
(`original_excess = excess`), so the package value equals the raw 11fe681 call.

**Inputs verified** (existence and SHA-256 in `AUDIT_SET.tsv`; no gap; nothing substituted or rebuilt):
- **C4.** Model `graphs/NA19650.loo.gfa` (the LOO graph without NA19650); truth = the NA19650 haplotype paths in
  `graphs/c4.gfa`; no curation; reference `HG00096#1#haplotype1-0000024:31848049-32074186`. Reads: simulated from the
  truth paths in the held-out coverage experiment and subsampled to 10x and 30x of 30x (`library.txt` header); each
  case's own `library.txt`.
- **CYP2D6.** Model `data/graphs/loo.NA19650.gfa`; truth = NA19650 paths in `data/graphs/cyp2d6.gfa`; curation
  `tests/genotype_cyp2d6_curation.tsv` (sha256 `25b8bc4f…`, equal to the pilot's curation of record); reference
  `grch38#1#chr22:42031864-42245566`. Reads: `genotype_probe simulate` seed 20260921 salted with NA19650 at 30x; 5x is a
  subsample of it.
- **Historical outputs.** C4: `loo/NA19650/call.{10,30}.*` with `.eval.{summary,blocks}.tsv`. CYP2D6:
  `data/runs_a5/loo.NA19650.{5,30}.original.*` (pre-dosage-repair outputs) with `data/eval_a5/…{5,30}.raw.*`.

## 2. Settings: fresh run vs historical vs frozen baseline

Fresh composition (`m1b_run.py`): the baseline cell's base profile, mode `full` with
`profiles/genotype_dev_performance_v1.tsv`, overrides `time_limit_seconds` = timeout − 20 s and
`footprint_limit_bytes` = 6 GiB, plus `early_incumbent experimental_v1`, `union_fragment_blocks
fine_core_4096_experimental`, `repair_round1_origin_admission reuse_only_experimental`; `--preset short-read`,
`--threads 4`, `--output-schema sequence_first_v2`, no `--dosage-mode`. Each case's own graph, reads, `library.txt` and
curation. `outcome.py` lists every difference between the fresh `call.settings.tsv` and the M1-A frozen baseline
settings of the matching baseline cell (C4 → `c4_30`; CYP2D6 → `cyp_5`/`cyp_30`).

**C4 historical (11fe681, `genotype_c4_v6_campaign_local.tsv`, no preset, 8 threads).** The v6 scientific values equal the
short-read preset and the code defaults now in force (k 31, s 15, min_variant_bp 50, switch 0.01, proposal_state_limit
256, max_mosaic_candidates 128, phase_audit_candidates 0, mismatch_band 6, repair_rounds 32, repair_local_pairs 32,
repair_max_neighbours 670, tie branching 1 with 16 branches, representation_guard 1). Its cache budgets and consumption
caps equal the overlay's (`p1_run.selfcheck`: overlay caches = `perf_bounded.C4V6_CACHE`). Differences:

| setting | historical v6 | fresh (= baseline) | kind |
|---|---|---|---|
| `max_physical_origins` | 700,000,000 | 500,000,000 (default) | work cap: can only refuse, never change a completed result |
| `footprint_limit_bytes` | 4 GiB | 6 GiB | resource guard |
| `time_limit_seconds` | 3,600 | timeout − 20 s | resource guard |
| threads | 8 | 4 | operational |
| scoring backends, union blocks, round-1 admission, S2-A, output schema | 11fe681 code | 233cb65 frozen route | computational route (S2-L: scientifically identical on every completed pair) |
| `library_sensitivity_audit` | absent (later feature) | 1 (default; declared library → NOT_APPLICABLE) | audit only |
| dosage | none in 11fe681; package = 26401b4 repair replay, not applied at sv:7 | not requested (sequence_first_v2: audit-only) | reporting |

The C4 baseline cell (`c4_30`, HG00096) and these cases share the same base profile; HG00096's inputs are not used.

**CYP2D6 historical (pilot a5 binary, `--preset short-read --profile resource_profile_a4.tsv --threads 4 --dosage-mode
repair`).** a4 sets `repair_max_neighbours 958` and `max_state_classes 600000000`, exactly as the baseline `cyp_*` base
profiles. Differences: footprint 4 → 6 GiB; time 2,700 s → timeout − 20 s; `--dosage-mode repair` → not requested (the
historical comparison uses the pre-repair `.original` outputs and the `raw` evaluation); code revision.

## 3. Outcome rules (`outcome.py`, declared before execution)

Target block: **PERSISTS** = same unordered called allele pair as the historical call; **DISAPPEARS** = fresh excess
≤ 100 bp; **CHANGES** = otherwise; **NOT_ASSESSED** = the caller or its evaluation did not complete. Allele indices are
compared only when the model fingerprints are equal (checked). A disappeared error is reported and not diagnosed from the
archived prediction. A time or resource stop is reported as such (S2-A checkpoint PROVISIONAL), never as a call, and no
cap is raised.

## 4. Structure (`block-sequences` + `structure.py`; truth only in this diagnostic evaluator)

`genotype_probe block-sequences` (new, diagnostic) writes every allele of the target block with its carrying templates,
and the donor's truth haplotypes projected onto the same block. `structure.py` counts repeat-unit copies from bytes, not
lengths: anchors are 31-mers sampled every 25 bp from the shortest template-carried allele ≥ 10 kb, unique in it and
conserved in ≥ 90 % of carried alleles. An anchor is informative when its copy count varies between carried alleles.
Copies = the median count over informative anchors, with support = the share at that count; fewer than 50 informative
anchors is LOW_INFORMATION. Totals separate **dosage** (copies hap1 + hap2) from **allocation** (the unordered pair):
same total with a different pair = ALLOCATION_ERROR; different total = DOSAGE_ERROR.

## 5. Alternatives and score competition (declared as an amendment after the outcome; ≤ 12 per persistent case)

Only for an error that PERSISTS or CHANGES. Order of sources:
1. **saved evaluated alternatives of the fresh run** (`mosaic_scores.tsv` / `mosaic_candidates.tsv` /
   `sequence_alternatives.tsv`): the reported candidate and the best-scoring saved candidate of each structural class at
   the target block (per-haplotype unit copies), with their saved totals as expectations;
2. then declared single-block substitutions of the reported diplotype at the target block (panel alleles with the true
   unit copies, chosen with truth **in this evaluator only**), each in both homologue assignments.

Never a whole-locus oracle assembled from independent block floors. Scoring: `genotype_probe score-set` (the unchanged
reference scorer) with `--run` = the fresh run (fingerprint checked); every saved expectation must reproduce within
|Δ| ≤ 1e-6 or the unit refuses. Saved totals are printed to 1e-6, so this is tolerance reproduction, not bit identity.
Decomposition per alternative: fragment term, exposure term, sequence-pair prior, total; numerical uncertainty = the
fresh run's `numeric_error_bound`; `--decompose` against the reported candidate gives per-fragment classes (block
interior / junction / outside / one-sided).

Interpretation: a near-equal total means weak observed discrimination under this model and these reads, not proof of
non-identifiability. Whether the fresh caller proposed and evaluated a structurally appropriate alternative is read from
its saved candidates (membership by unit-copy class), not inferred.

## 6. Budget and order

1,200 s cumulative (development runs included), one process at a time, ≤ 4 threads, ≤ 300 s per unit, 6 GiB envelope.
Order: input check → C4 10x, 30x → CYP2D6 5x, 30x (each caller then its evaluation) → block sequences → structure →
outcome. Manifests are frozen and committed before execution; amendments are declared before further execution; failed
units stay recorded.
