# M15-D: the basic short-read command (frozen before execution)

**What this is:** the ordinary `panvar genotype` command made usable without a profile. It covers the basic interface only.

**What this is not:**
- no filter redesign or advanced-help reorganisation;
- no model tuning and no new search algorithm;
- no public persisted-index command;
- no reserved-family run;
- no merge, push, tag or cleanup.

**Logistics:**
- **Branch:** `m15-d-basic-command`, from `m15-c-interface-review` `c0e2dea`.
- **Results:** `tests/results/m15d_basic_command_20261021/`.
- **Budget:** 2,400 s, including 300 s of disclosed development (see `m15d_run.py`).
- **Resources:** 4 threads, one process at a time, 6 GiB.

**Normal usage:**
```
panvar genotype --graph locus.gfa --reads1 sample.R1.fq.gz --reads2 sample.R2.fq.gz --threads 4 --out sample
```
`--reference-path` is needed only when the graph does not name a single reference path. That happens when there is no
`RS:Z` header, or when the `RS:Z` sample has several paths. The refusal lists the candidates.

## 1. Before / after defaults (an ordinary command with no profile)

"Before" is `c0e2dea` (M15-C). "After" is this branch. A flag or profile still overrides every row.

| setting | before (no profile) | after (no profile) | kind |
|---|---|---|---|
| union, adaptive, gap and census backends (`union_scoring_backend`, `adaptive_window_discovery`, `indexed_union_store_layout`, `indexed_union_evidence_layout`, `adaptive_fused_seed`, `combination_prepared_library_terms`, `gap_recursion_backend`, `dosage_census_backend`, `marker_proposal_reuse_table_bytes`) | reference routes | the release candidate's accelerated routes | **exact backend** (gated identical on completed runs, I/S batches and M10-A) |
| `union_fragment_blocks`, `repair_round1_origin_admission` | off / existing | `fine_core_4096_experimental` / `reuse_only_experimental` | **exact route** (S2-F/G, S2-J) |
| cache and store budgets (`emission_cache_bytes` 512 MiB, `indexed_union_store_bytes` 2 GiB, `adaptive_fine_store_bytes` 512 MiB, `repair_origin_cache_bytes` 768 MiB, `union_placement_cache_bytes` 512 MiB, `union_origin_cache_bytes` 1 GiB, `union_prefilter_cache_bytes` 128 MiB, `start_sort_scratch_bytes` 320 MiB) | mostly 0 (repair origin 384 MiB) | the release candidate's values | retained memory for reuse; not a work cap |
| consumption capacities (`max_union_cached_*_consumed`, `adaptive_fine_max_served_placements`) | built-in | 2e9 each (release candidate) | bounded served work of the reuse routes |
| `max_state_classes` | 2e8 | 6e8 (release candidate) | executed-work cap (typed refusal when it binds) |
| `footprint_limit_bytes` | 0 (none) | 6 GiB | total-footprint guard (sampled; exit 125, or 3 PROVISIONAL with a checkpoint) |
| `time_limit_seconds` | 0 (none) | 0 (none) | unchanged: no time limit is imposed by default. The tested profiles' 900 s was a campaign watchdog |
| `repair_max_neighbours` | 670 (CYP2D6 and ACOT needed profiles: 958, 1004) | `model`: the model's own structural neighbourhood maximum, `Σ_blocks 2·min(repair_local_pairs, a(a+1)/2)`, recorded beside the bound | **locus-independent policy** (see §2) |
| `repair_junction_proposals` (JP1.1) | off | `experimental_v1` | **candidate generation**: can change predictions (M10-A changed one candidate pool and one conditional GQ) |
| `early_incumbent` | off | `experimental_v1` | budget-stop behaviour: completed runs identical; a stop keeps a PROVISIONAL checkpoint |
| `output_schema` | legacy | `sequence_first_v2` | the report contract (exit 0 / 3 / 1) |
| `compound_annotation` | off | `contract_v2` | report |
| `timing_detail` | 0 | 1 | report |
| dosage | off (legacy); **audit implied by v2** | **off unless `--dosage-check`** (v2 no longer implies dosage); reported NOT_EVALUATED (`dosage_mode_off`) | diagnostic only; never changes the prediction |
| `--dosage-mode repair` | legacy only (refused under v2) | unchanged: legacy interface only (`--output-schema legacy --dosage-mode repair`) | |
| library | omitted = auto | unchanged: omitted = auto (identical) | estimated: insert mean/sd (decisive anchored pairs, ≥ 500), fragment-start intensity (marker depth fit); fixed: `support_sigmas` 4, `substitution_error` 0.002, `outlier_weight` 0.001, gap model. Refusal is `LIBRARY_REFUSED`, never a fallback. The 0.9×/1.1× intensity sensitivity audit runs under auto |
| profile | needed for the tested configuration | optional expert override (flag > profile > preset > built-in) | |

**Unchanged** (already in the preset): `repair_rounds` 32, tie branching, representation guard, `kmer_size` 31,
`syncmer_s` 15, `min_variant_bp` 50, `mismatch_band` 6, `switch_probability` 0.01, `proposal_state_limit` 256,
`phase_audit_candidates` 0.

**Not in the preset, so the user decides:** a time limit, the dosage stage, the reference, the panel curation, the
candidate count and diagnostics.

## 2. Policies

**No locus dispatch.**
- The preset is one list. Nothing reads a locus name, a repository profile or an archived result at run time.
- **Neighbourhood bound (`model`):** the locus-specific caps (C4 670, CYP2D6 958, ACOT 1004) existed only to admit each
  locus. Each equalled or exceeded the model's own maximum (670 / 958 / 998).
  - The model maximum is a closed-form function of the block allele counts and the preset's `repair_local_pairs` (32), so
    every round's neighbourhood stays bounded by it.
  - The bound is recorded (`repair_neighbour_bound`, `repair_max_neighbours model`).
  - A numeric cap from a flag or profile still applies, and a neighbourhood above it still no-calls `neighbourhood_bound`.
  - The other work caps and the footprint guard are unchanged: no safety cap is removed.

**Exact backends vs candidate generation:**
- **Exact backends:** routes gated identical on completed runs.
- **Candidate generation:** JP1.1 can change predictions. It is included because the release candidate (M15-B) was
  evaluated with it, and it is recorded as its own row.
- **Budget-stop behaviour:** `early_incumbent`.

**Preset features whose requirements a user removes:**
- With `--adaptive-repair off`, `--repair-rounds 0` or `--output-schema legacy`, the preset-sourced `output_schema` (→
  legacy), `early_incumbent`, `compound_annotation`, `repair_junction_proposals` and `union_fragment_blocks` (when its
  requirements fail) are switched off.
- They are recorded in `.settings.tsv` as `preset … (inactive: <reason>)`.
- The same features set explicitly by a flag or profile keep their refusals.

**Dosage:**
- `--dosage-check` is the existing audit (`dosage_mode audit`), recorded with source `flag --dosage-check`.
- It is consistent with an explicit `--dosage-mode audit`.
- It contradicts an explicit `--dosage-mode off|repair` or `--dosage-audit off`. Contradictions are refused (exit 2) before
  any output.
- Dosage is not a profile key.

## 3. Tests (frozen; development controls only)

**Inputs:**
- the M11-D development cells: C4 NA19705, CYP2D6 HG02769 and ACOT HG03209, at 5×, with the M15-C graphs, curations and
  reference paths;
- the declared libraries are the simulations' own.

**No reserved family is used.**

| unit | command | frozen expectations |
|---|---|---|
| `build`, `cli_contract` | every target; the CLI contract test | build passes; the contract test checks the preset keys and sources, the `model` bound, preset deactivation, `--dosage-check` contradictions (both orders), repair only under legacy, and the help text |
| `ctest.full` | the full enabled suite, serial, 900 s per test | every enabled test passes. Two tests were updated for the new defaults: `curation_contract` is pinned to the legacy interface it tests; the S2-G route test declares its unblocked control and treats the exact caches' served-work counters as volatile |
| `ex.invalid` | `--dosage-check` with `--dosage-mode off` (both orders), with `--dosage-audit off`, and with the legacy repair; `--dosage-mode repair` under the default v2; explicit v2 with `--adaptive-repair off`; no `--reference-path` on a graph without `RS:Z`; `--preset long-read`; a missing library file | exit 2; the message names the problem; a sentinel file already at the prefix is unchanged and nothing else is written. For the missing library file (an input error found when the run reads it), only exit 2 and the message are checked, and the prefix state is recorded |
| `ex.noprofile_{c4,cyp2d6,acot}` | no profile, no `--library`, no dosage flag | the exit matches the reported outcome (0 PREDICTION / 3 PROVISIONAL / 1 NO_PREDICTION); every accelerated key is recorded with source `preset short-read (default)`; no time limit; nothing from a profile; `repair_max_neighbours model` with its bound; the library is estimated (`auto: insert from …`); dosage off and NOT_EVALUATED (`dosage_mode_off`); for a prediction, the v2 FASTA rule (panel pair: no FASTA and extraction references verified; off-panel: one FASTA with both records matching the SHA-256; never the legacy FASTA). **Recorded, not checked:** the estimated parameters, the library sensitivity outcome, filters, runtime and peak footprint |
| `ex.auto_equivalence_c4` | `--library auto` | identical to the omitted-library run: prediction, library estimate and sensitivity tables, sequence and block tables, every effective setting; library status rows equal apart from timing |
| `ex.declared_{c4,cyp2d6,acot}` | no profile, `--library FILE --dosage-check` | exit 0; the prediction (rank-1 score, block paths, haplotype SHA-256 and origin, status, fingerprint), `.sequence_haplotypes.tsv`, `.equivalence_set.tsv` and FASTA are identical to the archived M11-D accelerated run; settings differ only in `time_limit_seconds`, `block_diagnostics`, `repair_max_neighbours` and the absent built-in `repair_joint_search off`; the dosage status is equal apart from timing, with source `flag --dosage-check`. Runtimes are recorded |
| `ex.dosage_absent_{c4,cyp2d6,acot}` | the same as `declared`, without `--dosage-check` | identical prediction, sequence tables, FASTA, scores and candidates; without the flag, `dosage_mode off` and DOSAGE_FILTER NOT_EVALUATED (`dosage_mode_off`); with it, the audit COMPLETE and the filter evaluated from it; GENOTYPE_FILTER identical; the audit adds only `.dosage*` diagnostics and no prediction FASTA |
| `ex.unrepaired_view` | C4 `--adaptive-repair off` | runs (exit 0 or 1) under legacy; `output_schema`, `early_incumbent`, `compound_annotation` and `repair_junction_proposals` are off and recorded inactive; the exact backends stay selected; no sequence report |
| `ex.auto_refused` | C4 with only the first 300 read pairs, auto library | `LIBRARY_REFUSED`, exit 1, no prediction, no FASTA, the estimate written |
| `ex.resource_stop` | ACOT with an expert override profile containing only `time_limit_seconds 8` | exit 3 (PROVISIONAL) or 124, never a completed prediction; the override is recorded from the profile and everything else from the preset |

**Failure handling:**
- A failed check is a result. It is recorded and never retuned.
- An amendment is declared before any rerun.
- If an added audit exhausts resources, it is reported separately.

## 4. Deliverable

- the implemented basic command;
- this before/after table;
- the tested examples (with the estimated parameters, outcomes and runtime differences);
- the remaining limitations.

The batch stops for review. The next steps are a filter review, then the advanced help.
