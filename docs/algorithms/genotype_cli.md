# `panvar genotype`: user-facing command contract and presets (normative design, 2026-09-28; partly implemented 2026-09-29)

**Implemented at `c758ae2`:**
- **`--preset short-read`:** §3 kind 2, only the production-algorithm and technology keys.
- **`--preset long-read`:** refused as not implemented.
- **`--library auto|FILE`:** §4, and `library_estimate.md` §5.
- **The normal-user flags:** `--graph`, `--out`, `--mosaic-candidates`, `--adaptive-repair on|off` and
  `--repair-rounds`.
- **Declared-reference resolution:** an `RS:Z` header tag, or a refusal that lists the candidates.

**Implemented since:**
- **`5b1f059`:** an omitted `--library` means `--library auto`.
- **`headline_filters_v2`:**
  - the advanced scientific flags (`--help-advanced`), each exactly its profile key;
  - per-setting source records in `<prefix>.settings.tsv`;
  - the three headline filters (`genotype_headline_filters.md`);
  - an omitted `--preset` means `--preset short-read` (the §2 default), recorded as `preset short-read (default)`;
    the unrepaired development view needs `--adaptive-repair off`;
  - `--phase-audit-candidates` is a development-only control, refused under adaptive repair.

**Not implemented:**
- **The resource predictor** (§3).

**Existing flags keep working unchanged.**
- `--library FILE` reproduces the frozen calls unchanged (gated).
- The block and locus tables gain trailing headline columns.
- The repair-mode `dosage_filter` component is renamed `dosage_stage_filter`. Frozen historical tables (the sealed
  campaign, replay and packages) keep the old name.
- **One intended behaviour change:** a command with neither `--preset` nor a profile used to give the unrepaired
  development view; it now runs the production algorithm. Profile runs and explicit `--preset short-read` runs are
  unchanged.

## 1. The basic command

    panvar genotype \
      --graph locus.gfa.gz \
      --reads1 sample_R1.fastq.gz \
      --reads2 sample_R2.fastq.gz \
      --out sample \
      --preset short-read \
      --threads 8

**The algorithm.** There is one production algorithm: panel-anchored adaptive completion.
- It scores every panel diplotype and up to K recombinant mosaics, repairs the incumbent by exact neighbour search, and
  branches exactly at improving-move ties.
- It adds the representation guard and the combined block quality.

**Not modes of the command:**
- **Leave-zero-out, leave-one-out and synthetic mosaic are evaluation scenarios.** They are made by preparing inputs:
  a graph with a donor removed (`genotype_probe loo-gfa`), or a staged graph with a curation manifest. There is no
  scheme flag, and the command cannot tell which scenario it runs.
- **No preset or flag names a locus.** There is no public C4 preset. The frozen C4 study is reproduced by its profile
  (§5), not by a preset.

## 2. Controls

**Normal controls, shown first in `--help`:**

| flag | values | default | meaning |
|---|---|---|---|
| `--graph FILE` | GFA / GFA.gz | required | the locus pangenome graph (today `--gfa`; kept as an alias) |
| `--reads1`, `--reads2` | FASTQ / FASTQ.gz | required | paired short reads |
| `--out PREFIX` | path | required | output prefix (today `--out-prefix`; kept as an alias) |
| `--preset NAME` | `short-read` | `short-read` | a technology preset (§3); `long-read` is reserved and refused: "long-read scoring is not implemented" |
| `--threads N` | ≥ 1 | 1 | workers; results are thread-invariant |
| `--mosaic-candidates K` | ≥ 0 | preset (128) | off-panel mosaics scored beyond the panel diplotypes (today `max_mosaic_candidates`) |
| `--adaptive-repair on\|off` | | `on` | `off` is the unrepaired development view, labelled as such in the status; never a production call |
| `--repair-rounds N` | ≥ 1 | preset (32) | the repair's round cap; reaching it is REPAIR_INCOMPLETE, never a truncated call |
| `--dosage-mode off\|audit\|repair` | | `off` | `genotype_dosage.md` §20–23. It stays `off` until loci other than C4 validate the dosage stage |
| `--library auto\|FILE` | | target `auto`, **not enabled** | `auto` estimates the sample quantities (§4) and takes the rest from the preset; `FILE` is the expert override. Until `auto` is reviewed, `FILE` is required as today |

**Advanced controls, in a separate help section.** Each is recorded in the status with its source (preset, profile or
flag).

| group | controls (current names) |
|---|---|
| markers | `kmer_size`, `syncmer_s`, syncmer offset, marker-frequency filters, declared `marker_background` |
| Li–Stephens / recombination | `switch_probability`, `proposal_state_limit` |
| read model | `mismatch_band`, `substitution_error`, `gap_open`, `gap_extend`, `outlier_weight`, `support_sigmas` |
| phase audit | `phase_audit_candidates` |
| repair search | `repair_local_pairs`, `repair_max_neighbours`, `repair_tie_branching`, `repair_tie_max_branches` |
| dosage | the dosage policy v2 constants (alpha, bootstrap replicates and seed, SUSPECT and UNDERPOWERED thresholds, minimum features), `--dosage-allocation-test` |
| resources | every `max_*` scoring cap, the `repair_max_*` caps, cache budgets (`*_cache_bytes`), `footprint_limit_bytes`, `time_limit_seconds` |
| inputs | `--reference-path` (the explicit coordinate-reference override), `--panel-curation`, `--profile FILE` |
| diagnostics | `--library-diagnostic`, the `test_*` and `diagnostic_*` hooks (never in `--help`'s normal section) |

**Coordinate reference.** By default the reference is the graph's declared reference sample: the rGFA stable-sequence
or reference-sense tags. With none declared, or several, the command refuses and lists the candidate paths. It never
picks one silently. `--reference-path` overrides.

**Precedence.** An explicit flag wins over `--profile`, which wins over the preset, which wins over the built-in
constants. Every value and its source is written to `mosaic_status.tsv`, as `overridden_settings` does today.

## 3. The `short-read` preset: technology, not locus

The preset names a technology's assumptions and nothing about a locus. Each setting is exactly one of four kinds:

| kind | settings | where it comes from |
|---|---|---|
| **1. Data-derived** (per sample) | read lengths; insert centre and spread where estimable; fragment-start intensity; marker depth per copy and background (already fitted from invariant anchors today); quality summaries | the reads, against the model's invariant anchors (§4), recorded with method and support |
| **2. Technology preset** (transfers between loci) | `substitution_error` 0.002; `gap_open` 1e-4; `gap_extend` 0.5; `outlier_weight` 0.001; `support_sigmas` 4; `kmer_size` 31; `syncmer_s` 15; `mismatch_band` 6; `switch_probability` 0.01; `proposal_state_limit` 256; `max_mosaic_candidates` 128; `repair_rounds` 32; `repair_local_pairs` 32; `repair_tie_branching` on; `repair_tie_max_branches` 16; `representation_guard` on; `phase_audit_candidates` 0 | declared constants of the preset |
| **3. Locus / model-derived** | blocks (`min_variant_bp` 50 is the block rule's threshold, a model parameter); graph paths; the marker census and anchors; panel size; the resource predictions | the graph and the model built from it |
| **4. Expert overrides** | every scientifically meaningful control above | flags, `--profile` |

**As implemented (`c758ae2`, `genotype_preset_settings`):**
- **Keys the preset sets:** `repair_rounds` 32, `repair_tie_branching` 1 and `representation_guard` 1, which switch on
  the production algorithm.
- **Keys it restates, equal to the built-in defaults:** `kmer_size` 31, `syncmer_s` 15, `min_variant_bp` 50,
  `mismatch_band` 6, `switch_probability` 0.01, `proposal_state_limit` 256, `phase_audit_candidates` 0.
- **Keys it leaves at the built-in defaults, never the C4 profile's values:** `max_mosaic_candidates`,
  `repair_local_pairs`, `repair_tie_max_branches`, every cap, every cache budget, and the footprint and time limits
  (no limit unless given). The status lists the preset's keys (`preset_keys`).
- **Library constants** (`genotype_preset_library`): support 4, substitution 0.002, outlier 0.001, the built-in gap
  model, and fragment-start intensity 0.05.
  - The intensity is a declared convention: λ = depth / (4 L) at 30× with 2 × 150 bp pairs. It is not inferred.

**Rules the preset keeps:**
- **No C4 resource cap enters `short-read`.** The frozen profile's caps, cache budgets, and footprint and time limits
  were sized on C4 (`genotype_mosaic_current.md`, v5/v6). A cap changes only when a run refuses, never what it calls.
  - The generic preset therefore declares none. Caps are model-derived (kind 3); a resource predictor from
    blocks, alleles, panel size and coverage is future work.
  - Until then a generic run is bounded by `footprint_limit_bytes`, defaulting to a fraction of physical memory, and
    by no time limit unless given.
- **No dosage assumption enters `short-read`.**
  - Dosage policy v2 (pooled features, α 0.001, bootstrap 200, thresholds) is the dosage stage's declared policy.
  - The dosage *default* is `off`. The frozen C4 short-read invocation sets `--dosage-mode repair` explicitly.
- **Which kind-2 values are C4-studied rather than transferable.** Several were chosen or confirmed on C4, and are
  marked so in the preset file until a second locus confirms them: `repair_rounds` 32, `max_mosaic_candidates` 128,
  `repair_tie_max_branches` 16, and the representation-guard bound ln 1000.
  - They are search budgets and conservative thresholds, not locus facts.
  - Reaching one is always reported. It never produces a silent truncation.

## 4. The library under `--library auto` (target; measurement in `library_estimate.md`)

**Estimated, with method, support and fallback status in `<out>.library_estimate.tsv`:**
- read lengths (reads);
- insert centre and spread (decisive paired placement against the model's invariant anchors);
- fragment-start intensity (anchor depth).

**Taken from the preset:** `substitution_error`, the gap model, `outlier_weight` and `support_sigmas`.

**Refusal instead of a guess.** The insert falls back, and the run refuses with a message naming `--library FILE`,
when any of these holds:
- accepted placements are below the minimum;
- the insert distribution is multimodal;
- the outlier fraction exceeds its bound.

A generic insert centre is never substituted. `--library FILE` replaces every value it declares, and the table records
the override.

## 5. The frozen C4 short-read invocation

This reproduces `final_26401b4`:

    panvar genotype --profile profiles/genotype_c4_v6_campaign_local.tsv --dosage-mode repair \
      --gfa <graph> --reference-path 'HG00096#1#haplotype1-0000024:31848049-32074186' \
      --reads1 … --reads2 … --library <declared library> --threads 8 --out-prefix <prefix> [--panel-curation <manifest>]

The profile is frozen (sha256 74a3723b…). The C4 freeze manifest names every hash.
