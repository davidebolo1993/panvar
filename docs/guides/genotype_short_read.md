# Short-read structural genotyping with `panvar genotype` (EXPERIMENTAL)

**This is the authoritative usage guide for the experimental short-read release candidate.**
- **Algorithm:** `docs/algorithms/genotype.md`.
- **Module overview:** `docs/modules/genotype.md`.
- **Status:** experimental. Nothing in the outputs is calibrated confidence.

**The ordinary command needs no profile:**

```
panvar genotype --graph locus.gfa --reads1 sample.R1.fq.gz --reads2 sample.R2.fq.gz --threads 4 --out sample
```

What this command uses:
- the short-read preset's accelerated configuration (§2);
- a library estimated from the reads;
- the sequence-first v2 report;
- no dosage audit unless `--dosage-check` is given.

Add `--reference-path NAME` when the graph does not name a single reference path, and `--panel-curation FILE` to remove
paths from the active panel. Profiles are optional expert overrides.

## 0. What is supported, and the evidence behind it

**The evidence was produced with one complete profile per locus.** The no-profile command reproduces those runs' predictions
exactly when given the same declared library and `--dosage-check` (M15-D `ex.declared_*`). The only settings that differ are
the profiles' 900 s watchdog, their locus-specific neighbourhood caps (now the model bound, §3) and diagnostics. The
profiles stay in the repository as an exact record:

| locus | profile (reproduces the evidence configuration) | evidence |
|---|---|---|
| C4 | `profiles/genotype_experimental_sr_v1_c4.tsv` | development only (M11-D: 12 cells at 5×) |
| CYP2D6 | `profiles/genotype_experimental_sr_v1_cyp2d6.tsv` | development only (M11-D: 6 pilot donors at 5×) |
| ACOT | `profiles/genotype_experimental_sr_v1_acot.tsv` | development (12 families) **and 18 reserved families** (M15-B: 18 of 18 completed predictions, median per-block excess 0, max 478 bp, 0 catastrophes ≥ 1 kb) |

**All of this evidence uses:**
- **simulated** paired 100 bp reads at 5×, with a **declared `--library FILE`** (the simulation's parameters), on graphs and
  curations described in each package;
- predictions made against a panel that excludes the sample (leave-one-out).

**Not established:**
- real reads;
- `--library auto` accuracy;
- other coverages;
- other loci;
- independent C4/CYP2D6 accuracy;
- calibrated confidence.

**A profile layer** (`profiles/genotype_dev_performance_v1.tsv`, `profiles/genotype_dev_m10a_v1.tsv`) **is not a complete
configuration.**

**The automatic library** is the default, but no accuracy evidence uses it. M15-D records its estimates and outcomes at the
three development cells; that is a descriptive record, not an accuracy result.

## 1. Required inputs and reference orientation

| input | flag | notes |
|---|---|---|
| locus graph | `--graph FILE` (alias `--gfa`, `-i`) | GFA with `P` paths. Every path defines the block frame, whether or not it is active |
| reads | `--reads1 FILE --reads2 FILE` | paired FASTQ, mate order preserved. Short reads only (`--preset long-read` is refused) |
| output prefix | `--out PREFIX` (alias `--out-prefix`, `-o`; default `genotype`) | every report is written at the prefix, and a new run removes run-owned reports there |
| coordinate reference | `--reference-path NAME` (alias `-r`; default: the single path of the graph's `RS:Z` sample; without `RS:Z`, or with several paths of that sample, the run is refused and the candidates are listed) | the **locus frame and orientation**. Predicted haplotypes are spelled in this frame. Extraction references give each panel path's orientation relative to it (`forward` / `reverse_complement`) |
| active panel | `--panel-curation FILE` (optional) | one path name per line, `#` comments. Named paths are **removed from the active panel** but still define blocks. An allele carried only by removed paths is uncallable. A block left with no active allele refuses the run (`MODEL_INCOMPLETE`) |
| library | `--library auto\|FILE` (omitted = `auto`, identical) | **`auto` estimates:** insert mean and sd, from at least 500 decisive anchored pairs; fragment-start intensity, from the marker depth fit. **`auto` keeps fixed:** `support_sigmas` 4, `substitution_error` 0.002, `outlier_weight` 0.001 and the gap model. **Refusals:** too few pairs, a multimodal, contaminated or discordant insert, or a non-converged marker fit give `LIBRARY_REFUSED` (exit 1, no prediction), never a fallback. `auto` also runs a 0.9×/1.1× intensity sensitivity audit (`LIBRARY_SENSITIVE`). **`FILE`** is the exact reproducibility route: exactly `insert_mean_bp`, `insert_sd_bp`, `support_sigmas`, `substitution_error`, `outlier_weight`, `fragment_start_intensity` (and optionally `marker_background`), each once; unknown, repeated or missing keys are refused. **The evidence above used `FILE`** |
| bubbles | `--bubbles FILE` (optional) | precomputed bubbles CSV; default: detected from the graph |

## 2. Preset, profile and precedence

**Precedence:** explicit flag > profile > preset > built-in default, whatever the argument order. Every effective setting and
its source is written to `<prefix>.settings.tsv`.

**`--preset short-read`** is always applied, whether or not it is given. It names no locus and reads no repository file.
It sets:
- the production algorithm and the short-read constants (marker k-mer and syncmer sizes, site span, copying prior,
  proposal width, placement band, refinement rounds, tie branching, representation guard);
- the exact computational backends and routes, their cache budgets and consumption capacities, and the state-class work cap;
- candidate generation: junction proposals (`repair_junction_proposals experimental_v1`, JP1.1), which can change
  predictions;
- budget-stop behaviour: `early_incumbent experimental_v1`;
- reporting: `output_schema sequence_first_v2`, `compound_annotation contract_v2`, `timing_detail`;
- resources: the 6 GiB footprint guard and the model-derived neighbourhood bound (`repair_max_neighbours model`); **no
  time limit**.

**Every value, and its built-in counterpart, is in the option inventory (§3).** That table is generated from the parser
and the preset, and the test suite checks it.

**Not set by the preset:** the dosage stage (`--dosage-check`), the reference, the panel curation and diagnostics.

**When a flag removes a preset feature's requirement,** the preset-sourced feature is switched off instead of refusing the
run. This happens with `--adaptive-repair off`, `--repair-rounds 0` or `--output-schema legacy`.
- `.settings.tsv` records it as `preset … (inactive: <reason>)`.
- Under the unrepaired view, `output_schema` becomes legacy.
- The same features set explicitly by a flag or profile keep their refusals.

**`--profile FILE`:** `key<TAB>value` lines.
- The file must name itself (`profile <name>`).
- A key set twice in one file is refused.
- **Repeated `--profile` files apply in order, a later file silently overriding an earlier one.**
  - **Limitation:** the settings source then names the profile that set each key, but there is no merged-profile digest.
  - Prefer one complete profile.

**A resource-adjusted run:** a one-key expert profile (for example `profile my-limits` + `time_limit_seconds 900`) on top of
the preset. Alternatively, copy a complete profile and change only resource keys (§4).

## 3. Options by purpose: the inventory

`panvar genotype --help` lists the everyday options. `--help-advanced` lists every option, grouped by purpose:
1. input and panel overrides;
2. the copying prior and the marker / placement model;
3. candidate generation and bounded refinement;
4. resource limits, work caps, cache budgets and exact backends;
5. reporting and diagnostics;
6. legacy and development compatibility.

**Columns** (each entry gives its accepted values, effective short-read default, units, what it can change, and its
status):
- **changes:**
  - **inference:** the prediction may change;
  - **computation:** time or memory only; results gated identical on completed runs;
  - **reporting:** outputs or annotations only;
  - **completion:** whether or how a run completes (caps, guards, budget stops).
- **status:**
  - **everyday:** in `--help`;
  - **supported:** an advanced flag;
  - **expert profile only:** no flag; `key<TAB>value` in a `--profile` file;
  - **legacy** or **development:** kept working, not for ordinary use.
- **Aliases** are listed once, beside the preferred spelling.

**The inventory** (generated by `genotype_cli_contract_test --print-inventory`; the test suite fails if this table differs
from the parser and the preset):

<!-- option-inventory:begin -->
| group | option | aliases | profile key | values | effective default (short-read) | units | changes | status |
|---|---|---|---|---|---|---|---|---|
| everyday | `--graph FILE` | --gfa, -i | — | a GFA path | required | - | input | everyday |
| everyday | `--reads1 FILE / --reads2 FILE` | — | — | FASTQ paths (.gz accepted) | required | - | input | everyday |
| everyday | `--out PREFIX` | --out-prefix, -o | — | a path prefix | genotype | - | reporting | everyday |
| everyday | `--reference-path NAME` | -r | — | a path name | the single path of the graph's RS:Z sample (refused, candidates listed, when absent or ambiguous) | - | inference | everyday |
| everyday | `--threads N` | — | threads | N >= 1 | 1 | worker threads | computation | everyday |
| everyday | `--dosage-check` | — | — | flag | off (dosage_check NOT_REQUESTED) | - | reporting | everyday |
| everyday | `--build-panel-index FILE` | — | — | a new file path | none | - | computation | everyday |
| everyday | `--panel-index FILE` | — | — | an index path | none (prepared in memory) | - | computation | everyday |
| everyday | `-q, --quiet` | — | — | flag | off | - | reporting | everyday |
| everyday | `-h, --help / --help-advanced` | — | — | flag | - | - | reporting | everyday |
| input | `--library auto\|FILE` | — | — | auto, or a library file | auto | - | inference | supported |
| input | `--panel-curation FILE` | — | — | a file of path names | none | - | inference | supported |
| input | `--bubbles FILE` | — | — | a bubbles CSV | computed from the graph | - | inference | supported |
| input | `--min-variant-bp N` | — | min_variant_bp | N >= 1 | 50 | bp | inference | supported |
| input | `--preset short-read` | — | — | short-read | short-read | - | inference+computation | supported |
| input | `--profile FILE` | — | — | key<TAB>value file naming itself | none | - | any | supported |
| model | `--switch-probability P` | — | switch_probability | 0 < P <= 0.5 | 0.01 | per block boundary | inference | supported |
| model | `--kmer-size K` | — | kmer_size | 1 <= K <= 31 | 31 | bp | inference | supported |
| model | `--syncmer-s S` | — | syncmer_s | 1 <= S <= K | 15 | bp | inference | supported |
| model | `--mismatch-band N` | — | mismatch_band | N >= 0 | 6 | edits per mate | inference | supported |
| model | `copying_transition` | — | copying_transition | constant_v1 \| distance_ref_core_v1_experimental | constant_v1 | - | inference | expert profile only |
| model | `copying_length_bp` | — | copying_length_bp | > 0 (hex double) | 0x0p+0 | bp | inference | expert profile only |
| search | `--mosaic-candidates K` | --max-mosaic-candidates | max_mosaic_candidates | K >= 0 | 128 | mosaics | inference+computation | supported |
| search | `--proposal-state-limit N` | — | proposal_state_limit | N >= 1 | 256 | beam states | inference | supported |
| search | `--adaptive-repair on\|off` | — | — | on \| off | on | - | inference | supported |
| search | `--repair-rounds N` | — | repair_rounds | N >= 0 | 32 (built-in 0) | rounds | inference+completion | supported |
| search | `--repair-local-pairs N` | — | repair_local_pairs | N >= 1 | 32 | allele pairs per block | inference | supported |
| search | `--repair-max-neighbours model\|N` | — | repair_max_neighbours | model, or N >= 0 (0: none) | model (the model-derived bound) (built-in 670) | candidates per round | completion | supported |
| search | `--repair-tie-branching on\|off` | — | repair_tie_branching | on \| off | on (built-in off) | - | inference | supported |
| search | `--repair-tie-max-branches N` | — | repair_tie_max_branches | N >= 1 | 16 | branches | completion | supported |
| search | `repair_junction_proposals` | — | repair_junction_proposals | off \| experimental_v1 | experimental_v1 (built-in off) | - | inference | expert profile only |
| search | `early_incumbent` | — | early_incumbent | off \| experimental_v1 | experimental_v1 (built-in off) | - | completion | expert profile only |
| search | `repair_joint_search` | — | repair_joint_search | off \| experimental_v1 | off | - | inference | expert profile only |
| resources | `footprint_limit_bytes` | — | footprint_limit_bytes | N >= 0 (0: none) | 6 GiB (built-in 0) | bytes | completion | expert profile only |
| resources | `time_limit_seconds` | — | time_limit_seconds | S >= 0 (0: none) | 0 | seconds | completion | expert profile only |
| resources | `max_physical_origins` | — | max_physical_origins | N >= 0 | 500000000 | origins | completion | expert profile only |
| resources | `max_state_classes` | — | max_state_classes | N >= 0 | 600000000 (built-in 200000000) | state classes | completion | expert profile only |
| resources | `max_piece_hits` | — | max_piece_hits | N >= 0 | 2000000000 | piece hits | completion | expert profile only |
| resources | `max_mass_cells` | — | max_mass_cells | N >= 0 | 200000000 | cells | completion | expert profile only |
| resources | `repair_max_physical_origins` | — | repair_max_physical_origins | N >= 0 | 4900000000 | origins | completion | expert profile only |
| resources | `repair_max_state_classes` | — | repair_max_state_classes | N >= 0 | 4200000000 | state classes | completion | expert profile only |
| resources | `repair_max_piece_hits` | — | repair_max_piece_hits | N >= 0 | 14000000000 | piece hits | completion | expert profile only |
| resources | `repair_max_latent_pairs` | — | repair_max_latent_pairs | N >= 0 | 14000000000 | pairs | completion | expert profile only |
| resources | `max_union_cached_placements_consumed` | — | max_union_cached_placements_consumed | N >= 0 | 2000000000 (built-in 0) | placements | completion | expert profile only |
| resources | `max_union_cached_origins_consumed` | — | max_union_cached_origins_consumed | N >= 0 | 2000000000 (built-in 0) | origins | completion | expert profile only |
| resources | `max_union_cached_states_consumed` | — | max_union_cached_states_consumed | N >= 0 | 2000000000 (built-in 0) | states | completion | expert profile only |
| resources | `max_union_cached_prefilters_consumed` | — | max_union_cached_prefilters_consumed | N >= 0 | 2000000000 (built-in 0) | prefilters | completion | expert profile only |
| resources | `adaptive_fine_max_served_placements` | — | adaptive_fine_max_served_placements | N >= 0 | 2000000000 (built-in 0) | placements | completion | expert profile only |
| resources | `emission_cache_bytes` | — | emission_cache_bytes | N >= 0 | 512 MiB (built-in 0) | bytes | computation | expert profile only |
| resources | `repair_origin_cache_bytes` | --repair-origin-cache-bytes | repair_origin_cache_bytes | N >= 0 | 768 MiB (built-in 384 MiB) | bytes | computation | expert profile only |
| resources | `repair_tie_seed_cache_bytes` | --repair-tie-seed-cache-bytes | repair_tie_seed_cache_bytes | N >= 0 | 0 | bytes | computation | expert profile only |
| resources | `union_placement_cache_bytes` | — | union_placement_cache_bytes | N >= 0 | 512 MiB (built-in 0) | bytes | computation | expert profile only |
| resources | `union_origin_cache_bytes` | — | union_origin_cache_bytes | N >= 0 | 1 GiB (built-in 0) | bytes | computation | expert profile only |
| resources | `union_prefilter_cache_bytes` | — | union_prefilter_cache_bytes | N >= 0 | 128 MiB (built-in 0) | bytes | computation | expert profile only |
| resources | `start_sort_scratch_bytes` | — | start_sort_scratch_bytes | N >= 0 | 320 MiB (built-in 0) | bytes | computation | expert profile only |
| resources | `indexed_union_store_bytes` | — | indexed_union_store_bytes | N >= 0 | 2 GiB (built-in 0) | bytes | computation+completion | expert profile only |
| resources | `adaptive_fine_store_bytes` | — | adaptive_fine_store_bytes | N >= 0 | 512 MiB (built-in 0) | bytes | computation+completion | expert profile only |
| resources | `marker_proposal_reuse_table_bytes` | — | marker_proposal_reuse_table_bytes | N >= 0 | 4 MiB (built-in 0) | bytes | computation | expert profile only |
| resources | `union_scoring_backend` | — | union_scoring_backend | reference \| indexed_fine_experimental | indexed_fine_experimental (built-in reference) | - | computation | expert profile only |
| resources | `adaptive_window_discovery` | — | adaptive_window_discovery | reference \| fine_experimental | fine_experimental (built-in reference) | - | computation | expert profile only |
| resources | `indexed_union_store_layout` | — | indexed_union_store_layout | exact \| legacy_i2a | exact | - | computation | expert profile only |
| resources | `indexed_union_evidence_layout` | — | indexed_union_evidence_layout | reconstructed \| fine_views_experimental | fine_views_experimental (built-in reconstructed) | - | computation | expert profile only |
| resources | `adaptive_fused_seed` | — | adaptive_fused_seed | 0 \| 1 | 1 (built-in 0) | - | computation | expert profile only |
| resources | `combination_prepared_library_terms` | — | combination_prepared_library_terms | 0 \| 1 | 1 (built-in 0) | - | computation | expert profile only |
| resources | `gap_recursion_backend` | — | gap_recursion_backend | scaled_reference \| fixed_geometry_experimental | fixed_geometry_experimental (built-in scaled_reference) | - | computation | expert profile only |
| resources | `dosage_census_backend` | — | dosage_census_backend | reference \| interned_experimental | interned_experimental (built-in reference) | - | computation | expert profile only |
| resources | `union_fragment_blocks` | — | union_fragment_blocks | off \| fine_core_4096_experimental | fine_core_4096_experimental (built-in off) | - | computation | expert profile only |
| resources | `repair_round1_origin_admission` | — | repair_round1_origin_admission | existing \| reuse_only_experimental | reuse_only_experimental (built-in existing) | - | computation | expert profile only |
| resources | `union_progressive_schedule` | — | union_progressive_schedule | off \| experimental_v1 | off | - | computation | expert profile only |
| resources | `adaptive_emission_handoff` | — | adaptive_emission_handoff | 0 \| 1 | 0 | - | computation | expert profile only |
| reporting | `--dosage-mode off\|audit\|repair` | — | — | off \| audit \| repair | off | - | reporting | supported |
| reporting | `--output-schema sequence_first_v2\|legacy\|sequence_first_v1` | — | output_schema | sequence_first_v2 \| legacy \| sequence_first_v1 | sequence_first_v2 (built-in legacy) | - | reporting+completion | supported |
| reporting | `--representation-guard on\|off` | — | representation_guard | on \| off | on (built-in off) | - | reporting | supported |
| reporting | `compound_annotation` | — | compound_annotation | off \| contract_v2 | contract_v2 (built-in off) | - | reporting | expert profile only |
| reporting | `compound_interval_declarations` | — | compound_interval_declarations | a file | — | - | reporting | expert profile only |
| reporting | `block_support_rule` | — | block_support_rule | legacy \| full_data_likelihood_v1_experimental \| full_data_likelihood_v2_experimental | legacy | - | reporting | expert profile only |
| reporting | `library_sensitivity_audit` | — | library_sensitivity_audit | 0 \| 1 | 1 | - | reporting+computation | expert profile only |
| reporting | `timing_detail` | — | timing_detail | 0 \| 1 | 1 (built-in 0) | - | reporting | expert profile only |
| reporting | `block_diagnostics` | — | block_diagnostics | off \| diagnostic_v1 | off | - | reporting | expert profile only |
| reporting | `emission_lookup_timing` | — | emission_lookup_timing | 0 \| 1 | 0 | - | reporting | expert profile only |
| legacy | `--max-candidates` | — | — | - | refused | - | - | legacy |
| legacy | `--dosage-audit off\|shadow` | — | — | off \| shadow | off | - | reporting | legacy |
| legacy | `--phase-audit-candidates N` | — | phase_audit_candidates | N >= 0 | 0 | candidates | reporting | development |
| legacy | `--dosage-allocation-test` | — | — | flag | off | - | reporting | development |
| legacy | `--library-diagnostic` | — | — | flag | off | - | reporting | development |
| legacy | `diagnostic_* and test_* profile keys` | — | — | - | off | - | - | development |
<!-- option-inventory:end -->

**Terms.** Read these before interpreting results:
- **Donor redraw:** `--switch-probability` is the copying prior's per-boundary donor-redraw probability, a regularisation of
  the panel model. It is **not a biological recombination probability**.
- **Mosaic candidates:** `--mosaic-candidates` sets the **initial** off-panel mosaics. The refinement may generate and score
  further mosaics beyond that number.
- **Bounded refinement:** local moves scored exactly within declared caps. It is **not a globally exhaustive search**, and a
  completed ranking is conditional on what was explored.
- **`min_variant_bp`** constructs the variant sites (blocks). It is **not a filter on final calls**.
- **Library estimation** learns the insert distribution and the fragment-start intensity only. The substitution, outlier
  and gap parameters stay fixed.
- **Support is conditional on the explored candidates:** READS means no explored rival projecting to a different genotype
  was within bounds of the winner's fragment term.
- **Prior-dependent** (`read_support PRIOR_DEPENDENT`): the reads alone did not separate the call from an explored rival.
  This **does not establish that the prior alone caused the selection**.
- **Phase** is assessed once per locus (the decision's phase support), **not independently at each boundary**.
- **Locus summaries** must be read with their coverage counts and `not_assessed` list. A positive value covers only the
  assessed blocks, and `[assessed/applicable]` marks partial coverage.
- **Nothing here is calibrated reliability,** and there is no quality threshold.

## 4. Resource limits, caches and work caps

These are four different things. None of them changes a completed prediction while it does not bind. The values are in
the inventory (§3, group 4).

| kind | keys | behaviour when it binds |
|---|---|---|
| **total-footprint guard** | `footprint_limit_bytes` | **sampled** every 100 ms, so it is not an allocation ceiling and a short spike can exceed it. Exit 125, or 3 PROVISIONAL with an early-incumbent checkpoint |
| **time guard** | `time_limit_seconds` (none under the preset) | exit 124 before any checkpoint (NO_PREDICTION), or 3 PROVISIONAL after one |
| **cache and store budgets** | `*_cache_bytes`, `*_store_bytes`, `start_sort_scratch_bytes`, `marker_proposal_reuse_table_bytes` | they bound **retained** memory for reuse. Work beyond a budget is recomputed. A store that cannot hold the required data (`indexed_union_store_bytes`, `adaptive_fine_store_bytes`) ends the run with a recorded resource refusal. They do not cap work |
| **executed-work caps** | `max_physical_origins`, `max_state_classes`, `max_piece_hits`, `max_mass_cells`, `repair_max_*`, `max_union_cached_*_consumed`, `adaptive_fine_max_served_placements` | a typed refusal, no-call or PROVISIONAL with the cap named. They do not degrade silently |

**More RAM alone does not raise any work cap,** and it does not guarantee completion. A one-key expert profile sets any of
these, for example `profile my-limits` + `time_limit_seconds 900`.

**Not implemented:** memory-per-thread budgeting and disk spilling.

## 5. Reporting and diagnostics

- **The report contract:** `--output-schema` (default `sequence_first_v2`). The primary view is `assessment_v1` (§7.0). The
  unrepaired view reports with `legacy`.
- **The dosage check:** `--dosage-check` (an optional audit of the reported alleles; it never changes the prediction).
  - The expert spelling is `--dosage-mode off|audit|repair`.
  - `repair` (may replace a block) needs `--output-schema legacy` and is refused under v2.
  - A contradiction with `--dosage-check` is refused.
- **Profile-only report keys** (`compound_annotation`, `compound_interval_declarations`, `block_support_rule`,
  `library_sensitivity_audit`, `timing_detail`, `block_diagnostics`, `emission_lookup_timing`) are in the inventory.
- **The older headline fields** (`genotype_filter`, `dosage_filter`, `phase_filter`, `call_summary`) are **superseded
  compatibility diagnostics** for this experimental release (§7).

## 6. Outputs: what you get

| outcome | exit (v2) | what exists |
|---|---|---|
| **completed prediction** | 0 | `.sequence_locus.tsv` `prediction PREDICTION`, `completion completed_ranking` (even when uncertain); full block tables |
| **PROVISIONAL** | 3 | the best **fully evaluated** candidate after a budget stop or declared cap (early-incumbent checkpoint). It is **not a completed search** and carries no decision set and no confidence |
| **no prediction** | 1 | `prediction NO_PREDICTION` with `no_call_reason` (for example `neighbourhood_bound` with no checkpoint to keep). A legacy `distinct_tie` is **not** a no-prediction when the tie is only over labellings of one sequence pair: that is a PREDICTION (exit 0) with `legacy_no_call_reason distinct_tie` |
| time / footprint stop before any checkpoint | 124 / 125 | minimal status (`TIME_LIMIT` / `RESOURCE_REFUSAL`) and a `NO_PREDICTION` `.sequence_locus.tsv` (`completion budget_stop_before_first_evaluated_candidate`) |
| argument or declared-library error | 2 | **refused before the prefix is touched**: every file already there is left as it was |
| input, model or output error found during the run | 2 | the run already owns the prefix: earlier run-owned outputs are gone and a `FAILED` status (or, after the sequence report was published, `locus_output FAILED: …`) says what happened |

**Exit agreement:** the process exit, `.sequence_locus.tsv` `exit_code` and the status (`sequence_exit_code`, and
`exit_code` for stops and failures) state the same value in every outcome above.

**What is refused before the prefix is touched:**
- contradictory or invalid arguments and profiles;
- a declared `--library FILE` that is missing, unreadable or invalid (bad keys, values or format).

**What is not:** errors that can only be found while running, such as unreadable reads, a graph or model refusal, a
resource stop or a failed write. Those happen after the run has cleared its prefix. The guarantee is not "every error is
transactional".

The legacy schema uses 0 call / 1 no-call / 2 / 124 / 125.

**Sequences under v2:**
- **Both haplotypes on the panel:** **no FASTA** (`fasta_reason PANEL_REPRESENTED_PAIR`). `.sequence_haplotypes.tsv` gives
  each haplotype's SHA-256, length, `origin PANEL`, every sequence-equivalent panel path and an extraction reference
  (`gfa_path=…;interval=0-<len>;orientation=…`), with the GFA SHA-256 and the model fingerprint.
  - **To recover the bytes:** `python3 tests/analysis/s1d/extract_panel_haplotype.py --check <graph.gfa> <prefix>.sequence_haplotypes.tsv`.
  - This is a repository script, not a `panvar` subcommand.
- **Either haplotype off-panel** (a mosaic): `.sequence_haplotypes.fa` holds **both** records, each annotated with its origin
  (`OFF_PANEL_PAIR`). The legacy `.mosaic_haplotypes.fa` is suppressed under v2.

## 7. Decision guide: filters, ambiguity, audit and support

### 7.0 The primary view: `assessment_v1`

Read this first. `.sequence_locus.tsv` (mirrored in `.locus_genotype.tsv`) and `.sequence_blocks.tsv` (mirrored in the
block tables) carry it. Every assessment concerns the **published prediction**, and a legacy block fact is used only where
it verifiably describes the same alleles. The full contract is `tests/analysis/m15g/CONTRACT.md`.

| field | values |
|---|---|
| `prediction` | PREDICTION · PROVISIONAL · NO_PREDICTION (never overridden by the legacy call class) |
| `sequence_resolution` | RESOLVED · GENOTYPE_RESOLVED_SEQUENCE_SET · SEQUENCE_TIED · NOT_DETERMINED · NA |
| `genotype_check` (+ coverage) | NONE_DETECTED · CONFLICT (representation conflict, library sensitivity) · NOT_ASSESSED |
| `dosage_check` (+ coverage) | AGREES · SUSPECT · DISAGREES · UNAVAILABLE · NOT_REQUESTED · NOT_APPLICABLE: always about the reported alleles. Corrections are listed separately in `dosage_alternatives` (diagnostic) |
| `read_support` (+ coverage) | READS · PRIOR_DEPENDENT · NOT_ASSESSED (unknown support is never READS) |
| `phase` | READ_SUPPORTED · PRIOR_ONLY · UNRESOLVED · NOT_ASSESSED · NOT_APPLICABLE. **Locus-level** (phase support is decided once per locus) and applicable when two or more blocks are heterozygous, adjacent or not |
| `findings` / `not_assessed` | every concrete finding and every check that could not run, each with its block ids; both always written |
| `assessment_summary` | the scan line, with `[assessed/applicable]` wherever coverage is partial |

**How to read it:**
- **Missing evidence is not a failure,** and it never hides a finding.
- **There is no overall PASS** and nothing calibrated.

**Example, C4 with `--dosage-check`:** `PREDICTION;SEQUENCE_RESOLVED;GENOTYPE_NONE_DETECTED;DOSAGE_AGREES[1/2];SUPPORT_PRIOR_DEPENDENT;PHASE_PRIOR_ONLY`.
- `findings`: `PRIOR_DEPENDENT@sv:1,ctx:2,ctx:4,sv:5,sv:7;PHASE_PRIOR_ONLY@locus`.
- `not_assessed`: `DOSAGE:UNDERPOWERED@sv:1`.

**The headline filters below** (GENOTYPE_FILTER, DOSAGE_FILTER, PHASE_FILTER, `call_summary`) are kept unchanged for this
experimental release and labelled `headline_filters_status superseded_by_assessment_v1`. Known misreadings:
- they mark a v2 prediction NO_CALL when the legacy decision tied;
- NOT_EVALUATED outranks LOW_QUALITY;
- DOSAGE_FILTER PASS can describe an unapplied correction.

**First:** is it a completed prediction (exit 0, `PREDICTION`, `completed_ranking`)? If not, report its status. Never treat a
PROVISIONAL or absent result as a completed one.

**Then read these six independently:**

| question | fields | values and meaning |
|---|---|---|
| **sequence ambiguity** | `sequence_status` | `RESOLVED`: one sequence pair in the decision set. `GENOTYPE_RESOLVED_SEQUENCE_SET`: the genotype agrees, the sequences do not. `SEQUENCE_TIED`. **Not a confidence claim** |
| **representation (labelling) ambiguity** | per block `block_status`; `.representation_intervals.tsv` | `REPRESENTATION_AMBIGUOUS`: explored labellings **of the same sequence pair** differ at the block, so the **sequence may still be unique**. `GENOTYPE_UNRESOLVED`: the decision set disagrees on the block genotype. Computed over **explored** labellings only |
| **phase** | `phase_support`, PHASE_FILTER | read-supported → PASS; prior-only → LOW_QUALITY `MODEL_ONLY_PHASE`; unresolved → NO_CALL; next to a homozygous block → NOT_APPLICABLE. A genotype-equivalent call has no phase and no FASTA |
| **headline filters** | GENOTYPE_FILTER, DOSAGE_FILTER, PHASE_FILTER (per block and per locus, with `REASON@ids`) | **PASS means the computational and representation checks passed.** It is **not** calibrated reliability, and **no combination of PASS values establishes an exact or biologically correct genotype.** LOW_QUALITY reasons include a representation conflict, `LIBRARY_SENSITIVE` (`auto` only), unresolved dosage content and `MODEL_ONLY_PHASE`. NOT_EVALUATED means a check could not complete (dosage off, audit failed or underpowered, library audit incomplete). NO_CALL means a no-call |
| **dosage audit** (only with `--dosage-check`) | `.dosage_status.tsv` (COMPLETE / FAILED / NOT_RUN, with note); block `dosage_status`. A requested audit that does not complete is also reported on stderr, even with `--quiet`, and DOSAGE_FILTER is NOT_EVALUATED (`audit_FAILED`), never consistent. Known case: the CYP2D6 development cell, where the audit's unit library refuses sv:9 (`alleles disagree on their flanks or anchor trims`) | **disagreement:** CONTENT_REPAIRED, TRIGGERED_UNREPAIRED, SUSPECT. **Unavailable:** audit FAILED or NOT_RUN, or block UNDERPOWERED / NOT_EVALUATED. The audit **never** replaces haplotypes under v2. Its representatives appear only in dosage diagnostics (`REPRESENTATIVE_ONLY_NOT_RESOLVED`) |
| **conditional support** | `genotype_support` (READS / PRIOR_ONLY), `error_mass_lower/upper`, `legacy_representative_GQ`, `candidate_scope CONDITIONAL`, `structure_confidence NA` | conditional on the **explored** candidate pool. PRIOR_ONLY: the reads alone do not separate the call from an explored rival. **Uncalibrated:** in development, error mass did not predict block accuracy better than a constant baseline, and catastrophic blocks carried GQ up to 443 |

**Descriptive warnings** (not rejection rules; see M15-B `CONTRACT.md`):

| warning | trigger |
|---|---|
| `REPEATED_UNIT_STRUCTURE` | YES / NO / **UNKNOWN**, from the census route; unavailable means UNKNOWN |
| `DOSAGE_AUDIT_DISAGREES` / `DOSAGE_AUDIT_UNAVAILABLE` | the dosage audit fields above |
| `PRIOR_ONLY_SUPPORT` | `genotype_support = PRIOR_ONLY` |

**These warnings do not identify every error.**

**Why PASS is not calibrated reliability — a worked example** (M15-B, reserved ACOT families GB21/HG02572 and NA19468):
- **The call:** block sv:9 called as 902 + 424 bp against a truth of 424 + 424. That is one extra ~478 bp segment: excess 478
  and 475 bp, normalized excess 0.36.
- **How both blocks were reported:**
  - GENOTYPE_FILTER **PASS**;
  - `REPEATED_UNIT_STRUCTURE` NO (the census route is not repeated-unit);
  - no dosage disagreement;
  - support READS for one, PRIOR_ONLY for the other.
- **What it shows:** no filter or warning flagged them. These were the two largest errors among the 18 reserved families.

## 8. Development-only interfaces (not part of the public CLI)

**Diagnostic and test settings** (refused or unsupported for production use):
- `diagnostic_*` keys and `test_*` hooks;
- `--dosage-allocation-test`, `--library-diagnostic`, `--dosage-audit off|shadow` (an older spelling of `--dosage-mode
  off|audit`);
- the cache flags `--repair-origin-cache-bytes` and `--repair-tie-seed-cache-bytes`;
- `phase_audit_candidates`.

**No CLI exists for these:**
- **In-process multi-sample reuse** (one process analysing several samples against one prepared locus): library/API-level,
  exercised only by development drivers (for example `genotype_prepared_locus_driver`). The persisted index has a public
  CLI (§8a).
- **Batch / multi-sample runs:** one locus and one sample per invocation.

**Long reads:** not implemented (`--preset long-read` is refused).

## 8a. Build once, genotype many: the panel index

The panel index is the persisted, **sample-independent** preparation of one locus: the structural model and the marker
panel (format: the I3-C/I3-D prepared-locus file). Building it once saves that preparation in every later run.
- **What it does not save:** reading and scoring the sample. Library estimation, markers, depth, candidates, placements,
  refinement and the dosage audit all still run per sample.
- **How much time it saves depends on the locus.** Measured times are in `tests/analysis/m15e/RESULT.md`. There is no
  universal speedup.

**Build** (no reads and no library):

```
panvar genotype --graph locus.gfa --reference-path NAME [--panel-curation FILE] [--bubbles FILE] \
  --build-panel-index locus.pvpl
```

- **Outputs:** the index, published atomically, and a build record on stdout. The record gives the schema, preparation,
  panel-extension and marker-panel versions; the identity (graph SHA-256, reference, bubbles, `min_variant_bp`, curation,
  k, s); the model fingerprint; the file size and SHA-256; and the preparation and write times.
- **An existing file is never overwritten** (refused, exit 2).
- **Refused in build mode:** sample arguments (`--reads1/2`, `--library`, `--out`, `--dosage-check`, `--dosage-mode`) and
  `--panel-index`.
- **Nothing else is written:** no output prefix is cleared, no sample report is created and no cache is left. The one
  temporary, a sorted graph next to the index, is removed.

**Genotype each sample:**

```
panvar genotype --graph locus.gfa --reference-path NAME [--panel-curation FILE] [--bubbles FILE] \
  --panel-index locus.pvpl --reads1 S1.R1.fq.gz --reads2 S1.R2.fq.gz --threads 4 --out S1
panvar genotype --graph locus.gfa --reference-path NAME [--panel-curation FILE] [--bubbles FILE] \
  --panel-index locus.pvpl --reads1 S2.R1.fq.gz --reads2 S2.R2.fq.gz --threads 4 --out S2
```

**What is verified before the output prefix is touched:**
- **The file:** integrity (header SHA-256), bounded decoding, schema and preparation versions, and reconstruction digests.
- **This sample's inputs against the index's identity:** graph content, reference, bubbles, curation, `--min-variant-bp`,
  `--kmer-size` and `--syncmer-s`.
- **Inputs are still read:** `--graph` and the other inputs are read and hashed for this check, and the extraction
  references point into the graph.

**Refusals:** a missing, truncated, corrupt or incompatible index is refused with a named reason (exit 2, for example
`io`, `truncated`, `checksum_mismatch`, `incompatible_preparation`, or the differing identity field). It is never replaced
by an in-memory preparation.

**What stays the same:**
- **Rebuilding is never automatic,** and a different git revision does not by itself make an index incompatible.
- **Output:** the status records `prepared_locus_origin loaded`, `prepared_locus_panel_source` and `panel_index*`. A panel
  that had to be rebuilt from the stored model (an index without a panel section, or with one this binary does not
  support) is reported on stderr and in the status, never as a loaded panel.

**Rebuild the index when** any of these change: the graph, the reference path, the panel curation, the bubbles (or
`--min-variant-bp` when bubbles are computed), `--kmer-size`, `--syncmer-s`, or a new panvar version that raises the
schema or preparation version (the load says so).

**These do not invalidate an index:** reads, library (`auto` or `FILE`), threads, `--dosage-check`, output schema,
resource limits and every other sample or search option.

**Never stored in an index:** library estimates, read counts or depth, candidate pools, placements, refinement state and
dosage results.

## 9. Tested commands (M15-D `runs/ex.*`; development controls only)

**The ordinary command** (C4; CYP2D6 and ACOT add their `--panel-curation`; these development graphs carry no `RS:Z`, so
they name the reference):

```
panvar genotype --graph C4.loo.gfa --reference-path 'HG00096#1#haplotype1-0000024:31848049-32074186' \
  --reads1 r1.fq --reads2 r2.fq --threads 4 --out sample.c4
```

- **Exit status:** it matches the reported outcome.
- **Settings:** every accelerated key is recorded from the preset, the library is estimated from the reads, and dosage is
  NOT_EVALUATED.
- **Sequences:** the v2 FASTA rule holds.
- The estimated library and the runtime are listed in `tests/analysis/m15d/RESULT.md`.

**Reproducing the evidence configuration** without a profile:

```
panvar genotype ... --library library.txt --dosage-check --out sample.c4
```

The prediction is identical to the archived M11-D accelerated run at all three cells.

**The other tested commands:**

| test | command | result |
|---|---|---|
| omitted vs explicit auto | `--library auto` | identical prediction, library estimate and settings |
| dosage absent vs `--dosage-check` | the declared-library command without `--dosage-check` | identical prediction; DOSAGE_FILTER NOT_EVALUATED (`dosage_mode_off`) without it; the audit adds only `.dosage*` diagnostics |
| unrepaired view | `--adaptive-repair off` | runs under the legacy schema; the preset's repair-only features are recorded inactive |
| automatic library refusal | 300 read pairs | `LIBRARY_REFUSED`, exit 1, no prediction, no FASTA |
| binding time limit | a one-key expert profile with `time_limit_seconds 8` (ACOT) | exit 3 PROVISIONAL or 124; never a completed prediction |
| `--dosage-check` with `--dosage-mode off` / `--dosage-audit off` / the legacy repair | | refused, exit 2, the prefix untouched |
| `--dosage-mode repair` under the default v2 | | refused, exit 2 (use `--output-schema legacy --dosage-mode repair`) |
| explicit `--output-schema sequence_first_v2` with `--adaptive-repair off` | | refused, exit 2 |
| no `--reference-path` on a graph without `RS:Z` | | refused, exit 2, the candidate paths listed |
| `--preset long-read` | | refused, exit 2 |

The M15-C profile-based commands (§0 profiles, resource-adjusted copies, refusals for duplicate, unknown or unnamed profile
keys and the removed `--max-candidates`) are unchanged.

## 10. Recommended future removals and aliases (not implemented; each needs separate approval)

**Remove** (after the experimental release):
- the superseded headline columns, or move them behind a flag;
- `--dosage-audit` (an older spelling of `--dosage-mode`);
- the refused `--max-candidates` message;
- the `legacy_i2a` store layout;
- the `--repair-origin-cache-bytes` / `--repair-tie-seed-cache-bytes` diagnostic flag spellings (the profile keys remain).

**Consider aliases or flags:**
- `--time-limit SECONDS` and `--memory-limit BYTES` for the two run guards, which are now profile-only;
- `--curation` as a short alias of `--panel-curation`.

**Keep as is:** `--mosaic-candidates` (preferred) with `--max-mosaic-candidates` as its alias; the `-i` / `-o` / `-r` short
aliases.

