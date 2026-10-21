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
- **the algorithm and constants:**
  - `kmer_size` 31, `syncmer_s` 15, `min_variant_bp` 50, `switch_probability` 0.01;
  - `proposal_state_limit` 256, `mismatch_band` 6, `phase_audit_candidates` 0;
  - `repair_rounds` 32, `repair_tie_branching` on, `representation_guard` on.
- **exact computational backends and routes** (§3, the performance table), plus `union_fragment_blocks` and
  `repair_round1_origin_admission`.
- **cache budgets and consumption capacities** (§4), and `max_state_classes` 6e8.
- **candidate generation:** `repair_junction_proposals experimental_v1` (JP1.1). It can change predictions (§3).
- **budget-stop behaviour:** `early_incumbent experimental_v1`.
- **reporting:** `output_schema sequence_first_v2`, `compound_annotation contract_v2`, `timing_detail` 1.
- **resources:**
  - `footprint_limit_bytes` 6 GiB;
  - `repair_max_neighbours model`;
  - **no time limit.**

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

## 3. Scientific and search controls

The **Affects predictions?** column means:
- **yes:** changing the value can change the reported haplotypes;
- **route:** a computation path gated identical on completed runs;
- **report:** outputs or annotations only.

**Core settings** (flag = profile key; same parser):

| flag / key | values | default (built-in → preset → tested profile) | affects predictions? | notes and refusals |
|---|---|---|---|---|
| `--kmer-size` / `kmer_size` | 1..31 | 31 → 31 | yes (markers) | |
| `--syncmer-s` / `syncmer_s` | 1..k | 15 → 15 | yes | |
| `--min-variant-bp` / `min_variant_bp` | ≥ 1 | 50 → 50 | yes (block decomposition) | |
| `--switch-probability` / `switch_probability` | (0, 0.5] | 0.01 → 0.01 | yes (copying prior) | |
| `--proposal-state-limit` / `proposal_state_limit` | ≥ 1 | 256 → 256 | yes (proposal beam) | |
| `--mismatch-band` / `mismatch_band` | ≥ 0 | 6 → 6 | yes (placement model) | |
| `--mosaic-candidates` / `max_mosaic_candidates` | ≥ 0 | 128 | yes | `--max-candidates` is removed and refused |
| `--adaptive-repair` | on, off | on | yes | **bounded coordinate refinement:** local moves scored exactly within declared caps, **not a global optimum**. `off` is the unrepaired development view |
| `--repair-rounds` / `repair_rounds` | ≥ 0 | 0 → 32 | yes | resolved with `--adaptive-repair` / tie branching. A contradiction (off with N > 0, on with 0, tie branching with 0 rounds) is refused |
| `--repair-local-pairs` / `repair_local_pairs` | ≥ 1 | 32 | yes | |
| `--repair-max-neighbours` / `repair_max_neighbours` | `model`, or ≥ 0 (0: none) | 670 → **`model`** (tested profiles: C4 670, CYP2D6 958, ACOT 1004) | yes, when numeric (a neighbourhood above it no-calls `neighbourhood_bound`) | **`model`:** the declared bound is the model's own structural maximum, `Σ_blocks 2·min(repair_local_pairs, a(a+1)/2)` (C4 670, CYP2D6 958, ACOT 998). It is recorded as `repair_neighbour_bound`. The locus-specific caps existed only to admit each locus, and each equalled or exceeded this maximum |
| `--repair-tie-branching` / `repair_tie_branching` | on, off | off → on | yes | |
| `--repair-tie-max-branches` / `repair_tie_max_branches` | ≥ 1 | 16 | yes (beyond it: no-call) | |
| `--representation-guard` / `representation_guard` | on, off | off → on | report (GENOTYPE_FILTER) | |
| `--phase-audit-candidates` | ≥ 0 | 0 | development only | refused under adaptive repair |

**Experimental profile keys** (profile only; "preset" = the short-read preset's value):

| key | values (tested value in **bold**) | default | affects predictions? |
|---|---|---|---|
| `early_incumbent` | off, **experimental_v1** | off → preset experimental_v1 | completed runs identical; on a budget stop it keeps a PROVISIONAL checkpoint. Needs `sequence_first_v2` and adaptive repair, otherwise refused |
| `union_fragment_blocks` | off, **fine_core_4096_experimental** | off → preset fine_core_4096 | route |
| `repair_round1_origin_admission` | existing, **reuse_only_experimental** | existing → preset reuse_only | route |
| `repair_junction_proposals` | off, **experimental_v1** | off → preset experimental_v1 | **yes:** adds junction proposals (it changed one candidate pool and one block's conditional GQ in M10-A) |
| `repair_joint_search` | **off**, experimental_v1 | off | yes (not in the tested configuration) |
| `copying_transition` / `copying_length_bp` | **constant_v1**, distance_ref_core_v1_experimental / > 0 | constant_v1 | yes (not in the tested configuration) |
| `block_support_rule` | **legacy**, full_data_likelihood_v1/v2_experimental | legacy | report (support semantics) |
| `union_progressive_schedule`, `adaptive_emission_handoff` | off / 0 (tested) | off / 0 | route (not in the tested configuration) |

**Performance backends** (route; identical results gated in the I/S batches and in M10-A; all set by the preset):

| key | tested and preset value |
|---|---|
| `union_scoring_backend` | indexed_fine_experimental |
| `adaptive_window_discovery` | fine_experimental |
| `indexed_union_store_layout` | exact |
| `indexed_union_evidence_layout` | fine_views_experimental |
| `adaptive_fused_seed` | 1 |
| `combination_prepared_library_terms` | 1 |
| `gap_recursion_backend` | fixed_geometry_experimental |
| `dosage_census_backend` | interned_experimental |
| `marker_proposal_reuse_table_bytes` | 4194304 |

## 4. Resource limits, caches and work caps

These are **four different things.** None of them changes a completed prediction while it does not bind.

| kind | keys | behaviour when it binds |
|---|---|---|
| **total footprint guard** | `footprint_limit_bytes` (built-in 0: none; preset and tested 6 GiB) | a **sampled** guard: the physical footprint is polled every 100 ms, so it is not an allocation ceiling, and a short spike can exceed it. Exit 125, or a PROVISIONAL checkpoint with exit 3 under `early_incumbent` |
| **time guard** | `time_limit_seconds` (0: none, also under the preset; tested profiles 900; set it with a one-key profile) | exit 124 before any checkpoint (NO_PREDICTION, as in the tested 8 s example), or 3 PROVISIONAL after one |
| **cache and store budgets** (the preset sets the tested values) | `emission_cache_bytes`, `repair_origin_cache_bytes`, `union_placement_cache_bytes`, `union_origin_cache_bytes`, `union_prefilter_cache_bytes`, `indexed_union_store_bytes`, `adaptive_fine_store_bytes`, `start_sort_scratch_bytes`, `repair_tie_seed_cache_bytes`, `marker_proposal_reuse_table_bytes` | they bound **retained** memory for reuse. Work beyond a budget is recomputed, or a store that cannot hold the required data ends the run with a recorded resource refusal (`.mosaic_status.tsv` names it). They do not cap work |
| **executed-work caps** | `max_physical_origins` (5e8), `max_state_classes` (2e8; preset and tested 6e8), `max_piece_hits`, `max_mass_cells`, `repair_max_physical_origins`, `repair_max_state_classes`, `repair_max_piece_hits`, `repair_max_latent_pairs`, `max_union_cached_*_consumed`, `adaptive_fine_max_served_placements` | a typed refusal, no-call or PROVISIONAL with the cap named. They do not degrade silently |

**More RAM alone does not raise any work cap,** and it does not guarantee completion. Large loci at high coverage reach work
caps such as the 500 M physical-origin cap independently of memory.

**Not implemented:** memory-per-thread budgeting and disk spilling. `--threads N` (default 1; tested 4) sets workers;
memory is not divided per thread.

## 5. Reporting and diagnostics

| key / flag | values | default | role |
|---|---|---|---|
| `--output-schema` / `output_schema` | legacy, sequence_first_v1, **sequence_first_v2** | legacy → preset v2 (legacy under the unrepaired view) | the report contract and the exit statuses (§6). It no longer implies a dosage stage |
| `--dosage-check` | flag | not given | runs the existing dosage audit and writes its diagnostics (`dosage_mode_source flag --dosage-check`). **It never changes the prediction** and adds run time. Without it, dosage is not run: DOSAGE_FILTER is NOT_EVALUATED (`DOSAGE_NOT_EVALUATED:dosage_mode_off`) |
| `--dosage-mode` (expert spelling) | off, audit, repair | off | `audit` = `--dosage-check`. A contradiction with `--dosage-check` (`off`, `repair`, or `--dosage-audit off`) is refused, exit 2. `repair` is refused under v2; the legacy repair, which may replace a block (`<prefix>.original.*` kept), needs `--output-schema legacy --dosage-mode repair`. Dosage is not a profile key |
| `compound_annotation` | off, **contract_v2** | off | report (`.sequence_compound.tsv`) |
| `compound_interval_declarations` | file | none | report; model-bound declarations, **not used in the tested configuration** |
| `block_diagnostics` | **off**, diagnostic_v1 | off | report (`.block_diagnostics.tsv`: census route, JP tie); identical science verified |
| `timing_detail` | 0, **1** | 0 | report (stage timers) |
| `library_sensitivity_audit` | 0, 1 | 1 | under `--library auto`, the 0.9×/1.1× audit behind `LIBRARY_SENSITIVE` |
| `emission_lookup_timing` | 0, 1 | 0 | measurement only |
| `-q`, `--quiet` | | | |

## 6. Outputs: what you get

| outcome | exit (v2) | what exists |
|---|---|---|
| **completed prediction** | 0 | `.sequence_locus.tsv` `prediction PREDICTION`, `completion completed_ranking` (even when uncertain); full block tables |
| **PROVISIONAL** | 3 | the best **fully evaluated** candidate after a budget stop or declared cap (early-incumbent checkpoint). It is **not a completed search** and carries no decision set and no confidence |
| **no prediction** | 1 | `prediction NO_PREDICTION` with `no_call_reason` (for example `distinct_tie`, `neighbourhood_bound`) |
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
- **The persisted prepared-locus index** (`.panel.pvpl`) and sequential multi-sample reuse: library/API-level, exercised only
  by development drivers (for example `genotype_prepared_locus_driver`). The tested configuration does not use an index.
- **Batch / multi-sample runs:** one locus and one sample per invocation.

**Long reads:** not implemented (`--preset long-read` is refused).

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
