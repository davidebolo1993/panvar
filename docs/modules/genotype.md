# Module `genotype`

CLI: `panvar genotype`. Experimental short-read genotyping; not a validated production caller.

## What it does

Genotypes one diploid sample at one locus using paired short reads and a pangenome graph.
The graph supplies a panel of haplotypes and a chain of context and variant blocks.
The prediction is two whole-locus sequences, with an unordered allele pair at each block.
Whole haplotypes can be panel paths or mosaics of panel block alleles; missing block content is not invented.

The short-read preset includes the accelerated scoring routes, bounded refinement, junction-aware proposals,
checkpointing and sequence-first reporting. An ordinary run needs neither a profile nor a declared library.
Optional `--dosage-check` adds diagnostics; it never changes the prediction.

Read the [algorithm](../algorithms/genotype.md) for scoring, search and numerical bounds.
Completion, uniqueness and read support are different questions. None is calibrated confidence.

## Required inputs

```bash
panvar genotype --graph locus.gfa \
  --reads1 sample.R1.fq.gz --reads2 sample.R2.fq.gz \
  --threads 4 --out results/sample
```

The examples use `panvar` on PATH; after a local build use `./build/panvar` instead.
Create the parent output directory before running.

| Input | Requirement |
|---|---|
| `--graph FILE` | A locus GFA with haplotype paths (`P` or `W` records). Paths must traverse the whole modeled locus; incomplete paths can cause a model refusal. |
| `--reads1 FILE --reads2 FILE` | Paired FASTQ, optionally gzip-compressed, with mates in matching order. Long reads are not supported. |
| `--reference-path NAME` | Sets the coordinate frame and orientation. Optional only when the graph's `RS:Z` sample identifies exactly one path. Otherwise the run refuses and lists candidate paths. |
| `--out PREFIX` | Output prefix, not a directory; default `genotype`. Choose a separate prefix for every sample and locus. |

### Library: automatic by default

Omitting `--library` is identical to `--library auto`. The program estimates insert mean and standard deviation,
and fragment-start intensity. At least 500 decisively anchored concordant pairs are required. Too few pairs,
discordant or contaminated inserts, multiple insert modes, or an invalid depth fit cause `LIBRARY_REFUSED`:
no prediction and no silent fallback.

Automatic estimation does not learn every model parameter. The short-read values of `support_sigmas` (4),
`substitution_error` (0.002), `outlier_weight` (0.001) and the gap model remain fixed.
With the default sensitivity audit, the caller also checks the effect of changing intensity by ±10%.

Use `--library FILE` when parameters are independently known or exact reproducibility is needed.
The file contains one whitespace-separated `key value` pair per line:

```text
insert_mean_bp 350
insert_sd_bp 40
support_sigmas 4
substitution_error 0.002
outlier_weight 0.001
fragment_start_intensity 0.0125
```

These are illustrative values, not a recommended library for arbitrary reads. All six keys are required exactly
once; `marker_background` is optional. Blank lines and `#` comments are allowed. Unknown, repeated or missing
keys, trailing content and invalid/non-finite values are refused. The gap model cannot be supplied in this file.

### Panel curation does not rebuild the graph

`--panel-curation FILE` removes the listed haplotype paths from the callable panel. Use one exact path name per
line; blank lines and `#` comments are accepted. This is an exclusion list, not a list of paths to retain.

All graph paths still define the block boundaries and trimming. Nodes and edges are not deleted. For example,
removing the only path carrying a deletion leaves its bubble in the graph, but the deletion allele is unavailable
unless a retained path also carries it. Shared allele sequences retain only their active sources. A model with no
usable active allele at a required block is refused. Thus leave-one-out curation does not hide the held-out path's
contribution to the graph frame.

## Key options

`panvar genotype --help` lists everyday options. `panvar genotype --help-advanced` lists all options and expert
profile keys, including accepted values, defaults and their effects. The complete generated table is below.

### Panel index: build once, reuse across samples

```bash
panvar genotype --graph locus.gfa --reference-path REF \
  --build-panel-index locus.pvpl

panvar genotype --graph locus.gfa --reference-path REF --panel-index locus.pvpl \
  --reads1 sample.R1.fq.gz --reads2 sample.R2.fq.gz --threads 4 --out results/sample
```

Use the same curation, bubbles and marker settings for the build and sample commands when overriding them.
Build mode takes no reads, library, sample output prefix or dosage check, and never overwrites an existing index.
It publishes the index atomically and prints its versions, identity, source provenance, size, hash and timings.

The index stores the sample-independent model and marker panel. It does not store a sample's reads, library,
depth, candidates or search state. Loading still reads and hashes the graph; it is not an index-only run.
Integrity and compatibility are checked before the output prefix is touched. A missing, corrupted or incompatible
index is an error, not a request to rebuild or fall back silently.

Rebuild after changing graph content, reference, curation, bubbles, minimum site size or marker `k`/`s`, or when
the format/preparation version becomes incompatible. Reads, library, threads, dosage checking and sample/search
settings do not invalidate it. A different Git revision alone does not invalidate it; producer provenance and
compatibility identity are separate. Loading is bounded to a 1 GiB file and 3.5 GiB accounted decoding.

Illustrative back-to-back development runs, four threads, automatic library:

| Locus | Index | Build / load | Cold / indexed wall time |
|---|---|---|---|
| C4 | 44.5 MB | 2.0 / 0.47 s | 12.7 / 11.6 s |
| CYP2D6 | 39.2 MB | 2.2 / 0.43 s | 22.9 / 21.6 s |
| ACOT | 160.2 MB | 7.0 / 1.71 s | 34.5 / 32.0 s |

These are example timings, not a promised speedup. Other reserved ACOT 5× calls took 80–198 s; inputs and search
work matter. In the index example ACOT peak memory fell from 3.23 to 1.88 GB.

### Model and search controls

| Option | Meaning |
|---|---|
| `--mismatch-band N` | Maximum edits per mate, counting substitutions and gap bases. It is not a block-size tolerance. |
| `--mosaic-candidates K` | Maximum initial off-panel mosaic diplotypes sent to full scoring, in addition to panel pairs. Refinement can evaluate further mosaics. |
| `--proposal-state-limit N` | Beam width while constructing candidates block by block: how many partial allele-prefix pairs survive each layer. These are not N finished proposals. |
| `--adaptive-repair on\|off` | Bounded single-block refinement, on by default. Each round examines proposals across blocks and makes the best certified improving move, rather than choosing the first weak block. |
| `--repair-local-pairs N` | Local allele-pair shortlist per block; both homologue assignments are tested. It is not exhaustive over every possible pair. |
| `--repair-rounds N` | Maximum refinement rounds. Zero, like repair off, reports the already-scored union using the same v2 output contract. |
| `--switch-probability P` | Copying model's donor-redraw probability per block boundary, not a biological recombination rate. |
| `--min-variant-bp N` | Site construction: retain a site when some path spans at least N bp inside it. Not a final-call size filter. |

For example, a beam width of 256 retains up to 256 partial candidates at each layer, then up to 128 complete
non-panel proposals are fully scored. Search breadth and the number of final proposals therefore need separate
controls. Panel diplotypes remain the baseline.

Unrefined runs keep the v2 reporting rules, but do not run the representation guard. Their `genotype_check`
is therefore `NOT_ASSESSED`; do not interpret disabling refinement as an equivalent faster preset.

### Optional expert profile

The preset is complete. `--profile FILE` is only for overrides; no locus-specific profile is required.
Precedence is explicit flag > profile > preset > built-in default, regardless of argument order.
Every effective setting and its source appears in the diagnostics settings table.

A profile uses literal TAB-separated `key<TAB>value` rows and requires one `profile<TAB>NAME` row.
For example, this creates a two-setting file with actual TABs:

```bash
printf 'profile\tmy_limits\ntime_limit_seconds\t900\nfootprint_limit_bytes\t4294967296\n' > limits.tsv
panvar genotype --graph locus.gfa --reads1 sample.R1.fq.gz --reads2 sample.R2.fq.gz \
  --threads 4 --profile limits.tsv --out results/sample
```

Blank lines and `#` comments are allowed. Repeated profiles, duplicate keys, unknown keys and malformed values
are refused. A profile cannot select an output contract or dosage mode. It may only restate
`output_schema sequence_first_v2` and the internal historical `phase_audit_candidates 0`.
When a flag removes a preset feature's prerequisite, preset-sourced dependents are recorded as inactive;
an explicitly incompatible flag/profile setting still refuses.

### Resources

| Control | Behavior |
|---|---|
| `footprint_limit_bytes` | Default 6 GiB. Sampled every 100 ms, not a hard allocation ceiling; short spikes can exceed it. |
| `time_limit_seconds` | Default 0: no time limit. A profile can set a bound. |
| Cache/store budgets | Bound retained data, not total process memory. Reuse caches can fall back to recomputation; mandatory stores can refuse when too small. |
| Work caps | Bound executed origins, state classes, placements and neighbours. A cap can leave a provisional result or no prediction; more RAM does not automatically raise it. |

No memory-per-thread partitioning or disk spilling is implemented. Exact computational routes preserve scores
when their runs complete; resource limits can still change completion and which checkpoint survives.

<details>
<summary>Complete option and profile-key inventory</summary>

Generated from the parser and preset; the regression suite checks this table. `inference` can change predictions;
`computation` changes implementation cost; `reporting` changes annotations/files; `completion` changes whether/how
the run finishes. Entries marked development or legacy are not ordinary user controls; removed spellings are
listed as refusals, not recommendations. The profile keys use their literal spelling.

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
| resources | `repair_origin_cache_bytes` | — | repair_origin_cache_bytes | N >= 0 | 768 MiB (built-in 384 MiB) | bytes | computation | expert profile only |
| resources | `repair_tie_seed_cache_bytes` | — | repair_tie_seed_cache_bytes | N >= 0 | 0 | bytes | computation | expert profile only |
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
| reporting | `reporting contract (output_schema)` | — | output_schema | sequence_first_v2 | sequence_first_v2 (built-in legacy) | - | reporting | fixed |
| reporting | `--representation-guard on\|off` | — | representation_guard | on \| off | on (built-in off) | - | reporting | supported |
| reporting | `compound_annotation` | — | compound_annotation | off \| contract_v2 | contract_v2 (built-in off) | - | reporting | expert profile only |
| reporting | `compound_interval_declarations` | — | compound_interval_declarations | a file | — | - | reporting | expert profile only |
| reporting | `block_support_rule` | — | block_support_rule | legacy \| full_data_likelihood_v1_experimental \| full_data_likelihood_v2_experimental | legacy | - | reporting | expert profile only |
| reporting | `library_sensitivity_audit` | — | library_sensitivity_audit | 0 \| 1 | 1 | - | reporting+computation | expert profile only |
| reporting | `timing_detail` | — | timing_detail | 0 \| 1 | 1 (built-in 0) | - | reporting | expert profile only |
| reporting | `block_diagnostics` | — | block_diagnostics | off \| diagnostic_v1 | off | - | reporting | expert profile only |
| reporting | `emission_lookup_timing` | — | emission_lookup_timing | 0 \| 1 | 0 | - | reporting | expert profile only |
| legacy | `--max-candidates` | — | — | - | refused | - | - | legacy |
| legacy | `phase_audit_candidates` | — | phase_audit_candidates | 0 | 0 | candidates | reporting | internal |
| legacy | `--output-schema` | — | — | - | refused | - | - | removed |
| legacy | `--dosage-mode, --dosage-audit` | — | — | - | refused | - | - | removed |
| legacy | `--phase-audit-candidates` | — | — | - | refused | - | - | removed |
| legacy | `--repair-origin-cache-bytes, --repair-tie-seed-cache-bytes` | — | — | - | refused | - | - | removed |
| legacy | `--dosage-allocation-test, --library-diagnostic` | — | — | - | refused | - | - | removed |
| legacy | `diagnostic_* and test_* profile keys` | — | — | - | off | - | - | development |
<!-- option-inventory:end -->

</details>

## Outputs

For `--out results/sample`, the prefix is `results/sample`, not `results/sample/`:

```text
results/
├── sample.sequence_locus.tsv
├── sample.sequence_blocks.tsv
├── sample.sequence_haplotypes.tsv
├── sample.sequence_haplotypes.fa       # only when either haplotype is off-panel
└── sample.diagnostics/
    ├── sample.settings.tsv
    ├── sample.mosaic_status.tsv
    ├── sample.mosaic_candidates.tsv
    ├── sample.mosaic_scores.tsv
    ├── sample.equivalence_set.tsv
    ├── sample.locus_genotype.tsv
    ├── sample.block_genotypes.tsv
    └── …                             # optional library, dosage and search diagnostics
```

Files depend on the outcome and enabled checks; the tree is not a guarantee that every diagnostic exists.
The primary report contract is fixed (`sequence_first_v2`). No `--output-schema` selection is needed or accepted.

### Primary files: start here

| File | Read it for |
|---|---|
| `sample.sequence_locus.tsv` | Key/value run outcome: `prediction`, `completion`, `exit_code`, sequence resolution, assessments, coverage and reasons. It is also the outcome record for provisional/failed runs whenever one can be written. |
| `sample.sequence_blocks.tsv` | One row per modeled block, contexts included: genotype, block resolution, possible genotypes, assessments and pool-conditional error-mass bounds. Allele ids are local to this model. |
| `sample.sequence_haplotypes.tsv` | Two representative haplotypes when a prediction/checkpoint exists: sequence SHA-256, length, panel origin, paths/orientations, extraction reference, graph hash and model fingerprint. `hap1`/`hap2` order is canonical, not parental identity. |
| `sample.sequence_haplotypes.fa` | Representative sequence pair, both records, only when at least one is off-panel. A representative of an ambiguous or provisional result does not resolve that uncertainty. |

An illustrative excerpt of the locus record (TAB-separated, not a single wide row):

```text
key	value
prediction	PREDICTION
completion	completed_ranking
exit_code	0
sequence_resolution	RESOLVED
genotype_check	NONE_DETECTED
dosage_check	NOT_REQUESTED
read_support	PRIOR_DEPENDENT
phase	PRIOR_ONLY
```

When both sequences exactly match active panel paths, no FASTA is written intentionally:
`fasta_reason PANEL_REPRESENTED_PAIR`. The haplotype table provides extraction references instead.
The repository utility checks every listed `P`-path extraction against its recorded hash and length:

```bash
python3 scripts/extract_genotype_haplotypes.py --check locus.gfa results/sample.sequence_haplotypes.tsv
```

To write one sequence, copy its extraction reference from the table:

```bash
python3 scripts/extract_genotype_haplotypes.py locus.gfa \
  'gfa_path=PATH;interval=0-LENGTH;orientation=forward;whole_locus_path' --fasta haplotype.fa
```

Replace PATH and LENGTH; use the recorded orientation, not a guessed one. The utility currently reads `P` records;
it does not extract directly from `W`-only graphs. Off-panel sequence FASTA contains both haplotypes, not an
assembly or a dosage correction. `NO_PREDICTION` has no prediction FASTA; a failure to write is an error, not the
intentional panel-represented case.

### Outcomes and exit codes

| `prediction` / outcome | Exit | Interpretation |
|---|---|---|
| `PREDICTION` | 0 | A completed ranking, possibly with unresolved sequence/genotype ties. Not a correctness guarantee or a global optimum. |
| `PROVISIONAL` | 3 | Best fully evaluated checkpoint retained after a cap/stop. Search is incomplete; no completed decision set or confidence. |
| `NO_PREDICTION` | 1 | No retained candidate, e.g. an unsupported model or refused automatic library. Read the reason. |
| Time/footprint stop without a checkpoint | 124 / 125 | No prediction; minimal status and a sequence-locus outcome record when writable. |
| Argument, input or output error | 2 | Refusal or runtime/reporting failure. A result file left by a reporting failure does not make the run successful. |

Only argument/profile errors, declared-library validation, and index compatibility preflight leave an existing
prefix untouched. Runtime graph/read/model/write failures are not transactional. A run removes its own previous
report names, including old flat diagnostic names, but not unrelated files. Diagnostics are not deleted automatically.
The diagnostic directory is created when the run first touches its prefix. Files cannot be guaranteed when the
filesystem itself is unwritable. Check process exit as well as the locus `exit_code` and status.

### Assessments: independent checks, not one PASS/FAIL

The primary view is `assessment_v1`. Locus fields are rows in `.sequence_locus.tsv`; block fields are columns
in `.sequence_blocks.tsv`. Read findings and missing checks together. A missing check is not evidence against the
call, but neither is it a passed check. With no prediction assessments are `NA`.

#### Sequence resolution

This describes agreement among sequence groups in the numerically inseparable decision set, conditional on the
explored candidates. It is separate from support and from how DNA is labelled across blocks.

| `sequence_resolution` | When it occurs |
|---|---|
| `RESOLVED` | One unordered sequence pair remains. Its DNA may still have several block labellings. |
| `GENOTYPE_RESOLVED_SEQUENCE_SET` | Several sequence pairs remain, but every reported block has the same unordered genotype. Phase/sequence assignment is not resolved. |
| `SEQUENCE_TIED` | Several pairs remain and at least one block genotype cannot be uniquely identified. Exit can still be 0. |
| `NOT_DETERMINED` | A provisional checkpoint, not a completed decision. |
| `NA` | No prediction. |

Block `block_resolution` is `IDENTIFIED`, `LABELS_AMBIGUOUS`, `GENOTYPE_UNRESOLVED`, or `PROVISIONAL_LABEL`.
`LABELS_AMBIGUOUS` means explored labellings of the same sequence pair disagree at the block. For example, bases
can be assigned to either of two neighbouring blocks while the concatenated DNA stays identical. This does not
by itself mean the sequence is wrong. Consult the representation interval diagnostic rather than treating each
label difference as an independent event.

#### Genotype check

`genotype_check` combines the representation guard and automatic-library sensitivity findings, only where the
underlying block evidence is verified to describe the published alleles.

| Value | When it occurs |
|---|---|
| `NONE_DETECTED` | The representation guard produced an assessed result without a representation conflict or recorded library sensitivity. Not proof of correctness; a missing library check is still listed separately. |
| `CONFLICT` | A representation conflict and/or library-sensitive genotype was recorded. Reasons name the finding. |
| `NOT_ASSESSED` | The guard did not run, its result is unavailable, or its block alleles cannot be verified against this prediction. Examples: unrefined runs, provisional results, legacy-tied predictions without matching block rows. |

A representation conflict is disagreement between evidence components: markers and the expected-count penalty
favor another block genotype, while fragment placement favors the reported one, above the declared guard margin.
It is not synonymous with ambiguous block labelling. Library sensitivity means that the ±10% intensity audit
changes the block genotype. Neither check independently estimates copy number; that is the dosage check below.

#### Dosage check

Enable with `--dosage-check`. It compares copy features implied by the reported alleles against read-derived
counts and their uncertainty. It is not simply a vote for whichever alternative haplotypes have the most markers.
Feature counts can depend on content, and matching counts need not establish allocation across homologues.
The [algorithm's dosage section](../algorithms/genotype.md#9-audit-dosage-without-changing-the-prediction)
gives the trigger, suspect and power criteria.

| `dosage_check` | When it occurs |
|---|---|
| `AGREES` | The assessed reported copy features are consistent with the count evidence. This does not prove the sequence or allocation. |
| `SUSPECT` | Count disagreement reaches the suspect criterion but is not reported as a triggered contradiction. |
| `DISAGREES` | A count contradiction triggers, whether or not a possible correction was found. The prediction is kept. |
| `UNAVAILABLE` | Requested but not usable: audit failure/refusal, underpowered features, missing evaluation, or unverified reported alleles. Read `dosage_check_reason` / `not_assessed`. |
| `NOT_REQUESTED` | The flag was not supplied. |
| `NOT_APPLICABLE` | The completed audit identifies no dosage-sensitive comparison for this block (or no applicable blocks at the locus). |

Alternative corrections are separate: block `dosage_alternative` and locus `dosage_alternatives`. Values can
describe a complete alternative, conditional content, a refused correction or no correction found. Even files
named `.dosage_repaired*` contain audit alternatives, not replacements for the primary prediction.
`--dosage-mode` and `--dosage-audit` are removed; no public dosage repair mode exists.

#### Read support

The current support rule compares the winner's fragment-placement term with explored rivals within numerical
bounds. It excludes the copying prior and does not use the exposure term for this contender comparison.
An informative fragment is a read pair, not two independent mates.

| `read_support` | When it occurs |
|---|---|
| `READS` | No explored fragment contender has a different unordered genotype at this block. This can be vacuous if the pool lacks useful rivals. |
| `PRIOR_DEPENDENT` | At least one explored rival is not clearly below the winner by fragment term and differs at this block. Does not mean no read informs the block, or that the prior alone caused selection. |
| `NOT_ASSESSED` | Support is missing or cannot be verified for the published genotype. Unknown is never converted to `READS`. |

#### Phase

`phase` is assessed once per locus, not independently for each block boundary. It concerns which alleles share
a homologue, not whether each unordered block genotype is correct. Canonical haplotype numbering is not parental phase.

| Value | When it occurs |
|---|---|
| `READ_SUPPORTED` | A unique, verified decision with read-supported phase under the explored-candidate rule. It does not establish that every boundary has a spanning molecule. |
| `PRIOR_ONLY` | The unique decision has fragment contenders, so the reads do not independently settle its phase under this rule. |
| `UNRESOLVED` | Several sequence pairs remain in the decision set. This takes precedence over the heterozygous-block applicability check. |
| `NOT_ASSESSED` | Provisional result, missing verification, or unavailable phase evidence. A legacy-tied but uniquely resolved v2 sequence can have this value. |
| `NOT_APPLICABLE` | A resolved prediction has fewer than two heterozygous blocks and no unknown genotypes. Two heterozygous blocks need not be adjacent. |

#### Coverage, findings and missing checks

Locus `genotype_check_coverage` and `read_support_coverage` count assessed / reported blocks.
`dosage_check_coverage` excludes `NOT_APPLICABLE` blocks from its denominator.
The summary gives precedence to an observed conflict, dosage disagreement/suspicion, or prior dependence;
unavailable blocks do not hide those findings. A positive result can cover only part of the locus.

For example, `DOSAGE_AGREES[1/2]` means one of two applicable blocks was assessed and agreed, not that both passed.
`not_assessed` could contain `DOSAGE:UNDERPOWERED@sv:1`. Findings such as `PRIOR_DEPENDENT@sv:3` remain visible too.
Always read the coverage fields themselves: the compact summary does not add brackets when zero checks ran.
Provisional records carry no completed support assessment; a requested dosage check may be `UNAVAILABLE`.

### Diagnostics: inspect when needed

Every name below is relative to `sample.diagnostics/` and prefixed with `sample`.
Diagnostics are stage-dependent; absence may mean a stage did not run. Check status/reasons, not just filenames.

| Suffix or family | Purpose |
|---|---|
| `.settings.tsv` | Effective `key`, `value`, `source`, including profile and preset provenance. |
| `.mosaic_status.tsv` | Detailed completion, refusal, resource, index, library and search record. |
| `.mosaic_candidates.tsv`, `.mosaic_scores.tsv` | Candidate labellings and score components for the evaluated pool. |
| `.equivalence_set.tsv`, `.sequence_alternatives.tsv` | Ranked sequence groups, decision-set membership, representative and alternatives. |
| `.representation_intervals.tsv` | Conservative enclosing intervals of explored synonymous block labellings. Not a new biological event call. |
| `.sequence_compound.tsv` | Compound interval annotations when enabled and declarations match; empty/header-only output can be legitimate. |
| `.locus_genotype.tsv`, `.block_genotypes.tsv`, `.genotypes.tsv` | Compatibility tables; the last contains variant blocks only. The assessment view is mirrored here, but these are not the primary interface. |
| `.library_estimate.tsv`, `.library_sensitivity*.tsv` | Automatic parameter estimates and sensitivity evidence where applicable. |
| `.dosage_status.tsv`, `.dosage_blocks.tsv`, `.dosage_census.tsv`, `.dosage_features.tsv`, `.dosage_feature_copies.tsv` | Audit completion, per-block conclusions, feature inventory and count evidence. |
| `.dosage_allocation*.tsv`, `.dosage_repair*.tsv`, `.dosage_provenance*.tsv`, `.dosage_*.fa` | Diagnostic alternative content/allocation, scoring, provenance and sequences when produced. Never substitute them for primary haplotypes automatically. |
| `.repair_rounds.tsv`, `.repair_memory.tsv`, `.tie_*.tsv`, `.block_diagnostics.tsv` | Search steps, resources, ties and optional diagnostic flags. |

Compatibility columns `genotype_filter`, `dosage_filter`, `phase_filter`, `call_summary` are superseded,
not a second recommended filtering system. They retain historical meanings: `NO_CALL` can sit beside a valid
v2 prediction, missing evidence can mask a finding, and dosage `PASS` can refer to an unapplied alternative.
Use `assessment_v1`; do not combine old PASS values into a reliable-call label. GQ/error mass remain conditional
on the explored pool and labellings; a printed zero error mass can be rounding, not certainty.

`scripts/plot_genotype_results.py` summarizes runs; `--help` describes inputs. Block error/GQ plotting utilities
also remain available. Truth-based summaries require separately supplied evaluation tables.

## Limitations

The best evaluated candidate is not a global optimum. Short reads can leave repeated-copy allocation and phase
unidentified; missing panel content and ungenerated competitors can produce confident-looking wrong calls.
No calibrated confidence, reliable-call threshold, automatic dropping rule or guarantee of always returning a
genotype is implemented. `NONE_DETECTED`, `AGREES` and `READS` describe checks, not validated accuracy.

The known CYP2D6 sv:9 audit refusal (alleles disagree on flanks or anchor trims) leaves dosage `UNAVAILABLE` and
warns even under `--quiet`; it does not change the prediction. Other loci may have analogous refusals.
Block checks for legacy-tied v2 predictions remain unassessed without verified matching block data.
Direct caller API runs without the command layer do not produce the primary assessment view. Concurrent
callers in one process are outside the validated contract; in-process sequential reuse is API-only.
The ACOT 30× reference/ordinary-route equivalence gap remains; the corresponding backend comparison was at 5×.
`--preset long-read` is refused. This release does not establish that long reads are the only possible remedy.

### Evidence and its scope

| Evaluation | What was observed |
|---|---|
| C4 and CYP2D6 | Development evidence only; no independent held-out accuracy validation. Large structural errors occurred. |
| Reserved ACOT | 18 prior-exposed families, simulated paired 100 bp reads at 5×, declared library, leave-one-out panel within the full-panel graph frame. All completed; 296 informative blocks: median/q90/q99/max excess 0/8/98/478 bp; 66.6% at the availability floor; zero ≥1 kb catastrophes. |
| Automatic library | Development checks at three loci reproduced declared-library predictions. Not independent accuracy validation. |

The ACOT families had been used for masked prior scoring and helped define the full graph frame, so they are not
pristine independent validation. Zero catastrophic families out of 18 still permits roughly a 15–17% upper
family-level rate bound. Two roughly 478 bp sv:9 errors were not flagged: both old genotype filters said PASS,
neither had a repeated-unit or dosage-disagreement warning. No current annotation combination establishes correctness.

Reported excess is an offline truth-based evaluation metric, not a runtime filter. For a block, sum edit distances
under the better whole-pair homologue assignment, then subtract the closest available active-model block-pair
distance (the availability floor). Normalized excess divides this by the larger of the truth and predicted diploid
block lengths. Attaining the floor means best available, not necessarily identical to truth. Whole-region distance
and QV are separate diagnostics; do not subtract a per-block floor from a whole-region distance.
All-empty blocks are not informative for accuracy, and normalization is undefined at zero diploid length.

Performance on real reads, other loci and broader coverages is not established. Computational checks ran on one
macOS arm64 machine (Apple clang 21, CMake 4.3.3, Python 3.14.7); portability is not validated.

## Example

Each invocation processes one sample and one locus; there is no multi-sample manifest CLI.
Reuse a panel index across separate invocations and keep separate output prefixes:

```bash
mkdir -p results
panvar genotype --graph locus.gfa --reference-path REF --build-panel-index locus.pvpl
for sample in S1 S2; do
  panvar genotype --graph locus.gfa --reference-path REF --panel-index locus.pvpl \
    --reads1 "reads/${sample}.R1.fq.gz" --reads2 "reads/${sample}.R2.fq.gz" \
    --threads 4 --dosage-check --out "results/${sample}"
  code=$?
  printf '%s\texit=%s\n' "$sample" "$code"
done
```

Replace the example graph, reference and reads with actual inputs. Inspect each process exit and its primary
locus record before combining results. A shell loop does not convert provisional or absent calls into completed ones.
For several loci, use a directory per locus and build a compatible index for each graph/curation combination.
