# Module `genotype`

> **Experimental short-read release candidate:** the authoritative usage guide (inputs, every flag and profile key, resources, outputs, filters, evidence and tested commands) is [docs/guides/genotype_short_read.md](../guides/genotype_short_read.md).

CLI: `panvar genotype`

## What it does

`panvar genotype` genotypes a diploid locus from paired short reads using a pangenome graph whose active paths form the haplotype panel. Each invocation genotypes one locus. The call is a diplotype: two whole-locus haplotypes, each a sequence of panel block alleles, reported per block as an allele pair.
- Candidates are every unordered pair of active panel haplotypes, plus a bounded number of mosaics that switch between panel haplotypes at block boundaries.
- Each evaluated candidate is scored by a paired-read likelihood under the declared read model, plus a Li–Stephens switching prior.
- A bounded repair search improves the best panel pair one block at a time and branches where moves tie. The call is decided over every candidate evaluated.
- Annotations then describe the call without changing it: a representation check, read support and phase support, and, under automatic library estimation, a library-sensitivity audit.
- An optional dosage stage tests the copy numbers the call implies. In repair mode it may replace a block with an eligible complete diplotype repair.

Algorithm and worked trace: [algorithms/genotype.md](../algorithms/genotype.md).

## Required inputs

- `--graph <graph.gfa>` (alias `--gfa`): the locus graph. Every path (`P` or `W` line) is a panel haplotype unless curated out, and must cross the whole locus; a path covering only part of it refuses the run.
- `--reads1 <r1.fq>`, `--reads2 <r2.fq>`: paired FASTQ, optionally gzipped, mates in the same order.
- `--out <prefix>` (alias `--out-prefix`; default `genotype`). Earlier `genotype` outputs under the prefix are removed first.
- The coordinate reference: `--reference-path`, or the graph's `RS:Z` sample, which must name exactly one path. Otherwise the run is refused and the candidate paths are listed.
- The read library:
  - **Omitted or `--library auto`:** the insert distribution and fragment-start intensity are estimated from the reads. An invalid estimate refuses the run (`LIBRARY_REFUSED`); there is no fallback.
  - **`--library FILE`:** the exactly reproducible route. One `key value` per line: `insert_mean_bp`, `insert_sd_bp`, `support_sigmas`, `substitution_error`, `outlier_weight`, `fragment_start_intensity`, and optionally `marker_background`, each once. Blank lines and `#` comments are allowed. Unknown, repeated or missing keys, trailing content and non-finite values are refused, with the file, line and key named. The gap model is fixed and cannot be declared.

```bash
panvar genotype --graph locus.gfa --reads1 R1.fq.gz --reads2 R2.fq.gz --threads 4 --out sample
```

This runs `--preset short-read`, the default preset. No profile is needed. The preset carries:
- the production algorithm and the short-read constants;
- the accelerated exact computational backends, with their cache budgets and capacities;
- junction proposals (`repair_junction_proposals experimental_v1`, a candidate-generation setting);
- the early evaluated incumbent;
- the `sequence_first_v2` report with the compound annotation;
- a 6 GiB footprint guard;
- the model's own repair-neighbourhood bound (`repair_max_neighbours model`).

It sets no time limit and no dosage stage. The usage guide lists every value. Settings are applied in a fixed order, whatever the order of the arguments: a flag over `--profile` over the preset over the built-in default. Every effective setting and its source go to `<prefix>.settings.tsv`.

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `--library auto\|FILE` | estimate the library or declare it | `auto` (omitted = auto) |
| `--dosage-check` | run the [dosage audit](#dosage-modes) of the reported alleles; never changes the prediction | not run (`dosage_check NOT_REQUESTED`) |
| `--panel-curation <file>` | remove named paths from the active panel ([below](#panel-curation)) | — |
| `--threads <N>` | worker threads; results do not depend on them | `1` |
| `--mosaic-candidates <K>` | the **initial** off-panel mosaics evaluated beside the panel pairs (the refinement may add more); `0` evaluates panel pairs only | `128` |
| `--bubbles <file>` | precomputed bubbles CSV instead of detected variant sites | detected |
| `--profile <file>` | optional expert overrides, never required: one `key<TAB>value` file with a `profile<TAB>name` row; given once (a repeated `--profile` is refused); a profile cannot select a reporting contract (`output_schema` may only restate `sequence_first_v2`) or set the dosage stage. `profiles/genotype_experimental_sr_v1_<locus>.tsv` reproduce the evidence configurations | — |
| `--adaptive-repair on\|off`, `--repair-rounds <N>` | the bounded coordinate refinement (local moves scored exactly within declared caps; not a global optimum) and its round cap; `off` or `0` rounds gives the unrefined view, with tie branching off: the already-scored union pool, reported with the v2 contract (its tie and phase outcomes follow the v2 rules; no representation guard). Contradictory combinations are refused | `on`, `32` |
| `--build-panel-index <file>` / `--panel-index <file>` | build the sample-independent panel index once / genotype against it ([below](#panel-index-build-once-genotype-many)) | — |
| `-q, --quiet` | no progress messages (warnings are still printed) | off |

`panvar genotype --help-advanced` lists every option, grouped by purpose, with its accepted values, effective short-read default, units and what it can change (inference, computation, reporting or completion). The guide's option inventory is the same table, generated from the parser and checked by the test suite.

### Panel index (build once, genotype many)

`--build-panel-index FILE` writes the locus's sample-independent preparation (structural model and marker panel) and
stops. It needs no reads or library, writes no sample output and never overwrites an existing file.

`--panel-index FILE` genotypes against it:
- the index's integrity, versions and identity are verified, and the sample's graph, reference, curation, bubbles and
  marker settings are checked against it, before the output prefix is touched;
- a missing, corrupt or incompatible index is a named refusal (exit 2), never an in-memory fallback;
- `--graph` is still read and hashed.

The index saves locus preparation only: everything that reads the sample still runs. Rebuild it when the graph,
reference, curation, bubbles, `--min-variant-bp`, `--kmer-size` or `--syncmer-s` change. Details:
[the usage guide, §8a](../guides/genotype_short_read.md).

### Panel curation

`--panel-curation` names paths to remove from the active panel, one per line. An allele carried only by removed paths is not callable, and no stage proposes, repairs or reports it. A sequence-identical allele on a retained path stays, with only retained provenance. The removed paths still define the graph's blocks. A panel that cannot spell the locus without them refuses the run (`MODEL_INCOMPLETE`).

### Dosage modes

The command has one dosage control, `--dosage-check`, which runs the audit. The modes below are the internal ones
(`dosage_mode` in `GenotypeCommandOptions`); `repair` is reachable only through the internal API, under the legacy
contract, for regression tests.

| mode | what it does |
|------|--------------|
| `off` | no dosage stage (the command's default; `dosage_check NOT_REQUESTED`) |
| `audit` | `--dosage-check`: the dosage census, tests and repair search; writes `<prefix>.dosage_*` and never changes the reported genotype. It is not necessarily cheaper than `repair` |
| `repair` | internal only: as `audit`, and applies each eligible complete diplotype repair; the caller's outputs are kept as `<prefix>.original.*` |

A dosage-sensitive block ends in one of four states:
- **Complete diplotype repair** (`COMPLETE_DIPLOTYPE_REPAIR`): the repaired content is complete and the diplotype is uniquely resolved within the declared model and search: the enumeration is complete and, for repeated units, the top allocation-equivalence set contains one unordered diplotype. "Complete" describes that scope, not biological truth. In repair mode it replaces the caller's block only when the run is a strict call, it is the only one at the locus, the reported haplotypes reproduce, and both repaired alleles belong to the active panel; otherwise the reason it was not applied is recorded.
- **Conditional dosage content** (`CONDITIONAL_DOSAGE_CONTENT`): the corrected copy numbers are determined but the diplotype is not. This includes an unresolved allocation-equivalence set (`dosage_reason=unresolved_allocation_equivalence_set`): several distinct arrangements of the same content fit equally, even when every one of them was enumerated, so the coverage scopes can all read `COMPLETE`. The reported alleles and FASTA stay the caller's. A representative arrangement (`serialization_representative=true`) is written only to the dosage diagnostics (`<prefix>.dosage_repaired.fa`, `.dosage_repaired_block_alleles.fa`, `.dosage_allocation_equivalence.tsv`), labelled `REPRESENTATIVE_ONLY_NOT_RESOLVED`.
- **Triggered but unrepaired:** the copy numbers are contradicted but no admissible correction was found.
- **Consistent, not evaluated, or not dosage-sensitive.**

## Outputs

**Layout (M15-L, a breaking change of output paths).**
- **Beside the prefix:** only the primary sequence reports. These are `<prefix>.sequence_locus.tsv`,
  `<prefix>.sequence_blocks.tsv` and `<prefix>.sequence_haplotypes.tsv`, plus `<prefix>.sequence_haplotypes.fa` when the v2
  FASTA rule writes it.
- **In `<prefix>.diagnostics/`:** every other report, under its usual name `<name>.<suffix>`. The tables below list them
  by suffix.
- **The outcome record** (`.sequence_locus.tsv`) is written beside the prefix whenever it is produced.
- **The internal API** (`GenotypeCommandOptions` without `diagnostics_directory`) keeps the flat layout.

| file | contents |
|------|----------|
| `<prefix>.locus_genotype.tsv` | the locus summary, one row: run status, call class, no-call reason, the haplotypes and unordered genotype when called, score and separation, dosage summary, the three headline filters with reasons and per-status counts |
| `<prefix>.block_genotypes.tsv` | per block, context included: allele pair, lengths, zygosity, read and phase support, conditional posterior and GQ, representation fields, block quality, block headline filters |
| `<prefix>.genotypes.tsv` | the variant blocks only, with each allele's panel source alleles and the same block headline filters |
| `<prefix>.mosaic_haplotypes.fa` | the two haplotypes of a strict call |
| `<prefix>.mosaic_status.tsv`, `<prefix>.settings.tsv` | the run record and every effective setting with its source |
| `.mosaic_candidates.tsv`, `.mosaic_scores.tsv`, `.equivalence_set.tsv`, `.repair_rounds.tsv`, `.tie_*.tsv` | the search record |
| `.library_estimate.tsv`, `.library_sensitivity*.tsv`, `.dosage_*`, `.original.*` | library estimate and audit; dosage stage |

What a run writes depends on its outcome:

| outcome | exit | locus summary and run record | block tables | FASTA |
|---|---|---|---|---|
| strict call | 0 | yes | yes, phased | yes |
| genotype-equivalent call | 0 | yes | yes, unordered | no |
| no-call (model, library, search or tie) | 1 | yes | no | no |
| time or footprint limit | 124 / 125 | minimal status (and, under `sequence_first_v2`, a `NO_PREDICTION` sequence record) | no | no |
| argument or declared-library error | 2 | refused before the prefix is touched: existing files are left as they were | no | no |
| input, model or output error during the run | 2 | a `FAILED` status; earlier run-owned outputs are already removed | no | no |

The process exit, `.sequence_locus.tsv` `exit_code` and the status's `sequence_exit_code` (and `exit_code` for stops and
failures) always agree. Only argument, profile and declared-library errors are refused before the prefix is touched;
errors found while running are not.

### Reading the call

**The primary view is `assessment_v1`** (sequence-first v2). It is written as rows in `.sequence_locus.tsv` and columns in
`.sequence_blocks.tsv`, mirrored in the locus and block tables:
- `prediction`, `sequence_resolution`;
- `genotype_check`, `dosage_check` (about the reported alleles; corrections are in `dosage_alternatives`), `read_support`
  and locus-level `phase`, each with coverage;
- `findings` and `not_assessed`.

It uses a legacy block fact only where that fact verifiably describes the published alleles. The headline filters below
are kept unchanged but superseded (`headline_filters_status`). Contract: `tests/analysis/m15g/CONTRACT.md`; usage: the
guide §7.0.

- **Strict call** (`COMPLETE`): one best pair of haplotypes. A *genotype-equivalent call* (`CONDITIONAL_GENOTYPE_EQUIVALENT`) has several best pairs that agree on every block's unordered genotype but not on phase, so it reports no FASTA and no phase.
- **No-call:** the reason is in `no_call_reason`, for example `distinct_tie`, `tie_branching_incomplete` or `neighbourhood_bound`.
- **Panel alleles:** every reported block allele is an active panel allele, and `allele*_sources` names its panel sources. A whole haplotype may still be off-panel as a mosaic.
- **FASTA headers:** they carry `phase_support` and `representation=`, which say how far the `hap1`/`hap2` order is read-supported. An applied complete repair adds `dosage_repair_role=RESOLVED`. A serialization representative never reaches this file: integration refuses any result with `serialization_representative=true` (`unresolved_allocation_equivalence_set`), whatever its class.
- **Dosage provenance:** sequences produced by the dosage stage are compared with the panel in `<prefix>.dosage_provenance.tsv`.

**Headline filters.** There is no single PASS/FAIL. `call_summary` is `GENOTYPE_<s>;DOSAGE_<s>;PHASE_<s>`. Each filter takes the most severe applicable status per unit, ranked NO_CALL > NOT_EVALUATED > LOW_QUALITY > PASS > NOT_APPLICABLE, and the locus takes the worst unit, with reasons as `REASON@ids`.

| filter | unit | PASS | other outcomes |
|---|---|---|---|
| GENOTYPE_FILTER | block | representation check passed, no library sensitivity, and the dosage stage, if run, left no exact diplotype unresolved | LOW_QUALITY for a representation conflict, a library-sensitive block, or an unresolved exact diplotype (conditional content, including an unresolved allocation-equivalence set; triggered-unrepaired; a complete repair not applied); NOT_EVALUATED when the representation check or the library audit could not complete |
| DOSAGE_FILTER | block | consistent; a complete repair applied, or available in audit mode; or conditional content whose content is complete (the reported alleles stay the caller's, so GENOTYPE_FILTER is LOW_QUALITY) | LOW_QUALITY for a complete repair not applied, suspect or triggered-unrepaired; NOT_EVALUATED for mode `off`, an incomplete stage, or underpowered; NOT_APPLICABLE when not dosage-sensitive |
| PHASE_FILTER | adjacent pair of heterozygous blocks | read-supported phase | LOW_QUALITY `MODEL_ONLY_PHASE` for prior-only phase; NO_CALL for unresolved phase; NOT_APPLICABLE next to a homozygous block |

- **A no-call** is NO_CALL in all three filters. A run that stopped before its blocks existed reports one `locus` unit per filter.
- **What PASS means:** the computational and representation checks passed. **PASS is not calibrated reliability,** and no
  combination of PASS values establishes an exact or biologically correct genotype. In the reserved ACOT evaluation (M15-B),
  two ~478 bp sv:9 length errors were GENOTYPE_FILTER PASS
  ([guide §7](../guides/genotype_short_read.md#7-decision-guide-filters-ambiguity-audit-and-support)).
- **Component fields:** `filter`, `genotype_support`, `quality_filter`, the library-sensitivity fields and `dosage_stage_filter` remain in the tables. `quality_filter` also requires read support, which GENOTYPE_FILTER does not.

`scripts/plot_genotype_results.py` summarises many runs as tables, figures and a short report. It is truth-free unless truth tables are given; see its `--help`. `scripts/plot_genotype_block_error.py` draws ordered-block raincloud figures of excess and normalised excess from a block-evaluation table; see its `--help`. `scripts/plot_genotype_block_gq.py` draws the same blocks as points coloured by conditional GQ; see its `--help`.

### Sequence-first output (experimental, opt-in)

**How to enable it.** Through the internal API only (`output_schema sequence_first_v1` in `GenotypeCommandOptions`): since M15-J no flag or profile selects a reporting contract. The command reports v2 (below).

**What it adds.** It reports the caller's existing decision by sequence pair instead of by block label. It changes no score, search, legacy output, dosage or phase decision. It needs the adaptive repair; other routes are refused.

**Outputs:**
- `<prefix>.sequence_locus.tsv`;
- `<prefix>.sequence_blocks.tsv`;
- `<prefix>.sequence_haplotypes.fa`, a representative marked as such;
- `<prefix>.representation_intervals.tsv`;
- `<prefix>.sequence_alternatives.tsv`.

**The prediction.** After a completed ranking it is a representative: the canonical member of the numerically inseparable decision set, chosen by SHA-256 of the sequences. It is not necessarily the strictly highest floating-point score. Labellings that spell the same pair are pooled. `representative_labelling_hap1` spells FASTA hap1, and interval `hapA` is the same sequence (by SHA-256). The status is one of:
- `RESOLVED`: one sequence pair in the decision set (not a confidence claim);
- `GENOTYPE_RESOLVED_SEQUENCE_SET`;
- `SEQUENCE_TIED`.

**Block status:** `IDENTIFIED`, `REPRESENTATION_AMBIGUOUS` or `GENOTYPE_UNRESOLVED`. It is computed over the *explored* labellings only (`projection_coverage=EXPLORED`), so an unexplored labelling of the same sequence may exist.

**Identified blocks** carry conditional error-mass bounds, with each sequence pair counted once. The existing GQ is kept as `legacy_representative_GQ`.

**Representation intervals** enclose every block where explored labellings of one sequence pair differ, including unchanged context between coupled alternatives. They list only observed label tuples, and they are not claimed minimal.

**Stops at a declared cap.** A round or cumulative cap after a committed round reports `PROVISIONAL`, from the last exactly scored incumbent, with no decision set and no confidence.

**When nothing is reported.** Watchdog stops, refusals and invalid input give no prediction. The five reports belong to the run: a later run at the same prefix, or a watchdog stop, removes them.

**Dosage.** The report and its FASTA describe the caller's decision *before* the dosage stage (`dosage_scope`, `post_dosage_sequence`, `dosage_mode`). When a dosage repair (internal API, legacy contract only) changes the legacy reported genotype, this FASTA is not the post-dosage sequence.

**Exit status under this schema:**
- 0: a prediction, even when uncertain;
- 3: provisional;
- 1: no prediction;
- 2: the report could not be written completely (nothing partial is left; `sequence_publication FAILED`);
- 2: errors.

**Confidence.** `structure_confidence` is `NA` (`UNCALIBRATED`).

### Sequence-first output v2 (experimental; the short-read preset's report)

**How it is selected.** The short-read preset sets `output_schema sequence_first_v2`, recorded in `.settings.tsv`; it is the
command's only reporting contract (M15-J: `--output-schema` removed; a profile may only restate `sequence_first_v2`).

**Unrefined runs** (`--adaptive-repair off` or `--repair-rounds 0`; M15-K) report v2 too. The pool is the union candidates
of the ranking already computed, with their scores, sequence digests and priors. The same `decide_pooled_genotype` decides
over it: nothing is scored, searched or refined beyond the union pass.
- **Outcomes follow v2, not the historical unrefined rules.** A genotype-equivalent tie is a prediction (formerly a block
  consensus, exit 1); a tie that disagrees at a block is `SEQUENCE_TIED`; phase is the pooled decision's.
- **No representation guard runs** on this route, so `genotype_check` is `NOT_ASSESSED`.
- **The legacy tables are written as before.**
- **Internal, legacy contract only:** the historical K-prefix phase audit and development prefix ladders.

`sequence_first_v1` and `legacy` remain internal, unchanged.

**What it adds to v1:**
- **Panel membership.** Each predicted haplotype, completed or provisional, is classified by exact sequence identity, in the locus frame, against the active (curated) panel.
- **Extraction references.** `<prefix>.sequence_haplotypes.tsv` lists each haplotype's SHA-256, length, origin (`PANEL`/`OFF_PANEL`), every sequence-equivalent panel path with its orientation, and an extraction reference. The reference has the form `gfa_path=…;interval=0-<length>;orientation=forward|reverse_complement`, applied to the whole locus path. The table also carries the GFA SHA-256, the model fingerprint and the selected block path.
- **The FASTA rule.**
  - When both haplotypes are panel-represented, no FASTA is written (`fasta_reason PANEL_REPRESENTED_PAIR`).
  - When either is off-panel, one FASTA holds both, each record annotated with its origin (`OFF_PANEL_PAIR`).
  - No prediction means no FASTA (`NO_PREDICTION`).

  An intentional absence is recorded as such, never as a failed publication.
- **The rule covers the whole prefix.** `.sequence_haplotypes.fa` is the only prediction FASTA. The legacy `.mosaic_haplotypes.fa` is removed once the dosage stage and the locus table have read it, and so is every empty dosage FASTA. A non-empty dosage FASTA holds audit diagnostics (repair-candidate representatives), so it stays. In `.locus_genotype.tsv`:
  - `main_fasta` links `.sequence_haplotypes.fa` when that file exists, and is `NA` otherwise;
  - `original_fasta` and the haplotype record ids are `NA`;
  - `dosage_representative_fasta` links only a dosage FASTA that stays.

  Every other legacy table and audit diagnostic is unchanged.

**Dosage is diagnostic under v2:**
- **default:** not run, and reported as not evaluated. v2 no longer implies the audit;
- **`--dosage-check`:** the audit, which never changes the prediction (the only public dosage control; `--dosage-mode`
  and `--dosage-audit` were removed in M15-J);
- **`repair`:** not reachable from the command line or a profile. Through the internal API it is refused under v2;
  historical repair results are reproduced only through the internal API with the legacy contract.

The command line, `execute_genotype_command` and the prepared-locus command route all resolve this through one function, `resolve_dosage_contract`. An API caller declares a mode with `dosage_mode_explicit` or a non-default `dosage_mode`. `dosage_mode_source` records the setting actually used. `run_genotype` and `run_genotype_prepared` call the caller alone, with no dosage stage, and record `dosage_mode not_run`.

The audit never changes the prediction or its sequences. It is not run when the legacy classification is a no-call: `dosage_status.tsv` then says NOT_RUN.

**Extraction:** `tests/analysis/s1d/extract_panel_haplotype.py --check <graph.gfa> <prefix>.sequence_haplotypes.tsv` extracts every listed panel path and checks its hash and length.

### Early evaluated incumbent (experimental; on under the short-read preset)

**How it is selected.** The short-read preset sets `early_incumbent experimental_v1`, or the profile key sets it explicitly.
- **Requirements:** it needs `sequence_first_v2` and the adaptive repair.
- **When they are missing:** a preset-sourced setting is switched off and recorded inactive; an explicit setting is refused.
- **When off:** it changes nothing.

**What happens.**
- **Selection.** After the markers, before the mosaic proposal and the union pass, one panel diplotype is selected. It is the argmax of the existing proposal marker emissions plus the exact path prior over the panel pairs, with ties going to panel order. This uses table lookups only and is not a winner claim.
- **Evaluation.** That pair is fully evaluated by the unchanged scorer, through the union's own route and caps, with the ranking's sequence prior.
- **Checkpoint.** It is checkpointed as a PROVISIONAL report that follows the v2 panel-reference and conditional-FASTA rules. The `checkpoint_*` rows give its stage, score components, evaluated count and semantics: fully evaluated; the best among the candidates evaluated so far; search incomplete. A checkpoint is never a global optimum, a calibrated claim or a confidence.
- **Replacement.** The checkpoint is replaced only by a better fully evaluated candidate: the union's best once the union pass completes, then any repair-scored candidate, committed at the round commit. The union check requires the union pass to reproduce the early candidate bit for bit.

**When the budget stops the run:**

| situation | outcome |
|---|---|
| time or footprint watchdog after a checkpoint | the checkpoint is kept as PROVISIONAL (exit 3); only its locus record is rewritten with the stop |
| watchdog before any checkpoint | typed NO_PREDICTION, `budget_stop_before_first_evaluated_candidate` (exit 124 or 125) |
| declared cap after a checkpoint (proposal budget, union cap, repair cumulative or round cap, neighbourhood bound) | the best fully evaluated candidate, as PROVISIONAL |
| scoring inconsistency, checkpoint publication failure or exception | explicit failure; no checkpoint is offered as a fallback |

**After an external SIGKILL.** The last locus record is the commit record. It names its haplotypes' SHA-256 and labellings. Content files that disagree with it belong to an uncommitted later checkpoint. The status then still reads RUNNING.

**Completed runs and accounting.** Scores, ranking and outputs are unchanged. The early evaluation is accounted in the status (`early_incumbent_*`, `checkpoint_*`), and the run-wide combination tally includes its work.

## Limitations

- Only paired short reads are supported.
- Reported block alleles are active panel alleles. A sequence absent from the active panel cannot be reported as an authoritative panel allele; the caller may select a representable alternative or make a no-call.
- GQ and the conditional posterior are conditional on the evaluated candidates and are not externally calibrated.
- Conditional dosage content corrects copy numbers, not the exact diplotype. When several arrangements of repeated units fit equally, which homologue carries each unit stays unresolved, and the caller's alleles are reported.
- Resource use grows with the active panel, the reads and the repair search. A numeric `--repair-max-neighbours` below a neighbourhood ends in a no-call (`neighbourhood_bound`). The preset's `model` bound records the model's own maximum and compares no external cap.
- No time limit is set by default. The 6 GiB footprint guard is sampled, not an allocation ceiling.

## Example

```bash
panvar genotype --graph locus.gfa --reads1 R1.fq.gz --reads2 R2.fq.gz --out sample --threads 4
panvar genotype --graph locus.gfa --reads1 R1.fq.gz --reads2 R2.fq.gz --out sample \
  --library library.txt --dosage-check --threads 4
```
