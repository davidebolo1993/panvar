# Module `genotype`

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
- The coordinate reference: `--reference-path`, or the graph's `RS:Z` sample, which must name exactly one path.
- The read library:
  - **Omitted or `--library auto`:** the insert distribution and fragment-start intensity are estimated from the reads. An invalid estimate refuses the run (`LIBRARY_REFUSED`); there is no fallback.
  - **`--library FILE`:** the exactly reproducible route. One `key value` per line: `insert_mean_bp`, `insert_sd_bp`, `support_sigmas`, `substitution_error`, `outlier_weight`, `fragment_start_intensity`, and optionally `marker_background`, each once. Blank lines and `#` comments are allowed. Unknown, repeated or missing keys, trailing content and non-finite values are refused, with the file, line and key named. The gap model is fixed and cannot be declared.

```bash
panvar genotype --graph locus.gfa --reads1 R1.fq.gz --reads2 R2.fq.gz --out sample
```

This runs `--preset short-read`, the default preset. Settings are applied in a fixed order, whatever the order of the arguments: a flag over `--profile` over the preset over the built-in default. Every effective setting and its source go to `<prefix>.settings.tsv`.

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `--library auto\|FILE` | estimate the library or declare it | `auto` |
| `--dosage-mode off\|audit\|repair` | the [dosage stage](#dosage-modes) | `off` |
| `--panel-curation <file>` | remove named paths from the active panel ([below](#panel-curation)) | — |
| `--threads <N>` | worker threads; results do not depend on them | `1` |
| `--mosaic-candidates <K>` | off-panel mosaics evaluated beside the panel pairs; `0` evaluates panel pairs only | `128` |
| `--bubbles <file>` | precomputed bubbles CSV instead of detected variant sites | detected |
| `--profile <file>` | a validated `key<TAB>value` settings file | — |
| `--adaptive-repair on\|off`, `--repair-rounds <N>` | the repair search and its round cap; `off` or `0` rounds gives the unrepaired view, with tie branching off. Contradictory combinations are refused | `on`, `32` |
| `-q, --quiet` | no progress messages | off |

`panvar genotype --help-advanced` lists the scientific flags (marker k-mer and syncmer sizes, variant-site span, switch probability, proposal width, mismatch band, repair-neighbourhood bounds, tie branching, representation check), resource keys and their defaults.

### Panel curation

`--panel-curation` names paths to remove from the active panel, one per line. An allele carried only by removed paths is not callable, and no stage proposes, repairs or reports it. A sequence-identical allele on a retained path stays, with only retained provenance. The removed paths still define the graph's blocks. A panel that cannot spell the locus without them refuses the run (`MODEL_INCOMPLETE`).

### Dosage modes

| mode | what it does |
|------|--------------|
| `off` | no dosage stage |
| `audit` | the dosage census, tests and repair search; writes `<prefix>.dosage_*` and never changes the reported genotype. It is not necessarily cheaper than `repair` |
| `repair` | as `audit`, and applies each eligible complete diplotype repair; the caller's outputs are kept as `<prefix>.original.*` |

A dosage-sensitive block ends in one of four states:
- **Complete diplotype repair** (`COMPLETE_DIPLOTYPE_REPAIR`): the repaired content and diplotype are complete within the declared model and search scope. "Complete" describes that scope, not biological truth. In repair mode it is applied when the run is a strict call, it is the only one at the locus, the reported haplotypes reproduce, and both repaired alleles belong to the active panel; otherwise the reason it was not applied is recorded. When several repeated-unit arrangements fit equally, the result keeps this class with `serialization_representative=true`, and its allocation is not resolved.
- **Conditional dosage content** (`CONDITIONAL_DOSAGE_CONTENT`): the corrected copy numbers are determined but the arrangement is not. The reported alleles stay the caller's. A representative arrangement is written only to `<prefix>.dosage_repaired.fa`, labelled `REPRESENTATIVE_ONLY_NOT_RESOLVED`.
- **Triggered but unrepaired:** the copy numbers are contradicted but no admissible correction was found.
- **Consistent, not evaluated, or not dosage-sensitive.**

## Outputs

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
| time or footprint limit | 124 / 125 | minimal status only | no | no |
| argument, input or output error | 2 | no locus summary | no | no |

### Reading the call

- **Strict call** (`COMPLETE`): one best pair of haplotypes. A *genotype-equivalent call* (`CONDITIONAL_GENOTYPE_EQUIVALENT`) has several best pairs that agree on every block's unordered genotype but not on phase, so it reports no FASTA and no phase.
- **No-call:** the reason is in `no_call_reason`, for example `distinct_tie`, `tie_branching_incomplete` or `neighbourhood_bound`.
- **Panel alleles:** every reported block allele is an active panel allele, and `allele*_sources` names its panel sources. A whole haplotype may still be off-panel as a mosaic.
- **FASTA headers:** they carry `phase_support` and `representation=`, which say how far the `hap1`/`hap2` order is read-supported. An applied repair adds `dosage_repair_role=RESOLVED`.
- **Dosage provenance:** sequences produced by the dosage stage are compared with the panel in `<prefix>.dosage_provenance.tsv`.

**Headline filters.** There is no single PASS/FAIL. `call_summary` is `GENOTYPE_<s>;DOSAGE_<s>;PHASE_<s>`. Each filter takes the most severe applicable status per unit, ranked NO_CALL > NOT_EVALUATED > LOW_QUALITY > PASS > NOT_APPLICABLE, and the locus takes the worst unit, with reasons as `REASON@ids`.

| filter | unit | PASS | other outcomes |
|---|---|---|---|
| GENOTYPE_FILTER | block | representation check passed, no library sensitivity, and the dosage stage, if run, left no exact diplotype unresolved | LOW_QUALITY for a representation conflict, a library-sensitive block, or an unresolved exact diplotype; NOT_EVALUATED when the representation check or the library audit could not complete |
| DOSAGE_FILTER | block | consistent; a complete repair applied, or available in audit mode; or a content correction whose content is complete (the reported alleles stay the caller's, so GENOTYPE_FILTER is LOW_QUALITY) | LOW_QUALITY for a complete repair not applied, suspect or triggered-unrepaired; NOT_EVALUATED for mode `off`, an incomplete stage, or underpowered; NOT_APPLICABLE when not dosage-sensitive |
| PHASE_FILTER | adjacent pair of heterozygous blocks | read-supported phase | LOW_QUALITY `MODEL_ONLY_PHASE` for prior-only phase; NO_CALL for unresolved phase; NOT_APPLICABLE next to a homozygous block |

- **A no-call** is NO_CALL in all three filters. A run that stopped before its blocks existed reports one `locus` unit per filter.
- **Usage policies, not filters:**
  - an exact genotype needs GENOTYPE PASS and DOSAGE PASS or NOT_APPLICABLE, and therefore a dosage stage;
  - dosage content needs DOSAGE PASS;
  - a phased genotype also needs PHASE PASS.
- **Component fields:** `filter`, `genotype_support`, `quality_filter`, the library-sensitivity fields and `dosage_stage_filter` remain in the tables. `quality_filter` also requires read support, which GENOTYPE_FILTER does not.

`scripts/plot_genotype_results.py` summarises many runs as tables, figures and a short report. It is truth-free unless truth tables are given; see its `--help`.

## Limitations

- Only paired short reads are supported.
- Reported block alleles are active panel alleles. A sequence absent from the active panel cannot be reported as an authoritative panel allele; the caller may select a representable alternative or make a no-call.
- GQ and the conditional posterior are conditional on the evaluated candidates and are not externally calibrated.
- Conditional dosage content corrects copy numbers, not the exact diplotype. A complete diplotype repair with `serialization_representative=true` does not resolve which homologue carries each repeated unit.
- Resource use grows with the active panel, the reads and the repair search. A repair neighbourhood larger than `--repair-max-neighbours` ends in a no-call (`neighbourhood_bound`).

## Example

```bash
panvar genotype --graph locus.gfa --reads1 R1.fq.gz --reads2 R2.fq.gz --out sample --threads 8
panvar genotype --graph locus.gfa --reads1 R1.fq.gz --reads2 R2.fq.gz --out sample \
  --library library.txt --dosage-mode repair --threads 8
```
