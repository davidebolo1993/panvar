# Module `genotype`

CLI: `panvar genotype`

## What it does

Genotypes one structural locus from paired short reads. It is a general graph genotyper, not specific to any locus. The locus is described by a pangenome graph whose paths are the panel haplotypes. The call is a diplotype: two whole-locus haplotypes, built from the panel's block alleles, then projected to one allele pair per block.

- Every unordered pair of panel haplotypes is scored, together with up to `--mosaic-candidates` off-panel mosaics. A mosaic switches between panel haplotypes at block boundaries.
- Read pairs are placed on each candidate, and the candidates are ranked by an exact paired-read likelihood plus a Li–Stephens recombination prior.
- The best panel pair is then repaired by exact search over single-block changes, branching exactly where moves tie. The decision is made over every candidate scored.
- The call is annotated, never changed, by a representation check, read support and phase support, and (under automatic library estimation) a library-sensitivity audit.
- An optional dosage stage asks whether the reads contradict the copy numbers the call implies. With `--dosage-mode repair` it replaces a block's genotype when a complete correction is found.

Algorithm and worked trace: [algorithms/genotype.md](../algorithms/genotype.md).

`genotype` reads a locus graph directly, for example one built by `rebuild` or a graph builder. It does not need `bubble`, `panphorte` or `call` output: without `--bubbles`, it finds the variant sites itself.

## Required inputs

- `--graph <graph.gfa>` (alias `--gfa`): the locus graph. Every path (`P` or `W` line) is a panel haplotype and must cross the whole locus; a haplotype that covers only part of it refuses the run.
- `--reads1 <r1.fq>` and `--reads2 <r2.fq>`: paired reads, FASTQ, optionally gzipped, mates in the same order. With the default `--mismatch-band 6`, reads shorter than 84 bp cannot be seeded.
- `--out <prefix>` (alias `--out-prefix`; default `genotype`). A run first removes every earlier `genotype` output under the prefix.
- The coordinate reference. Without `--reference-path`, the graph's declared reference sample (`RS:Z` header tag) is used. It must name exactly one path; otherwise the run refuses and lists the candidates.
- The read library: the insert-size distribution and the fragment-start intensity.
  - **`--library` omitted, or `--library auto`:** both are estimated from the reads. Inserts come from read pairs placed identically on every panel haplotype: at most 20,000 pairs examined, at least 500 required. The intensity comes from the marker depth fit. The other library constants come from the short-read preset. An estimate that fails validation refuses the run (`LIBRARY_REFUSED`), and there is no fallback. The validation fails on: fewer than 500 pairs; a multimodal insert distribution; more than 5% insert outliers; more than 5% non-FR pairs; or a depth fit that does not converge.
  - **`--library FILE`:** declares the library, one `key value` per line. Required keys: `insert_mean_bp`, `insert_sd_bp`, `support_sigmas`, `substitution_error`, `outlier_weight`, `fragment_start_intensity`. `marker_background` is optional. Nothing is estimated and no sensitivity audit runs, so this is the exactly reproducible route.

The simplest production command is:

```bash
panvar genotype --graph locus.gfa --reads1 R1.fq.gz --reads2 R2.fq.gz --out sample
```

It runs the default `--preset short-read`, which is the production algorithm with short-read constants. The preset sets the scientific defaults below; it sets no locus-specific, resource or dosage value. A profile overrides it key by key. The precedence, whatever the argument order, is: flag, then `--profile`, then preset, then built-in default. Every effective setting and its source is written to `<prefix>.settings.tsv`. Invalid values refuse the run before any work. `--preset long-read` is refused as not implemented.

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `--library auto\|FILE` | estimate the library from the reads, or declare it ([above](#required-inputs)) | `auto` |
| `--dosage-mode off\|audit\|repair` | the [dosage stage](#dosage-modes) | `off` |
| `--threads <N>` | worker threads; results do not depend on the thread count | `1` |
| `--mosaic-candidates <K>` | off-panel mosaics scored beside every panel pair; `0` scores panel pairs only | `128` |
| `--panel-curation <file>` | panel haplotypes to exclude, one path name per line (`#` starts a comment) | — |
| `--bubbles <file>` | a precomputed bubbles CSV instead of detecting variant sites | detected |
| `--reference-path <name>` | the coordinate reference path | graph `RS:Z` |
| `--profile <file>` | a validated `key<TAB>value` settings file that names itself with a `profile` line | — |
| `--adaptive-repair on\|off` / `--repair-rounds <N>` | the exact repair and its round cap; `off` gives the unrepaired development view | `on` / `32` |
| `-q, --quiet` | no progress messages | off |

`panvar genotype --help-advanced` lists the scientific flags, each identical to its profile key:

| flag | short-read preset |
|------|-------------------|
| `--kmer-size`, `--syncmer-s` | `31`, `15` |
| `--min-variant-bp` | `50` |
| `--switch-probability` | `0.01` |
| `--proposal-state-limit` | `256` |
| `--mismatch-band` | `6` |
| `--repair-local-pairs`, `--repair-max-neighbours` | `32`, `670` |
| `--repair-tie-branching`, `--repair-tie-max-branches` | `on`, `16` |
| `--representation-guard` | `on` |

Resource caps, caches, and time and memory limits are profile keys only. `--phase-audit-candidates` is a development-only control for the unrepaired view; adaptive repair refuses a non-zero value.

### Dosage modes

| mode | what it does |
|------|--------------|
| `off` | no dosage work; DOSAGE_FILTER is NOT_EVALUATED |
| `audit` | runs the dosage stage and writes `<prefix>.dosage_*`; the reported genotype is unchanged |
| `repair` | as `audit`, and a block with an eligible complete repair gets the repaired genotype; the caller's outputs are kept as `<prefix>.original.*` |

The stage separates three kinds of result, which are never mixed:
- **Ordinary call:** the caller's block alleles and haplotypes.
- **Complete repair** (`repair_class=COMPLETE_DIPLOTYPE_REPAIR`): the reads determine both repaired block alleles. In repair mode it replaces the block's genotype when the run is a strict call, it is the only complete repair at the locus, the reported haplotypes are reproduced, and both repaired sequences are panel alleles. Otherwise it is recorded with the reason it was not applied.
- **Conditional dosage content** (`CONDITIONAL_DOSAGE_CONTENT`): the reads determine the corrected copy numbers (`corrected_content`) but not which arrangement carries them. The block keeps the caller's alleles; a representative arrangement is written only to `<prefix>.dosage_repaired.fa`, marked `REPRESENTATIVE_ONLY_NOT_RESOLVED`, and is never a call.

## Outputs

| file | contents |
|------|----------|
| `<prefix>.locus_genotype.tsv` | the locus-level call, one row, written in every mode: run status, call class and no-call reason; the two haplotypes (IDs, SHA-256, lengths, block paths); the unordered genotype and its score and separation; panel-template switch counts; dosage summary; the three headline filters with reasons and per-status counts; links to the other files |
| `<prefix>.block_genotypes.tsv` | one row per block, context blocks included (`kind`, `reportable`): the allele pair, lengths, zygosity, read and phase support, conditional posterior and GQ, the representation fields, block quality, and the per-block headline filters |
| `<prefix>.genotypes.tsv` | the same per-block call for variant blocks only, with each allele's panel source alleles |
| `<prefix>.mosaic_haplotypes.fa` | the two whole-locus haplotypes of a strict call |
| `<prefix>.mosaic_status.tsv` | `key value` run record: status, provenance digests, model, markers, search, ranking and repair details |
| `<prefix>.settings.tsv` | every effective setting and its source (`default`, `preset short-read`, `profile <path>`, `flag --…`) |
| `<prefix>.mosaic_candidates.tsv`, `.mosaic_scores.tsv`, `.equivalence_set.tsv`, `.repair_rounds.tsv`, `.tie_*.tsv` | the search record: every candidate scored, its terms, the top set and the repair trajectory |
| `<prefix>.library_estimate.tsv`, `.library_sensitivity.tsv`, `.library_sensitivity_blocks.tsv` | under `--library auto`: the estimate and the 0.9×/1.1× intensity audit |
| `<prefix>.dosage_*`, `<prefix>.original.*` | with `--dosage-mode audit\|repair` / `repair` |

Exit status is 0 when a genotype was emitted (strict or genotype-equivalent call) and 1 for a no-call. Invalid arguments and unreadable inputs exit 2; the optional time and memory limits (profile keys `time_limit_seconds`, `footprint_limit_bytes`) exit 124 and 125.

### Reading the call

**Call class.** A *strict call* (`run_status=COMPLETE`) has one best pair of haplotypes. A *genotype-equivalent call* (`CONDITIONAL_GENOTYPE_EQUIVALENT`) has several best pairs that agree on every block's unordered genotype but not on phase: it reports the unordered genotype only, with no FASTA and no phase. A *no-call* reports its status and `no_call_reason`, for example `distinct_tie` (best candidates disagree at a block) or `tie_branching_incomplete`, and writes no genotype tables.

**Haplotype sequences.** Each block's two alleles are panel block alleles: `allele1_sources` and `allele2_sources` name the panel alleles, and the locus table counts the fewest switches between panel haplotypes needed to spell each haplotype. A haplotype can therefore be off-panel as a whole, while every block of it is a panel allele. The FASTA header's `phase_support` and `representation=` say how far its `hap1`/`hap2` order is supported. After an applied complete repair the header gains `dosage_repair_role=RESOLVED`. Sequences generated by the dosage stage are compared with the panel in `<prefix>.dosage_provenance.tsv`, which reports the nearest panel haplotype and allele by edit distance.

**Headline filters.** There is no single PASS/FAIL. Three independent verdicts summarise the call, rendered together as `call_summary=GENOTYPE_<s>;DOSAGE_<s>;PHASE_<s>`:

| filter | unit | PASS means |
|---|---|---|
| GENOTYPE_FILTER | block | the representation check passed; under `--library auto`, the genotype was stable at both audit intensities; and, if the dosage stage ran, it did not leave the exact diplotype unresolved |
| DOSAGE_FILTER | block | the dosage stage ran and found the copy numbers consistent, or a complete repair resolved them |
| PHASE_FILTER | adjacent pair of heterozygous blocks | the reads, not only the prior, support the phase |

Each takes the most severe applicable status: NO_CALL, then NOT_EVALUATED, then LOW_QUALITY, then PASS, then NOT_APPLICABLE. The locus value is the worst unit, with reasons listed as `REASON@ids`.
- **LOW_QUALITY:** a call was made, but a check flagged it, e.g. `REPRESENTATION_CONFLICT`, `LIBRARY_SENSITIVE`, or `MODEL_ONLY_PHASE` when the phase rests on the prior.
- **NOT_EVALUATED:** the check could not be made, e.g. `DOSAGE_NOT_EVALUATED:dosage_mode_off` (the default), or an incomplete library audit.
- **NO_CALL:** nothing was called.
- **NOT_APPLICABLE:** the check does not apply, e.g. a block that is not dosage-sensitive, or a boundary next to a homozygous block.

The usage policies are selections, not filters. An exact genotype is usable when GENOTYPE_FILTER is PASS and DOSAGE_FILTER is PASS or NOT_APPLICABLE; with `--dosage-mode off` it therefore needs the dosage stage first. Dosage content is usable when DOSAGE_FILTER is PASS. A phased genotype is usable when it is an exact genotype and PHASE_FILTER is PASS.

The headline filters summarise component annotations that stay in the tables: `filter` (the representation check), `genotype_support` (`READS` or `PRIOR_ONLY`), `quality_filter`, the library-sensitivity fields and, in repair mode, `dosage_stage_filter`. A component can be stricter than the headline. For example, `quality_filter` is LOW_QUALITY for `PRIOR_ONLY` read support, while GENOTYPE_FILTER does not use read support. After an applied repair the components describe the reported genotype; the caller's values are kept in `original_*` columns.

### Reporting

`scripts/plot_genotype_results.py` summarises many runs as tables, figures and a short report. It reads a directory of `panvar genotype` prefixes or a two-table package. It never changes a call. Without truth tables it is truth-free QC: call classes, headline filters, support, dosage states, runtime. Accuracy figures need explicit truth tables (`--evaluation-blocks`, `--dosage-content-evaluation`, `--phase-boundaries`). Reruns are byte-identical. Its options are documented under `--help`.

## Limitations

- **The panel bounds the answer.** Every reported block allele is a panel allele. A sequence absent from the panel cannot be called; the closest representable allele is called instead. A switch between panel haplotypes can only occur at a block boundary.
- **Allocation is unresolved.** Short reads may establish how many copies of a repeated unit a sample carries without showing which homologue carries each copy. A complete repair with an unresolved allocation is reported as such, and conditional dosage content never becomes a call.
- **Some sequence is structurally unobservable.** Where the alternatives differ only beyond the reach of a read pair, the reads cannot separate them. The call is then decided by the prior (`PRIOR_ONLY`, `MODEL_ONLY_PHASE`), or the run no-calls on a tie. GQ is conditional on the scored candidates and uncalibrated.
- **No long-read support.** Only paired short reads are modelled, and `--preset long-read` is refused. Single-end reads are not supported.
- **Validation scope.** The method and the resource caps have been studied in depth on one complex locus with simulated reads, which carry constant base qualities. The dosage stage is off by default for that reason.
- **Resources.** These figures are an example, not a guarantee: on that locus, the median 30× sample took 1–2.5 minutes, and the peak memory was up to 3.6 GiB, at 8 threads. Larger panels and more alleles per block cost more. A repair neighbourhood larger than `--repair-max-neighbours` refuses the run (`neighbourhood_bound`).

## Example

A minimal run, with the library estimated from the reads:

```bash
panvar genotype --graph locus.gfa --reads1 R1.fq.gz --reads2 R2.fq.gz --out sample --threads 8
```

The same run with a declared library and dosage repair:

```bash
cat > library.txt <<'EOF'
insert_mean_bp 400
insert_sd_bp 50
support_sigmas 4
substitution_error 0.002
outlier_weight 0.001
fragment_start_intensity 0.05
EOF
panvar genotype --graph locus.gfa --reads1 R1.fq.gz --reads2 R2.fq.gz --out sample \
  --library library.txt --dosage-mode repair --threads 8
cut -f1-5 sample.block_genotypes.tsv
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i} NR==2{print $c["call_class"], $c["call_summary"]}' sample.locus_genotype.tsv
```
