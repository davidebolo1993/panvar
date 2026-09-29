# Reporting `panvar genotype` results

`scripts/plot_genotype_results.py` turns many `panvar genotype` outputs into summary tables, figures and a short
report. It only reads its inputs. It never genotypes, never changes a call and never infers truth.

**Dependencies:** Python ≥ 3.10 with pandas and matplotlib (numpy comes with pandas). It was tested with Python 3.14,
pandas 3.0, matplotlib 3.11 and numpy 2.5. Nothing else is needed.

## Input

`--input` is either of two forms, which normalise to the same two tables (`data/runs.tsv`, `data/blocks.tsv`):

- **A directory of output prefixes.** Every `<prefix>.locus_genotype.tsv` below it is a run (`<prefix>.original.*`
  dosage-repair originals are not). `run_id` is the prefix path relative to the directory. The tool reads the locus
  table, `mosaic_status.tsv`, `block_genotypes.tsv`, `dosage_blocks.tsv`, `dosage_status.tsv`, `settings.tsv` and
  `library_sensitivity.tsv`.
- **A package:** a directory with `runs.tsv` and `blocks.tsv` in that normalised schema. Runs are keyed by `run_id`, or by
  `dataset` and `job`, which become `dataset/job`.

**Block order** is read from the outputs: the model order of each block table, or `block_index`. Runs whose block
sequences differ are refused, because one report covers one locus model.

**Outputs from before `headline_filters_v2`:**
- Their `dosage_filter` column is the dosage-stage component, since renamed `dosage_stage_filter`. It is read as that,
  never as the headline filter.
- They have no headline figures.
- Mixing them with newer outputs is refused.

## Two modes

**Truth-free QC** needs only the outputs and draws:
- call class and no-call reasons;
- the three headline filters and their reasons;
- conditional GQ;
- genotype and phase support;
- dosage states;
- repair rounds;
- runtime and memory;
- library source and sensitivity;
- settings provenance.

Nothing in this mode is labelled accurate or erroneous.

**Evaluation** is enabled only by explicit truth tables. Without them every accuracy figure is refused, and the reason
is recorded.

| option | content |
|---|---|
| `--evaluation-blocks` | The evaluator's per-block fields (`unordered_distance`, `availability_floor`, `excess`, and for normalised excess `truth1_bp`, `truth2_bp`, `called1_bp`, `called2_bp`). Either one table with a run key and `block_id`, or a manifest (run key, `path`) of per-run evaluator block tables. |
| `--dosage-content-evaluation` | `content_correct` (and `content_truth_source`) per block. |
| `--representative-diagnostics` | Conditional representatives' sequence distances. Non-authoritative. |
| `--switch-metadata` | Mosaic plans: `subject`, `switch_class`, and `switch_before_block_hap1` and `switch_before_block_hap2`. Each lists the 0-based block indices i at which that homologue switches source template, between blocks i−1 and i. `source_a` and `source_b` are optional. |
| `--phase-boundaries`, `--phase-calls`, `--phase-projected` | The phase evaluation's boundary table, per-call table and secondary sequence-projected table. |

**Run metadata** (`--run-metadata`) gives each run an `experiment` and a `coverage`, and optionally a `subject`, a
`modality` and an `execution_group`. Without it every run is in one group.

**Modality:**
- When run metadata supplies a modality, it leads every aggregate's grouping.
- With several modalities, every figure is drawn once per modality (`<figure>.<modality>.png`), and runtime is never
  pooled across them.
- With one modality, every figure is byte-identical to a run without the column.
- A missing modality is never assumed, and an empty one is refused.

## Definitions

- **Raw distance, floor and excess** are the evaluator's `unordered_distance`, `availability_floor` and `excess`, used
  as they are. Excess is the leave-one-out accuracy metric. Floor attained means excess 0.
- **Normalised excess** is `excess / max(truth1_bp + truth2_bp, called1_bp + called2_bp)`, from the evaluator's recorded
  allele lengths: the fewest alignment columns any alignment of the pair can have. The figures are refused when those
  lengths are absent.
- **A truth-undefined block**, where the evaluator records no truth length, is excluded from every accuracy denominator
  and counted.
- **EXACT_GENOTYPE_USABLE** is exactly `GENOTYPE_FILTER == PASS and DOSAGE_FILTER in {PASS, NOT_APPLICABLE}`.
- **DOSAGE_CONTENT_USABLE** is `DOSAGE_FILTER == PASS`, judged on dosage content.
- **PHASE_FILTER** is judged only against phase truth. It is never folded into exact-genotype filtering.
- **Catastrophic** means `excess >= --catastrophic-bp`, 1,000 by default.
- **Intervals** are 95% Wilson intervals.
- **GQ** is conditional on the scored candidate set. It is shown only as a diagnostic, and no threshold is recommended.
- **Runtime** is pooled only within a compatibility group: settings digest, threads, execution group and modality. The
  outputs record no host or binary identity, so supply it as `execution_group`.

## Output

| path | content |
|---|---|
| `data/*.tsv` | every aggregate behind a figure, including point-level tables; numerators beside denominators |
| `figures/*.png`, `figures/*.pdf` | the figures; the PNG and PDF metadata carry no timestamps or software versions |
| `report.md` | inputs, denominators and exclusions, definitions, a linked figure index with captions, and skipped figures with reasons |
| `manifest.tsv` | the command, tool and library versions, the script commit and file hashes, and every input and output file with its SHA-256 |

**Reproducibility:** a rerun on the same inputs gives byte-identical `data/`, `figures/` and `report.md`. Only the
manifest's command line differs, when `--out-dir` differs.

**Safety:** an `--out-dir` that overlaps an input location, or that is not empty, is refused. A malformed schema,
mixed block models, metadata that does not cover every run, and evaluation rows for unknown runs or blocks are
refused with exit status 2.

## Examples

**QC only** (no truth):

    python3 scripts/plot_genotype_results.py --input results/genotype_runs --out-dir report_qc \
      --run-metadata runs_metadata.tsv --title "Cohort QC"

**Truth evaluation** (the C4 report; the inputs are joined from sealed tables by a separate adapter):

    python3 scripts/plot_genotype_results.py --input c4_inputs/package --out-dir c4_report \
      --run-metadata c4_inputs/run_metadata.tsv --evaluation-blocks c4_inputs/evaluation_manifest.tsv \
      --switch-metadata c4_inputs/switch_metadata.tsv \
      --dosage-content-evaluation dosage_content_blocks.tsv \
      --phase-boundaries boundaries_primary.tsv --phase-calls calls_primary.tsv \
      --phase-projected loo_sequence_projected_boundaries.tsv \
      --panel-groups "LZO=lzo,fresh_lzo;mosaic LZO=mosaic;LOO=loo,fresh_loo"

`--panel-groups` sets the rows of the ordered-block figures and places every experiment exactly once.

## Tests

`tests/genotype_report_test.py`, registered as `panvar_genotype_report`, runs on small synthetic tables. It checks:
- counts and denominators computed by hand;
- the normalisation formula;
- truth-undefined exclusion;
- every refusal;
- the modality contract;
- byte-identical reruns;
- unmodified inputs;
- directory-mode loading of outputs from before the headline filters.
