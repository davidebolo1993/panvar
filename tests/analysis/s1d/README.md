# S1-D: `output_schema sequence_first_v2` — panel references, conditional FASTA, dosage audit only

**What it is:** a new, explicitly versioned experimental contract. `legacy` (the default) and `sequence_first_v1` are unchanged: v1 output is byte-identical to S1-C. The one exception is the unknown-schema error message, which now lists `sequence_first_v2`.

**Branch:** `s1-d-panel-references`, from S1-C `c744035`.

| commit | content |
|---|---|
| `cf477ba` | implementation |
| `c947bfc` | amendment 2: the repair-refusal ordering fix |
| `ac8b794` | validation tools and frozen manifest |
| docs | this README and the module page |

**What did not change:** the model, priors, search, caps, the persisted prepared-locus format, and every historical result.

## 1. Panel membership decides the FASTA

**The rule.** Each predicted haplotype, completed or provisional, is classified against the **active** panel: the model's templates after curation. The test is **exact sequence identity in the locus frame**.

**How it is computed:**
- Template sequences are spelled in memory from the existing model, so no panel FASTA is written or rebuilt per sample.
- Lookup is by SHA-256, and every hit is verified byte for byte.
- No truth, label, copying history, switch count or similarity is used.
- A copying history that spells an existing panel sequence is panel-represented.
- Sequence-equivalent aliases all match. They are listed in byte order, and the first is the reporting alias. Choosing it changes nothing in the model.

| outcome | FASTA | `fasta_reason` |
|---|---|---|
| both haplotypes panel-represented | **not written**; extraction references reported | `PANEL_REPRESENTED_PAIR` |
| at least one off-panel | one FASTA with **both** haplotypes; each record carries `origin=PANEL panel_path=…` or `origin=OFF_PANEL` | `OFF_PANEL_PAIR` |
| no prediction | none | `NO_PREDICTION` |

**Provisional predictions** follow the same rule and keep `provisional=yes`.

**Intentional absence is never a failure.** The status says `sequence_fasta not_emitted` and `sequence_fasta_reason`, with `sequence_publication complete`.

## 2. Extraction references

**`<prefix>.sequence_haplotypes.tsv`**, one row per haplotype, with these columns:
- `sha256`, `length_bp`;
- `origin`;
- `panel_paths` (every equivalent template) and `path_orientations` (source path relative to the locus chain);
- `extraction`;
- `graph_sha256` (the GFA) and `model_fingerprint`;
- `block_path` (the selected labelling);
- `representation_intervals`;
- `provisional`.

**The locus record** repeats the origins and paths, and adds `fasta_emitted`, `fasta_reason`, the membership rule, and the graph and model identity. The prepared-index identity stays in the status (`prepared_locus_file_sha256`, loaded route only), so fresh and loaded outputs are identical.

**The extraction reference** is `gfa_path=<template>;interval=0-<length>;orientation=<forward|reverse_complement>;whole_locus_path`. Each template path spans the locus graph whole: concatenate its oriented segments, then reverse-complement when the orientation says so.

```
python3 tests/analysis/s1d/extract_panel_haplotype.py <graph.gfa> "gfa_path=…;interval=0-…;orientation=forward" [--fasta out.fa]
python3 tests/analysis/s1d/extract_panel_haplotype.py --check <graph.gfa> <prefix>.sequence_haplotypes.tsv   # every panel path of every PANEL row
```

This was tested on every fixture, and on the retained C4 graph for both HG00096 haplotypes.

**Off-panel haplotypes** are reported as spelled candidate sequences, not assemblies: `extraction` = `not_extractable:…`. They keep their block-path provenance and representation intervals.

## 3. Dosage under v2: diagnostic only

| `--dosage-mode` | behaviour |
|---|---|
| not given | **audit**; `dosage_mode_source contract_default_sequence_first_v2` |
| `off` | kept; the legacy rows record `DOSAGE_NOT_EVALUATED:dosage_mode_off` |
| `repair` | **refused**, at argument resolution and, for API callers, before the run touches its prefix. The message directs historical reproduction to the legacy contract |

**What the audit never does:** synthesize a sequence, choose a genotype, change the prediction or its hashes, or force labels onto ambiguous blocks.

**When it is not run** (a legacy no-call), `dosage_status.tsv` says NOT_RUN explicitly. The only difference audit makes to legacy outputs is the dosage filter pair of the legacy block table, which records the mode.

Legacy repair behaviour is unchanged.

## 4. Lifecycle

`.sequence_haplotypes.tsv` and its staging name are run-owned, like the other reports:
- removed at run start and by the watchdog;
- published before the locus record, which is the commit record.

An off-panel run followed by a panel-only run at the same prefix removes the earlier FASTA. A later failed run removes everything it owns.

## 5. Limitations

- **Whole-path frame.** Extraction assumes each template path spans the locus graph whole, which is how locus graphs are built here. It was verified on every fixture template matched and on the C4 graph. A graph whose paths extend beyond the locus frame would need path offsets, which the persisted format does not carry.
- **Reverse paths.** Panel membership is in the locus frame. A template reverse to the chain is matched in chain orientation and reported `reverse_complement`.
- **Unexplored labellings.** As in S1-B/S1-C, block statuses rest on explored labellings. Membership itself is sequence-exact and independent of exploration.
- **Dosage scope.** The audit is not run on runs whose legacy classification is a no-call. That includes label-ambiguity predictions; the audit says so rather than labelling ambiguous blocks.
