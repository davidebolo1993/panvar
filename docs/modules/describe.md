# Describe Module (Module 4)

Date: 2026-05-06

CLI entrypoint:

- `panvar describe`

## What it does

`describe` consumes the region-level VCF from `panvar call` and writes per-bubble tables that summarize haplotypes with a consistent feature schema.

For each bubble, every haplotype/path gets a feature vector with:

- presence/count of event classes (`INS`, `DEL`, `INV`, `DUP` evidence)
- size-binned event counts
- optional gene-region overlap counts from GTF
- a `feature_signature` label (same vector -> same signature)

## Required inputs

- `--vcf-in <call.region.vcf>`

## Optional inputs

- `--gtf <path>` for gene/exon/CDS/UTR/intron overlaps
- `--gene-match <expr>` repeatable, case-insensitive regex filter on gene names
  - if omitted: intersections are counted against all genes in GTF
  - if provided: intersections are counted only for matching genes
- `--size-bins <csv>` thresholds for length bins (default `100,1000`)

## Outputs

Inside `--out-dir` (default `describe_out`):

(only bubbles that have at least one VCF event record are emitted)

- `describe.index.tsv`
- one event table per bubble:
  - `bubble_<id>.events.tsv`
- one haplotype feature table per bubble:
  - `bubble_<id>.haplotype_features.tsv`

`bubble_<id>.events.tsv` includes per-event metadata, carriers, and overlap flags.

`bubble_<id>.haplotype_features.tsv` includes one row per haplotype with:

- `has_ins`, `has_del`, `has_inv`, `has_dup`
- `n_ins`, `n_del`, `n_inv`, `n_dup`
- `n_events_gene`, `n_events_exon`, `n_events_cds`, `n_events_utr`, `n_events_intron`
- size-bin counts for `ins/del/inv/dup`
- `feature_signature`, `signature_group_size`

## Algorithm overview

For each VCF record:

1. Read event info (`BUBBLE_ID`, `EVENT`, `SVTYPE`, `SVLEN`, `INS_SUBTYPE`, etc.).
2. Determine carrier haplotypes from per-sample `GT`.
3. Optionally annotate event interval overlap with GTF-derived classes (`gene`, `exon`, `CDS`, `UTR`, derived `intron`).
4. Aggregate per-haplotype features within each bubble.
5. Assign a stable signature ID to identical feature vectors.

## Key options

- `--vcf-in <path>`
- `--out-dir <dir>`
- `--gtf <path>`
- `--gene-match <expr>` (repeatable)
- `--size-bins <csv>`
- `--quiet`

## Example (C4)

```bash
./build/panvar describe \
  --vcf-in /tmp/panvar_call_bounded.region.vcf \
  --out-dir /tmp/panvar_describe_c4 \
  --gtf tests/real_data/gencode.v49.annotation.gtf.gz \
  --gene-match C4 \
  --size-bins 100,1000
```
