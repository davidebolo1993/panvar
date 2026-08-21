# Module `benchmark`

CLI: `panvar benchmark`

## What it does

Answers two different questions about the caller's own output, and keeps them apart.

What was there to be found, and what became of it. Independently of any call, each haplotype's walk over a bubble is compared with the reference's, and every divergent block becomes a truth event classified by size and by whether a record covers it:

| class | definition |
|-------|------------|
| `called` | at least `--min-sv-bp` and some emitted record shares a node with the block |
| `missed` | at least `--min-sv-bp` and no emitted record covers it |
| `below_threshold` | under `--min-sv-bp` — variation the caller never set out to emit |

How much of the sequence comes back. Four reconstructions over the same bubbles against the same truth, differing only in what each is allowed to use:

| level | reconstruction built from | answers |
|-------|---------------------------|---------|
| `graph` | reference walk plus the haplotype's own steps at every block sharing a node with any call | optimistic node-discovery ceiling: how far the union of called nodes reaches |
| `called` | the same, restricted to blocks a specific record is attributed to and that reach `--min-sv-bp` | what the retained records would reach if each reproduced its block exactly |
| `carrier` | `called`, plus the requirement that this haplotype's genotype names a record overlapping the block | the same ceiling once genotypes are applied |
| `region_vcf` | the reference plus only the edits this haplotype's genotype names | what a consumer reconstructing this sample from the VCF actually gets |

The first three implant the haplotype's own true sequence and are therefore ceilings, not VCF reconstructions. Only `region_vcf` reconstructs what the records themselves say. No level uses reads: this measures the fidelity of the compact representation, not genotyping accuracy.

The region-VCF residual is then partitioned into five terms that sum to it exactly, so the shortfall can be attributed rather than just measured:

| term | bases lost to |
|------|---------------|
| `out_of_scope` | truth events below the size threshold, which the caller was never asked to emit |
| `discovery_or_attribution` | eligible truth events no retained record covers |
| `carrier_missed` | covered by a record, but this haplotype is not genotyped as carrying it |
| `representation` | right record, right carrier, but the encoding does not reproduce the sequence |
| `false_positive_damage` | an edit applied where the haplotype has no eligible truth event |

A do-nothing baseline, the plain reference with no edits, is computed alongside and is the denominator of the headline figures:

```
gap_closed          = (baseline_delta - region_vcf_delta) / (baseline_delta - graph_delta)
variation_recovered = (baseline_delta - <level>_delta)  /  baseline_delta
```

`gap_closed` measures against the achievable bound, `variation_recovered` against all the variation there was. The baseline also makes the metric self-checking: a haplotype genotyped as carrying nothing must reconstruct byte-identically to it.

Algorithm and worked trace: [algorithms/benchmark.md](../algorithms/benchmark.md).

## Required inputs

- `-i, --gfa <graph.gfa>` — the same graph used by `call`; paths, links, overlaps and names are validated before spelling.
- one of `-b, --bubble-prefix-in <prefix>` or `-c, --bubbles-csv-in <path>`; bubble ids and graph nodes are cross-checked.
- `--variant-nodes <path>` — `call`'s matching `<prefix>.variant_nodes.tsv`. Its schema, variant ids, bubble ids and node membership are validated because these node sets define truth-event attribution.
- `-r, --reference-path <name>` — the diff baseline, resolved through the shared resolver (exact name, else a unique case-insensitive substring; ambiguity is an error, never settled by file order).
- `-o, --out-prefix <prefix>`.

Optional: `--vcf <path>` accepts either VCF from `call`; `vcf_mode` reports which contract was detected. Region-VCF ids must match `variant_nodes.tsv` completely and agree on `BUBBLE_ID`. An allele VCF instead has one `bubbleN_ALLELES` record per bubble, so per-call carrier and loss-attribution terms are not applicable; allele mode measures explicit-allele serialization.

`<prefix>.region.vcf` scores the merged, interpreted output; `<prefix>.alleles.vcf` (from `call --allele-vcf`) scores the explicit per-bubble alleles. Sample columns join graph paths by exact name and must contain haploid `GT`; record ids, `BUBBLE_ID`, row widths and duplicate names are validated. A partial sample/path join is accepted but reported because its quality values describe only the joined subset.

All outputs are staged and committed only on success, and no output may name an input.

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `--vcf <path>` | add the `region_vcf` level and the baseline it is measured against; accepts the region VCF or the allele VCF | off |
| `--called-bubbles-only` | score only bubbles carrying `≥ 1` call. Ascertainment-biased by construction: the bubble the caller said nothing about is exactly where a miss lives | off — every reference-traversed bubble is scored |
| `--min-sv-bp <N>` | event-size threshold for the ledger and the `called` reconstruction; set it to what `call` ran with | `50` |
| `--dup-model <cn\|cnbp>` | lay down `CN` copies of `RU_LEN` in place of `REF_CN`, or apply the per-sample `CNBP` delta. Both tile an inferred reference span, so a DUP is reconstructed at the right length out of approximately right sequence and is counted `heuristic` | `cnbp` |
| `--no-truth-events` | skip the per-event ledger table | off |
| `--threads <N>` | worker threads over haplotypes | `0` (auto) |
| `-q, --quiet` | disable progress logs | — |

## Outputs

| file | contents |
|------|----------|
| `<prefix>.truth_events.tsv` | one row per (haplotype, bubble, truth event): `ref_bp`, `hap_bp`, `size_bp`, `class`, and the `variant_id` / `svtype` it maps to. The ledger everything else aggregates |
| `<prefix>.qv.tsv` | per (haplotype, bubble): `svtypes`, `is_carrier`, `δ`, `aln_len`, `qv`, the run split `resid_run_lt_bp` / `resid_run_ge_bp`, the ledger counts `truth_events` / `truth_called` / `truth_missed` / `truth_below` / `truth_missed_bp`, `called_delta` / `called_qv`, and `carrier_delta` / `carrier_qv`. With `--vcf`, also `gt_delta` / `gt_aln_len` / `gt_qv`, the baseline `ref_delta`, and `gt_called_carrier` / `gt_true_carrier` |
| `<prefix>.qv_by_haplotype.tsv` | the same per haplotype: `identity`, the run split, the ledger, `called_sum_delta` / `called_qv` / `called_quintile`, `carrier_sum_delta` / `carrier_qv`, plus `qv`, cosigt `band`, `qv_max`, `qv_ratio`, `quintile`. With `--vcf`, the `gt_*` mirror plus `ref_sum_delta`, `ref_qv` and `gap_closed` (literal `NA` where undefined) |
| `<prefix>.qv_summary.tsv` | `truth_event` / `truth_bp` (the ledger, `ALL` and per svtype), `quintile` / `haplotype` / `called_quintile` distributions, `residual_run`, `excluded` (bubbles with no reference walk, haplotype-bubble pairs not traversed, pairs decomposed coarsely, haplotypes with no VCF column, VCF samples with no path), `called_recon`. With `--vcf`: `gt_quintile` / `gt_haplotype`, `gt_gap` (pooled deltas, `gap_closed_pooled`, `gap_closed_mean`, `gap_closed_undefined`, `worse_than_baseline`), `variation_recovered` for every level, `carrier_quintile` / `carrier_recon` and the `loss_bp` partition, `gt_carrier` (TP/FP/FN/TN with precision and recall) and `gt_records` (`applied` / `unplaceable` / `clamped` / `unhandled` / `ref_mismatch` / `heuristic`) |

Nothing is dropped silently: every haplotype-bubble pair excluded from a denominator is counted in the `excluded` scope, so no rate is read as though it covered everything.

Results can be plotted via `scripts/plot_benchmark.R`.

## Limitations

- The first three levels implant the haplotype's own true block, so they bound what the graph and the retained records could achieve rather than measuring what the VCF says. `region_vcf` is not bounded by them, since it applies every record the haplotype carries rather than only attributed eligible blocks, which is why the last two partition terms are signed.
- `called` means a record shares a node with the block, not that it covers or reproduces it. It is an upper bound on discovery; `missed` is the firm negative.
- Duplication reconstruction is heuristic: the inferred reference span is tiled, so the length is right and the sequence only approximately.
- Allele-VCF mode has one record per bubble while the node sidecar has one per call, so they share no identifier. The `carrier` level and the per-call loss terms are reported as not applicable in that mode rather than as zero.
- There is no class for an eligible event removed by a frequency, support or resource policy: it is indistinguishable from one never discovered until the caller emits a decision audit.
- False-positive damage is attributed per haplotype and bubble, so a second spurious record in a bubble that also holds a valid truth event is not isolated.
- `--min-sv-bp` reclassifies the ledger and gates the `called` reconstruction but does not re-run discovery, so lowering it here finds nothing new. Use the threshold sweep, which re-runs both.

## Example

See the [walkthrough](../walkthrough.md) for this module in a full end-to-end run.
