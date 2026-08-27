# Module `benchmark`

CLI: `panvar benchmark`

## What it does

`benchmark` compares every haplotype walk with the reference inside each bubble. It answers two questions:

1. Was the haplotype variation large enough to call, and did any call find it?
2. How much of the true haplotype sequence can be reconstructed from the calls?

Truth differences are classified before the VCF is inspected:

| class | meaning |
|-------|---------|
| `called` | size is at least `--min-sv-bp` and an emitted record shares a node with the difference |
| `missed` | size is at least `--min-sv-bp`, but no emitted record covers it |
| `below_threshold` | smaller than `--min-sv-bp`; outside the requested call set |

The module then builds progressively more realistic reconstructions:

| level | what is reconstructed | interpretation |
|-------|-----------------------|----------------|
| `baseline` | unchanged reference | variation present before applying calls |
| `graph` | paste the true haplotype sequence wherever any call touches the same graph neighbourhood | loose discovery ceiling |
| `called` | paste truth only for differences attributed to a retained record | call-attribution ceiling |
| `carrier` | as `called`, but only when this haplotype's `GT` carries the record | genotyping ceiling |
| `region_vcf` | apply the VCF record's own edit to the reference | sequence a VCF consumer receives |

The `graph`, `called` and `carrier` levels paste the known truth, so they are diagnostic ceilings. Only `region_vcf` tests whether the VCF itself represents the haplotype correctly.

Two headline summaries are available:

```text
variation_recovered = (baseline_delta - level_delta) / baseline_delta
gap_closed           = (baseline_delta - region_vcf_delta)
                       / (baseline_delta - graph_delta)
```

`variation_recovered` uses all true variation as its denominator. `gap_closed` compares the VCF with the graph-level ceiling.

Algorithm and worked trace: [algorithms/benchmark.md](../algorithms/benchmark.md).

## Required inputs

- `-i, --gfa <graph.gfa>` — the graph used by `call`.
- `-b, --bubble-prefix-in <prefix>` or `-c, --bubbles-csv-in <path>` — matching bubbles.
- `--variant-nodes <path>` — `call`'s matching `variant_nodes.tsv`.
- `-r, --reference-path <name>` — reference path.
- `-o, --out-prefix <prefix>`.

Add `--vcf <path>` to evaluate genotypes and reconstruct sequence from either the region VCF or the per-bubble allele VCF. Sample columns are matched to graph paths by exact name. A partial match is allowed but reported.

## Key options

| flag | purpose | default |
|------|---------|---------|
| `--vcf <path>` | enable VCF-based reconstruction and carrier metrics | off |
| `--called-bubbles-only` | score only bubbles containing a call; useful diagnostically, but excludes completely missed bubbles | off |
| `--min-sv-bp <N>` | minimum truth-event size; normally match the value used by `call` | `50` |
| `--dup-model <cn\|cnbp>` | reconstruct duplications from copy count or copy-number base-pair change | `cnbp` |
| `--no-truth-events` | do not write the event ledger | off |
| `--threads <N>` | worker threads over haplotypes | `0` (auto) |

## Outputs

| file | contents |
|------|----------|
| `<prefix>.truth_events.tsv` | one row per haplotype, bubble and truth difference |
| `<prefix>.qv.tsv` | reconstruction metrics per haplotype and bubble |
| `<prefix>.qv_by_haplotype.tsv` | the same metrics summed over each haplotype |
| `<prefix>.qv_summary.tsv` | pooled event, reconstruction and loss summaries |

### Reading the reconstruction fields

The prefixes encode the reconstruction level. The unprefixed fields are the graph ceiling for historical compatibility.

| fields | level |
|--------|-------|
| `ref_delta`, `ref_qv` | unchanged-reference baseline |
| `delta`, `aln_len`, `qv`, `identity` | graph ceiling |
| `called_delta`, `called_qv` | call-attribution ceiling |
| `carrier_delta`, `carrier_qv` | genotype-aware ceiling |
| `gt_delta`, `gt_aln_len`, `gt_qv`, `gt_identity` | reconstruction from the VCF itself |
| `gap_closed` | fraction of the baseline-to-graph gap recovered by the VCF |

`delta` is edit distance to the true haplotype; lower is better. `identity` and `qv` are alternative presentations of the same alignment result. Fields requiring a VCF are written as `.` when `--vcf` is absent.

The event ledger's main fields are `ref_bp`, `hap_bp`, `size_bp`, `class`, `variant_id` and `svtype`. The summary additionally reports:

- event counts and missed bases;
- `variation_recovered` at each level;
- carrier TP/FP/FN/TN;
- the five-part residual loss: `out_of_scope`, `discovery_or_attribution`, `carrier_missed`, `representation` and `false_positive_damage`;
- VCF records that were applied, heuristic, unplaceable, clamped, unhandled or reference-mismatched;
- observations excluded from a denominator.

### Plotting

```bash
Rscript scripts/plot_benchmark.R --table <prefix>.qv_by_haplotype.tsv --locus <name> --out <prefix>
```

## Limitations

- `called` means node overlap with an emitted record, not exact sequence representation. Use `region_vcf` for the latter.
- DUP reconstruction is heuristic because `CN` and `CNBP` provide copy amount or length change, not necessarily the inserted sequence.
- Allele-VCF records are per bubble rather than per merged call, so per-call carrier and loss-attribution terms are not applicable in allele mode.
- Changing `--min-sv-bp` only reclassifies an existing run. To measure recall at a lower calling threshold, rerun `call` and `benchmark` at that threshold.

## Example

See the [walkthrough](../walkthrough.md) for a complete run.
