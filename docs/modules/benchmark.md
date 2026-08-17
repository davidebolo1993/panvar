# Module `benchmark`

CLI: `panvar benchmark`

## What it does

Answers two different questions about the caller's own output, and keeps them apart.

**What was there to be found, and what became of it** — the truth event ledger. Independently of any call, the haplotype's walk over a bubble is node-aligned to the reference walk, and each maximal divergent block between shared anchors is one truth event, sized `max(ref_bp, hap_bp)` from the walks. Every event is then classified:

| class | definition |
|-------|------------|
| `called` | size `≥ --min-sv-bp` and the block contains a node of a specific emitted record |
| `missed` | size `≥ --min-sv-bp` and no emitted record covers it |
| `below_threshold` | size `< --min-sv-bp` — variation the caller never set out to emit |

Event size is a property of the two walks. It does not move with the alignment, with the calls, or with anything else the run does, which is what makes it usable as truth. `called` is a statement about **discovery** — a record covering this block exists — not about whether this haplotype was genotyped as carrying it. That second question is the carrier table below.

**How much of the sequence comes back** — three reconstructions over the same bubbles, against the same truth, on the same QV scale, differing only in what each is allowed to use.

| level | reconstruction built from | answers |
|-------|---------------------------|---------|
| `graph` | reference walk + the haplotype's own steps at every block sharing a node with **any** call at that bubble | optimistic upper bound: can the graph hold this haplotype, and did the caller flag the divergent blocks |
| `called` | the same substitution restricted to blocks that map to a **specific** record and reach `--min-sv-bp` | the reference with exactly the retained calls implanted, and nothing else |
| `genotype` | reference sequence + only the edits this haplotype's `GT` names, from the VCF alone (`--vcf`) | what a consumer reconstructing this sample from the VCF actually gets |

`graph` is deliberately optimistic and is not a genotyping score: sharing a node with some call is not the same as matching one, so a called block can authorise copying a neighbouring uncalled allele, and no genotype is read at all. `called` is the honest version of the same idea. The gap between `called` and `genotype` is what the VCF's *representation* costs, separated from what discovery costs.

A fourth **baseline** — plain reference, no edits — is always computed alongside, and is the denominator of the headline:

```
gap_closed          = (baseline_delta - genotype_delta) / (baseline_delta - graph_delta)
variation_recovered = (baseline_delta - <level>_delta)  /  baseline_delta
```

`gap_closed` measures against the *achievable* bound; `variation_recovered` against all of the variation there was. Reported for all three levels, the pair says which of representation or graph is the limit. **When `baseline_delta == graph_delta` there is no distance to close and `gap_closed` is `NA`, not `1.0`** — the two are equal precisely when a bubble was missed entirely, so reporting `1.0` scored a total miss as a perfect result.

The baseline also makes the metric self-checking: a haplotype genotyped as carrying nothing must score byte-identically to it, so a coordinate or orientation error in the projection shows up immediately instead of hiding inside a plausible QV.

### Carrier calls

Per (haplotype, bubble): does the VCF genotype this haplotype as carrying an alt, against whether its walk really diverges by at least `--min-sv-bp`. Truth comes from the ledger, so it is a property of the walks. Reported for `ALL` only — a bubble-level judgement has no single SV type, and fanning it out over every type present in the bubble reported one observation several times. The per-type partition is the ledger, where each event maps to at most one record.

`FN` here means the record exists but this haplotype's `GT` is `0` at every record of the bubble: a **genotyping** miss, distinct from a `missed` event, which is a **discovery** miss.

### On the residual split

`resid_run_lt_bp` / `resid_run_ge_bp` divide the residual by contiguous alignment-run length. **This is a property of the edit path edlib returns, not of any event.** A clean 60 bp deletion of ordinary non-repetitive sequence comes back as fourteen runs of 1–10 bases, because co-optimal paths distribute the gap over chance matches; every real locus therefore reported `over_threshold_bp = 0`. The columns are kept as a description of the residual's shape. The called/missed question belongs to the ledger.

### On `--min-sv-bp`

It sets the event-size threshold for the ledger and for the `called` reconstruction. It does **not** re-run discovery, so lowering it here does not find anything new — it reclassifies. A real threshold experiment runs `call` at each threshold and `benchmark` on that run's output at the same threshold, over a fixed bubble set.

Algorithm and worked trace: [algorithms/benchmark.md](../algorithms/benchmark.md).

## Required inputs

- `-i, --gfa <graph.gfa>` — the same passed graph the calls were made on. Validated: every path step must name a node the graph has, consecutive steps must be joined by a link, path names must be unique and link overlaps must be zero, because every walk is spelled by concatenating whole segments.
- one of `-b, --bubble-prefix-in <prefix>` or `-c, --bubbles-csv-in <path>`. Bubble ids must be unique and every node named must exist in the graph.
- `--variant-nodes <path>` — `call`'s `<prefix>.variant_nodes.tsv`. If none of its bubble ids appear in the CSV the two are from different runs and the run is refused.
- `-r, --reference-path <name>` — the diff baseline, resolved through the shared resolver (exact name, else a unique case-insensitive substring; ambiguity is an error, never settled by file order).
- `-o, --out-prefix <prefix>`.

Optional: `--vcf <path>` — either VCF `call` produces. `<prefix>.region.vcf` scores the merged, interpreted output; `<prefix>.alleles.vcf` (from `call --allele-vcf`) scores the lossless one, which is a serialization ceiling rather than a call-sensitivity score. Sample columns are joined to graph paths by exact name; each column is one haplotype, so genotypes are haploid and nothing has to be phased. Duplicate sample columns are refused. A partial join is accepted and **prominently reported**, because a QV over a subset is not the QV of the run.

All outputs are staged and committed only on success, and no output may name an input.

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `--vcf <path>` | add the genotype level and the baseline it is measured against | off |
| `--called-bubbles-only` | score only bubbles carrying `≥ 1` call. Ascertainment-biased by construction: the bubble the caller said nothing about is exactly where a miss lives | off — every reference-traversed bubble is scored |
| `--min-sv-bp <N>` | event-size threshold for the ledger and the `called` reconstruction; set it to what `call` ran with | `50` |
| `--dup-model <cn\|cnbp>` | lay down `CN` copies of `RU_LEN` in place of `REF_CN`, or apply the per-sample `CNBP` delta. Both tile an inferred reference span, so a DUP is reconstructed at the right *length* out of approximately right sequence and is counted `heuristic` | `cnbp` |
| `--no-truth-events` | skip the per-event ledger table | off |
| `--threads <N>` | worker threads over haplotypes | `0` (auto) |
| `-q, --quiet` | disable progress logs | — |

## Outputs

| file | contents |
|------|----------|
| `<prefix>.truth_events.tsv` | one row per (haplotype, bubble, truth event): `ref_bp`, `hap_bp`, `size_bp`, `class`, and the `variant_id` / `svtype` it maps to. The ledger everything else aggregates |
| `<prefix>.qv.tsv` | per (haplotype, bubble): `svtypes`, `is_carrier`, `δ`, `aln_len`, `qv`, the run split `resid_run_lt_bp` / `resid_run_ge_bp`, the ledger counts `truth_events` / `truth_called` / `truth_missed` / `truth_below` / `truth_missed_bp`, and `called_delta` / `called_qv`. With `--vcf`, also `gt_delta` / `gt_aln_len` / `gt_qv`, the baseline `ref_delta`, and `gt_called_carrier` / `gt_true_carrier` |
| `<prefix>.qv_by_haplotype.tsv` | the same per haplotype: `identity`, the run split, the ledger, `called_sum_delta` / `called_qv` / `called_quintile`, plus `qv`, cosigt `band`, `qv_max`, `qv_ratio`, `quintile`. With `--vcf`, the `gt_*` mirror plus `ref_sum_delta`, `ref_qv` and `gap_closed` (literal `NA` where undefined) |
| `<prefix>.qv_summary.tsv` | `truth_event` / `truth_bp` (the ledger, `ALL` and per svtype), `quintile` / `haplotype` / `called_quintile` distributions, `residual_run`, `excluded` (bubbles with no reference walk, haplotype-bubble pairs not traversed, pairs decomposed coarsely, haplotypes with no VCF column, VCF samples with no path), `called_recon`. With `--vcf`: `gt_quintile` / `gt_haplotype`, `gt_gap` (pooled deltas, `gap_closed_pooled`, `gap_closed_mean`, `gap_closed_undefined`, `worse_than_baseline`), `variation_recovered` for all three levels, `gt_carrier` (TP/FP/FN/TN with precision and recall) and `gt_records` (`applied` / `unplaceable` / `clamped` / `unhandled` / `ref_mismatch` / `heuristic`) |

Nothing is dropped silently: every haplotype-bubble pair excluded from a denominator is counted in the `excluded` scope, so no rate is read as though it covered everything.

Results can be plotted via `scripts/plot_benchmark.R`.
