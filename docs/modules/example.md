# Worked Examples (Clustering + Calling)

Date: 2026-05-15

This page gives small concrete examples for the core `allele` and `call` behaviors.

## Example 1: Allele Clustering (walk default vs sequence mode)

Bubble has 4 unique alleles:

- `A1`: `ACCTGAAATTTGG`
- `A2`: `ACCTGAAATTTGG` (identical to `A1`)
- `A3`: `ACCTGAAACTTGG` (1 mismatch vs `A1`)
- `A4`: `ACCTGCCCCCTTGG` (larger change)

With default walk mode (`--cluster-mode walk`, now default), step mismatches are weighted by node bp length.
For sequence-based behavior, use `--cluster-mode sequence --min-similarity 0.90`:

1. `A1` and `A2` deduplicate first.
2. Pairwise distances are computed on sequence tokens.
3. `A3` is usually grouped with `A1`/`A2` (high similarity).
4. `A4` typically forms a separate cluster.

If `unique_alleles > --max-upgma-alleles`, the algorithm switches from UPGMA tree-cut to threshold-graph connected components for that bubble.

## Example 2: Predefined Clusters (`--clusters-json`)

If JSON contains:

```json
{
  "sampleA#hap1": "C4_long",
  "sampleB#hap2": "C4_long",
  "sampleC#hap1": "C4_short"
}
```

Then `allele` uses these labels directly for cluster assignment (no distance clustering for those paths).  
`--cluster-sequences-csv` still exports one representative sequence per resulting bubble+cluster.

## Example 3: Calling from CIGAR vs Split Evidence

Reference:

- `R`: `AAAACCCCGGGGTTTT`

Cluster representative:

- `Q`: `AAAACCCCNNNNGGGGTTTT`

If best minimap2 CIGAR includes `4I` between `CCCC` and `GGGG`, `call` emits an `INS` candidate from CIGAR.

If there are multiple split chains with a large query/reference span deviation, `call` also emits split-derived candidates.

Final per-cluster event set is the union:

- CIGAR-derived events
- split-derived events
- optional fallback net `INS`/`DEL` if both above miss a clear large shift

Then nearby same-type fragments can be compacted by within-cluster merge rules.

## Example 4: Cross-Cluster Merge in Final VCF

Two clusters may carry the same event:

- Cluster 2: `INS` at ref position ~10000, inserted seq length 220
- Cluster 5: `INS` at ref position ~10012, inserted seq length 224

They can merge into one VCF record when:

- event type matches
- coordinate windows pass (`--vcf-merge-window-bp` or lenient fallback window)
- sequence similarity passes (`--vcf-merge-min-seq-sim`)
- bounded edit fraction passes (effective cap `min(--vcf-merge-max-edit-frac, 1 - --vcf-merge-min-seq-sim)`)

`--vcf-merge-lenient-min-ref-jaccard` is only used in lenient fallback for non-INS events to require interval overlap on reference.

## Example 5: Debug Status Outputs

With `panvar call --debug`:

- global status: `<debug>/debug_summary.tsv`
- bubble status: `<debug>/bubble_<id>/bubble_status.tsv`
- cluster status table: `<debug>/bubble_<id>/cluster_status.tsv`
- per-cluster quick note: `<debug>/bubble_<id>/cluster_<id>/status.txt`

So missing `dotplot.svg` or absent events are traceable to explicit statuses such as:

- `skipped:no-precomputed-assignments`
- `skipped:no-minimap-best-hit`
- `ok:no-sv-events`
- `ok:sv-events`

## Example 6: Why One Allele Can Have Many Paths

In `allele_assignments.csv`, multiple paths can share the same `allele_id` within one bubble.
This is expected: `allele_id` identifies a deduplicated unique allele walk, while each CSV row is a path assignment.

In `allele_clusters.csv`:

- `member_alleles` / `member_allele_count` refer to unique alleles in that cluster
- `total_path_support` counts how many paths map to those alleles

So `member_allele_count=1` with `total_path_support=220` means:
one unique allele sequence/walk is shared by 220 paths.
