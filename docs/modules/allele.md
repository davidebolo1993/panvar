# Allele Module (Module 2)

Date: 2026-04-20

CLI entrypoint:

- `panvar allele`

## What it does

Given module-1 sites (`--bubble-prefix-in` or `--bubbles-csv-in`), this module:

1. extracts one canonical source-to-sink allele per path per bubble
2. deduplicates identical alleles
3. clusters alleles by similarity (or by predefined path groups via JSON)
4. picks one representative allele per cluster
5. writes cluster and per-path assignment outputs

## Inputs and outputs

Required inputs:

1. `--gfa <graph.gfa>`
2. one of:
   - `--bubble-prefix-in <module1-prefix>` (auto uses `<module1-prefix>.bubbles.csv`)
   - `--bubbles-csv-in <module1.bubbles.csv>`

`allele` accepts both the current simplified module-1 CSV schema and older legacy module-1 schemas.

Core outputs:

- `*.allele_clusters.csv`
- `*.allele_assignments.csv`

Optional outputs:

- representative sequences: `--cluster-sequences-csv`
- similarity diagnostics: `--similarity-out-dir`
- ODGI viz inputs: `--odgi-viz-out-dir`

Current module-2 handoff schemas:

- `*.allele_clusters.csv`:
  - `bubble_id,source,sink,cluster_id,representative_allele_id,total_path_support,member_alleles`
- `*.allele_assignments.csv`:
  - `bubble_id,source,sink,path_name,cluster_id,allele_id,allele_length,interval_start,interval_end,source_to_sink`

Notes:

- only fields consumed by downstream modules are retained in these CSVs
- older module-2 CSVs with extra columns remain acceptable to `panvar call` as long as the required columns above are present

Optional clustering input:

- `--clusters-json <path>` path->cluster-label JSON map
  - expected format: `{"path_name":"cluster_label", ...}`
  - if provided, module-2 skips distance-based clustering and uses these labels
  - paths missing from JSON are put in synthetic singleton labels (warning emitted)

`--cluster-sequences-csv` and `--clusters-json` are independent:

- `--clusters-json` controls how clusters are assigned
- `--cluster-sequences-csv` only exports representative sequences for whatever clusters were produced

## Algorithm overview

For each bubble:

1. scan all P/W paths and find the best source-to-sink interval crossing the bubble
2. canonicalize interval orientation to source->sink
3. build allele signature from oriented node steps and deduplicate
4. if `--clusters-json` is provided, assign clusters by predefined path label; otherwise:
5. build comparison tokens based on `--cluster-mode`
6. compute pairwise distances based on `--distance-mode`
7. cluster alleles
8. select representative allele (medoid-like tie breaking)
9. write bubble-level and path-level results

Best interval means:

- maximize number of bubble-internal nodes
- tie-break with shorter span
- tie-break again with earlier interval

## Cluster mode behavior

### `--cluster-mode sequence`

Token type:

- spelled DNA sequence for each allele

Clustering strategy:

- default path: UPGMA tree cut at `1 - min_similarity`
- if `unique_alleles > max_upgma_alleles`, switch to threshold-graph connected components
- in `--distance-mode auto`, very large sequence bubbles may switch to a sequence fast path:
  greedy threshold assignment using sketch-estimated distances

### `--cluster-mode walk`

Token type:

- oriented node-step token list (graph walk tokens)

Clustering strategy:

- default path: UPGMA tree cut at `1 - min_similarity`
- if `unique_alleles > max_upgma_alleles`, switch to threshold-graph connected components

## Distance mode behavior

### `sequence + exact`

- full edit distance on sequence tokens for all compared pairs
- strictest and usually slowest

### `sequence + auto`

- MinHash sketch prefilter (fast reject/fast accept where clear)
- bounded edit distance for uncertain near-threshold pairs
- very large bubbles can auto-switch to sequence fast path (sketch-estimated distances only)

### `walk + exact`

- full edit distance on walk tokens for all compared pairs

### `walk + auto`

- MinHash sketch rejection for clearly dissimilar pairs
- bounded edit distance for uncertain pairs
- for very large bubbles in threshold-graph path, can use sketch-estimated walk distances

## Representative sequence policy

- no consensus is built in module 2
- representative sequence is the sequence of the selected representative allele
- representative selection prefers:
  1. smaller max intra-cluster distance
  2. smaller mean intra-cluster distance
  3. higher path support
  4. smaller allele ID

## Key options

- `--min-similarity <X>`: threshold in `(0,1]` or percent (default `0.90`)
- `--cluster-mode sequence|walk` (default `sequence`)
- `--distance-mode auto|exact` (default `auto`)
- `--clusters-json <path>`: predefined path-group JSON (optional)
- `--threads <N>`: workers for distance calculations (`0` = auto)
- `--max-upgma-alleles <N>`: clustering-mode switch threshold (`0` disables)
  - if unique alleles in one bubble are `<= N`, clustering uses UPGMA tree cutting
  - if unique alleles are `> N`, clustering switches to threshold-graph connected components
  - this is not just for tree visualization; it changes the clustering algorithm for that bubble
- `--quiet`: disable progress logs

## Performance guidance

- start with `sequence + auto` for balanced sensitivity/speed
- use `walk + auto` when locus scale is very large and topology-first grouping is acceptable
- use `exact` modes for strict validation on selected loci
- similarity reports (`--similarity-out-dir`) are useful for debugging but add extra I/O

## Example

Default sequence clustering:

```bash
./build/panvar allele \
  -i tests/real_data/c4.gfa \
  -o tests/results/c4/allele \
  --bubble-prefix-in tests/results/c4/bubble
```

Walk clustering:

```bash
./build/panvar allele \
  -i tests/real_data/c4.gfa \
  -o tests/results/c4/allele_walk \
  --bubble-prefix-in tests/results/c4/bubble \
  --cluster-mode walk
```

Strict sequence distances:

```bash
./build/panvar allele \
  -i tests/real_data/c4.gfa \
  -o tests/results/c4/allele_exact \
  --bubble-prefix-in tests/results/c4/bubble \
  --cluster-mode sequence \
  --distance-mode exact
```

Predefined clustering from JSON:

```bash
./build/panvar allele \
  -i tests/real_data/c4.gfa \
  -o tests/results/c4/allele_json \
  --bubble-prefix-in tests/results/c4/bubble \
  --clusters-json tests/real_data/c4.clusters.json
```

## Module handoff

Use module-2 outputs in module-3:

- recommended shortcut: `panvar call --allele-prefix-in <prefix>`
- explicit files:
  - `panvar call --clusters-csv-in <prefix>.allele_clusters.csv`
  - `panvar call --assignments-csv-in <prefix>.allele_assignments.csv`
