# Describe Module (Module 4)

Date: 2026-06-07

CLI entrypoint:

- `panvar describe`

## What it does

`describe` converts each called bubble into sequence-derived haplotype features for downstream association analyses.

For each bubble, it:

1. extracts the canonical source-to-sink path sequence for every path crossing the bubble
2. counts canonical k-mers, or an optional sampled subset of k-mers, in each path sequence
3. removes non-discriminative k-mers
4. writes a compact feature map with graph-node provenance, a GWAS-style count matrix, and a sparse JSONL representation

A feature k-mer is removed as non-discriminative only when all paths have exactly the same count for that k-mer.
A feature k-mer is retained when it is absent from some paths, present only in some paths, or present in all paths with different copy counts.

## Required inputs

- `--gfa <graph.gfa>` / `-i <graph.gfa>`
- one of:
  - `--bubble-prefix-in <module1-prefix>` (auto uses `<module1-prefix>.bubbles.csv`)
  - `--bubbles-csv <module1.bubbles.csv>`

## Key options

- `--out-dir <dir>`: output directory (default `describe_out`)
- `--bubble-id <N>`: restrict to a bubble ID; repeatable
- `--kmer-size <K>`: k-mer size, `1..31` (default `31`)
- `--feature-mode <all|minimizer|syncmer>`: feature sampling mode (default `all`)
- `--minimizer-window <W>`: window of consecutive k-mers for minimizer mode (default `15`)
- `--syncmer-s <S>`: internal s-mer size for closed syncmer mode (default auto)
- `--max-wide-features <N>`: skip wide matrix above this number of features (default `250000`; `0` disables cap)
- `--force-wide`: write the wide matrix even above the cap
- `--no-wide-matrix`: write only feature map and sparse JSONL counts
- `--quiet`: disable progress logs

## K-mer Encoding

K-mers are encoded as 2-bit integers, with `A=0`, `C=1`, `G=2`, `T=3`.

By default, k-mers are canonicalized against their reverse complement. This means a k-mer and its reverse complement are counted as the same feature, using the smaller 2-bit encoding. K-mers containing non-ACGT bases are skipped.

The matrix uses compact feature names (`K1`, `K2`, ...). The feature map resolves each feature ID back to the encoded integer and DNA string.

Because of canonicalization, `describe` does not keep separate forward and reverse-complement features. If a forward k-mer is discriminative, its reverse-complement representation maps to the same feature ID. This avoids duplicating equivalent markers.

## Feature Sampling Modes

`--feature-mode all` is the default and counts every canonical k-mer. It is the most complete representation, but large bubbles can produce many features.

`--feature-mode minimizer` keeps only minimizer k-mers: for each window of `--minimizer-window` consecutive k-mers, the smallest canonical encoded k-mer is selected. The same selected k-mer position is counted once even if overlapping windows choose it repeatedly. This can reduce the feature space roughly by the window size while preserving local anchors.

`--feature-mode syncmer` keeps closed syncmers: a k-mer is retained when the smallest internal s-mer is located at either end of the canonical k-mer. This usually gives a more even sequence-space sampling than minimizers and avoids the window-boundary dependence of minimizers. Use `--syncmer-s` to tune density; smaller `s` generally keeps more k-mers.


## Outputs

Inside `--out-dir`:

- `describe.index.tsv`
- `describe.params.json`
- `bubble_<id>/kmer_features.tsv.gz`
- `bubble_<id>/kmer_matrix.tsv.gz` when enabled and below cap
- `bubble_<id>/kmer_counts.jsonl.gz`

### `describe.index.tsv`

One row per processed bubble:

- `bubble_id`
- `status`
- `paths`
- `discriminative_features`
- `feature_map_tsv_gz`
- `matrix_tsv_gz`
- `counts_jsonl_gz`
- `matrix_written`
- `matrix_reason`

### `bubble_<id>/kmer_features.tsv.gz`

Feature map columns:

- `feature_id`: integer feature ID
- `feature_name`: matrix column name (`K<feature_id>`)
- `encoded_kmer`: 2-bit encoded canonical k-mer
- `kmer`: decoded canonical DNA k-mer
- `paths_present`: number of paths with count > 0
- `min_count`: minimum count across all paths, including zero for absent paths
- `max_count`
- `total_count`
- `node_count`: number of bubble-local graph nodes touched by this k-mer across all paths
- `nodes`: semicolon-separated node IDs touched by this k-mer across all paths

### `bubble_<id>/kmer_matrix.tsv.gz`

Wide matrix columns:

- `bubble_id`
- `sample`
- `haplotype`
- `path_name`
- one column per retained k-mer feature (`K1`, `K2`, ...)

Values are k-mer counts. They can be used as binary markers (`count > 0`) or dosage-like numeric features.

### `bubble_<id>/kmer_counts.jsonl.gz`

Sparse representation, one JSON object per path. The `kmers` array stores compact `[feature_id, count]` tuples:

```json
{"bubble_id":1,"sample":"HG00096","haplotype":"1","path_name":"HG00096#1#...","path_length_bp":7941,"kmers":[[1,7],[2,2]]}
```

Use `bubble_<id>/kmer_features.tsv.gz` to map feature IDs back to sequences and to recover the bubble-local node union for each feature. Keeping node provenance in the feature map rather than every JSONL tuple keeps the sparse JSONL much smaller.

This JSONL is the sparse representation: absent/zero-count k-mers are omitted from each path row. The wide matrix is still emitted by default because many GWAS and association tools expect dense sample-by-feature tables.

## Association Interpretation

A k-mer-based association treats each retained k-mer as a marker.

- A significant binary association means presence/absence of that sequence word tracks the phenotype.
- A significant count association means copy number or repeated occurrence of that sequence word tracks the phenotype.
- The signal may tag an allele, insertion, deletion boundary, duplicated unit, gene-copy state, or another sequence context correlated with the causal change.

The k-mer itself is not automatically the causal variant; it is a sequence marker that should be interpreted together with the bubble, allele clusters, dotplots, and gene annotations.

## Algorithm overview

For each bubble:

1. Find each path's best source/sink interval crossing the bubble.
2. Canonicalize sink-to-source intervals into source-to-sink orientation.
3. Spell the path sequence from graph node sequences.
4. Count canonical 2-bit k-mers per path, applying minimizer/syncmer sampling if requested.
5. Accumulate per-k-mer counts across paths.
6. Drop k-mers whose count is identical in every path.
7. Re-scan retained k-mers, collect bubble-local node provenance for the feature map, and write sparse JSONL counts.
8. Write the feature map with sequence and node-union metadata.
9. Write the wide matrix unless disabled or above the feature cap.

The implementation streams per bubble: one discovery pass finds retained features, the sparse JSONL pass records only retained k-mer counts while collecting node provenance for the feature map, and the optional wide-matrix pass is written only when requested and below the feature cap. This avoids storing the full path-by-kmer matrix in memory.

## Example (C4)

```bash
./build/panvar describe \
  -i tests/real_data/c4.gfa \
  --bubble-prefix-in tests/results/c4/bubble \
  --bubble-id 1 \
  --out-dir tests/results/c4/describe \
  --kmer-size 31
```

For a smaller quick check:

```bash
./build/panvar describe \
  -i tests/real_data/c4.gfa \
  --bubble-prefix-in tests/results/c4/bubble \
  --bubble-id 1 \
  --out-dir /tmp/panvar_describe_c4_b1 \
  --kmer-size 21 \
  --max-wide-features 0
```

For a reduced feature set on large bubbles:

```bash
./build/panvar describe \
  -i tests/real_data/c4.gfa \
  --bubble-prefix-in tests/results/c4/bubble \
  --bubble-id 1 \
  --out-dir /tmp/panvar_describe_c4_b1_syncmer \
  --kmer-size 31 \
  --feature-mode syncmer \
  --syncmer-s 11
```
