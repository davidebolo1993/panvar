# Describe Module (Module 4)

Date: 2026-06-07

CLI entrypoint:

- `panvar describe`

## What it does

`describe` converts each called bubble into haplotype features for downstream association analyses,
in two parallel layers that share the same graph coordinates:

For each bubble, it:

1. extracts the canonical source-to-sink walk (and its sequence) for every path crossing the bubble
2. **node + edge dosage** (primary association substrate): counts how many times each path traverses
   each node and each oriented step-to-step edge
3. **k-mer features** (complementary): counts canonical k-mers in each path sequence, sampled with
   closed syncmers by default (or all k-mers / minimizers on request)
4. removes non-discriminative features from both layers
5. writes, per bubble, a node/edge dosage feature map + matrix, and a k-mer feature map (with
   graph-node provenance) + matrix + sparse JSONL

Both layers share the same two-part keep rule:

1. **Copy-number features are always kept.** A feature whose count *varies across the paths that
   carry it* (`min_nonzero_count != max_count`) is informative — it tracks copy number — and is
   retained even when present in every path.
2. **Otherwise a symmetric minor-presence (MAF-style) cut applies.** A feature is dropped when
   `min(present_paths, absent_paths) <= N`, where `N = --min-paths` (default `1`). This removes
   features with no contrast: those present in every path with one identical count, singletons
   (one path), and near-ubiquitous ones (all-but-one path). `--min-paths 0` disables the cut and
   keeps every discriminative feature (the previous behavior).

Because both layers are keyed by node IDs, a significant k-mer maps to nodes via its provenance, and
both map to the same graph location — and to the future graph-native variant call at that bubble.

The per-bubble `kept`/`candidates`/`discarded` counts in `describe.index.tsv` (and the run summary)
show how much each bubble was reduced, so you can tune `--min-paths` with the effect visible.

> **Scope of "discriminative".** Specificity here is evaluated *within a bubble*: a k-mer is
> kept when it distinguishes the bubble's own paths. `describe` does **not** check whether a
> k-mer is unique elsewhere in the genome. Enforcing genome-wide k-mer specificity — needed when
> querying a sequencing sample and attributing a hit back to a node/variant — is a planned
> follow-up tied to the variant-calling rework.

## Required inputs

- `--gfa <graph.gfa>` / `-i <graph.gfa>`
- one of:
  - `--bubble-prefix-in <module1-prefix>` (auto uses `<module1-prefix>.bubbles.csv`)
  - `--bubbles-csv <module1.bubbles.csv>`

## Key options

- `--out-dir <dir>`: output directory (default `describe_out`)
- `--bubble-id <N>`: restrict to a bubble ID; repeatable
- `--kmer-size <K>`: k-mer size, `1..31` (default `31`)
- `--feature-mode <all|minimizer|syncmer>`: feature sampling mode (default `syncmer`)
- `--minimizer-window <W>`: window of consecutive k-mers for minimizer mode (default `15`)
- `--syncmer-s <S>`: internal s-mer size for closed syncmer mode (default auto)
- `--min-paths <N>`: drop features with `min(present, absent)` paths `<= N`, keeping copy-number features (default `1`; `0` keeps all discriminative features)
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

`--feature-mode syncmer` is the default and keeps closed syncmers: a k-mer is retained when the smallest internal s-mer is located at either end of the canonical k-mer. This gives an even sequence-space sampling, avoids the window-boundary dependence of minimizers, and degrades gracefully under substitutions (a single base change disturbs only the syncmers overlapping it). Use `--syncmer-s` to tune density; smaller `s` generally keeps more k-mers. When `--syncmer-s` is unset, `s` is chosen automatically from `k` as `max(1, min(11, (k+2)/3))`.

`--feature-mode all` counts every canonical k-mer. It is the most complete representation, but large bubbles can produce many features; use it when you want the exhaustive marker set rather than a sampled one.

`--feature-mode minimizer` keeps only minimizer k-mers: for each window of `--minimizer-window` consecutive k-mers, the smallest canonical encoded k-mer is selected. The same selected k-mer position is counted once even if overlapping windows choose it repeatedly. This can reduce the feature space roughly by the window size while preserving local anchors.

## Choosing k and feature mode

- **k-mer size.** The default `k=31` is deliberate: 31 bp is long enough to be (near) unique
  genome-wide, which matters downstream when a sequencing sample is queried for these markers and
  each hit must be attributable to this locus rather than a paralog elsewhere. The graph paths
  `describe` reads are assembled haplotypes (no sequencing error), so there is no error-tolerance
  reason to shorten `k` for the description itself. Lower `k` only when a bubble is shorter than
  31 bp or when you intentionally want markers shared across more paths.
- **Feature mode.** `syncmer` (default) is the right balance for a compact, robust marker set.
  Switch to `all` when you need every distinguishing k-mer (e.g. exhaustive fine-mapping within a
  small bubble) and accept a wider matrix. `minimizer` is offered for compatibility but is
  generally dominated by `syncmer` for this use case.


## Node and edge dosage (primary association substrate)

For each bubble, `describe` emits a node/edge dosage table built from the same canonical walks used
for the k-mer features:

- **node dosage** — a **descriptive traversal count**: how many times each path traverses each graph
  node in the bubble. It is **not** a copy-number call. A count > 1 can mean a genuine tandem
  duplication (`…A,A…`) *or* the same node simply revisited elsewhere in the walk (`…A,B,A…`), which
  is just sequence structure; and for short nodes (e.g. length 1) a repeated traversal carries no CN
  meaning at all. Treat node dosage as presence/abundance evidence, not as duplication.
- **edge dosage**: how many times each path uses each oriented transition between consecutive nodes,
  written as `from_node±>to_node±`. Edges are adjacency-aware, so they are the better tandem signal:
  a tandem block `…A,B,A,B…` yields edge `A+>B+` with dosage 2 and a self-loop array `…X,X…` yields
  `X+>X+`, whereas a node revisited non-adjacently produces different flanking edges and no repeated
  edge. Because orientation is part of the edge identity, an inverted traversal also appears as
  distinct edge features (the inversion signal lives here, not in the node columns).

The same discriminative filter applies: constant features (e.g. the bubble's source/sink nodes, which
every path crosses exactly once) are dropped. These tables are the **primary substrate for
phenotype association** — a node or edge maps directly to a graph location and, downstream, to a
graph-native variant call at that bubble. The k-mer tables are the complementary sequence-level
layer (read-queryable, with node provenance). Node dosages match `inspect`'s `node_counts.tsv`
totals for the same path and node, by construction.

> **Copy number is deferred to variant calling, not decided here.** Real CN/duplication needs tandem
> detection — finding *adjacent* repeated blocks above a minimum motif and region size, and reporting
> the run length rather than a raw count. That judgement belongs to the (planned) graph-native calling
> step. These dosage tables only describe traversals; the edge layer is where adjacency is visible.

## Outputs

Inside `--out-dir`:

- `describe.index.tsv`
- `describe.params.json`
- `bubble_<id>/graph_features.tsv.gz` — node/edge dosage feature map
- `bubble_<id>/graph_matrix.tsv.gz` — node/edge dosage matrix (when the wide matrix is enabled)
- `bubble_<id>/kmer_features.tsv.gz`
- `bubble_<id>/kmer_matrix.tsv.gz` when enabled and below cap
- `bubble_<id>/kmer_counts.jsonl.gz`

### `describe.index.tsv`

One row per processed bubble. For both layers it reports how many features were seen
(`candidates`), how many survived the discriminative filter (`kept`), and how many were dropped as
non-informative (`discarded = candidates - kept`), so you can see at a glance how aggressively each
bubble was reduced:

- `bubble_id`
- `status`
- `paths`
- `kmer_candidates`, `kmer_kept`, `kmer_discarded`: k-mer features before/after the filter
- `feature_map_tsv_gz`
- `matrix_tsv_gz`
- `counts_jsonl_gz`
- `matrix_written`
- `matrix_reason`
- `node_candidates`, `node_kept`, `node_discarded`: node dosage features before/after the filter
- `edge_candidates`, `edge_kept`, `edge_discarded`: edge dosage features before/after the filter
- `graph_features_tsv_gz`, `graph_matrix_tsv_gz`, `graph_matrix_written`

The same kept/candidates/discarded totals (summed over bubbles) are printed in the `describe` run
summary.

### `bubble_<id>/graph_features.tsv.gz`

Node/edge dosage feature map columns:

- `feature_id`: integer feature ID (per type)
- `feature_name`: matrix column name (`N<id>` for nodes, `E<id>` for edges)
- `feature_type`: `node` or `edge`
- `label`: node ID, or `from_node±>to_node±` for an edge
- `paths_present`: number of paths with dosage > 0
- `min_count`: minimum dosage across all paths, including zero for absent paths
- `max_count`
- `total_count`

### `bubble_<id>/graph_matrix.tsv.gz`

Columns: `bubble_id`, `sample`, `haplotype`, `path_name`, then one column per retained node feature
(`N1..`) and edge feature (`E1..`); values are dosage counts.

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

## Using the node/edge matrix for association

`graph_matrix.tsv.gz` is a **graph-path × feature** dosage table (one row per panel haplotype that
crosses the bubble; columns `N1..`/`E1..`). To associate features with a phenotype across a cohort of
**genotyped samples**, you go from panel haplotypes to samples and then test each feature:

1. **Project to per-sample dosage.** For each genotyped sample, COSIGT assigns two panel haplotypes
   (e.g. `sampleX → {HG…#1, HG…#2}`). The sample's dosage for a feature is the **sum of the two
   assigned haplotypes' rows** for that column (so a feature on both haplotypes scores 2). This gives
   a sample × feature matrix. *(This projection step is not built into panvar yet — it is the next
   round-2 task, pending the COSIGT genotype-table format. For now the matrix is per-panel-haplotype.)*
2. **Join the phenotype** (sample → trait) and **test each feature** against it: logistic regression
   for a binary trait, linear for a quantitative one, with covariates as needed; correct for multiple
   testing across features.
3. **Map hits back:** a significant `N*`/`E*` column → its node/edge in `graph_features.tsv.gz` → a
   graph location in the bubble → (round 2) the variant called there. Edge hits are the adjacency-aware
   ones (tandem/junction structure); node hits are presence/abundance.

Illustrative per-feature test in Python (after the projection of step 1 has produced
`sample_x_feature` and a `phenotype` series):

```python
import pandas as pd, statsmodels.api as sm

# X: samples (rows) x features (N*/E* columns) of dosages; y: binary phenotype per sample
X = sample_x_feature            # from step 1 (sum of the 2 COSIGT haplotype rows per sample)
y = phenotype.loc[X.index]

results = []
for feat in X.columns:
    Xf = sm.add_constant(X[[feat]])
    fit = sm.Logit(y, Xf).fit(disp=0)        # or sm.OLS for a quantitative trait
    results.append((feat, fit.params[feat], fit.pvalues[feat]))

assoc = pd.DataFrame(results, columns=["feature", "beta", "pval"]).sort_values("pval")
# then map assoc.feature -> graph_features.tsv.gz (label = node id or from±>to±)
```

For large cohorts, write the sample × feature matrix in PLINK/BGEN-style dosage format and run a
proper mixed-model tool (GEMMA, regenie, PLINK) instead of per-feature loops — `describe` deliberately
emits tables for these external tools rather than computing association itself.

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
