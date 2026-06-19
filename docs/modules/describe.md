# Describe Module (Module 4)

Date: 2026-06-18

CLI entrypoint:

- `panvar describe`

## What it does

`describe` converts each called bubble into haplotype features for downstream association analyses,
in two parallel layers that share the same graph coordinates:

For each bubble, it:

1. extracts the canonical source-to-sink walk (and its sequence) for every path crossing the bubble
2. **k-mer features** (primary association substrate): counts canonical k-mers in each path sequence, sampled with closed syncmers by default (or all k-mers); read-queryable and the tested GWAS substrate
3. **node + edge dosage** (complementary): counts how many times each path traverses each node and each oriented step-to-step edge; graph-local, mapping straight to a node/edge and the variant `call` makes there
4. removes non-discriminative features from both layers
5. writes, per bubble, a k-mer feature map (with graph-node provenance) + matrix + sparse JSONL, and a node/edge dosage feature map + matrix

It is the **GWAS feature-extraction** step: besides the per-bubble tables it writes a single pooled
k-mer file (`<out-dir>/fsm_kmers.txt.gz`, the `<feature> | strain:count` fsm-lite text format that
count-based and presence/absence-based association tools read directly). For a full, worked sample-level
GWAS (copy number → quantitative trait, multiplicity vs presence/absence, Manhattan/QQ, traceback) see
the **[GWAS example](../gwas_example.md)**.

**Sample-level output (`--samples <cosigt.tsv>`):** a GWAS tests **samples**, not haplotypes. Pass a
cosigt-style table (one sample per line: `sample <tab> hap1 <tab> hap2 …`, haplotype names = graph path
names) and `describe` additionally writes a per-**sample** file for **each** substrate —
`<out-dir>/fsm_kmers.samples.txt.gz` (k-mers) and `<out-dir>/fsm_graph.samples.txt.gz` (node/edge
dosage) — where each strain is a **sample** and the value is the **summed dosage over its assigned
haplotypes** (a haplotype listed twice counts twice = homozygous). The per-haplotype tables are still
written.

**Two scopes for the k-mer markers:**
- **whole-graph (default)** — every called bubble is described; markers span all the variable sites.
- **variant-restricted (`--variant-nodes <call.variant_nodes.tsv>`)** — only the bubbles in that sidecar
  are processed, and within each, k-mer/syncmer generation is confined to the **called-variation nodes**
  (bases from other nodes are masked, so no k-mer spans a non-variant node). This is the
  genotyping/association scope: markers only where `call` found variation. `call` writes
  `<prefix>.variant_nodes.tsv` (`variant_id, bubble_id, svtype, node_ids`) automatically — its node set
  per variant is the **candidate pool**: every node any merged carrier's allele uses for that variant.
  Masking is then applied **per path**: each haplotype is sketched on `(its own walk) ∩ (the pool)`, i.e.
  **its own** variant nodes — the pool only supplies candidates, so a node a haplotype does not traverse
  has no effect on it. (This is also what makes deletions testable: the deleted nodes sit in the pool, the
  reference-like haplotypes keep them and the deletion carriers mask them, giving the contrast.)
  - **`--variant-flank-bp <N>`** widens the unmasked window: a node is also kept when its path-distance
    to the nearest variant node is `≤ N` bp on either side, so k-mers spanning the variant's **flanking
    sequence** (e.g. nearby SNPs) survive the masking. `0` (default) is strict variant-node-only.
    Running the association with `--variant-flank-bp 0` vs a positive value shows how much flanking
    context adds.

Both layers share the same two-part keep rule:

1. **Copy-number features are always kept.** A feature whose count *varies across the paths that carry it*
   is informative — it tracks copy number — and is retained even when present in every path.
2. **Otherwise a symmetric minor-presence (MAF-style) cut applies.** A feature is dropped when the smaller
   of {paths carrying it, paths lacking it} is `≤ --min-paths` (default `1`). This removes features with no
   contrast: those present in every path with one identical count, singletons (one path), and
   near-ubiquitous ones (all-but-one path). `--min-paths 0` disables the cut and keeps every discriminative
   feature.

Because both layers are keyed by node IDs, a significant k-mer maps to nodes via its provenance, and
both map to the same graph location — and to the future graph-native variant call at that bubble.

The per-bubble `kept`/`candidates`/`discarded` counts in `describe.index.tsv` (and the run summary)
show how much each bubble was reduced, so you can tune `--min-paths` with the effect visible.

## Required inputs

`describe` runs after `panphorte`, so its graph is the **panphorte-normalized/sorted GFA** and its
bubbles come from the **panphorte** prefix:

- `--gfa <panphorte_prefix.normalized.sorted.gfa>` / `-i`
- one of:
  - `--bubble-prefix-in <panphorte_prefix>` (auto uses `<panphorte_prefix>.bubbles.csv`)
  - `--bubbles-csv <panphorte_prefix.bubbles.csv>`

## Key options

```text
panvar describe -i <panphorte.normalized.sorted.gfa> (-b <panphorte_prefix> | -c <bubbles.csv>) [-o <dir>] [options]
```

- `-i, --gfa <path>`: panphorte-normalized/sorted GFA (required)
- `-b, --bubble-prefix-in <prefix>`: panphorte output prefix (auto-uses `<prefix>.bubbles.csv`)
- `-c, --bubbles-csv <path>`: panphorte bubbles CSV (required if no prefix)
- `-o, --out-dir <dir>`: output directory (default `describe_out`)
- `--bubble-id <N>`: restrict to a bubble ID; repeatable
- `-k, --kmer-size <K>`: k-mer size, `1..31` (default `31`)
- `--feature-mode <all|syncmer>`: feature sampling mode (default `syncmer`)
- `--syncmer-s <S>`: internal s-mer size for closed syncmer mode (default auto)
- `--min-paths <N>`: drop features with `min(present, absent)` paths `<= N`, keeping copy-number features (default `1`; `0` keeps all discriminative features)
- `--max-wide-features <N>`: skip wide matrix above this number of features (default `250000`; `0` disables cap)
- `--force-wide`: write the wide matrix even above the cap
- `--no-wide-matrix`: write only feature map and sparse JSONL counts
- `--variant-nodes <tsv>`: restrict k-mer generation to `call`'s `<prefix>.variant_nodes.tsv` (only those
  bubbles' variant nodes contribute k-mers; see "Two scopes" above)
- `--variant-flank-bp <N>`: with `--variant-nodes`, also keep nodes within `N` bp of a variant node so
  flanking-SNP k-mers are retained (default `0` = strict variant-node-only)
- `--samples <cosigt.tsv>`: write per-sample dosage files for both substrates (see "Sample-level output")
- `--no-pyseer`: do not write the pooled `fsm_kmers.txt.gz`
- `--quiet`: disable the progress bar

## K-mer Encoding

K-mers are encoded as 2-bit integers, with `A=0`, `C=1`, `G=2`, `T=3`.

By default, k-mers are canonicalized against their reverse complement. This means a k-mer and its reverse complement are counted as the same feature, using the smaller 2-bit encoding. K-mers containing non-ACGT bases are skipped.

The matrix uses compact feature names (`K1`, `K2`, ...). The feature map resolves each feature ID back to the encoded integer and DNA string.

Because of canonicalization, `describe` does not keep separate forward and reverse-complement features. If a forward k-mer is discriminative, its reverse-complement representation maps to the same feature ID. This avoids duplicating equivalent markers.

## Feature mode and k

Two sampling modes, set with `--feature-mode`:

- **`syncmer` (default)** keeps closed syncmers: a k-mer is retained when its smallest internal s-mer
  sits at either end of the canonical k-mer. This samples sequence space evenly and degrades gracefully
  under substitutions (a single base change disturbs only the syncmers overlapping it). `--syncmer-s`
  tunes density (smaller `s` keeps more k-mers); unset, `s` is chosen from `k` as `max(1, min(11, (k+2)/3))`.
- **`all`** counts every canonical k-mer — the exhaustive marker set. Use it for fine-mapping a small
  bubble; large bubbles can produce many features and a wide matrix.

The default **`k=31`** is a good fit here: the called bubbles are well over 31 bp, so every variable site
yields k-mers, and the graph paths are assembled haplotypes (low sequencing error), so there is no
error-tolerance reason to shorten `k`. Each marker is tested locally at its own bubble, so genome-wide
k-mer uniqueness is not a concern. Lower `k` only for bubbles shorter than 31 bp, or when you
intentionally want markers shared across more paths.

## Node and edge dosage (complementary association substrate)

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
every path crosses exactly once) are dropped. These tables are a **complementary substrate for phenotype
association** — a node or edge maps directly to a graph location and, downstream, to the variant `call`
makes at that bubble — alongside the k-mer tables (the primary, read-queryable layer with node provenance).
With `--samples` they are aggregated to per-sample dosage exactly like the k-mers
(`fsm_graph.samples.txt.gz`), so the same association tests run on either substrate. Node dosages match
`inspect`'s `node_counts.tsv` totals for the same path and node, by construction.

## Outputs

Inside `--out-dir`:

- `describe.index.tsv`
- `describe.params.json`
- `bubble_<id>/graph_features.tsv.gz` — node/edge dosage feature map
- `bubble_<id>/graph_matrix.tsv.gz` — node/edge dosage matrix (when the wide matrix is enabled)
- `bubble_<id>/kmer_features.tsv.gz`
- `bubble_<id>/kmer_matrix.tsv.gz` when enabled and below cap
- `bubble_<id>/kmer_counts.jsonl.gz`
- `fsm_kmers.txt.gz` — pooled per-haplotype k-mer file (unless `--no-pyseer`); see below
- with `--samples`: `fsm_kmers.samples.txt.gz` (per-sample k-mer dosage) and
  `fsm_graph.samples.txt.gz` (per-sample node/edge dosage) — same fsm-lite layout, keyed by sample

### `fsm_kmers.txt.gz`

The discriminative k-mers from every processed bubble, pooled into one **fsm-lite-format** file — the
`<feature> | strain:count` text format count-based and presence/absence-based association tools consume
directly. One line per k-mer:

```text
<kmer_sequence> | <path>:<count> <path>:<count> ...
```

- The **strain** id is the **path name** (one haplotype = one strain); map your phenotypes to these names.
  The `--samples` files use the same layout but the strain id is a **sample** and the count is its
  diploid dosage (node/edge features in `fsm_graph.samples.txt.gz` use the node id, or `from±>to±` for
  an edge, as the line label).
- The value is the **true per-strain count** (multiplicity), so copy-number expansions stay faithful: a
  count-based tool reads the dosage, a presence/absence tool binarizes any count `> 0` to "present".
- Only **carriers** (count `> 0`) are listed; a path missing from a line — including one that does not
  traverse that bubble at all — is implicit **absence (0)**. (fsm-lite has no missing-data state, so
  "reference allele" and "locus absent" are both encoded as 0.)
- The k-mer set is exactly the per-bubble **discriminative** features (count varies, or
  `min(present,absent) > --min-paths`), pooled — not the raw exhaustive k-mer set. The same canonical
  k-mer seen in two bubbles has its counts summed.

See the [Association](#association) section for how this is tested, and
[gwas_example.md](../gwas_example.md) for a worked run.

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
{"bubble_id":1,"sample":"sampleA","haplotype":"1","path_name":"sampleA#1#...","path_length_bp":7941,"kmers":[[1,7],[2,2]]}
```

Use `bubble_<id>/kmer_features.tsv.gz` to map feature IDs back to sequences and to recover the bubble-local node union for each feature. Keeping node provenance in the feature map rather than every JSONL tuple keeps the sparse JSONL much smaller.

This JSONL is the sparse representation: absent/zero-count k-mers are omitted from each path row. The wide matrix is still emitted by default because many GWAS and association tools expect dense sample-by-feature tables.

## Association

**The unit is a sample.** A GWAS tests genotyped **samples**, not panel haplotypes. In `panvar`'s model a
sample's genotype is the **two haplotypes COSIGT assigns it**, and a marker's value for that sample is the
**sum over those two haplotypes** (a feature on both scores twice = homozygous). For k-mers, `describe
--samples <cosigt.tsv>` produces exactly this per-sample file (`fsm_kmers.samples.txt.gz`) directly.

**Primary substrate — k-mers (tested).** The pooled k-mer file is the substrate we test: sequence
markers, read-queryable, with node provenance. A marker can be read two ways — **presence/absence**
(binarize count `> 0`) or **count/multiplicity** (the value is local copy number). The distinction is
decisive at a **copy-number locus**: the repeat-unit k-mer is present in *every* sample, so presence/absence
has no contrast and misses it, while the **count** test recovers it — the whole point of carrying true
counts. `scripts/gwas_demo.py` reads the fsm file + a phenotype and runs **both** tests (linear regression
for a quantitative trait; Fisher exact / Mann–Whitney U for a binary one; Benjamini-Hochberg `q` +
Bonferroni), as a stand-in for a count-based or presence/absence-based association method. On a
literature-plausible copy-number phenotype the count test concentrates on the copy-number bubble while
presence/absence finds nothing there — worked run, Manhattan/QQ and traceback in
[gwas_example.md](../gwas_example.md). The same fsm-lite file also feeds external association tools;
`panvar` does not run the association itself.

**Complementary substrate — node/edge dosage.** `graph_matrix.tsv.gz` (`N*`/`E*` columns) is a per-**panel
haplotype** dosage table that maps directly to a node/edge and the variant `call` makes there (edges carry
the adjacency/tandem signal, nodes carry presence/abundance). With `--samples` it is aggregated to the same
COSIGT per-sample dosage as the k-mers (`fsm_graph.samples.txt.gz`), so `scripts/gwas_demo.py
--substrate graph` runs the identical presence/absence-vs-count tests on it. At a copy-number locus the
count test localizes the signal to the repeat node and its self-loop edge.

**Interpretation.** A significant marker is not automatically the causal variant — it is a sequence/graph
marker. Map a hit back through `kmer_features.tsv.gz` / `graph_features.tsv.gz` to its nodes, the bubble,
and the variant called there, and read it together with the allele clusters (`inspect`) and annotations.

## Algorithm overview

The encoding, syncmer sampling, and discriminative filter are traced on a tiny worked dataset
in [algorithm_example.md](../algorithm_example.md). For each bubble:

1. Find each path's best source/sink interval crossing the bubble.
2. Canonicalize sink-to-source intervals into source-to-sink orientation.
3. Spell the path sequence from graph node sequences.
4. Count canonical 2-bit k-mers per path, applying syncmer sampling unless `--feature-mode all`.
5. Accumulate per-k-mer counts across paths.
6. Drop k-mers whose count is identical in every path.
7. Re-scan retained k-mers, collect bubble-local node provenance for the feature map, and write sparse JSONL counts.
8. Write the feature map with sequence and node-union metadata.
9. Write the wide matrix unless disabled or above the feature cap.

The implementation streams per bubble: one discovery pass finds retained features, the sparse JSONL pass records only retained k-mer counts while collecting node provenance for the feature map, and the optional wide-matrix pass is written only when requested and below the feature cap. This avoids storing the full path-by-kmer matrix in memory.

## Example

```bash
# whole graph, default closed-syncmer markers (k = 31), on the panphorte-normalized graph
./build/panvar describe \
  -i results/real_data/lpa/panphorte/panphorte.normalized.sorted.gfa \
  --bubble-prefix-in results/real_data/lpa/panphorte/panphorte \
  --out-dir results/real_data/lpa/describe \
  --kmer-size 31
```

Restrict markers to `call`'s variant nodes (genotyping/association scope), keeping 50 bp of flanking
context, and aggregate to per-sample dosage for both substrates:

```bash
./build/panvar describe \
  -i results/real_data/lpa/panphorte/panphorte.normalized.sorted.gfa \
  --bubble-prefix-in results/real_data/lpa/panphorte/panphorte \
  --variant-nodes results/real_data/lpa/call/call.variant_nodes.tsv \
  --variant-flank-bp 50 \
  --samples results/real_data/lpa/gwas/samples.tsv \
  --out-dir results/real_data/lpa/describe
```

Exhaustive markers for fine-mapping one bubble:

```bash
./build/panvar describe \
  -i results/real_data/lpa/panphorte/panphorte.normalized.sorted.gfa \
  --bubble-prefix-in results/real_data/lpa/panphorte/panphorte \
  --bubble-id 7 \
  --out-dir results/real_data/lpa/describe \
  --feature-mode all \
  --kmer-size 21 \
  --max-wide-features 0
```
