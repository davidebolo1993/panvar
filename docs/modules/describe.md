# Module `describe`

CLI: `panvar describe`

## What it does

Converts the called bubbles into per-haplotype genotype features for association, on three substrates that share the same graph coordinates and are each exported as a BIMBAM mean-genotype dosage matrix — the canonical input to [`associate`](associate.md) and tools like [GEMMA](https://github.com/genetics-statistics/GEMMA):

- variant — one dosage row per VCF record (or per ALT, multiallelic) from `call`'s VCF (Variant Call Format) (`--variant-vcf`): copy number for a `DUP` — which must carry a usable `CN`, since falling back to `GT` presence would silently substitute a different quantity — and presence for `DEL`/`INS`/`INV`
- k-mers — canonical k-mers per path, [closed-syncmer](../algorithms/describe.md) sampled by default, each carrying its node provenance;
- graph — per-node and per-oriented-edge traversal counts, a graph-local dosage;

The k-mer and graph substrates are built from the graph and emitted by default; the variant substrate is emitted when its VCF is supplied. Any one can be produced on its own (`--only-kmers`/`--only-graph`/`--only-variant`). The k-mer and graph substrates pass through a discriminative filter that drops features which do not vary across haplotypes — those present in all (or all but `--min-paths`) paths, unless their copy number varies — so only informative markers reach the matrix; the variant substrate emits every call and leaves frequency filtering to `associate`.
Pass a [cosigt](https://github.com/davidebolo1993/cosigt) table (`sample <tab> hap1 <tab> hap2 …`) to `--samples <cosigt.tsv>` and `describe` also writes a per-sample BIMBAM for each substrate (under `sample/<substrate>/`) whose value sums the dosage over the sample's haplotypes, so a diploid genotype is the sum of the two. See [GWAS example](../gwas.md).

Outputs are transactional: everything `describe` owns is built in a sibling staging directory and moved into place only once every substrate pass has succeeded. A rerun with a different substrate selection, `--bubble-id` or `--no-wide-matrix` setting removes the previous run's now-stale families rather than leaving them beside the new ones looking current; unrelated files in the output directory are left alone.

Algorithm and worked trace: [algorithms/describe.md](../algorithms/describe.md).



## Required inputs

`describe` must run on the same graph and bubble prefix as `call` so feature node ids line up with the VCF: the `panphorte` `.normalized.sorted.gfa` or the `bubble` `.sorted.gfa`.

- `-i, --gfa <path>` — the call substrate from `panphorte`/`refine`/`bubble` (node ids should match the CSV to `-b`);
- one of `-b, --bubble-prefix-in <prefix>` or `-c, --bubbles-csv <path>`.
- `--variant-nodes <call.variant_nodes.tsv>` to restrict features to `call`'s variant scope (recommended; this is what keeps `describe` and the VCF in lockstep).

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `-o, --out-dir <dir>` | output directory | `describe_out` |
| `-k, --kmer-size <K>` | k-mer size, 1–31 | `31` |
| `--feature-mode <all\|syncmer>` | keep all canonical k-mers vs sampled [closed syncmers](../algorithms/describe.md#k-mer-encoding--syncmer-selection) | `syncmer` |
| `--min-paths <N>` | discriminative-filter cut: drop features with `min(present,absent) ≤ N` (copy-number features always kept; `0` keeps all) | `1` |
| `--variant-nodes <tsv>` | restrict both substrates to `call`'s variant nodes (the genotyping scope; [masking detail](../algorithms/describe.md#node--edge-dosage-complementary-substrate)) | — |
| `--variant-flank-bp <N>` | with `--variant-nodes`, widen the scope by N bp. Base-granular for k-mers — exactly N bases at the node end facing the variant — and node-granular for graph dosage, since a node dosage is a property of the whole node. It therefore selects more nodes than bases; that asymmetry is deliberate. Requires `--variant-nodes` | `k-1` |
| `--samples <cosigt.tsv>` | also write per-sample BIMBAM (diploid summed dosage) | — |
| `--variant-vcf <vcf>` | also emit a VARIANT-level BIMBAM — one row per VCF record, or per ALT for a multiallelic one — from `call`'s region VCF (uncompressed; `.vcf.gz` is refused). A row is not necessarily one independent biological event: a paired DEL+INS representation, or correlated CN and sequence records, produce several rows for one event. It is a coarser and more interpretable unit than a k-mer, not an independent one; correlation is `associate`'s and `Meff`'s problem | — |
| `--only-kmers` / `--only-graph` / `--only-variant` | emit exactly one substrate; passing two is an error rather than an intersection that selects nothing (`--only-variant` needs `--variant-vcf` and no GFA — Graphical Fragment Assembly) | k-mer + graph |
| `--scale-dosage` | rescale each feature's dosage to the `0..2` range (per-feature min-max) | off (raw counts) |
| `--no-bimbam` | skip the BIMBAM dosage matrices + their `feature_annot` sidecars | off |
| `--no-wide-matrix` | write only the feature map and the sparse JSONL (JSON Lines; skip the dense per-bubble matrix) | off |
| `--max-wide-features <N>` / `--force-wide` | wide-matrix cap / override | `250000` |
| `--threads <N>` / `--quiet` | workers / quiet | `0` / off |

## Outputs

Each BIMBAM substrate is written to its own self-contained folder, `<out-dir>/<level>/<substrate>/`, where `<level>` is `haplotype` (per-haplotype, always) or `sample` (per-sample diploid, only with `--samples`) and `<substrate>` is `kmers`, `graph`, or `variant`. A folder holds three files: the matrix, its sidecar, and its column order. So an `associate` run points at one folder and takes all three from it.

| file (in each `<level>/<substrate>/`) | contents |
|------|----------|
| `bimbam_<substrate>.bimbam.gz` | BIMBAM mean-genotype dosage (`id, A, B, dose…`; missing = `NA`); the `associate`/GEMMA input |
| `feature_annot.<substrate>.tsv.gz` | `feature_id, layer, encoding, bubbles, nodes` sidecar for the BIMBAM rows (variant adds `svtype, gene, AF, AN`) |
| `samples.txt.gz` | column order for that matrix — haplotype names under `haplotype/`, sample names under `sample/` |

`kmers`/`graph` come from the graph; `variant` needs `--variant-vcf` and carries the variant-level dosage (`CN` for `DUP`, allele indicator for multiallelic, `GT` 0/1 otherwise). The `sample/` level appears only with `--samples`, summing each sample's haplotype dosages into a diploid genotype. Alongside the substrate folders, at the top of `<out-dir>`:

| file | contents |
|------|----------|
| `describe.index.tsv` | per-bubble kept/candidates/discarded for both substrates |
| `describe.params.json` | the run's resolved parameters (k, feature mode, filter, wide-matrix flags) |
| `bubble_<id>/graph_features.tsv.gz` | node and edge dosage map (`feature_type` = node/edge) |
| `bubble_<id>/kmer_features.tsv.gz` + `kmer_counts.jsonl.gz` | k-mer map (with `nodes` provenance) + per-path sparse counts |
| `bubble_<id>/{kmer,graph}_matrix.tsv.gz` | dense feature × path matrices — only without `--no-wide-matrix` (the gene drivers pass `--no-wide-matrix`, so these are normally absent; the feature maps + JSONL carry the same information) |

By default dosages are raw counts (not rescaled to 0–2), so a haplotype carrying 50 copies shows 50; `NA` = a haplotype that doesn't traverse the feature's bubble (distinct from `0` = traverses but reference). A path taking the direct source-to-sink edge, that is a pure deletion, is a traverser and reads `0`, not `NA`.

Missingness is all-or-nothing on purpose: a pooled feature is finite only when every bubble contributing it is observable on that path, and a diploid sample only when every assigned haplotype traverses. Taking "any" instead reported a partial sum as though it were whole, biasing dosage downward on the substrate `associate` tests. `panvar associate` tests these raw counts directly. `--scale-dosage` maps each feature to the `0..2` range instead — a per-feature linear map, so it doesn't change `associate`'s linear-model result; it's there for external tools (e.g. GEMMA) that assume a diploid `0..2` dosage and would otherwise drop copy-number markers.

Each `feature_annot.<substrate>.tsv.gz` (one row per BIMBAM feature, in row order) carries the provenance of each genotype:

| column | meaning |
|--------|---------|
| `feature_id` | the k-mer sequence (k-mer substrate) or the node id/oriented-edge key (graph substrate). Node and edge ids share one namespace; `encoding` is what distinguishes them |
| `layer` | `kmer` or `graph` |
| `encoding` | how the feature is built: `syncmer`/`all` for k-mers, `node`/`edge` for graph dosage, `dosage` for variant rows |
| `bubbles` | the bubble id(s) the feature belongs to |
| `nodes` | the graph node(s) the feature localizes to — the traceback into `call`'s `variant_nodes.tsv` |

The BIMBAM matrices themselves hold `feature_id`, two allele-label columns `A`,`B`, then one dosage column per haplotype (or per sample, with `--samples`) in the folder's `samples.txt.gz` order. The per-bubble `*_features.tsv.gz` maps name each feature to its nodes/encoding, and `kmer_counts.jsonl.gz` holds the per-path sparse counts (`[feature_id, count]` tuples); the dense `*_matrix.tsv.gz` is the same data materialised as a feature × path table and is written only when the wide matrix is enabled (i.e. not under `--no-wide-matrix`).

The BIMBAM exports feed [`panvar associate`](associate.md) directly (it tests the dosage, so copy-number loci are first-class). A significant marker traces back through its `nodes`/`bubbles` to `call`'s `variant_nodes.tsv` and from there to the variant `call` typed — and, with `call --gtf` plus `associate --node-genes`, on to the gene it sits in. A worked end-to-end run with the concepts in context is in [GWAS example](../gwas.md).

## Limitations

- The discriminative filter drops features that do not vary across the panel, so a feature carried by every haplotype does not reach the matrix even if it is biologically interesting.
- Missingness is all-or-nothing: a feature spanning several bubbles is finite only when every contributing bubble is observable on that haplotype, and a diploid sample only when every assigned haplotype traverses. Taking any-instead-of-all would report a partial sum as though it were whole.
- Dosages are raw counts by default rather than a diploid 0-to-2 scale, so an external tool assuming the latter needs `--scale-dosage`.
- `--variant-flank-bp` admits whole nodes for graph dosage and individual bases for k-mers, so the two substrates see slightly different neighbourhoods of the same variant.
- Node and edge features share one identifier namespace. A node id containing the oriented-edge separator is refused rather than renamed, since renaming would break the join back to `call`'s node sets.
- Carrier maps for all bubbles are held in memory at once; `--max-wide-features` bounds the dense matrices, not this.

## Example

See the [walkthrough](../walkthrough.md) for this module in a full end-to-end run.
