# Module `describe`

CLI: `panvar describe`

## What it does

Converts the called bubbles into per-haplotype genotype features for association, on three substrates that share the same graph coordinates and are each exported as a BIMBAM mean-genotype dosage matrix — the canonical input to [`associate`](associate.md) and tools like [GEMMA](https://github.com/genetics-statistics/GEMMA):

- **variant** — one dosage per structural-variant (SV) call from `call`'s VCF (Variant Call Format) (`--variant-vcf`): copy number for a `DUP`, presence for `DEL` / `INS` / `INV`
- **k-mers** — canonical k-mers per path, [closed-syncmer](../algorithms/describe.md#terms) sampled by default, each carrying its node provenance;
- **graph** — per-node and per-oriented-edge traversal counts, a graph-local dosage;

The k-mer and graph substrates are built from the graph and emitted by default; the variant substrate is emitted when its VCF is supplied. Any one can be produced on its own (`--only-kmers` / `--only-graph` / `--only-variant`). The k-mer and graph substrates pass through a discriminative filter that drops features which do not vary across haplotypes — those present in all (or all but `--min-paths`) paths, unless their copy number varies — so only informative markers reach the matrix; the variant substrate emits every call and leaves frequency filtering to `associate`.
Pass a [cosigt](https://github.com/davidebolo1993/cosigt) table (`sample <tab> hap1 <tab> hap2 …`) to `--samples <cosigt.tsv>` and `describe` also writes a per-sample BIMBAM for each substrate (`bimbam_{kmers,graph,variant}.samples.bimbam.gz`) whose value sums the dosage over the sample's haplotypes (diploid `CN` = `CN_A` + `CN_B`). See [gwas/example.md](../gwas/example.md).

Algorithm and worked trace: [algorithms/describe.md](../algorithms/describe.md).



## Required inputs

`describe` must run on the same graph and bubble prefix as `call` so feature node ids line up with the VCF: the `panphorte` `.normalized.sorted.gfa` or the `bubble` `.sorted.gfa`.

- `-i, --gfa <path>`; 
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
| `--variant-flank-bp <N>` | with `--variant-nodes`, also keep nodes within N bp of a variant node | `k-1` |
| `--samples <cosigt.tsv>` | also write per-sample BIMBAM (diploid summed dosage) | — |
| `--variant-vcf <vcf>` | also emit a VARIANT-level BIMBAM (one dosage row per SV call) from `call`'s region VCF — the honest GWAS unit for `associate --unit variant` | — |
| `--only-kmers` / `--only-graph` / `--only-variant` | emit only one substrate (`--only-variant` needs `--variant-vcf` and no GFA — Graphical Fragment Assembly) | k-mer + graph |
| `--scale-dosage` | rescale each feature's dosage to the `0..2` range (per-feature min-max) | off (raw counts) |
| `--no-bimbam` | skip the BIMBAM dosage matrices + `feature_annot.tsv.gz` | off |
| `--no-wide-matrix` | write only the feature map and the sparse JSONL (JSON Lines; skip the dense per-bubble matrix) | off |
| `--max-wide-features <N>` / `--force-wide` | wide-matrix cap / override | `250000` |
| `--threads <N>` / `--quiet` | workers / quiet | `0` / off |

## Outputs

| file | contents |
|------|----------|
| `bimbam_{kmers,graph}.bimbam.gz` | canonical BIMBAM mean-genotype dosage (`id, A, B, dose…`; missing = `NA`); the `associate`/GEMMA input |
| `feature_annot.tsv.gz` | `feature_id, layer, encoding, bubbles, nodes` sidecar for the BIMBAM rows |
| `bimbam.samples.txt.gz` | sample (column) order shared by both BIMBAM files |
| `bimbam_{kmers,graph}.samples.bimbam.gz` + `bimbam.samples.samples.txt.gz` + `feature_annot.samples.tsv.gz` | with `--samples`: per-sample versions |
| `bimbam_variant.bimbam.gz` + `bimbam_variant.samples.txt.gz` + `feature_annot.variant.tsv.gz` | with `--variant-vcf`: the variant-level dosage (`CN` for `DUP`, allele indicator for multiallelic, `GT` 0/1 otherwise), its column order, and a sidecar adding `svtype, gene, AF, AN`. With `--samples`, also the diploid `bimbam_variant.samples.bimbam.gz` (+ `.samples.samples.txt.gz`) |
| `describe.index.tsv` | per-bubble kept/candidates/discarded for both substrates |
| `describe.params.json` | the run's resolved parameters (k, feature mode, filter, wide-matrix flags) |
| `bubble_<id>/graph_features.tsv.gz` | node and edge dosage map (`feature_type` = node/edge) |
| `bubble_<id>/kmer_features.tsv.gz` + `kmer_counts.jsonl.gz` | k-mer map (with `nodes` provenance) + per-path sparse counts |
| `bubble_<id>/{kmer,graph}_matrix.tsv.gz` | dense feature × path matrices — only without `--no-wide-matrix` (the gene drivers pass `--no-wide-matrix`, so these are normally absent; the feature maps + JSONL carry the same information) |

By default dosages are raw counts (not rescaled to 0–2), so a haplotype carrying 50 copies shows 50; `NA` = a haplotype that doesn't traverse the feature's bubble (distinct from `0` = traverses but reference). `panvar associate` tests these raw counts directly. `--scale-dosage` maps each feature to the `0..2` range instead — a per-feature linear map, so it doesn't change `associate`'s linear-model result; it's there for external tools (e.g. GEMMA's minor-allele-frequency (MAF) model) that assume a diploid `0..2` dosage and would otherwise drop copy-number markers.

`feature_annot.tsv.gz` (one row per BIMBAM feature, in row order) carries the provenance of each genotype:

| column | meaning |
|--------|---------|
| `feature_id` | the k-mer sequence (k-mer substrate) or the node id / oriented-edge key (graph substrate) |
| `layer` | `kmer` or `graph` |
| `encoding` | how the feature is built (`syncmer`/`all` for k-mers; `node`/`edge` for graph dosage) |
| `bubbles` | the bubble id(s) the feature belongs to |
| `nodes` | the graph node(s) the feature localizes to — the traceback into `call`'s `variant_nodes.tsv` |

The BIMBAM matrices themselves hold `feature_id`, two allele-label columns `A`,`B`, then one dosage column per haplotype (or per sample, with `--samples`) in `bimbam.samples.txt.gz` order. The per-bubble `*_features.tsv.gz` maps name each feature to its nodes/encoding, and `kmer_counts.jsonl.gz` holds the per-path sparse counts (`[feature_id, count]` tuples); the dense `*_matrix.tsv.gz` is the same data materialised as a feature × path table and is written only when the wide matrix is enabled (i.e. not under `--no-wide-matrix`).

## Association

The BIMBAM exports feed [`panvar associate`](associate.md) directly (it tests the dosage, so copy-number loci are first-class). A significant marker traces back through its `nodes`/`bubbles` to `call`'s `variant_nodes.tsv` and from there to the variant `call` typed — and, with `call --gtf` plus `associate --node-genes`, on to the gene it sits in. A worked end-to-end run with the concepts in context is in [gwas/example.md](../gwas/example.md).

## Example

See the [LPA walkthrough](../walkthrough.md) for this module in a full end-to-end run.
