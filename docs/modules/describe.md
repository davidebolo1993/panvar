# Describe Module (Module 4)

CLI: `panvar describe`

## What it does

Converts each called bubble into haplotype features for association on two mirrored substrates that share the
same graph coordinates. The primary, read-queryable substrate is k-mers: canonical k-mers per path,
sampled with [closed syncmers](../algorithms/describe.md#terms) by default, each carrying its node
provenance. The complementary, graph-local substrate is node/edge dosage: traversal counts per node and
per oriented edge. Features that do not discriminate haplotypes are dropped from both substrates. The pooled
cohort genotypes are exported as BIMBAM mean-genotype dosage (`bimbam_{kmers,graph}.bimbam.gz`), the
canonical export read by [`associate`](associate.md) and by GEMMA. Encoding, the discriminative filter, and a
worked trace: [algorithms/describe.md](../algorithms/describe.md).

Sample-level (`--samples <cosigt.tsv>`): a GWAS tests samples, not haplotypes. Pass a cosigt table
(`sample <tab> hap1 <tab> hap2 …`) and describe also writes per-sample BIMBAM
(`bimbam_{kmers,graph}.samples.bimbam.gz`) whose value is the summed dosage over the sample's haplotypes
(diploid CN = CN_A + CN_B). See [gwas/example.md](../gwas/example.md).

## Required inputs

`describe` runs after `panphorte` (graph = `.normalized.sorted.gfa`, bubbles = the `panphorte` prefix):

- `-i, --gfa <path>`; one of `-b, --bubble-prefix-in <prefix>` or `-c, --bubbles-csv <path>`.

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
| `--only-kmers` / `--only-graph` | emit only one substrate | both |
| `--no-bimbam` | skip the BIMBAM dosage matrices + `feature_annot.tsv.gz` | off |
| `--no-wide-matrix` | write only the feature map + sparse JSONL (skip the dense per-bubble matrix) | off |
| `--max-wide-features <N>` / `--force-wide` | wide-matrix cap / override | `250000` |
| `--threads <N>` / `--quiet` | workers / quiet | `0` / off |

## Outputs

| file | contents |
|------|----------|
| `bimbam_{kmers,graph}.bimbam.gz` | canonical BIMBAM mean-genotype dosage (`id, A, B, dose…`; missing = `NA`); the `associate`/GEMMA input |
| `feature_annot.tsv.gz` | `feature_id, layer, encoding, bubbles, nodes` sidecar for the BIMBAM rows |
| `bimbam.samples.txt.gz` | sample (column) order shared by both BIMBAM files |
| `bimbam_{kmers,graph}.samples.bimbam.gz` + `bimbam.samples.samples.txt.gz` + `feature_annot.samples.tsv.gz` | with `--samples`: per-sample versions |
| `describe.index.tsv` | per-bubble kept/candidates/discarded for both substrates |
| `bubble_<id>/graph_features.tsv.gz` + `graph_matrix.tsv.gz` | node and edge dosage map + matrix (`feature_type` = node/edge) |
| `bubble_<id>/kmer_features.tsv.gz` + `kmer_matrix.tsv.gz` + `kmer_counts.jsonl.gz` | k-mer map (with `nodes` provenance) + matrix + sparse counts |

Dosages are raw counts (not rescaled to 0–2), so a CN-50 KIV-2 haplotype shows 50; `NA` = a haplotype
that doesn't traverse the feature's bubble (distinct from `0` = traverses but reference).

`feature_annot.tsv.gz` (one row per BIMBAM feature, in row order) carries the provenance of each genotype:

| column | meaning |
|--------|---------|
| `feature_id` | the k-mer sequence (k-mer substrate) or the node id / oriented-edge key (graph substrate) |
| `layer` | `kmer` or `graph` |
| `encoding` | how the feature is built (`syncmer`/`all` for k-mers; `node`/`edge` for graph dosage) |
| `bubbles` | the bubble id(s) the feature belongs to |
| `nodes` | the graph node(s) the feature localizes to — the traceback into `call`'s `variant_nodes.tsv` |

The BIMBAM matrices themselves hold `feature_id`, two allele-label columns `A`,`B`, then one dosage column
per haplotype (or per sample, with `--samples`) in `bimbam.samples.txt.gz` order. The per-bubble
`*_features.tsv.gz` maps name each feature to its nodes/encoding; the paired `*_matrix.tsv.gz` is the dense
feature × path count table.

## Association

The BIMBAM exports feed [`panvar associate`](associate.md) directly (it tests the dosage, so copy-number
loci are first-class) and GEMMA unchanged. A significant marker traces back via `nodes`/`bubbles` →
`call`'s `variant_nodes.tsv` → the variant, and (with `call --gtf` + `associate --node-genes`) to a gene.
Worked end-to-end run: [gwas/example.md](../gwas/example.md); concepts: [gwas/primer.md](../gwas/primer.md).

## Example

Matches `scripts/genes/lpa.sh` / `tests/gwas/run_lpa_real.sh`:

```bash
./build/panvar describe \
  -i results/real_data/lpa/panphorte/panphorte.normalized.sorted.gfa \
  --bubble-prefix-in results/real_data/lpa/panphorte/panphorte \
  --out-dir results/real_data/lpa/describe \
  --kmer-size 31 --no-wide-matrix \
  --variant-nodes results/real_data/lpa/call/call.variant_nodes.tsv
# add --samples <cosigt.tsv> for per-sample (diploid) BIMBAM, as the GWAS driver does.
```

Algorithm & worked example: [algorithms/describe.md](../algorithms/describe.md). References:
[references.md](../references.md#describe).
