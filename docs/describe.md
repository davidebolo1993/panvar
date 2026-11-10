# Module `describe`

CLI: `panvar describe`

## What it does

Turns each haplotype's walks through the sites into genotype features for association, and writes them as BIMBAM dosage matrices, the input of `associate`. There are three substrates:

- k-mers: canonical k-mers spelled from each walk, closed syncmers by default, each traced to the nodes it lies on;
- graph: how often each walk visits each node and each oriented edge of the site;
- variant: one row per record of `call`'s VCF (`--variant-vcf`)


The k-mer and graph substrates keep only features that vary across the haplotypes (`--min-paths`), and at the sample level across the samples. The variant substrate keeps every record and leaves frequency filtering to `associate`. With a [cosigt](https://github.com/davidebolo1993/cosigt) table (`--samples`), each substrate is also written per sample, the sum over the sample's haplotypes.

## How it works

Closed syncmers are those from the original [syncmer paper](https://doi.org/10.7717/peerj.10805); the matrices are in [BIMBAM](https://github.com/genetics-statistics/GEMMA) mean-genotype format.

### Variant rows

One row per record of `call`'s VCF, one value per haplotype:

- a sequence record (`DEL`, `INS`, `INV`, `SUB`) gives its `GT`: 1 for a carrier, 0 otherwise;
- a `CNV` record gives its `CN`, so a haplotype like the reference reads its copy number, not 0;
- a record with several ALTs gives one row per ALT (`<id>_a1`, `<id>_a2`, ...);
- `.` reads `NA`.

Every record is written; frequency filtering is left to `associate`.

```text
call's VCF
record  SVTYPE  FORMAT  h1   h2   h3   h4
V1      CNV     GT:CN   1:3  0:2  1:4  .:.
V2      DEL     GT      1    0    1    0
V3      2 ALTs  GT      1    2    0    .

variant matrix
row     h1  h2  h3  h4
V1      3   2   4   NA
V2      1   0   1   0
V3_a1   1   0   0   NA
V3_a2   0   1   0   NA
```

### k-mer and graph features

1. Each site is read crossing by crossing, from source to sink, as `inspect` reads it. A path crossing twice gives both crossings; a deletion crosses like any other path and reads 0 at the nodes it skips. A node two sites share is counted once per visit, and nested sites are refused, so no visit counts twice.
2. k-mers: each walk is scanned for canonical k-mers, by default only closed syncmers, which sample the sequence evenly so a substitution changes only the few that overlap it. Each is counted per path over all its crossings: three copies read 3.
3. Graph: each path's count of every node and oriented edge. Edges carry adjacency, so a duplication or an inversion shows as edges no other path has.
4. A feature is kept when its count varies among the paths carrying it, or when more than `--min-paths` paths carry it and more than `--min-paths` lack it; otherwise there is nothing to contrast.
5. With `--variant-nodes`, only sites with a record are read, and only the records' nodes, the site's boundaries and a flank (`--variant-flank-bp`), so the features describe the variants `call` found.
6. A feature kept at several sites, typically a k-mer, becomes one row, the sum of its counts. A path is observed in it only where it crosses every one of those sites, since a partial sum would read as a low dosage.

With `--samples`, every row is also summed over each sample's haplotypes; a sample is observed only where each of its haplotypes is.

## Required inputs

- `-i, --gfa <path>`: the sorted GFA `call` ran on, so node ids match its outputs. Not needed with `--only-variant` alone.
- `-b, --bubble-prefix <prefix>` (reads `<prefix>.bubbles.csv`) or `-c, --bubbles-csv <path>`: the site table `call` read, with its boundary orientation columns.
- `--variant-vcf <vcf>`: `call`'s region VCF, uncompressed and per haplotype (not `--diploid`: the per-sample level comes from `--samples`), for the variant substrate.

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `-o, --out-dir <dir>` | output directory | `describe_out` |
| `-r, --reference-path <name>` | the reference path, exact name or unique part: k-mer and graph rows get its coordinates (`chrom`, `pos`, `end`) | — |
| `--bubble-id <N>` | only this site; repeatable. Restricts every substrate | all |
| `--variant-nodes <tsv>` | `call`'s `variant_nodes.tsv`: only sites with a record are described, and of them only the records' nodes and the site boundaries give k-mer and graph features. The table and `--variant-vcf` must come from one `call` run. With `--only-variant`, it gives only the variant rows' nodes | — |
| `--variant-flank-bp <N>` | with `--variant-nodes`, also the N bases along the walk on either side of a record's nodes, or of a junction that skips them, for k-mers, and every node those bases reach for graph features, kept on every walk. With 0, no k-mer and no edge crosses a deletion's junction, except one between the boundaries | `k-1` |
| `--samples <tsv>` | cosigt table (`sample`, then one column per haplotype; a haplotype in several contigs lists them comma-separated): also write the sample level | — |
| `--only-kmers`, `--only-graph`, `--only-variant` | only these substrates; the flags combine. `--only-variant` needs `--variant-vcf` and no GFA | all three, the variant one with `--variant-vcf` |
| `-k, --kmer-size <K>` | k-mer size, 1 to 31; 2 or more with syncmers | `31` |
| `--feature-mode <syncmer\|all>` | closed syncmers only, or every k-mer | `syncmer` |
| `--syncmer-s <S>` | s-mer size of the closed syncmers | `(k+2)/3`, at most 11 |
| `--min-paths <N>` | drop a feature present in, or absent from, N paths or fewer, unless its count varies; 0 drops only features constant across the paths | `1` |
| `--scale-dosage` | rescale each BIMBAM feature to 0..2 (per-feature min-max), for tools that expect a diploid dosage; a constant feature reads 0 | off |
| `--no-bimbam` | no BIMBAM matrices, only the per-site files | off |
| `--no-wide-matrix` | per site, only the feature map and the sparse counts | off |
| `--max-wide-features <N>`, `--force-wide` | skip a site's dense matrix above N features (0 = no cap), or write it anyway | `250000` |
| `-t, --threads <N>`, `-q` | worker threads (`0` = all cores); no progress logs | `0`, off |

## Outputs

Each BIMBAM matrix sits in its own folder, `<out-dir>/<level>/<substrate>/`, with `<level>` either `haplotype` or `sample` (with `--samples`) and `<substrate>` one of `kmers`, `graph`, `variant`:

| file | contents |
|------|----------|
| `bimbam_<substrate>.bimbam.gz` | one row per feature: `id, A, B`, then one dosage per column, `NA` when missing |
| `feature_annot.<substrate>.tsv.gz` | one row per matrix row, in the same order (columns below) |
| `samples.txt.gz` | the matrix's columns: haplotype paths, or samples |

At the top of `<out-dir>`, `describe.index.tsv` lists the sites: per site its counts, its files in `bubble_<id>/`, and why a site or its dense matrix was skipped (`status`, `matrix_reason`).

Each `bubble_<id>/` holds the site's own features, before pooling:

| file | contents |
|------|----------|
| `kmer_features.tsv.gz` | the k-mers kept: `feature_id`, `feature_name` (`K1`, `K2`, ...), `encoded_kmer`, `kmer`, `paths_present`, `min_count`, `max_count`, `total_count`, `node_count`, `nodes` |
| `kmer_counts.jsonl.gz` | per path, its k-mer counts as `[feature_id, count]` pairs |
| `graph_features.tsv.gz` | the nodes and edges kept: `feature_id`, `feature_name` (`N1`..., `E1`...), `feature_type` (`node` or `edge`), `label` (node id, or `from±>to±`), `paths_present`, `min_count`, `max_count`, `total_count` |
| `kmer_matrix.tsv.gz`, `graph_matrix.tsv.gz` | dense tables, one row per path and one column per feature; not written with `--no-wide-matrix` or above `--max-wide-features` |

Sidecar columns (`feature_annot.<substrate>.tsv.gz`):

| column | meaning |
|--------|---------|
| `feature_id` | the k-mer's sequence, the node id or edge key, or the record ID (`<id>_a1`, `<id>_a2`... for each ALT of a multiallelic record) |
| `layer` | `kmer`, `graph` or `variant` |
| `encoding` | `syncmer` or `all` for k-mers, `node` or `edge` for graph features, `dosage` for variant rows |
| `bubbles` | the sites the feature was kept at, `;`-separated |
| `nodes` | the nodes it lies on: `;`-separated for k-mers, the id itself for graph features, and for a variant row the record's `node_ids` from `--variant-nodes` (else its `EVENT_NODES`), `,`-separated |
| `chrom`, `pos`, `end` | where the row lies on the reference, 1-based: a record's `POS` and `END` (else its last `REF` base); a k-mer, node or edge over its nodes on the reference (a node the reference walks more than once spans all its visits, so at a duplication it covers the copies), or, where it touches a node only other haplotypes walk, over its site's boundaries; `.` for k-mers and graph rows without `-r` |
| `svtype`, `gene`, `AF`, `AN` | variant rows only: the record's type, its genes (`GENES` decoded, `;`-separated), and its frequency and allele number (recomputed over the samples at the sample level) |

Node and edge ids share one namespace, so a node id containing `>` is refused.

## Limitations

- Canonical k-mers do not see an inversion between inverted repeats: across each breakpoint the k-mers are the reverse complements of the reference's. An inversion site therefore has graph features only, described by its oriented edges, and its `INV` record in the variant substrate.
- A feature carried by every haplotype in the same count is dropped, however interesting.
- A pooled k-mer is summed only over the sites where it varies, so its dosage is not always its copy count across the locus.
- Under `--variant-nodes`, the flank is counted in bases for k-mers and in whole nodes for graph features, so the two substrates see slightly different neighbourhoods.
- A variant row is not always one biological event: a replacement written as a DEL and an INS, or a copy-number record beside a sequence record, gives correlated rows. `associate` accounts for that correlation.
- Sites nested in one another are refused, except inside an inversion site. A table must carry the boundary orientation columns.

## Example

`scripts/describe_examples.sh` describes every graph `scripts/call_examples.sh` called, with its VCF, its `variant_nodes.tsv` and a samples table made from the path names, one folder per graph in `tests/examples/describe`.
