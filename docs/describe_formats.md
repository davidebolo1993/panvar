# `describe` output formats

A reader's guide to what `panvar describe` writes and how to read it. For what the module is for and
its full option list see [modules/describe.md](modules/describe.md); for how the features are built
see [algorithms/describe.md](algorithms/describe.md).

Every example below is synthetic. The shapes, column names and separators are real.

## The three substrates

`describe` turns called bubbles into per-haplotype genotype features on three substrates. They share
the same graph coordinates, so a hit on one traces back to the others.

| substrate | one feature is | its dosage means |
|-----------|----------------|------------------|
| `kmers` | a canonical k-mer, closed-syncmer sampled by default | how many times it occurs on that haplotype's walk |
| `graph` | a node, or an oriented edge | traversal count: node is presence/abundance, edge is adjacency-aware and the better tandem signal |
| `variant` | one VCF record, or one ALT of a multiallelic record | `CN` for a `DUP`, otherwise presence from `GT` |

`kmers` and `graph` are built from the graph and emitted by default. `variant` needs
`--variant-vcf`. Only `kmers` and `graph` pass the discriminative filter that drops features not
varying across haplotypes; `variant` emits every call and leaves frequency filtering to `associate`.

## What part of the graph the features cover

By default the features cover the whole bubble: every node, edge and k-mer on the haplotype walks
through it. Passing `--variant-nodes <call.variant_nodes.tsv>` restricts the `kmers` and `graph`
substrates to the nodes `call` actually typed, widened by `--variant-flank-bp` (default `k-1`).
Without the flag nothing is masked.

The flank has two granularities on purpose. For k-mers it admits exactly N bases at the neighbouring
node's facing end. For node and edge dosage it admits the whole node, because a node dosage is a
property of the entire node and there is no partial-node count to give. The same flag therefore
selects more nodes than bases.

How much this changes depends on the locus. Where most bubble-interior nodes are already variant
nodes the two runs are nearly identical; where a bubble carries a lot of uncalled interior, the
restriction is what keeps the matrix to the called variation.

## Directory layout

```
<out-dir>/
  describe.index.tsv                 per-bubble kept/candidate/discarded counts
  describe.params.json               the run's resolved parameters
  bubble_<id>/
    kmer_features.tsv.gz             the k-mer map for this bubble
    kmer_counts.jsonl.gz             per-path sparse counts
    kmer_matrix.tsv.gz               the same, dense (omitted under --no-wide-matrix)
    graph_features.tsv.gz            the node/edge map for this bubble
    graph_matrix.tsv.gz              the same, dense (omitted under --no-wide-matrix)
  haplotype/<substrate>/             one folder per substrate: kmers, graph, variant
    bimbam_<substrate>.bimbam.gz
    feature_annot.<substrate>.tsv.gz
    samples.txt.gz
  sample/<substrate>/                same three files, only with --samples
```

The `haplotype` level is always written. The `sample` level appears only with `--samples` and holds
diploid values: each sample's dosage is the sum over its assigned haplotypes.

## The BIMBAM matrices

`bimbam_<substrate>.bimbam.gz` is the file `associate` and GEMMA consume. It is comma separated:
the feature id, two allele-label columns `A` and `B` that exist to satisfy the format, then one
dosage per column in the order given by `samples.txt.gz` in the same folder.

```
GCTAAAGACAATTACATAACATACACGTCAG, A, B, 0, 1, 1, 0
CACGAAACTTGTTGGCCCAGTGTGAATCGCT, A, B, 2, 2, NA, 2
```

The feature id is the feature itself, and its form tells you the substrate:

```
TAAGGGTTAAGTAAGTGTGATGCATACGCCT, A, B, 1, 0, 1, 1     a k-mer: the sequence is the id
4021, A, B, 2, 1, 2, 2                                 a graph node
4021+>4023+, A, B, 2, 1, 2, 2                          a graph edge, oriented
bubble3_DEL_512, A, B, 0, 1, 0, 0                      a variant: the VCF record id
```

Nodes and edges share one identifier namespace. The `encoding` column of the sidecar is what
separates them, so do not try to tell them apart by parsing the id.

## The sidecar

`feature_annot.<substrate>.tsv.gz` has one row per BIMBAM row, in the same order, and carries the
provenance. This is the traceback into `call`'s `variant_nodes.tsv`.

```
feature_id                       layer  encoding  bubbles  nodes
GCTAAAGACAATTACATAACATACACGTCAG  kmer   syncmer   3        811;812;813;814
4021                             graph  node      7        4021
4021+>4023+                      graph  edge      7        4021+>4023+
```

The `variant` sidecar adds four columns:

```
feature_id       layer    encoding  bubbles  nodes            svtype  gene   AF        AN
bubble3_DEL_512  variant  dosage    3        512,514,515,517  DEL     .      0.030534  131
```

`samples.txt.gz` is one name per line: haplotype names under `haplotype/`, sample names under
`sample/`.

```
HG0AAAA#1#haplotype1-0000024:31848049-32074186
HG0AAAA#2#haplotype2-0000117:31856177-32049825
```

## The per-bubble tables

These are the working tables behind the matrices, one set per bubble. They carry count summaries the
BIMBAM does not.

`kmer_features.tsv.gz` maps each k-mer to its numeric id, its column name in the dense matrix, and
the nodes it localizes to:

```
feature_id  feature_name  encoded_kmer  kmer                             paths_present  min_count  max_count  total_count  node_count  nodes
1           K1            594           GCTAAAGACAATTACATAACATACACGTCAG  2              0          1          2            4           811;812;813;814
2           K2            747           CACGAAACTTGTTGGCCCAGTGTGAATCGCT  9              0          2          11           3           822;823;824
```

`graph_features.tsv.gz` is the same idea for nodes and edges, where `label` is the graph's own node
id and `feature_type` separates the two kinds:

```
feature_id  feature_name  feature_type  label        paths_present  min_count  max_count  total_count
1           N1            node          4021         131            1          2          201
2           N2            node          4023         3              0          1          3
3           E1            edge          4021+>4023+  128            1          2          194
```

`kmer_counts.jsonl.gz` is one JSON object per path per bubble, holding the counts sparsely. The
tuples are `[feature_id, count]`, where `feature_id` is the numeric id from `kmer_features.tsv.gz`,
not the sequence:

```json
{"bubble_id":3,"sample":"HG0AAAA","haplotype":"1","path_name":"HG0AAAA#1#haplotype1-0000024:31848049-32074186","path_length_bp":7941,"kmers":[[5,1],[7,2],[12,1]]}
```

`kmer_matrix.tsv.gz` and `graph_matrix.tsv.gz` are the same data materialised densely, one row per
path, with the `feature_name` values as columns:

```
bubble_id  sample   haplotype  path_name                                        K1  K2  K3
3          HG0AAAA  1          HG0AAAA#1#haplotype1-0000024:31848049-32074186   0   1   2
```

They are written only when the wide matrix is enabled, which is the default, and are capped by
`--max-wide-features`. Skipping them with `--no-wide-matrix` loses nothing: the feature map plus the
sparse JSONL carry the same information.

## The index and the parameters

`describe.index.tsv` is one row per bubble with candidate, kept and discarded counts for k-mers,
nodes and edges separately, the path to each file it wrote, and whether the dense matrix was
written or skipped and why. It is the first place to look when a bubble produced fewer features than
expected.

`describe.params.json` records the resolved run: `kmer_size`, `feature_mode`, `syncmer_s`,
`min_feature_paths`, the emit and matrix flags, `variant_nodes` and `variant_flank_bp`, and the
input sizes. It is what makes a matrix reproducible months later.

## Reading the values

`NA` is not `0`. `NA` means the haplotype does not traverse that feature's bubble at all, so the
feature is unobservable there. `0` means it traverses and does not carry the feature. A path taking
the direct source-to-sink edge, that is a pure deletion, is a traverser and reads `0`.

Missingness is all-or-nothing. A feature pooled across several bubbles is finite only when every
contributing bubble is observable on that path, and a sample only when every assigned haplotype
traverses. This is deliberate: taking any-instead-of-all would report a partial sum as though it
were a complete low dosage.

Dosages are raw counts by default, not a diploid 0-to-2 scale, so a haplotype carrying 50 copies
reads `50`. `associate` tests raw counts directly. `--scale-dosage` applies a per-feature min-max
rescale to `0..2` for external tools that assume diploid dosage.
