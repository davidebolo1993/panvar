# Module `inspect`

CLI: `panvar inspect`

## What it does

Shows, site by site, what every haplotype does there: the sequence it spells, the nodes it walks, and how many times it passes each node and edge. Use it to check a site before calling it, or to draw it as a heatmap.

- Works on the sorted graph and site table of `bubble`, `panphorte` or `refine`.
- `--cluster` groups crossings whose walks are alike.
- `--genotype-table` writes an allele table: a row per haplotype (or per sample, with `--samples`), a column per distinct allele.

## How it works

1. A crossing is one pass of a haplotype through the site, from one boundary to the other through the interior. A haplotype crossing twice (the two copies of a duplication) gets two rows, `<path>@1` and `<path>@2`; what it walks between them is not counted.
2. Each crossing is read from source to sink: its sequence, its walk, and its passes through each node (forward and reverse) and each link. A folded repeat walked once per copy shows its copy number here.
3. With `--genotype-table`, a haplotype's allele is the sequence of its crossings, so one base makes a new allele. Alleles are numbered by how many haplotypes carry them; with `-r` they are written on the reference's strand.
4. With `--cluster`, crossings with the same walk are joined (a haplotype crossing twice is two members, which can land in different groups), and distinct walks are compared by the bases they share between the boundaries ([Dice](https://doi.org/10.2307/1932409) similarity weighted by node length: twice the shared bases over the bases of both walks). Each walk, most carried first, joins the group whose representative it is most similar to, at `--cluster-similarity` or above, or starts a new group. A single-base change costs one base and an extra copy the whole copy, so groups follow structure more than SNPs.

## Required inputs

- `-i, --gfa <graph.gfa>`: the sorted GFA from `bubble`, `panphorte` or `refine`
- `-b, --bubble-prefix <prefix>` (reads `<prefix>.bubbles.csv`) or `-c, --bubbles-csv <path>`: the sites

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `-o, --out-prefix <p>` | the run's folder; a trailing `/` is dropped, and `.` takes the folder's own name | `inspect` |
| `--bubble-id <N>` | only this site; repeatable | all |
| `--cluster` | group crossings with similar walks | off |
| `--cluster-similarity <f>` | share of bases a walk needs in common with its group's representative; above 1 a percentage | `0.90` |
| `--genotype-table` | write the allele table, the alleles and a node dictionary | off |
| `--samples <cosigt.tsv>` | rows are samples, from a [cosigt](https://github.com/davidebolo1993/cosigt) `sample <tab> hap1 <tab> hap2` table; implies `--genotype-table` | — |
| `-r, --reference-path <name>` | with `--genotype-table`: alleles on this path's strand, node positions on it | — |
| `--fasta-out`, `--table-out`, `--edge-table-out <path>` | name an output file, for one `--bubble-id` | from `-o` |
| `-q, --quiet` | no progress logs | off |

## Outputs

Per site `N`, files are in `<prefix>/bubble_<N>/`, named `<name>.bubble_<N>.<file>` (`<name>` the folder's name):

| file | contents |
|------|----------|
| `paths.fa.gz`, `walks.fa.gz` | each crossing's sequence, and its node walk (`1+,7-,8+`) |
| `node_counts.tsv` | a row per crossing (`path_name`, `path_length_bp`, `crossing`, `crossings`), a `node.<id>` column per interior node: `total:forward:reverse` passes |
| `edge_counts.tsv` | the same rows, an `edge.<from±>to±>` column per link walked: passes |
| `node_lengths.tsv` | `node_id`, `length_bp`, in the columns' order |
| `clusters.tsv` | with `--cluster`: `cluster_id`, `n_paths` (the crossings in the group), `representative_path`, `members` (`;`-separated row labels) |
| `genotypes.tsv.gz` | with `--genotype-table`: `haplotype` (or `sample`), then an `alleleK` column per allele: 0/1 per haplotype, carriers per sample; `NA` with no allele |
| `alleles.fa.gz`, `alleles.walks.fa.gz` | with `--genotype-table`: each allele's sequence and walks, with its carriers |
| `nodes.json.gz` | with `--genotype-table`: every node on a crossing, with its sequence, length, carriers, rank and, with `-r`, its position |

A haplotype is `NA` at a site when it has no crossing there: its assembly starts or ends inside the site, it skips both boundaries, or it only touches them from outside. A path going straight from one boundary to the other is not `NA`: it carries a deletion. A sample is `NA` when one of its haplotypes is.

With `--genotype-table`, `<prefix>/<name>.sites.tsv` has a row per site: `bubble_id`, `source`, `sink`, `haplotypes`, `crossing`, `na`, `touching_only`, `reference_crosses`, `alleles`, `samples` and `samples_na` (`NA` without `--samples`).

Columns are in graph order, so on a sorted graph a heatmap reads along the reference. Files an earlier run left at the prefix and this run did not write are removed, so they are not read as its own: those in its folder, site files beside it, flat or in a folder (`<prefix>.bubble_<N>.<file>`, `<prefix>.bubble_<N>/`), and a `sites.tsv` in the folder or beside it. A `--bubble-id` run removes only its own sites' files, since the other sites' files are not stale, and any `sites.tsv` it does not write, since that table would not match this run's files.

Path names are refused when they would break an output: whitespace always (a FASTA id ends at the first space), a comma with `--genotype-table` (the allele FASTA lists its carriers comma-separated) and a `;` with `--cluster` (the members are `;`-separated).

### Visualization

`scripts/plot_inspect.R` draws a site as a heatmap: a row per crossing, a column per node (`--features nodes`, the default) or link (`--features edges`), columns wider for longer nodes, and rows grouped by cluster after `--cluster` (`--representatives-only` keeps one row per group). It lists its options under `--help`:

```
Rscript scripts/plot_inspect.R --table <prefix>/bubble_<N>/
Rscript scripts/plot_inspect.R --table <prefix>/bubble_<N>/ --features edges
```

## Limitations

- Alleles are distinct by exact sequence: a hypervariable site has close to one allele per haplotype. `--cluster` groups them by structure instead.
- Clusters follow bases: at `0.90`, two forms differing by under a tenth of their length share a group; a higher `--cluster-similarity` separates them.
- An inversion site is shown on its own nodes, without the sites inside it.

## Example

```bash
# every site: each haplotype's allele, and its node and link counts
panvar inspect -i locus_r.normalized.sorted.gfa -b locus_r -o locus_inspect
# one site, crossings clustered, and the allele table by sample
panvar inspect -i locus_r.normalized.sorted.gfa -b locus_r -o locus_inspect --bubble-id 3 --cluster \
  --samples samples.tsv -r GRCh38
```

`scripts/inspect_examples.sh` inspects the copy-number site of four loci from `tests/examples/bubble`, with clusters, allele tables and heatmaps, one folder per locus in `tests/examples/inspect`.
