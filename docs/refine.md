# Module `refine`

CLI: `panvar refine`

## What it does

Realigns the haplotypes inside each site, so a variant the graph drew in pieces becomes one record. A graph builder can split one indel over several nodes, or spell the same bases through different nodes in different haplotypes; `call` would then report pieces, or a difference where there is none.

- Every haplotype still spells exactly what it did.
- Folded repeats from `panphorte` are kept as they are, so copy numbers never change.
- The graph is then sorted and its sites found again, ready for `call`.

## How it works

The alignment is [abPOA](https://github.com/yangao07/abPOA)'s partial-order alignment.

1. The folded repeats are read from panphorte's unit-graph table (and any node linked to itself), and kept out of every change.
2. Sites that touch (one's sink is the next one's source) are joined into one region, so an indel split across neighbouring sites comes back in one piece.
3. A region is skipped when realigning it could lose something: a haplotype walks it without crossing it (a contig starting inside), the reference does not cross it, a crossing walks a node twice outside a folded repeat (an unfolded copy, whose copy number would become an indel), or crossings walk different folded repeats from the reference.
4. In each region, the distinct haplotype sequences between folded repeats are aligned with abPOA, reference first, and the alignment becomes nodes: a run of columns shared by the same sequences is one node. Identical haplotypes end on identical routes, and each keeps its own bases.
5. Twin nodes, the same sequence entered from and leaving to the same neighbours, are merged across the graph, so haplotypes differing only in which twin they walk are no longer different.
6. Every haplotype must spell what it spelled, over links the graph has, or nothing is written.
7. The graph is sorted and its sites found again as `bubble` does, without merging nearby sites, so a folded site keeps its own records. Site ids change, and panphorte's table is written again with the new node ids.

## Required inputs

- `-i, --gfa <graph.gfa>`: the sorted GFA from `panphorte` (a regular file, not a pipe)
- `-b, --bubble-prefix <prefix>` (reads `<prefix>.bubbles.csv`) or `-c, --bubbles-csv <path>`: the sites
- `-r, --reference-path <name>`: the reference path; a unique part of the name is enough
- `-o, --out-prefix <prefix>`: output prefix

panphorte's unit-graph table is read from the site table's prefix, or from `--panphorte-prefix`.

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `--panphorte-prefix <p>` | read `<p>.panphorte.unit_graph.tsv` | the site table's prefix |
| `--bubble-id <id[,id...]>` | only these sites; repeatable. An inversion site is refused: refine realigns the sites inside it, not the inversion | all |
| `--max-poa-bp <N>` | skip a region whose longest sequence to align is over N bp | `5000` |
| `--max-walks <N>` | skip a region with more than N distinct sequences to align | `500` |
| `--max-poa-work <N>` | skip a region whose longest × total bp is over N (`0` = no limit) | `0` |
| `--min-dup-bp <N>` | realign anyway when a crossing's repeated node visits total under N bp, turning those copies into an indel | `0` |
| `--min-bubbles <N>` | only realign regions of at least N sites | `1` |
| `--resnarl-min-variant-bp <N>`, `--merge-nearby-bp <N>`, `--no-fuse-ambiguous`, `--no-flip` | `bubble`'s options for the sites found again | `50`, `0`, off, off |
| `--gtf <path>` | genes per node of the refined graph | — |
| `-q, --quiet` | no progress logs | off |

The size limits count distinct sequences, so adding a haplotype identical to another never changes a decision.

## Outputs

| file | contents |
|------|----------|
| `<prefix>.normalized.sorted.gfa`, `<prefix>.bubbles.csv` | the refined graph, sorted, and its sites, for `call` |
| `<prefix>.panphorte.unit_graph.tsv` | panphorte's table with the refined graph's node ids; `call -b <prefix>` reads it |
| `<prefix>.refine.report.tsv` | one row per region: realigned or skipped, and why |
| `<prefix>.bandage_nodes.csv`, `<prefix>.bandage_genes.csv` | Bandage colours, and genes with `--gtf` |

`refine.report.tsv`, one row per region (numbers are `.` for a skipped region, names `.` where they do not apply):

| column | meaning |
|--------|---------|
| `region` | the region's number in this run |
| `n_bubbles` | sites joined in it |
| `source`, `sink` | its outer boundaries, as ids of refine's input |
| `decision` | `rebuilt` or `skipped` |
| `reason` | why (below) |
| `crossings` | haplotype crossings rewritten |
| `realigned_bp` | bases realigned, summed over every crossing |
| `distinct_sequences`, `aligned_bp` | the distinct sequences abPOA aligned, and their bases |
| `nodes_removed`, `nodes_added` | old interior nodes removed, new nodes made |
| `short_revisits` | crossings realigned although they walk a node twice, under `--min-dup-bp` |
| `blocking_path`, `blocking_node` | the haplotype and the node behind a skip |

Reasons:

- rebuilt: `plain`; `unit_residual` (realigned around folded repeats).
- skipped:
  - `min_bubbles`: fewer sites than `--min-bubbles`;
  - `not_a_clean_chain`: the joined sites do not reduce to two outer boundaries;
  - `overlapping_sites`: a node is also in another region;
  - `partial_paths`: a haplotype walks the region without crossing it;
  - `reference_not_traversing`: the reference does not cross it;
  - `unfolded_dup`: a crossing walks a node twice outside a folded repeat;
  - `unit_skeleton_mismatch`: a crossing walks different folded repeats from the reference;
  - `segment_too_long`, `poa_work`, `too_diverse`: over `--max-poa-bp`, `--max-poa-work` or `--max-walks`.

### Visualization

The refined graph opens in [BandageNG](https://github.com/asl/BandageNG) coloured by site, as `bubble`'s (`<prefix>.bandage_nodes.csv`), and the graph before and after compares with the scripts in [rebuild's Visualization](rebuild.md#visualization).

## Limitations

- A region with an unfolded copy is skipped for every haplotype, even if only one carries it, so refine does not reach the sites where duplications are called.
- Around a folded repeat, a region is realigned only when every crossing walks the same folded repeats as the reference; a site where some haplotypes have no copy is skipped (`panphorte --refine` realigns those flanks).
- The size limits decide per region, so a locus can be partly refined; the report says which limit fired.

## Example

`scripts/refine_examples.sh` refines every graph in `tests/examples/panphorte` and compares each with its input, in `tests/examples/refine`.
