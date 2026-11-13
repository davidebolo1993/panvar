# Module `bubble`

CLI: `panvar bubble`

## What it does

Finds the variant sites of a locus graph, in reference order, for the later modules. A site is a stretch of the graph between two boundary nodes every haplotype passes through; the part between them is the interior, and each way through it is an allele.

- The graph is sorted along the reference first, so every site reads on the reference strand.
- Sites may hold loops and reversed nodes, where copy-number variation and inversions sit.
- One event the graph split into pieces is joined back into one site.

## How it works

Sites are [snarls](https://github.com/vgteam/vg/wiki/Snarls-and-chains) as [vg](https://github.com/vgteam/vg) defines them, from the cactus decomposition of [Paten et al. (2018)](https://doi.org/10.1089/cmb.2017.0251). `--snarls` can also read `vg snarls` output.

1. Nodes are renumbered in the order the reference walks them, and nodes it reads backwards are turned around (unless `--no-flip`). Nodes walked round in a loop are placed together. This is `<prefix>.sorted.gfa`.
2. Snarls are found and only the outermost are kept, so each region is described once.
3. A tandem copy the graph stores once is a loop. A site inside it would miss the extra copies, so it is widened to the nearest nodes outside the loop (up to 100 kb each side).
4. Each site is scored from the haplotypes' walks: support, distinct alleles, and interior length. Going straight across is a deletion allele. A haplotype that crosses a site more than once has one allele, from its first crossing through the interior to its last, measured at full length; a straight-across crossing outside that stretch is left out. A walk crosses only when it enters through one boundary and leaves through the other: one that turns back, or passes a boundary on its outer side, does not.
5. Sites are filtered by interior length (`--min-variant-bp`, `--max-variant-bp`) and support (`--min-path-support`, `--min-alt-support`). A site with a node walked backwards passes the size floor whatever its length, as reversing a node is a change in itself; the other filters still apply.
6. Kept sites at most `--merge-nearby-bp` apart on the reference are merged. Sites whose insertion or deletion could be placed in each other are shifted along a repeat (as VCF normalisation does) and fused.
7. Sites must not share interior nodes; of two that do, the smaller is kept.
8. Pieces of one event each under the size floor, all filtered out, are joined so the event is not lost. Only pieces whose own alleles differ in length take part, and they are joined when they lie at most `--merge-nearby-bp` apart (`--min-variant-bp` when merging is off) with no kept site between them. The joined site is kept when the reference crosses it, one of its alleles differs in length from the reference's by `--min-variant-bp` or more, it passes the filters of step 5, and it shares no node with a kept site.
9. With `-r`, an inversion whose two ends lie in a repeat the graph stores once has no site, as the haplotype leaves and comes back through the same nodes. Such inversions are found by reading each haplotype's walk against the reference: a stretch that runs backwards for `--min-variant-bp` (50 bp) or more gets an inversion site, and the sites inside it name it as their parent. Haplotypes stored on the other strand give the same carriers. Inversion sites pass the support filters and `--max-variant-bp` like any site, but not the size floor. Two may share a boundary; of two that overlap more, the one more haplotypes carry is kept and the other reported as a warning, as is an inversion with no room for a boundary outside it.

## Required inputs

- `-i, --gfa <graph.gfa>`: the locus graph ([GFA](https://github.com/GFA-spec/GFA-spec), from [pggb](https://github.com/pangenome/pggb), [minigraph](https://github.com/lh3/minigraph) or `rebuild`)
- `-r, --reference-path <name>`: the reference path; a unique part of the name is enough. Not needed with `--snarls`

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `-o, --out-prefix <p>` | output prefix | `bubble_calls` |
| `--min-variant-bp <N>` | keep a site where some haplotype has N bp or more in the interior (`0` = off) | `50` |
| `--max-variant-bp <N>` | drop a site whose longest interior is over N bp (`0` = no limit) | `0` |
| `--min-path-support <N>` | haplotypes that must pass through a site | `0` |
| `--min-alt-support <N>` | haplotypes that must carry its most common non-reference allele | `0` |
| `--merge-nearby-bp <N>` | merge sites at most N bp apart (`0` = off) | `100` |
| `--no-fuse-ambiguous` | do not fuse sites by placement | off |
| `--no-flip` | sort without turning nodes around | off |
| `--snarls <path>` | take the boundaries from `vg snarls` output (JSONL); the graph is not sorted. Without `-r` the sites have no reference order and are not merged: for viewing only | — |
| `--emit-snarls-jsonl <path>` | also write the boundaries in `vg snarls` format | — |
| `--sorted-gfa-out`, `--bubbles-csv`, `--bandage-csv <path>` | name an output file | from `-o` |
| `--snarl-debug-tsv <path>` | every candidate site, and whether it was kept | — |
| `--gtf <path>` | mark the genes of a reference-coordinate GTF on the reference's nodes; needs a PanSN reference name with its interval (`sample#hap#contig:start-end`) to place the nodes on the genome, and is skipped with a warning otherwise | — |
| `-q, --quiet` | no progress logs | off |

## Outputs

| file | contents |
|------|----------|
| `<prefix>.sorted.gfa` | the sorted graph, input to `panphorte` and `call` |
| `<prefix>.bubbles.csv` | the site table, read by every later module |
| `<prefix>.bandage_nodes.csv` | [Bandage](https://github.com/asl/BandageNG) colours of site interiors (boundaries are left uncoloured): red for kept sites, blue for the other candidates longer than 1 bp or with a node walked backwards, so the many SNP candidates do not hide the larger ones; inversion sites are left to the sites inside them |
| `<prefix>.bandage_genes.csv` | with `--gtf`: genes per node |

A rerun at a prefix removes the default-named files an earlier run left there that this run does not write, so they are not read as this run's: `<prefix>.sorted.gfa` under `--snarls`, and `<prefix>.bandage_genes.csv` when no gene table is written. An input is never removed, and nothing is removed when the graph is read through a pipe or a `/dev` path such as `/dev/stdin`, since such a graph cannot be told apart from the file it came from.

`bubbles.csv` columns:

| column | meaning |
|--------|---------|
| `bubble_id` | ordinary sites in reference order, then inversion sites |
| `source`, `sink`, `source_orient`, `sink_orient` | the boundaries, in the order the reference reaches them, and the strand a crossing reads each on |
| `inside_node_count`, `total_node_count` | interior nodes, and with the two boundaries |
| `path_support` | haplotypes passing through, reference carriers included, so near the panel size at most sites; `--min-alt-support` filters on the other alleles |
| `distinct_alleles` | different alleles among them |
| `ref_allele_support` | haplotypes carrying the reference's allele |
| `alt_allele_support_max`, `alt_allele_support_min` | haplotypes carrying the most and the least common other allele |
| `min_inside_bp`, `max_inside_bp` | shortest and longest interior, counting every copy |
| `inside_nodes` | interior node ids, `;`-separated |
| `kind`, `parent_id` | `site` or `inversion`; a site inside an inversion site names it in `parent_id` |

### Visualization

The sorted graph opens in [BandageNG](https://github.com/asl/BandageNG); File › Load CSV data with `<prefix>.bandage_nodes.csv` colours the sites as above, and with `<prefix>.bandage_genes.csv` labels the nodes with their genes.

## Limitations

- The size floor measures interior length, not change: a small change inside a long allele passes it.
- A large complex region is one site; the variants inside it are not reported separately.

## Example

```bash
# the sites, along the reference path
panvar bubble -i locus.gfa.gz -r GRCh38 -o locus
# and the genes on the reference's nodes, for Bandage
panvar bubble -i locus.gfa.gz -r GRCh38 -o locus --gtf genes.gtf.gz
```

`scripts/bubble_examples.sh` finds the sites of the paralog and tandem-repeat loci in `tests/real_data` and of a rebuilt tangled locus, one folder per graph in `tests/examples/bubble`.
