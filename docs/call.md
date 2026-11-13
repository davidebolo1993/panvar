# Module `call`

CLI: `panvar call`

## What it does

Calls structural variants and copy number from a locus graph into a multi-sample VCF, one column per haplotype (or per sample with `--diploid`). At each site, every haplotype's walk is compared with the reference's.

- A difference is a deletion, an insertion, an inversion or a replacement. Pieces of one event are joined, and the same event in different haplotypes becomes one record.
- A folded repeat is called as a copy number: each haplotype's passes through the unit are counted. With `--cn`, a duplicated stretch the graph did not fold is counted too.
- A copy gained or lost is reported once, as a copy number.
- With `--gtf`, records name the genes they touch, and copy-number sites get a copy number per gene inferred.

| `SVTYPE` | meaning |
|----------|---------|
| `DEL`, `INS` | sequence lost or gained; a replacement with a net change of `--min-sv-bp` or more is a `DEL` or `INS` of that change, flagged `SUBSTITUTION` |
| `SUB` | a replacement with a smaller net change, written as `REF` and `ALT` |
| `INV` | reference sequence the haplotype walks reverse-complemented |
| `CNV` | a copy-number change: `<CNV:TR>` for a folded repeat (`CN_METHOD=UNIT`), `<CNV>` for a duplicated stretch (`CN_METHOD=MODULE`); `FORMAT/CN` is each haplotype's copy number |

## How it works

Base-level alignment of two alleles uses [ksw2](https://github.com/lh3/ksw2), and [minimap2](https://github.com/lh3/minimap2) for long alleles and `--classify-ins`.

1. Haplotypes walking the same nodes across a site share an allele, compared once. A haplotype crossing an assembly gap there is `.`; a site where the reference has a gap is skipped.
2. Folded repeats (from `panphorte`'s unit-graph table) are counted: a haplotype's copy number is its passes through the unit's anchor. `CNBP` is its change in bases, and `CNRESID` the part whole copies do not explain. A copy's orientation is that of its step on the anchor.
3. With `--cn`, a duplicated stretch the graph collapsed onto shared nodes (a module: nodes some path walks more than once) is counted the same way. With no reference copy, one copy's length is the median of the haplotypes' bp per copy (`CN_UNIT_BP`). A module is not counted where its passes are not copies, and a warning suggests `rebuild`: where the site is tangled (ten or more nodes with fifty neighbours or more on one side), or where the module would span more than 80% of the reference (scattered revisits, as a real duplication spans a few percent of a locus).
4. Each allele is compared with the reference's node by node, which finds the nodes the two walks share on either side of a difference. By default the bases between are then aligned, which places each indel exactly where node boundaries would cut it arbitrarily; a `MODULE` site stays node by node, so an extra copy and a novel insertion beside it are not read as one gap. Pieces of one event are joined: a replacement cut by short chance matches, and deletions or inversions closer than `--merge-distance-bp`. A record is written as its minimal change: the bases a replacement's two sides share at either end (different nodes can spell the same bases) are trimmed, and a deletion or insertion that leaves is left-aligned, as VCF normalisation does.
5. The same event in different haplotypes becomes one record: the same type, anchors within `--merge-distance-bp` plus the smaller event's size, length changes less than `--min-sv-bp` apart, and the same sequence (`--merge-seq-identity`); events from the node comparison may instead share their nodes (`--merge-jaccard`). A replacement also joins a deletion or insertion of the same net change. The record takes the event most carriers hold; each carrier's own event is in `carrier_events.tsv`.
6. Non-carriers whose event matches a record are added to it. A haplotype whose larger deletion removes the record's bases is `.`. Records with fewer carriers than `--min-haplotypes`, or an `AF` under `--min-alt-af`, are dropped.
7. A copy change already counted as a `CNV` is removed from its carriers' insertions and deletions, so it is reported once. A haplotype whose copies differ but whose module changes by less than `--min-sv-bp` is no carrier, and keeps its own events.
8. With `--gtf`, each gene at a copy-number site gets its own copy number (see [Per-gene copy number](#per-gene-copy-number) below).
9. At an inversion site (from `bubble`), the inversion is read from the walks: a run walked backwards between two runs walked forwards.

### Per-gene copy number

A copy-number site can hold several genes, often paralogs, so its count says how many copies there are but not of which gene. Each haplotype's sequence across the site is read for each gene:

- A gene belongs to the copies when most of its coding sequence lies on the reference's copies, so a gene between two copies of a repeat does not take the repeat's count. A site's only gene in the copies takes the site's count; the rules below split it between several.
- A gene's private 31-mers are those of its coding sequence found in no other gene of the site. Each copy carries them once, so the haplotype's count of them over their number is the gene's copy number, rounded.
- A copy that differs from the reference's loses private 31-mers (one base removes up to 31) and can round down. So the hits are also read where they lie: a run of them where the whole gene aligns at 98% identity or more, closer to it than to any other gene of the site, is a copy, and the count is raised to the number of such copies.
- Two near-identical paralogs have few private 31-mers, so they are told apart where they differ: the two genes are aligned to each other once, and at each differing position the haplotype's 31-mers say how many of its copies carry each gene's version. The median over positions splits the site's count, so one position converted from one gene to the other does not move it. Not done where both genes lie in every copy, or for three genes that close.
- A copy of a paralog pair is read exon by exon. When the exons at one end match one gene and the rest the other, the copy is a hybrid (`hybrid_copies`); this needs five exons or more that tell the two genes apart, so a pair with a pseudogene is not assessed. A hybrid would read as a copy of both genes, so on a haplotype holding one the copies found by alignment do not raise the count.
- A gene counted past the site's total (its private 31-mers also lie outside the copies), or with nothing to tell it apart (no private 31-mers, or a pair with no usable position on that haplotype), is given the site's total and marked `reliable` 0.

## Required inputs

- `-i, --gfa <graph.gfa>`: the sorted graph from `bubble`, `panphorte` or `refine`
- `-b, --bubble-prefix <prefix>` (reads `<prefix>.bubbles.csv`) or `-c, --bubbles-csv <path>`: the sites
- `-r, --reference-path <name>`: the reference; positions come from the interval in its name (`sample#hap#contig:start-end`)
- `-o, --out-prefix <prefix>`: output prefix

panphorte's unit-graph table is read from the site table's prefix (or `--panphorte-prefix`). Without it, a unit of one node (linked to itself) is still counted, but a unit of several nodes is counted only with `--cn`, as a duplicated stretch (`MODULE`). A carrier's events of whole copies of the unit beside it are part of that copy change (a gain of 4 can be spelled as 6 copies inserted and 2 deleted), so they are dropped and the change is reported once; a part of a copy stays an event of its own. Without `--cn`, its copy changes read as insertions and deletions. With the table, none of this applies.

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `--panphorte-prefix <prefix>` | read the unit-graph table from this prefix | the site table's |
| `--cn` | also count duplicated stretches the graph did not fold (`MODULE`) | off |
| `--gtf <path>` | genes per record, and copy number per gene (needs a PanSN reference name) | — |
| `--bubble-id <N>` | only this site; repeatable | all |
| `--min-sv-bp <N>` | smallest event reported | `50` |
| `--rescue-min-bp <N>` | smallest piece kept to join into a larger event | half `--min-sv-bp` |
| `--no-align-alleles` | compare node by node only | off |
| `--merge-distance-bp <N>`, `--merge-jaccard <X>`, `--merge-seq-identity <X>`, `--merge-size-ratio <X>` | when two events are the same; a size ratio of `0` uses `--merge-seq-identity` | `100`, `0.80`, `0.80`, `0` |
| `--min-haplotypes <N>`, `--min-alt-af <X>` | fewest carriers, as a count and as the record's `AF` (`AC/AN`; with `--diploid`, over haplotype slots) | `1`, `0` |
| `--classify-ins` | tell new insertions (`NOVEL`) from duplications (`DUP`) with minimap2 (`--minimap-preset`, `--minimap-best-n`, `--ins-dup-min-identity`), for every `INS` record by the bases it inserts | off |
| `--diploid` | one column per sample, with phased `GT` (`/` for a sample with unnumbered haplotypes); `CN` is the sample's total | off |
| `--no-per-bubble-vcf`, `--no-variant-nodes` | skip the per-site VCFs, or the node and carrier tables | off |
| `--max-seq-bp <N>` | leave `INSSEQ`, `DELSEQ` and `INVSEQ` out above N bases, so VCF tools read the lines; `0` writes them all | `20000` |
| `-t, --threads <N>` | worker threads (`0` = all cores) | `0` |
| `-q, --quiet` | no progress logs | off |

## Outputs

| file | contents |
|------|----------|
| `<prefix>.region.vcf` | every record, sorted by position |
| `<prefix>.bubble_<id>.vcf` | one site's records |
| `<prefix>.variant_nodes.tsv` | per record, the nodes it involves; read by `describe --variant-nodes` |
| `<prefix>.carrier_events.tsv` | per record and carrier, the event that carrier holds; read by `benchmark` |
| `<prefix>.node_genes.tsv`, `<prefix>.dup_gene_cn.tsv` | with `--gtf`: genes per node, and copy number per gene and haplotype |

VCF 4.2. `POS` is the base before the event, `END` its last reference base. `FILTER` is `PASS`: a record that fails a filter is not written. `GT` is `1` for a carrier, `0` for reference-like, `.` for no allele there. Main `INFO` fields:

| field | meaning |
|-------|---------|
| `SVTYPE`, `SVLEN`, `END` | type, length change (an inversion: its length), last base |
| `INSSEQ`, `DELSEQ`, `INVSEQ`, `SUBSTITUTION` | the event's sequence, left out above `--max-seq-bp` (20,000) so VCF tools can read the lines; a replacement's two sides |
| `IMPRECISE`, `CIPOS`, `CIEND`, `HOMLEN` | uncertain or sliding breakpoints: `IMPRECISE` marks a cluster of pieces (`END`-`POS` wider than \|`SVLEN`\|), a breakpoint at an N run of the reference, or a `PATH_INVERSION` whose breakpoints lie in repeats |
| `AN`, `AC`, `AF` | haplotypes with an allele at the site (not `.`), carriers, ratio |
| `BUBBLE_ID`, `START_NODE`, `END_NODE`, `EVENT_NODES` | the site and the event's nodes; a `PATH_INVERSION`, typed from the walks, has no nodes |
| `CN_METHOD`, `REF_CN`, `RU_LEN` | how copies were counted, the reference's count, a unit copy's length |
| `NMERGED`, `SVLEN_RANGE`, `MERGE_*` | how many carriers were merged and how much they differ |
| `GENES` | with `--gtf`: genes over the record; inside a copy of a `MODULE` site holding several genes, the paralogs its carriers hold |

`FORMAT`: `GT`, `CN` (`.` except on a `CNV`; on a `MODULE` record a haplotype that is no carrier, but whose change a written record of its own holds, reads `REF_CN`, as a dosage read off `CN` would otherwise count it twice); on a `CNV` also `CNBP` and `CNRESID` (from the copies measured), on a `MODULE` record `MODBP`, and on a `UNIT` record `CNINV`. `MODBP` is the change over the module's own nodes, which decides `GT` there: a carrier's `CN` differs from `REF_CN` and its `MODBP` reaches `--min-sv-bp`, unless that change lies inside a larger event of its own, whose record then holds it. `CNINV` counts the copies in the opposite orientation to the reference's copy (with no reference copy, to the orientation most haplotypes' copies take). A gain with more such copies than the reference has is flagged `INVDUP`; with no reference copy, only a haplotype holding a copy and its reverse complement is. The full list is in the VCF header.

`carrier_events.tsv`, one row per record and event; carriers holding the same event share a row:

| column | meaning |
|--------|---------|
| `variant_id`, `bubble_id` | the record and its site |
| `kind` | how the event was found: `seq` (base-level alignment), `del`, `ins`, `sub`, `inv` (node comparison), `cn` (a copy-number record); `unresolved` when it could not be placed on the reference |
| `ref_start`, `ref_end` | the reference bases it replaces, 1-based; an insertion has `ref_end` = `ref_start` - 1 and goes in before `ref_start` |
| `cn`, `walk_cn` | on a `CNV`: the copy number reported, and the copies the carrier's walk holds |
| `event_bp` | the event's size as call measured it: net change, or an inversion's length |
| `alt_bp`, `alt` | the bases put in, forward strand |
| `haplotypes` | the carriers |

`dup_gene_cn.tsv`, one row per gene and haplotype at a copy-number site:

| column | meaning |
|--------|---------|
| `bubble_id`, `variant_id` | the site and its `CNV` record (`.` for a module with no record) |
| `sample`, `genes`, `gene_id` | the haplotype, the gene's name and its id |
| `cn` | the gene's copy number on this haplotype, from the copies it holds (as is the module's total): for a haplotype that is no carrier of the site's record it can differ from the VCF's `CN`, which is then `REF_CN` |
| `reliable` | 1 when backed by the haplotype's k-mers, or for a site's only gene; 0 when it is counted past the module's total or cannot be told apart from its paralogs, and `cn` is then the module's total |
| `dosage`, `hits`, `priv_kmers` | `hits` / `priv_kmers`: the gene's private 31-mers found in the haplotype, over how many it has |
| `copies` | copies found by aligning the whole gene, where they could raise `cn` |
| `hybrid_copies` | copies that switch from the gene to its paralog (`.` where not assessed) |

### Visualization

`scripts/build_variant_node_data.R` bundles the graph, the site table and call's outputs into one `.rds` file, and `scripts/variant_node_heatmap_app.R` opens it in the browser (R packages data.table, and shiny and plotly for the viewer). For a site or the whole locus, it shows each haplotype's passes through the nodes, each record's nodes and each haplotype's genotype; clicking a record selects its carriers:

```
Rscript scripts/build_variant_node_data.R --gfa <graph.gfa> --prefix <call prefix> --bubbles <bubbles.csv>
VN_RDS=<call prefix>.viewer.rds Rscript scripts/variant_node_heatmap_app.R
```

## Limitations

- Only structural variants: no SNP or small-indel records.
- A merged record carries one representative event; each carrier's own is in `carrier_events.tsv`.
- On a `MODULE` record, copies of different lengths show up in `CNRESID`, and a partial copy that keeps more than half of the anchor counts as a copy. The record spans the site's interior, which can hold more than the module.
- A unit-graph table of several nodes is checked by its node ids, anchors and lengths only, so a table from another run whose ids all exist can pass; give the table panphorte wrote with this graph.
- A repeat copy below panphorte's `--min-similarity` is not folded, so it is not counted: the haplotype reads as one copy fewer plus an insertion, which `CNRESID` shows.
- A path crossing a site more than once is compared over one walk; a crossing outside it is flagged (`MULTI_CROSSING_UNSCORED`, with a warning) only when its length differs from every reference crossing by `--min-sv-bp` or more.
- Per-gene copy numbers come from 31-mers. A hybrid copy counts toward the gene it mostly matches, or in part toward both; a partial copy counts by the share of the gene it holds, rounded. A related sequence carrying part of a gene's alleles adds to its count. Copies closer than 10 kb, or than the gene's length, form one run of hits, so the alignment finds only one of them and cannot raise the count past it.

## Example

`scripts/call_examples.sh` calls every graph in `tests/examples/refine` (with a GTF, also per gene) and builds a viewer bundle for each, in `tests/examples/call`.
