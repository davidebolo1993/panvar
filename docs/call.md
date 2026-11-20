# Module `call`

CLI: `panvar call`

## What it does

Calls structural variants and copy number from a locus graph into a multi-sample VCF. The VCF has one column per haplotype, or one per sample with `--diploid`. At each site, every haplotype's walk is compared with the reference's.

- A difference is a deletion, an insertion, an inversion or a replacement. Pieces of one event are joined. The same event in different haplotypes becomes one record.
- A folded repeat is called as a copy number: each haplotype's passes through the repeat's unit are counted. With `--cn`, a duplicated stretch the graph did not fold is counted too.
- A copy gained or lost is reported once, as a copy number. At a duplicated stretch (`--cn`), the copy-number record names the copies each carrier gained or lost. The carrier's other changes there are records of their own.
- With `--gtf`, records name the genes they touch, and each gene at a copy-number site gets its own copy number.

| `SVTYPE` | meaning |
|----------|---------|
| `DEL`, `INS` | sequence lost or gained. A replacement whose net change is `--min-sv-bp` or more is a `DEL` or `INS` of that change, flagged `SUBSTITUTION` |
| `SUB` | a replacement (bases replaced by other bases) whose net change is under `--min-sv-bp`, spelled in `REF` and `ALT` |
| `INV` | reference sequence the haplotype walks reverse-complemented |
| `CNV` | a copy-number change. `<CNV:TR>` is a folded repeat (`CN_METHOD=UNIT`). `<CNV>` is a duplicated stretch (`CN_METHOD=MODULE`). `FORMAT/CN` is each haplotype's copy number |

## How it works

Two alleles are aligned base by base with [ksw2](https://github.com/lh3/ksw2), or with [minimap2](https://github.com/lh3/minimap2) when they are long. `--classify-ins` also uses minimap2.

1. Haplotypes that walk the same nodes across a site share an allele, which is compared once. A haplotype that crosses an assembly gap at the site is `.`, and a site where the reference has a gap is skipped.

   Where the graph collapsed the copies of a duplication onto the same nodes, a path crosses a site once per copy it holds. One walk from its first crossing to its last would take in the stretch between the copies, at every site of the copies, so a copy gained or lost would be written once per site. Such sites are read one copy at a time instead:
   - Each path's crossings of the sites, in its own order, are matched one to one and in order to the reference's. A pair is one site crossed in the same direction, so copies inverted against each other are kept apart. The sites crossed once on either side anchor the matching, and so do the reference's first and last nodes where a path walks them once, so the matching is decided once per path for all the sites of the copies.
   - A path is read forwards or backwards as it walks the reference's first and last nodes, as stored on either strand. An inversion between inverted copies would read almost as well backwards, with the copies swapped.
   - A path holding fewer copies than the reference is matched to the reference copies it is most like.
   - Each reference crossing of such a site is a site of its own. A path is read from the crossing matched to it, and is `.` where it has none, as it lacks that copy.
   - Between two matched crossings, a path holding copies the reference does not, or lacking some it holds, or walking another stretch than the reference's, is read from one to the other on both walks. A gained or lost copy is then one insertion or deletion, written once, that carries the path's own bases: the copy's changes are not written apart. A path holding the reference's copies there is `0`. Where the two crossings meet at one boundary node, the stretch starts at the first crossing, which it then holds.
   - Crossings that read the reference's backwards, in reverse order, inside an inversion site are an inversion, written at that site (step 9) and not again between the copies.
   - A crossing the matching cannot place, one outside the anchors, is read in one walk with the nearest crossing of its site and direction, and the records say so (`MULTI_CROSSING`).

   These sites and stretches keep the site's `BUBBLE_ID` (a stretch, the site's before it), with ids numbered apart.

2. Folded repeats come from `panphorte`'s unit-graph table. A haplotype's copy number is how many times it passes the unit's anchor, and a copy's orientation is that of its step on the anchor. `CNBP` is its change in bases, and `CNRESID` the part of that change whole copies do not explain.

3. With `--cn`, duplicated stretches the graph collapsed onto shared nodes, the nodes some path walks more than once, are counted too, as modules. Their copies are counted as in step 2. One copy's length is `CN_UNIT_BP`; with no reference copy, it is the median of the haplotypes' bases per copy. `CNBP` runs from a walk's first module node to its last and takes in sequence that replaces the reference's first or last module nodes, so it is the walk's whole change across its copies.

   A module is not counted where its passes are not copies, and a warning suggests `rebuild`. Nor is one counted on a stretch between copies (step 1): what it walks twice are pieces of the copies on either side, and the copy gained or lost there is an insertion or deletion. That is a tangled site, where ten or more nodes have fifty neighbours or more on one side, or a module that would span more than 80% of the reference: a real duplication spans a few percent of a locus, while such revisits are scattered.

4. Each allele is first compared with the reference's node by node, which finds the nodes both walks share on either side of a difference. By default the bases between those nodes are then aligned, because node boundaries can cut an indel at an arbitrary place. At a `MODULE` site the comparison stays node by node, so that an extra copy and a novel insertion beside it are not read as one gap.

   The alignment is banded by how far the node comparison drifts. Walks matched out of place at short nodes can understate that drift, so where the band cannot rule out a better alignment outside it, the alignment is redone with a band that can. Pieces of one event are joined: a replacement cut by short chance matches, and deletions or inversions closer than `--merge-distance-bp`.

   A record is written as its minimal change. Different nodes can spell the same bases, so the bases both sides of a replacement share at either end are trimmed, and a deletion or insertion that remains is left-aligned, as VCF normalisation does. Outside a `MODULE` site, a node replacement that is a plain insertion or deletion is written as one, so it merges with the same event found another way. It slides left no further than where the walks' reference starts, because the copies a `UNIT` record counts are not in the walks.

   An event is reported when its size is `--min-sv-bp` or more. Its size is its largest single change: an insertion or deletion is its length, and a replacement counts its larger side. So 300 bases replaced by 290 are reported, though the net change is 10, as a `SUB`. `benchmark` sizes differences the same way.

5. At a duplicated stretch a haplotype can differ in two ways at once: in its number of copies, and by other changes in or next to them. Before anything is merged, `call` separates the two for each haplotype, so each goes to the right record.

   A haplotype carries the copy-number record when its number of copies differs from the reference's and its sequence across the stretch grows or shrinks by `--min-sv-bp` or more (`MODBP`). This is read from its walk alone, so how events are merged cannot change it. The exception is a larger event that takes the copies with it: a deletion that removes a copy and more sequence beside it is reported as that deletion, and the haplotype reads `GT` 0 on the copy-number record. The same holds for an insertion that brings a copy and more.

   The haplotype's changes that overlap its copies are cut where the copies end, and what lies outside them is an ordinary event with a record of its own. A piece that passes through the module's anchor, the node each copy holds once, is a gained or lost copy, which is what `CN` counts. A piece that misses the anchor is a copy only when its sequence is a version of a reference copy, such as an inserted copy that has diverged, or, next to a lost copy, another version of what was lost with it. A deletion that misses the anchor is a change inside a copy that remains. At most twice the larger of `MODBP` and the changed copies' length is taken as copies, so a large unrelated event beside them is not absorbed.

   The copies go on the copy-number record as edits of the reference (`CN_FORMS`, and per carrier `CNFORM`), and with the haplotype's other records they rebuild its walk. Everything else is merged with other haplotypes' events in step 6.

   Where the walks share short nodes inside the copies by chance, the node comparison can cut one change into pieces. These are joined again among the parts that are not copies, so joining never changes the copies: deletions, or insertions, inside the copies and closer than `--merge-distance-bp`, and pieces of any kind under `--min-sv-bp`, at most 10 bases apart and one touching the copies, when together they reach `--min-sv-bp`.

   A copy-number record that `--min-haplotypes` or `--min-alt-af` would drop takes no copies from its carriers: their changes stay ordinary events, so dropping the record never loses them.

6. The same event in different haplotypes becomes one record. Two events are the same when all of these hold:
   - the same type;
   - anchors within `--merge-distance-bp` plus the smaller event's size;
   - length changes less than `--min-sv-bp` apart;
   - the same sequence (`--merge-seq-identity`), or, for events from the node comparison, the same nodes (`--merge-jaccard`).

   A replacement also joins a deletion or insertion of the same net change. The record writes one event: one of `--min-sv-bp` or more before a smaller one, then the one most carriers hold as written (its minimal change, so members the walks placed differently count together), then an aligned event before a node one, then the larger. At a `MODULE` site every copy shares the same nodes, so events merge around their representative and a record holds one event per haplotype; linked pairwise instead, a haplotype's events in different copies would chain into one record that writes only one of them. Each carrier's own event is in `carrier_events.tsv`.

7. Non-carriers whose event matches a record are added to it; at a `MODULE` site, each event joins one record only. A haplotype whose own larger deletion removes the record's bases is `.`. Records with fewer carriers than `--min-haplotypes`, or an `AF` under `--min-alt-af`, are dropped. A record that repeats another, the same position, alleles, sequence and genotypes, is written once: two sites whose walks both span one event would otherwise write it under two ids and a reader would apply it twice. It is written at the site holding most of the event's nodes, so it stays with the site it lies in, and at the first of them in reference order when they hold as many. Copy-number records are not compared, as their counts belong to their site.

8. With `--gtf`, each gene at a copy-number site gets its own copy number (see [Per-gene copy number](#per-gene-copy-number) below).

9. At an inversion site from `bubble`, the inversion is read from the walks: a run walked backwards between two runs walked forwards.

10. Records that are alternative alleles of one repeat site share `INFO/SITE` and `SITE_N`, so a reader can see them as one site; genotypes stay per record, and `bcftools query -f '%INFO/SITE\t%ID\n'` groups them. A site is two or more records that are not `CNV`s and span at most 1 kb each, so larger records do not chain unrelated sites. Each record is widened to the right by how far it can slide (`HOMLEN`, `CIPOS`, `CIEND`), and records within 10 bp of each other belong to one site. A site is tagged only when no haplotype carries two of its records, since alleles of one site exclude each other; a mixed group is left untagged and counted in a warning. A haplotype is read as `--diploid` reads it, so its contigs count together.

### Per-gene copy number

A copy-number site can hold several genes, often paralogs. Its count says how many copies there are, but not of which gene, so each haplotype's sequence across the site is read for each gene.

A gene belongs to the copies when most of its coding sequence lies on the reference's copies, so a gene between two copies of a repeat does not take the repeat's count. A site's only gene in the copies takes the site's count; with several, the count is split as follows.

A gene's private 31-mers are those of its coding sequence found in no other gene of the site. Each copy carries them once, so the gene's copy number is the haplotype's count of them over their number, rounded. A copy that differs from the reference's loses private 31-mers, up to 31 per differing base, and can round down. So the hits are also read where they lie: a run of hits where the whole gene aligns at 98% identity or more, closer to it than to any other gene of the site, is a copy, and the count is raised to the number of such copies.

Two near-identical paralogs have few private 31-mers, so they are told apart where they differ. They are aligned to each other once, and at each position where they differ the haplotype's 31-mers say how many of its copies carry each gene's version. The median over positions splits the site's count, so one position converted from one gene to the other does not move it. This is not done where both genes lie in every copy, or for three genes this close.

A copy of a paralog pair is also read exon by exon. When the exons at one end match one gene and the rest match the other, the copy is a hybrid (`hybrid_copies`). A hybrid would read as a copy of both genes, so on a haplotype holding one, the copies found by alignment do not raise the count. This needs five exons or more that tell the two genes apart, so a pair with a pseudogene is not assessed.

A gene counted past the site's total, because its private 31-mers also lie outside the copies, is given the site's total and marked `reliable` 0. So is a gene nothing tells apart: one with no private 31-mers, or a pair with no usable position on that haplotype.

## Required inputs

- `-i, --gfa <graph.gfa>`: the sorted graph from `bubble`, `panphorte` or `refine`
- `-b, --bubble-prefix <prefix>` (reads `<prefix>.bubbles.csv`) or `-c, --bubbles-csv <path>`: the sites
- `-r, --reference-path <name>`: the reference; positions come from the interval in its name (`sample#hap#contig:start-end`)
- `-o, --out-prefix <prefix>`: output prefix

`panphorte`'s unit-graph table is read from the site table's prefix, or from `--panphorte-prefix`. Without the table:

- A unit of one node, linked to itself, is still counted.
- A unit of several nodes is counted only with `--cn`, as a duplicated stretch (`MODULE`). A carrier's events made of whole copies of the unit beside it are part of that copy change, so they are dropped and the change is reported once. A gain of 4 can be spelled as 6 copies inserted and 2 deleted. A part of a copy stays an event of its own.
- Without `--cn`, its copy changes read as insertions and deletions.

With the table, none of this applies.

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
| `--max-seq-bp <N>` | write `INSSEQ`, `DELSEQ`, `INVSEQ` and each `CN_FORMS` sequence over N bases to `<prefix>.alt_seqs.fa` instead of the line, so VCF tools can read the lines; `0` writes them all in the line | `20000` |
| `--max-fa-seq-bp <N>` | leave sequences over N bases out of `<prefix>.alt_seqs.fa` too, to cap the output's size; the VCF then cannot rebuild those records | `0` (none) |
| `-t, --threads <N>` | worker threads (`0` = all cores) | `0` |
| `-q, --quiet` | no progress logs | off |

## Outputs

| file | contents |
|------|----------|
| `<prefix>.region.vcf` | every record, sorted by position |
| `<prefix>.bubble_<id>.vcf` | one site's records |
| `<prefix>.variant_nodes.tsv` | per record, the nodes it involves; read by `describe --variant-nodes` |
| `<prefix>.carrier_events.tsv` | per record and carrier, the event that carrier holds; read by `benchmark` |
| `<prefix>.alt_seqs.fa` | sequences too long for their record's line (`SEQ_FA`), each as `<ID>.<field>`; written when there are any. A run removes one an earlier run left at the prefix, and its index (`.fai`, `.gzi`) either way, whose offsets would no longer match |
| `<prefix>.node_genes.tsv`, `<prefix>.dup_gene_cn.tsv` | with `--gtf`: genes per node, and copy number per gene and haplotype |

The VCF is version 4.2.

- `POS` is the base before the event, and `END` its last reference base. A deletion removes bases `POS`+1..`END`. An insertion spans no reference bases: `END` = `POS`, and its bases go in after `POS`.
- `FILTER` is `PASS`: a record that fails a filter is not written.
- `GT` is `1` for a carrier, `0` for reference-like, `.` for no allele there.
- Without `--diploid`, each path is a column, the reference's too. The reference reads `0` and is left out of `AN`.

Main `INFO` fields, by purpose. The VCF header describes every field.

Position and size:

| field | meaning |
|-------|---------|
| `SVTYPE` | the record's type |
| `END` | the last reference base. On a `CNV`, `POS`+1..`END` covers the reference's copies |
| `SVLEN` | the length change; on an `INV`, the inverted length. Left out on a `MODULE` record, and on a `UNIT` record whose carriers both gain and lose copies: no one size fits all carriers, so each carrier's size is its `CNBP` |
| `AN`, `AC`, `AF` | haplotypes with an allele at the site (not `.`); carriers; `AC`/`AN`. The reference path is left out. A haplotype in several contigs counts once, and is a carrier where one of its contigs carries the record, as with `--diploid` |
| `BUBBLE_ID` | the site |
| `EVENT_NODES` | the event's nodes; a `PATH_INVERSION` has none |
| `SITE`, `SITE_N` | the repeat site the record is one allele of, and how many records the site has (step 10). The id is `bubble<B>_SITE_<k>`: `B` is the `BUBBLE_ID` of the site's first record, `k` counts the sites starting there from the left. Absent outside a site |

Sequence:

| field | meaning |
|-------|---------|
| `DELSEQ`, `INSSEQ`, `INVSEQ` | the bases deleted, inserted or inverted |
| `SUBSTITUTION` | the `DEL` or `INS` replaces bases with other bases; `SVLEN` is the net change. Where both sides fit `--max-seq-bp` and the breakpoints are exact, `REF` and `ALT` spell it after the base at `POS`, as on a `SUB` |
| `SEQ_FA` | which sequences are in `<prefix>.alt_seqs.fa` instead of the line: `INSSEQ`, `DELSEQ`, `INVSEQ`, or `CN_FORMS.<k>` for the k-th form. Each is the FASTA record `<ID>.<field>`, so one can be fetched by name after indexing, as `samtools faidx <prefix>.alt_seqs.fa <ID>.INSSEQ` |
| `INS_SUBTYPE` | with `--classify-ins`: `DUP` or `NOVEL` |

A sequence goes to the side FASTA when it is over `--max-seq-bp` (20,000). So every line stays short, and a novel insertion keeps its sequence as much as a copy does.

Breakpoints:

| field | meaning |
|-------|---------|
| `HOMLEN` | on an exact `DEL` or `INS`: how many bases the event can slide right and stay the same. `POS` is its leftmost place |
| `CIPOS`, `CIEND` | how far each breakpoint can lie from `POS` and `END` |
| `IMPRECISE` | the breakpoints are not exact. The record is a cluster of pieces (`END`-`POS` wider than \|`SVLEN`\|), or a breakpoint is at an N run of the reference, or it is a `PATH_INVERSION` with breakpoints in repeats |
| `PATH_INVERSION` | an inversion read from the walks at an inversion site |

Copy number:

| field | meaning |
|-------|---------|
| `CN_METHOD` | how copies were counted: `UNIT` (a folded repeat) or `MODULE` (a duplicated stretch, `--cn`) |
| `REF_CN` | the reference's copies |
| `RU_LEN` | on `UNIT`: one copy's length |
| `CN_UNIT_BP` | on `MODULE`: one copy's length |
| `CN_FORMS` | on `MODULE`: the copies the carriers gained or lost, each as an edit of the reference (below). Distinct copies are listed once, in reference order, numbered from 1 |

A `CN_FORMS` entry is `START-END:SEQ`: reference bases `START`..`END` are replaced by `SEQ`, in reference orientation.

- An insertion replaces no bases. Its stretch is empty, `END` = `START`-1, and `SEQ` goes in before base `START`.
- A deletion has `SEQ` empty.
- `SEQ` names the reference stretches it repeats where that is shorter: `@A-B` for bases `A`..`B`, `~A-B` for their reverse complement. These are joined with literal bases by `+`. A gained copy repeats the reference's copies, so its form stays short however long the copy.
- `SEQ` is `*` where its literal bases are over `--max-seq-bp`. The sequence is then in `<prefix>.alt_seqs.fa`.

Merging:

| field | meaning |
|-------|---------|
| `SVLEN_RANGE` | the smallest and largest size among the merged events, when they differ |
| `MERGE_MIN_SIM` | the weakest similarity of a merged event to the record's event, in the test that admitted it. Under that test's threshold, the event joined through others |
| `MERGE_TEST` | that test: `spliced`, `jaccard` or `seq` |

Genes:

| field | meaning |
|-------|---------|
| `GENES` | with `--gtf`: the genes over the record. Inside a copy of a `MODULE` site holding several genes, the paralogs its carriers hold |

`FORMAT` fields. The copy-number fields are `.` on records that are not `CNV`s.

| field | on | meaning |
|-------|----|---------|
| `GT` | every record | `1` carrier, `0` reference-like, `.` no allele there. A haplotype is `.` where it does not cross the site, crosses an assembly gap there, could not be aligned there, or its own larger deletion removes the record's bases. With `--diploid`, the sample's haplotypes are joined by `\|` (phased) or `/` |
| `CN` | `CNV` | the haplotype's copies of the unit (`UNIT`) or of the module (`MODULE`). On a `MODULE` record, a non-carrier whose change is in a record of its own reads `REF_CN`. Otherwise a dosage read off `CN` would count that change twice |
| `CNBP` | `CNV` | bases gained (+) or lost (−) against the reference, from the haplotype's walk. On `UNIT`, over the site. On `MODULE`, from its first module node to its last, with everything walked between |
| `CNRESID` | `CNV` | the part of `CNBP` that whole copies do not explain: `CNBP` − (`CN` − `REF_CN`) × one copy's length. `0` means the whole change is copy number |
| `MODBP` | `MODULE` | the change over the module's own nodes, each visit counted. It decides `GT` (step 5) |
| `CNFORM` | `MODULE` with forms | the carrier's copies, as `CN_FORMS` numbers joined by `+`. `.` for a non-carrier. Also `.` for a carrier with no copy written: none of its changes is a copy, or none can be written as an edit of one reference stretch. Such a change stays a record of its own. With `--diploid`, joined per haplotype as `GT` is |
| `CNINV` | `UNIT` | how many of the haplotype's copies run opposite to the reference's copy. With no reference copy, opposite to the orientation most haplotypes' copies take. A gain with more such copies than the reference has is flagged `INVDUP`. With no reference copy, only a haplotype holding a copy and its reverse complement is |
| `HCN` | `CNV`, `--diploid` | each haplotype's `CN`, in `GT` order, joined as `GT` is. `CN` is their sum, and so are `CNBP`, `CNRESID`, `MODBP` and `CNINV` |

`carrier_events.tsv` holds each carrier's own event, as found. Carriers holding the same event share a row.

| column | meaning |
|--------|---------|
| `variant_id`, `bubble_id` | the record and its site |
| `kind` | how the event was found: `seq` (aligned bases); `del`, `ins`, `sub`, `inv` (node comparison); `cn` (a copy-number record, below); `unresolved` (could not be placed on the reference) |
| `ref_start`, `ref_end` | the reference bases the event replaces, 1-based. An insertion has `ref_end` = `ref_start` - 1 and goes in before `ref_start` |
| `cn`, `walk_cn` | on a `CNV`: the copy number reported, and the copies the carrier's walk holds |
| `event_bp` | the event's size as `call` measured it: the net change, or an inversion's length |
| `alt_bp`, `alt` | the bases put in, forward strand |
| `haplotypes` | the carriers |

A `cn` row is the carrier's copy stretch, so no base is both in the stretch and in another event of the carrier.

- Where the allele was aligned, the stretch runs from its first copy to its last on either walk. It is cut where the alignment cuts the two walks. An insertion where the copies meet the rest keeps its side of them.
- Where it was compared node by node, the stretch runs up to the nodes the walk shares with the reference.

`dup_gene_cn.tsv` has one row per gene and haplotype at a copy-number site.

| column | meaning |
|--------|---------|
| `bubble_id`, `variant_id` | the site and its `CNV` record (`.` for a module with no record) |
| `sample`, `genes`, `gene_id` | the haplotype, the gene's name and its id |
| `cn` | the gene's copies on this haplotype, from the copies it holds. For a non-carrier of the site's record it can differ from the VCF's `CN`, which then reads `REF_CN` |
| `reliable` | `1` when backed by the haplotype's k-mers, or for a site's only gene. `0` when the gene is counted past the site's total or cannot be told apart from its paralogs; `cn` is then the site's total |
| `dosage`, `hits`, `priv_kmers` | `hits` / `priv_kmers`: the gene's private 31-mers found in the haplotype, over how many it has |
| `copies` | copies found by aligning the whole gene, which can raise `cn` |
| `hybrid_copies` | copies that switch from the gene to its paralog (`.` where not assessed) |

### Visualization

`scripts/build_variant_node_data.R` bundles the graph, the site table and call's outputs into one `.rds` file. `scripts/variant_node_heatmap_app.R` opens it in the browser (R packages data.table, and shiny and plotly for the viewer). For a site or the whole locus, it shows each haplotype's passes through the nodes, each record's nodes and each haplotype's genotype. Clicking a record selects its carriers.

```
Rscript scripts/build_variant_node_data.R --gfa <graph.gfa> --prefix <call prefix> --bubbles <bubbles.csv>
VN_RDS=<call prefix>.viewer.rds Rscript scripts/variant_node_heatmap_app.R
```

## Limitations

- Only structural variants: no SNP or small-indel records.
- A merged record carries one representative event; each carrier's own is in `carrier_events.tsv`.
- On a `MODULE` record, copies of different lengths show up in `CNRESID`, and a partial copy that keeps more than half of the anchor counts as a copy. The record spans the site's interior, which can hold more than the module.
- A unit-graph table of several nodes is checked by its node ids, anchors and lengths only, so a table from another run whose ids all exist can pass. Give the table panphorte wrote with this graph.
- A repeat copy below panphorte's `--min-similarity` is not folded, so it is not counted: the haplotype reads as one copy fewer plus an insertion, which `CNRESID` shows.
- A gained or lost copy of a collapsed duplication is an insertion or deletion, not a copy number.
- A crossing outside the anchors of step 1 is read in one walk with the nearest crossing of its site and direction (`MULTI_CROSSING`); one with no such crossing is not read.
- A site table without boundary orientations gives no crossings to match: each path is read in one walk from its first crossing to its last, and a crossing outside it is flagged (`MULTI_CROSSING_UNSCORED`, with a warning) only when its length differs from every reference crossing by `--min-sv-bp` or more.
- Per-gene copy numbers come from 31-mers. A hybrid copy counts toward the gene it mostly matches, or in part toward both. A partial copy counts by the share of the gene it holds, rounded. A related sequence carrying part of a gene's alleles adds to its count. Copies closer than 10 kb, or than the gene's length, form one run of hits, so the alignment finds only one of them and cannot raise the count past it.

## Example

```bash
# one VCF column per haplotype
panvar call -i locus_r.normalized.sorted.gfa -b locus_r -r GRCh38 -o locus_call
# also the copies of duplicated stretches the graph did not fold
panvar call -i locus_r.normalized.sorted.gfa -b locus_r -r GRCh38 -o locus_call_cn --cn
# the genes each record touches, and copy number per gene
panvar call -i locus_r.normalized.sorted.gfa -b locus_r -r GRCh38 -o locus_call_genes --gtf genes.gtf.gz
# one VCF column per sample
panvar call -i locus_r.normalized.sorted.gfa -b locus_r -r GRCh38 -o locus_call_diploid --diploid
```

`scripts/call_examples.sh` calls every graph in `tests/examples/refine` (with a GTF, also per gene) and builds a viewer bundle for each, in `tests/examples/call`.
