# Module `call`

CLI: `panvar call`

## What it does

Calls structural variants and copy number from a locus graph into a multi-sample VCF, one column per haplotype (or per sample with `--diploid`). At each site, every haplotype's walk is compared with the reference's.

- A difference is a deletion, an insertion, an inversion or a replacement. Pieces of one event are joined, and the same event in different haplotypes becomes one record.
- A folded repeat is called as a copy number: each haplotype's passes through the unit are counted. With `--cn`, a duplicated stretch the graph did not fold is counted too.
- A copy gained or lost is reported once, as a copy number. At a duplicated stretch (`--cn`), the copy-number record names the copies each carrier gained or lost, and the carrier's other changes there are records of their own.
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
3. With `--cn`, a duplicated stretch the graph collapsed onto shared nodes (a module: nodes some path walks more than once) is counted the same way. With no reference copy, one copy's length is the median of the haplotypes' bp per copy (`CN_UNIT_BP`). `CNBP` runs from each walk's first module node to its last, and takes in sequence that takes the place of the reference's first or last module nodes, so it is the walk's change across its copies. A module is not counted where its passes are not copies, and a warning suggests `rebuild`: where the site is tangled (ten or more nodes with fifty neighbours or more on one side), or where the module would span more than 80% of the reference (scattered revisits, as a real duplication spans a few percent of a locus).
4. Each allele is compared with the reference's node by node, which finds the nodes the two walks share on either side of a difference. By default the bases between are then aligned, which places each indel exactly where node boundaries would cut it arbitrarily; a `MODULE` site stays node by node, so an extra copy and a novel insertion beside it are not read as one gap. The alignment is banded by how far the node comparison drifts; where the band's score cannot rule out a better alignment outside it, it is redone with a band that can, since walks matched at short nodes out of place understate the drift. Pieces of one event are joined: a replacement cut by short chance matches, and deletions or inversions closer than `--merge-distance-bp`. A record is written as its minimal change: the bases a replacement's two sides share at either end (different nodes can spell the same bases) are trimmed, and a deletion or insertion that leaves is left-aligned, as VCF normalisation does. Outside a `MODULE` site, a node replacement that is a plain insertion or deletion is made one, so it merges with the same event found another way; it slides no further left than where the walks' reference starts, as copies a `UNIT` record counts are out of the walks.
5. At a `MODULE` site, each haplotype's copies are told apart from its other changes before anything is merged, from its walk and its own events:
   - A haplotype whose copy number differs and whose module changes by `--min-sv-bp` or more (`MODBP`) is a carrier, unless the change beside its copies goes their way and is more than the copies: a deletion that takes a copy with it is a deletion, and its carrier reads `GT` 0. No event enters this, so how events are merged cannot change it.
   - Its events touching its copies (from the first to the last module node on each walk, widened over sequence that takes the place of the reference's first or last module nodes) are cut at the copies' ends, on node boundaries; a part outside is an event of its own. A replacement of a copy by a diverged version of it is not cut.
   - A part holding passes through the module's anchor, what `CN` counts, is a copy, either way. A part with none is a copy only by its bases: a gain's insertion that is a version of the reference's copies, or, beside a loss, an insertion that is another version of what the loss took beyond its copies. A deletion with no anchor pass is a change inside a remaining copy. The copies are at most twice the larger of `MODBP` and the copies changed.
   - The copies are written on the record as edits of the reference (`CN_FORMS`, `CNFORM`), so the record and the carrier's other records rebuild its walk. The rest goes on to the merge.
   - Then pieces of one change are joined, on the parts left, so joining never changes the copies: deletions, or insertions, inside the copies closer than `--merge-distance-bp`; and pieces of any kind each under `--min-sv-bp`, 10 bases apart or less, one of them touching the copies, when together they reach it. The node comparison cuts one change where the walks meet at short nodes by chance inside the copies.
   - A record whose carriers would not pass `--min-haplotypes` and `--min-alt-af` at the lowest `AF` it can be written with takes no copies from them, so a filter at writing cannot drop copies with it.
6. The same event in different haplotypes becomes one record: the same type, anchors within `--merge-distance-bp` plus the smaller event's size, length changes less than `--min-sv-bp` apart, and the same sequence (`--merge-seq-identity`); events from the node comparison may instead share their nodes (`--merge-jaccard`). A replacement also joins a deletion or insertion of the same net change. The record takes the event most of its carriers hold as written (its minimal change, so members the walks placed differently but that write the same count together): one of `--min-sv-bp` or more before a smaller one, then the most carriers, then an aligned event before a node one, then the larger. At a `MODULE` site every copy shares the same nodes, so events merge around their representative and a record holds one event per haplotype: linked pairwise, a haplotype's events in different copies would chain into one record, which writes one of them. Each carrier's own event is in `carrier_events.tsv`.
7. Non-carriers whose event matches a record are added to it; at a `MODULE` site each event to one record only. A haplotype whose larger deletion removes the record's bases is `.`. Records with fewer carriers than `--min-haplotypes`, or an `AF` under `--min-alt-af`, are dropped.
8. With `--gtf`, each gene at a copy-number site gets its own copy number (see [Per-gene copy number](#per-gene-copy-number) below).
9. At an inversion site (from `bubble`), the inversion is read from the walks: a run walked backwards between two runs walked forwards.

### Event size

One rule sizes an event, here and in `benchmark`, so an event `call` leaves under `--min-sv-bp` is a small difference there too, and one it reports is a large one:

- An insertion or a deletion is its length.
- A change with bases on both sides is aligned (step 4's scoring). Gaps with a stretch between them that matches well enough to pay for a new gap (twelve bases) are separate changes, so scattered edits in a repeat do not add up to one large event. An alignment with no gap at all (an equal-length replacement, or an inversion the alignment spells base for base) has its mismatches as changes, grouped as `benchmark` groups a difference with no gap: mismatches at most ten columns apart, in a window holding one in five columns or more, are one change from the first differing base to the last, so it keeps its span; any other mismatch is 1 bp. With a gap, mismatches beside it are small edits, as in a repeat. Each change counts its larger side, so an unrelated replacement counts whole; the event's size is its largest change. A side over 20 kb is not aligned: the event's size is its larger side, the bases both sides share at either end left out.
- An event is reported when its size is `--min-sv-bp` or more. A replacement is then typed by its net change (see the table above).

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
| `--max-seq-bp <N>` | write `INSSEQ`, `DELSEQ`, `INVSEQ` and each `CN_FORMS` sequence (its literal bases) over N bases to `<prefix>.alt_seqs.fa` instead of the line, so VCF tools read the lines; `0` writes them all in the line | `20000` |
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

VCF 4.2. `POS` is the base before the event, `END` its last reference base. `FILTER` is `PASS`: a record that fails a filter is not written. `GT` is `1` for a carrier, `0` for reference-like, `.` for no allele there. Without `--diploid` each path is a column, the reference's too: it reads `0` and is left out of `AN`. Main `INFO` fields:

| field | meaning |
|-------|---------|
| `SVTYPE`, `SVLEN`, `END` | type, length change (an inversion: its length), last base. `SVLEN` is left out on a `MODULE` record, and on a `UNIT` record whose carriers both gain and lose copies: no one size fits their carriers, whose own is `CNBP` |
| `INSSEQ`, `DELSEQ`, `INVSEQ`, `SUBSTITUTION` | the event's sequence; above `--max-seq-bp` (20,000) it is in `<prefix>.alt_seqs.fa` instead, so VCF tools can read the lines. A replacement (`SUB`, or a `DEL`/`INS` flagged `SUBSTITUTION`) is spelled in `REF` and `ALT` where both sides fit and its breakpoints are exact |
| `IMPRECISE`, `CIPOS`, `CIEND`, `HOMLEN` | uncertain or sliding breakpoints: `IMPRECISE` marks a cluster of pieces (`END`-`POS` wider than \|`SVLEN`\|), a breakpoint at an N run of the reference, or a `PATH_INVERSION` whose breakpoints lie in repeats; its `CIPOS`/`CIEND` give how far each breakpoint can lie |
| `AN`, `AC`, `AF` | haplotypes with an allele at the site (not `.`), carriers, ratio. The reference path is left out, and a haplotype in several contigs counts once, a carrier where one of its contigs carries the record, as with `--diploid` |
| `BUBBLE_ID`, `EVENT_NODES` | the site and the event's nodes; a `PATH_INVERSION`, typed from the walks, has no nodes |
| `CN_METHOD`, `REF_CN`, `RU_LEN` | how copies were counted, the reference's count, a unit copy's length; a `CNV`'s `POS`+1..`END` covers the reference's copies |
| `CN_FORMS` | on a `MODULE` record, the copies its carriers gained or lost, each an edit of the reference: `START-END:SEQ`, reference bases `START`..`END` replaced by `SEQ` in reference orientation. An insertion replaces no bases: its stretch is empty, `END` = `START`-1, and `SEQ` goes in before base `START`. A deletion has `SEQ` empty. `SEQ` names the reference stretches it repeats where that is shorter, `@A-B` (bases `A`..`B`) or `~A-B` (their reverse complement), joined with literal bases by `+`: a gained copy repeats the reference's copies, so its form stays short however long the copy. It is `*` where its literal bases exceed `--max-seq-bp`, the sequence then in `<prefix>.alt_seqs.fa`. Distinct copies are listed once, in reference order |
| `SEQ_FA` | the record's sequences in `<prefix>.alt_seqs.fa`, by field: `INSSEQ`, `DELSEQ`, `INVSEQ`, or `CN_FORMS.<k>` for the k-th form. Each is the FASTA record `<ID>.<field>`, so the file can be indexed (`samtools faidx`) and one allele fetched by name. A side file keeps every sequence, a novel insertion as much as a copy, while each line stays short |
| `SVLEN_RANGE`, `MERGE_MIN_SIM`, `MERGE_TEST` | how much the merged events differ: their sizes, and the weakest similarity of a member to the representative in the test that admits it (under the test's threshold: the member joined through others) |
| `GENES` | with `--gtf`: genes over the record; inside a copy of a `MODULE` site holding several genes, the paralogs its carriers hold |

`FORMAT` fields; the copy-number ones are `.` on records that are not `CNV`s:

| field | on | meaning |
|-------|----|---------|
| `GT` | every record | `1` carrier, `0` reference-like, `.` no allele there. With `--diploid`, the sample's haplotypes joined by `\|` (phased) or `/` |
| `CN` | `CNV` | the haplotype's copies of the unit (`UNIT`) or of the module (`MODULE`). On a `MODULE` record a haplotype that is no carrier, but whose change a written record of its own holds, reads `REF_CN`, as a dosage read off `CN` would otherwise count it twice |
| `CNBP` | `CNV` | bases gained (+) or lost (−) against the reference, from the haplotype's walk: over the site on a `UNIT` record; on a `MODULE` record, from its first module node to its last, with everything walked between |
| `CNRESID` | `CNV` | `CNBP` minus the copies' worth of bases, (`CN` − `REF_CN`) × the unit's length: the change whole copies do not explain. `0` means the whole size difference is copy number |
| `MODBP` | `MODULE` | the change over the module's own nodes, each visit counted. It decides `GT`: a carrier's `CN` differs from `REF_CN` and its `MODBP` reaches `--min-sv-bp`, unless the change beside its copies goes their way and is more than the copies (step 5) |
| `CNFORM` | `MODULE` naming copies | the carrier's copies as indices into `CN_FORMS`, joined by `+`; `.` for a non-carrier, or for a carrier with no copy written: none of its changes is a copy, or none could be written as an edit of one reference stretch, and such a change stays a record of its own. With `--diploid`, joined per haplotype as `GT` is |
| `CNINV` | `UNIT` | how many of the haplotype's copies run opposite to the reference's copy (with no reference copy, to the orientation most haplotypes' copies take). A gain with more such copies than the reference has is flagged `INVDUP`; with no reference copy, only a haplotype holding a copy and its reverse complement is |
| `HCN` | `CNV`, `--diploid` | each haplotype's `CN`, in `GT` order, joined as `GT` is; `CN` is their sum, and `CNBP`, `CNRESID`, `MODBP` and `CNINV` are sums too |

The VCF header describes every field.

`carrier_events.tsv`, one row per record and event; carriers holding the same event share a row:

| column | meaning |
|--------|---------|
| `variant_id`, `bubble_id` | the record and its site |
| `kind` | how the event was found: `seq` (base-level alignment), `del`, `ins`, `sub`, `inv` (node comparison), `cn` (a copy-number record); `unresolved` when it could not be placed on the reference. A `cn` event is the carrier's copy stretch: where the allele was aligned, its first copy to its last on either walk, cut where the alignment cuts the two walks (an insertion where the copies meet the rest keeps its side of them), so no base is both in the stretch and in another event of the carrier; node by node, up to the nodes the walk shares with the reference |
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

### Reading the VCF

One record of each kind, trimmed: `…` stands for left-out columns, fields and bases, and sample columns, named by their paths, are shortened to `sample#haplotype`.

**A deletion.**

```
#CHROM POS      ID               REF ALT   … INFO                                                                                       FORMAT HG02109#1 HG00097#1 HG02984#1
chr14  73618322 bubble9_DEL_9275 C   <DEL> … END=73618631;SVTYPE=DEL;SVLEN=-309;HOMLEN=1;CIPOS=0,1;CIEND=0,1;…;AN=465;AC=11;…;DELSEQ=AGAATTCTCT…AGAAAATAAT GT:CN  1:.       0:.       .:.
```

- `POS` is the base before the deletion: bases `POS`+1..`END` (73618323..73618631) are removed, and `DELSEQ` spells them.
- `HOMLEN=1`: the same deletion can be placed one base further right. `POS` is the leftmost placement, and `CIPOS`/`CIEND` give the range.

**Genotypes.** In the record above, HG02109#1 carries the deletion (`1`) and HG00097#1 is reference-like (`0`). HG02984#1 is `.`, no allele: its own 2,435-base deletion (`bubble9_DEL_9199`) removes these bases. A haplotype is also `.` where it does not cross the site, crosses an assembly gap there, or could not be aligned there. `.` haplotypes are left out of `AN`.

**A long insertion.**

```
#CHROM POS      ID                 REF ALT   … INFO                                                                      FORMAT HG03452#1 HG00096#1 HG00096#2
chr22  42123231 bubble4_INS_2229_2 G   <INS> … END=42123231;SVTYPE=INS;SVLEN=24316;HOMLEN=4;CIPOS=0,4;…;AN=126;AC=1;…;SEQ_FA=INSSEQ GT:CN  1:.       0:.       0:.
```

```
samtools faidx <prefix>.alt_seqs.fa bubble4_INS_2229_2.INSSEQ
```

- An insertion spans no reference bases: `END` = `POS`, and the `SVLEN` bases go in after `POS`.
- Its 24,316 bases are over `--max-seq-bp`, so `INSSEQ` is not in the line: `SEQ_FA=INSSEQ` says it is the FASTA record `<ID>.INSSEQ` in `<prefix>.alt_seqs.fa`, fetched by name as above.

**A replacement.**

```
#CHROM POS       ID               REF                    ALT   … INFO                                                               FORMAT HG00097#1 HG00097#2 CHM13#0
chr1   109644537 bubble3_DEL_2564 ATCCACCCGC…TCGTGATCCG AGCCA … END=109648773;SVTYPE=DEL;SUBSTITUTION;SVLEN=-4232;…;AN=465;AC=278;… GT:CN  1:.       0:.       1:.
```

- `REF` and `ALT` both start with the base at `POS`: the 4,236 bases after it (`POS`+1..`END`) are replaced by `GCCA`.
- `SVLEN` is the net change. It reaches `--min-sv-bp`, so the record is a `DEL` flagged `SUBSTITUTION`; under it, the record would be a `SUB`.
- Both sides fit `--max-seq-bp` and the breakpoints are exact, so `REF` and `ALT` spell them. Otherwise `ALT` is symbolic, with `DELSEQ` and `INSSEQ` in the line or in `SEQ_FA`.

**A copy number (`--cn`).**

```
#CHROM POS      ID               REF ALT   … INFO                                                                                                                                                                                                                                                  FORMAT                          HG00097#1   CHM13#0           HG005#1          HG00272#1
chr14  73515671 bubble5_CNV_4014 C   <CNV> … END=73516622;SVTYPE=CNV;…;AC=332;…;REF_CN=2;CN_UNIT_BP=475;CN_METHOD=MODULE;…;CN_FORMS=73515672-73516147:,73515672-73516148:,73515672-73516153:,73515672-73516156:,73515672-73516160:,73515673-73515672:@73515673-73516147,73515673-73516147:,73515895-73516370:A,73516118-73516593:T GT:CN:CNBP:CNRESID:MODBP:CNFORM 0:2:0:0:0:. 1:1:-476:-1:-476:1 1:1:-475:0:-475:7 1:3:475:0:475:6
```

- `POS`+1..`END` holds the reference's `REF_CN` copies, and `CN_UNIT_BP` is the length of one. `CN` is each haplotype's count: HG00097#1 holds 2 (`0`), CHM13#0 and HG005#1 lost one, HG00272#1 gained one. A `MODULE` record has no `SVLEN`, as its carriers gain and lose copies of different lengths: a carrier's change is its `CNBP`, and `CNRESID` is the part whole copies do not explain.
- `CN_FORMS` lists the copies gained or lost as edits of the reference, numbered from 1. Entry 1, `73515672-73516147:`, deletes those 476 bases. Entry 6, `73515673-73515672:@73515673-73516147`, is an insertion before base 73515673 of a copy of bases 73515673..73516147, written as their coordinates instead of 475 bases.
- `CNFORM` gives the carrier's entries, joined by `+` when it holds several: CHM13#0 lost entry 1, a base more than a copy (`CNRESID` -1); HG005#1 lost entry 7, exactly one copy; HG00272#1 gained entry 6. A non-carrier is `.`. A carrier's entries and its other records rebuild its walk.

**Records and `carrier_events.tsv`.**

```
#CHROM POS       ID               REF ALT   … INFO
chr1   109654087 bubble4_DEL_2916 A   <DEL> … END=109654190;SVTYPE=DEL;SVLEN=-103;HOMLEN=8;CIPOS=0,8;CIEND=0,8;SVLEN_RANGE=-103,-86;…;AN=465;AC=182;…
```

```
variant_id        bubble_id kind ref_start ref_end   cn walk_cn event_bp alt_bp alt       haplotypes
bubble4_DEL_2916  4         seq  109654096 109654198 .  .       -103     0      .         HG00099#1,…
bubble4_DEL_2916  4         seq  109654096 109654194 .  .       -99      0      .         HG00097#2,…
bubble4_DEL_2916  4         seq  109654096 109654190 .  .       -86      9      GAAAAGAAA HG02572#2,…
…
```

- A record holds one representative event, the one the most carriers hold: here a 103-base deletion carried by 114 of 182, written left-aligned; its row has it 8 bases further right, as `HOMLEN` allows.
- `carrier_events.tsv` holds each carrier's own event as found, one row per distinct event with its carriers: 8 rows here, whose sizes `SVLEN_RANGE` spans. An insertion's row has `ref_end` = `ref_start`-1, as for the long insertion above (42123232, 42123231).
- A copy-number carrier's row has `kind` `cn` and covers its whole copy stretch: HG00272#1 at the `--cn` record above replaces 73515672..73516622 (951 bases) with its 1,426, the change `CN_FORMS` entry 6 writes as one insertion.

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
