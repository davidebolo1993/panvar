# Module `transcripts`

CLI: `panvar transcripts`

## What it does

Projects a reference annotation onto every haplotype walk of a locus graph. Each transcript becomes one haplotype-specific transcript (HST) per haplotype and copy of its gene, spliced from that haplotype's own sequence.

- A haplotype without a copy of a gene gets no transcript of it. One with two copies gets two, numbered along the walk.
- A gene whose alleles the graph keeps on nodes of their own, as in hyperpolymorphic exons, is placed by its sequence instead of its nodes.
- A gene a deletion cuts below a copy is reported as a partial copy. Its transcripts with every exon left are still written.
- Exons inside a repeat unit of the gene follow the haplotype's own number of units.
- A transcript whose exon the haplotype lacks is reported and not written. The gene's other transcripts still are.
- An exon boundary inside an indel moves to the nearest base the haplotype shares with the reference, and is flagged.
- A haplotype's own annotation (`--hap-annotations`, such as HPRC's CAT or Liftoff) adds its transcripts, cut straight from the walk. Each HST carries its source: `projected`, `annotated`, or both.

## How it works

1. The GTF's exons are placed on the reference path through the `contig:start-end` in its name. Each exon boundary becomes a reference node and an offset inside it, since boundaries rarely fall between nodes. A transcript reaching past an end of the path keeps its exons on the path, cut there. A transcript with no exon on the path is left out.

2. A gene's copies are found with one transcript, its marker: the one tagged Ensembl_canonical, else the one with the most exonic bases. Its coding exons are used where it has them, because paralogs are best told apart by what they encode. A tandem stretch of the marker (step 15) is left out, as haplotypes carry its unit a different number of times.

3. A gene's stretch of the reference walk is its span with at least 100 bases of flank on each side. A boundary inside a variant across the gene's first or last base needs the bases the walk shares past that variant.

4. On each haplotype walk, read in both orientations, the nodes that occur once in the gene's stretch are chained in order. A chain pays for every base of indel between its links, so it stays within one copy instead of jumping to the next. Across a tandem stretch it pays nothing, as haplotypes carry any number of units there. Scores are exact integers, so equal chains tie the same way on every platform.

5. Where a walk holds a stretch of the gene twice, moving from one copy to the next costs the same at every node. A link with an indel inside an exon costs a little more, so the move falls in an intron or outside the gene, and no exon is made of two copies.

6. Each chain is completed into a pass. The gaps between its nodes are aligned node by node, and so is the walk past its ends, over a little more than the reference stretch left. That search stops at a node of a better chain, which is another copy of the stretch. Genes with a tandem stretch take no such stop, as their pass must reach every unit. A chain is completed only when its anchored bases inside the marker's span reach a tenth of the best chain's; the others are stray matches of single nodes. Bases outside the span are not counted, as long flanks would hide a short gene's other copies.

7. A pass that walks at least `--min-exon-share` of the marker's exonic bases is a candidate copy. Its stretch of the walk is the part matching the marker's span. Candidates are accepted best first, by that share. A candidate loses to an accepted pass of its own gene that it overlaps. It also loses to an accepted pass of a gene apart on the reference when they share more than half of the longer stretch. Two genes folded onto the same nodes, as the copies of a duplicated gene are, thus compete for each stretch, and the one the stretch resembles more takes it. Passes that tie on every measure keep the order they were found in.

8. A short candidate lying inside a longer accepted pass of a gene apart on the reference is dropped when that pass matches the stretch at least as well. The rest compete again. Such a stretch is the part of the other gene's copy that resembles this gene, such as a pseudogene fragment. A candidate that matches better is kept, as a copy left by a deletion between two paralogs lies inside whatever gene surrounds them.

9. A gene left without a copy on a walk, and not lost to another gene, is looked for by sequence. Node sharing fails where exons lie mostly on allele-specific nodes: every allele is then a node of its own, and a whole copy shares too few of the marker's bases. The gene's best chain below the share points to a window of the walk, with a tenth of the stretch and 100 bases more on each side. The best chain is the one matching most of the marker's bases, then most of the stretch, as introns and flanks alone still point to the gene. A window over three times the stretch is given up: the chain's nodes then disagree on where the gene lies.

   The gene's span, with 100 bases on each side, is aligned inside that window within a band of a quarter of its bases. Failing that, the marker's span with the same flanks is aligned at any distance. Such a copy is noted `placed_by_sequence:marker_span`, and transcripts with an exon past the aligned stretch are not written for it.

10. Unit-cost alignment can set a stretch the walk lacks against whatever follows it, with many short gaps. So an aligned base counts only where the 15 bases on each side match to 0.7 of their matches, mismatches and gaps. Gaps with fewer than 12 aligned bases between them are one indel and are set again in one piece. Among the marker's exonic bases, those aligned are the coverage, and those equal among them the identity.

    Alleles of one gene stay within a few percent of each other over their exons, and paralogs mostly do not, so the identity must reach `--min-identity`. Coverage of `--min-exon-share` or more is a copy, its transcripts noted `placed_by_sequence`. Coverage from 0.2 up to the share is a partial copy, reported and not written. Less is no copy. Placements are accepted best first, whole copies before partial ones. One over an accepted copy of a gene apart on the reference is that gene's, as for passes (step 8). Genes with a tandem stretch are placed by nodes only.

11. In a copy placed by sequence, an exon the whole alignment leaves without bases is looked for alone. An intron that differs by kilobases between alleles can pull the alignment off a short exon at its end, though the exon is there. The exon is aligned between the exons placed around it, or within three times the reference's distance from its one placed neighbour, at `--min-identity`. One found that way is noted `aligned_alone`.

12. A gene still without a copy can hold whole transcripts: a deletion that takes the marker's exons leaves those that do not use them. Each transcript is projected through the gene's best chain, and one with every exon in place is written. The walk then has one copy of the gene below the share, noted `below_share`, and its other transcripts are `missing_exon`. An exon counts there only when `--min-exon-share` of its bases lie on nodes the reference walks nowhere outside the gene's span, as the chain may lie on a paralog. Below that, the exon is missing, or noted `shared_nodes` where other nodes bring it to the share.

13. Within each copy, an exon boundary is carried to the walk through the node holding it, at the same offset. Where the walk misses that node, the boundary lies inside a variant, and the walk's unmatched stretch there is compared with the reference's. Of the same length, it is a substitution (clustered SNPs the graph folds into one node): every base keeps its place, so the boundary keeps its offset and is noted `in_variant`. Otherwise the boundary moves to the nearest position the walk shares with the reference, and a tie keeps the variant inside the exon. In a sequence placement, a boundary at an aligned base keeps its place, and one at a base the walk lacks moves to the nearer end of that deletion.

    An exon sharing no base with the walk is missing. So is one whose matched nodes lie further apart on the walk than on the reference by more than the exon's length, with nodes of the gene between them: the match moved from one copy of a duplicated stretch to another inside it (`not_contiguous`).

14. One stretch of a walk can be two stretches of the reference. Where a segment lies twice inside a gene and the graph folds its copies onto the same nodes, a deletion from one copy into the next leaves one copy with alleles of both. The pass follows one reference copy, so the other copy's exons would be missing though the walk holds them. A second chain of the gene's nodes, on reference steps the pass leaves out and inside the pass on the walk, is taken as that other copy. An exon the pass lacks between that chain's nodes is carried through it (`merged_copy`) when it then lies in its place among the transcript's exons.

15. Exons that repeat with a tandem unit start at the same node and offset in each unit the reference walks. A transcript's exons that share a start mark its tandem stretch, when the units between those starts are of comparable length (the longest at most twice the shortest) and each holds its exons whole. Otherwise the equal starts are a node walked twice for another reason, and the exons are left to their nodes.

    The stretch is placed unit by unit, by sequence. On the reference, the window between the exons around the stretch is cut into a head, its units and a tail. On a haplotype, the head and the tail are aligned at the two ends of its window, and units tile the rest. Each next unit starts where the walk best holds a reference unit's first bases, or ends where it best holds a unit's last bases, so an indel near one end leaves the boundary found from the other.

    Each unit then takes the reference unit sharing most 15-mers with it. Ties go to the unit of the same rank when the counts agree, else to the nearest in relative position. That unit's exons, perhaps none, are lifted through their alignment, which keeps each indel in one piece, so an exon with an insertion holds it.

16. The exons are cut from the walk with the haplotype's alleles and joined in transcript order; a minus-strand transcript is reverse-complemented. The two splice-site bases at each end of every intron are compared with the reference's. An intron the haplotype shortens below those four bases is lost, and its exons are joined.

17. The reference walk is projected first and must give back every transcript as annotated. A transcript it does not give back is projected wrongly there, so its records on every walk are `not_as_annotated`, and a warning names it. A record with a more severe status (`missing_exon`, `partial_copy`, `no_copy`) keeps it, and `detail` notes `not_as_annotated`.

18. A haplotype's own annotation is in its assembly's coordinates, and the haplotype's path names its contig and interval (`sample#hap#contig:start-end`, the start 0-based as pggb writes it). Its transcripts with an exon on the path are cut from the walk at those coordinates, strand-aware. A model with coding rows and no exon rows takes its coding rows as exons. Exon rows naming a transcript the annotation lacks are skipped, with a warning counting them. A transcript reaching past the path is cut at its ends (`clipped_at_edge`). Each gene record of the annotation is one copy. The annotations of different haplotypes are read in parallel.

19. An annotated transcript's gene is found by its gene id, the source gene's where the annotation names one. It is found by its name only where the GTF does not know that id, so a gene the GTF has elsewhere is not taken for a projected gene of the same name. Gene and transcript ids are compared without an Ensembl version, so `ENST00000000001.2` in a GENCODE GTF is `ENST00000000001` in an annotation. Each annotated copy is matched to the projected copy of its gene it overlaps most on the walk. An annotated transcript with the projected one's id and sequence on that copy is that record, of both sources (`projected+annotated`). Any other is written beside the projected records.

20. Where the two sources differ, these rules apply.

    - A clean projected record (`ok`, `splice_site_changed`, `boundary_moved`, `clipped_at_edge`) wins, and the annotated one is `discordant`. The projection carries the reference's exons base by base, while an annotation lifted by alignment can shift a boundary or drop an exon.
    - Where the projection fails on a transcript (a missing exon, a partial copy, or no copy), the annotation fills it (`annotation_fill`) when it has the reference transcript's exon count. With another count it is another model, and is `discordant`.
    - An annotated copy the projection lacks that lies on another gene's projected copy is a `paralog_lift`, and is no copy. The two genes must be apart on the reference, and over half of its exonic bases must lie in that copy's stretch. A gene the GTF has off the path is checked against every gene's copy, as its place on the path is not known.
    - Where the projected record is itself in doubt (`intron_lost`, `not_as_annotated`), both stay.
    - A gene the GTF lacks entirely, by id and by name, is never a paralog lift, as its place on the reference cannot be known. Its records are written and indexed as annotated, as are novel isoforms of known genes. `recruit` does not count such a gene as the locus's own, so its transcripts elsewhere in a transcriptome still label its junctions shared.

    `discordant` and `paralog_lift` records are written but left out of a recruit index. With `--prefer-annotation` the annotation wins instead: an annotated transcript differing from a projected one makes that record `discordant`, and no annotated record is `discordant` or a `paralog_lift` (a lift is still noted). A gene's copies are the projected ones, followed along the walk by the annotated copies the projection lacks that hold a kept record.

## Required inputs

- `-i, --gfa <graph.gfa>`: a locus GFA (plain or gzipped) with one P or W line per haplotype
- `--gtf <path>`: the annotation in reference coordinates (Ensembl or GENCODE, plain or gzipped); its `exon` rows give the exons, its `CDS` rows the coding parts
- `-r, --reference-path <name>`, with `--gtf`: the reference path; a unique part of the name is enough. Its name must carry `contig:start-end`, as for `bubble --gtf`
- or `--hap-annotations <table.tsv>` alone, for annotated transcripts only; both together give both sources
- `-o, --out-prefix <prefix>`: output prefix

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `--hap-annotations <tsv>` | haplotypes' own annotations: one row per haplotype, `haplotype` (`sample#hap`, as in the path names) and `path` (its GFF3 or GTF, plain or gzipped); a header row naming `haplotype` is skipped | — |
| `--prefer-annotation` | with `--hap-annotations`: annotated transcripts win where they differ from projected ones, and none is left out of a recruit index | off |
| `--contig-alias <tsv>` | contig names of the annotations matched to the paths' contigs: each row lists the names of one contig (a chromAlias table works); without it a name matches as it is, or once a `sample#hap#` prefix is removed | — |
| `--min-exon-share <F>` | the share of the marker's exonic bases a pass must walk, or a sequence placement must cover, to be a copy; a lower value also counts partial copies. On a copy below the share, the share of each exon that must lie on the gene's own nodes | `0.8` |
| `--min-identity <F>` | the identity of the marker's exonic bases a stretch must reach to be placed by sequence, as a copy or a partial copy; also that of an exon looked for alone | `0.93` |
| `-t, --threads <N>` | worker threads, one haplotype or annotation each (`0` = all cores) | `0` |
| `-q, --quiet` | no progress logs | off |

## Outputs

| file | contents |
|------|----------|
| `<prefix>.hst.fa` | one sequence per HST, 5' to 3', named `transcript|gene|haplotype|copy`, with `source=` and `status=` after the name. An annotated transcript written beside the projected ones adds `|annotated:<record>`, so two records of one source transcript keep apart |
| `<prefix>.hst.tsv` | one row per transcript, haplotype and copy; a row with copy 0 where the haplotype has no copy |
| `<prefix>.hst.summary.tsv` | one row per gene: copies per haplotype and rows per status |

`hst.tsv`, in the file's order (`.` where a value does not apply):

| column | meaning |
|--------|---------|
| `transcript_id`, `transcript_name`, `gene_id`, `gene_name`, `biotype`, `strand` | the transcript, from the GTF or the annotation |
| `haplotype` | the path |
| `copy`, `copies` | this copy, numbered along the walk as written, and the gene's copies on the walk (with haplotypes' annotations, the copies they fill too) |
| `source` | `projected`, `annotated` or `projected+annotated` (both, the same sequence) |
| `record` | the annotation's own id of the transcript (GFF3 `ID`, GTF `transcript_id`; several where annotated transcripts join one projected record); `.` for a projected one |
| `status` | below |
| `orientation` | `+` where the walk passes the gene in the reference's orientation, `-` where it passes it reversed; `.` for an annotated record or no copy |
| `path_start`, `path_end` | the HST's first and last exonic base on the walk as written, 0-based, end excluded |
| `length`, `exons` | the HST's length and exons |
| `ref_length`, `ref_exons` | the same on the reference path (the part on the path, for a transcript reaching past it); for an annotated record, those of the projected transcript of its id |
| `same_as_ref` | whether the HST is the reference transcript's sequence |
| `moved_boundaries` | exon boundaries moved off a variant |
| `changed_splice_sites` | intron ends whose two splice-site bases differ from the reference's |
| `exon_starts`, `exon_lengths` | each exon's first base on the walk (as `path_start`) and its length, in transcript order |
| `detail` | notes, `;`-joined (below) |

Exons and introns are numbered in transcript order, 1-based, over the whole annotated transcript. A transcript whose tandem exons follow the haplotype is numbered over its own exons.

| note in `detail` | meaning |
|------------------|---------|
| `e3.start:-3`, `e3.end:+2` | exon 3's 5' or 3' boundary moved off a variant; the exon is that many reference bases longer (+) or shorter (−) |
| `e5.end:in_variant` | exon 5's boundary lies in a substitution and kept its offset |
| `e2:missing` | the copy lacks exon 2 |
| `e2:not_contiguous` | exon 2's matched nodes lie in two copies of a duplicated stretch (step 13) |
| `e1:walk_end` | the walk ends before exon 1 |
| `e1:not_aligned` | exon 1 lies outside the stretch aligned for a copy placed by sequence |
| `e1:aligned_alone` | exon 1 was found by aligning it alone (step 11) |
| `e1:merged_copy` | exon 1 was found on the reference's other copy of a segment (step 14) |
| `e4:overlaps_neighbour` | exon 4 would overlap the exon before it on the walk |
| `i1.donor:GT>AT`, `i1.acceptor:AG>AC` | a changed splice site, on the transcript's strand |
| `i2:intron_lost` | intron 2 is shortened below its splice sites |
| `clipped:3':200` | annotated exonic bases past the path at the transcript's 3' end (or 5') |
| `tandem:14>46` | the tandem stretch's exons on the reference and on the haplotype |
| `tandem_resized:1` | tandem exons holding an indel |
| `tandem_unplaced:1` | tandem exons a unit's alignment leaves without bases |
| `tandem:missing`, `tandem:shared_nodes` | no tandem exon placed, or one placed on nodes the reference walks outside the gene |
| `rival:<gene>:<share>><share>` | the paralog this copy was taken from, and both shares |
| `placed_by_sequence`, `placed_by_sequence:marker_span` | a copy placed by sequence, over the gene's span or the marker's |
| `best_share:<F>` | the gene's best share by nodes on this walk |
| `identity:<F>`, `coverage:<F>` | the marker's identity and coverage in a sequence placement |
| `below_share` | a copy counted for its whole transcripts, the marker below the share (step 12) |
| `e2:shared_nodes` | on such a copy, exon 2 walked on nodes the reference walks outside the gene too |
| `taken_by:<gene>`, `inside:<gene>` | no copy: the stretch went to that gene's copy, or lies inside it |
| `annotation_fill` | an annotated transcript filling what the projection fails on |
| `on:<gene>` | an annotated transcript lying on that gene's copy |
| `exons:<n>/<m>` | an annotated transcript's exons against the reference transcript's |
| `annotated:<record>` | a projected record the annotation overrules |
| `not_as_annotated` | the reference walk does not give the transcript back, on a record whose status is more severe |

Statuses, most severe first:

| status | meaning |
|--------|---------|
| `no_copy` | the walk has no copy of the gene, nor a whole transcript of it |
| `partial_copy` | the walk holds part of the gene: its marker covered below the share, at `--min-identity` or more, with no whole transcript. No sequence is written and the copy count stays 0 |
| `missing_exon` | the copy lacks an exon of this transcript; no sequence is written |
| `paralog_lift` | an annotated transcript on another gene's copy (with `--hap-annotations`) |
| `discordant` | a record the other source contradicts (with `--hap-annotations`) |
| `not_as_annotated` | the reference walk does not give the transcript back as annotated, so its projection is not trusted on any walk; a record with a more severe status keeps it and notes `not_as_annotated` in `detail` |
| `intron_lost` | an intron shortened below its two splice sites, its exons joined |
| `clipped_at_edge` | the transcript reaches past the reference path, and only its part on the path is projected |
| `boundary_moved` | an exon boundary lies inside an indel and was moved |
| `splice_site_changed` | a splice site differs from the reference's |
| `ok` | every boundary carried exactly |

`hst.summary.tsv` has one row per gene, in reference order. `transcripts` counts the projected ones, `clipped` those reaching past the reference path, and `outside` those with no exon on the path. `haplotypes_by_copies` reads `0:2,1:120,2:9` for 2 haplotypes with no copy, 120 with one and 9 with two. Then come `hsts` written, `same_as_ref`, `annotated` (HSTs from the haplotypes' annotations) and the rows of each status (`discordant` and `paralog_lift` only with `--hap-annotations`). `copies_differ` lists, for a gene some haplotype's annotation has, each annotated haplotype whose own copies differ from the run's, as `sample#hap:run/annotated`.

## Limitations

- A copy is found through the nodes it shares with the reference copy, or by sequence where a chain of shared nodes points to it. A copy no chain reaches is not found, nor is a second copy on a walk where the first was found by nodes.
- Placement by sequence takes the window of the gene's best chain only. Where that chain lies on a paralog, the gene's own copy elsewhere on the walk is not looked at.
- No identity separates every allele from every paralog. A paralog the reference does not annotate can reach `--min-identity` and be placed as the gene: a gene the reference haplotype lacks, but another lineage carries in its place, is one such case. `identity` in `detail` shows how close each placement was. A higher `--min-identity` keeps such stretches out, at the cost of the gene's most distant alleles.
- An exon looked for alone takes the first match in its window at `--min-identity`, so a short exon can match a stretch of similar sequence there.
- Whole transcripts of a gene without a copy are looked for on its best chain only, and counted as one copy. A transcript with an exon on nodes a paralog shares is not written there.
- An exon is found on a merged copy only between two nodes of its own reference copy that the walk passes in order.
- A short exon between diverged introns can fall where the alignment does not hold, and be missed.
- A haplotype's tandem unit cut short still takes the layout of the reference unit it resembles most; exons outside what it holds are counted in `tandem_unplaced`. Among identical units, the layout follows their order.
- An exon whose bases lie in two copies of a duplicated stretch is not written.
- A copy mixing the differences of two near-identical paralogs goes to the one it shares more coding bases with; `rival` gives the margin.
- Only the reference annotation is projected, so isoforms and genes the reference lacks come only from a haplotype's own annotation. A changed splice site is reported, not judged.
- An annotated transcript joins a projected one only by the same transcript id (for CAT, its source transcript) on the copy it overlaps. An annotation with other ids, such as Liftoff's RefSeq ids, stays beside the projected transcripts, unjudged.
- A filling transcript is judged by its exon count and where it lies, not by its boundaries.
- A gene the GTF lacks is never a `paralog_lift`, so one an annotation lifted onto a paralog's copy is kept as a copy of its own. A GTF of the whole genome, rather than of the loci, lets lifts of genes from elsewhere be judged.
- With `--prefer-annotation`, a truncated annotated model replaces the projected transcript it differs from.

## Example

```bash
# every transcript on every haplotype
panvar transcripts -i locus.gfa.gz -r GRCh38 --gtf genes.gtf.gz -o locus
# count partial copies too
panvar transcripts -i locus.gfa.gz -r GRCh38 --gtf genes.gtf.gz -o locus --min-exon-share 0.5 -t 2
# with the haplotypes' own annotations beside the projected ones
panvar transcripts -i locus.gfa.gz -r GRCh38 --gtf genes.gtf.gz --hap-annotations hap_annotations.tsv \
  --contig-alias chrom_alias.tsv -o locus
# the same, the annotations winning where they differ
panvar transcripts -i locus.gfa.gz -r GRCh38 --gtf genes.gtf.gz --hap-annotations hap_annotations.tsv \
  --prefer-annotation -o locus
# the haplotypes' own annotations alone
panvar transcripts -i locus.gfa.gz --hap-annotations hap_annotations.tsv -o locus
```

`scripts/transcripts_examples.sh` projects the loci GTF onto the paralog graphs in `tests/real_data` and prints each gene's copies per haplotype, in `tests/examples/transcripts`.
