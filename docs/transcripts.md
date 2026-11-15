# Module `transcripts`

CLI: `panvar transcripts`

## What it does

Projects a reference annotation onto every haplotype walk of a locus graph. Each transcript becomes one haplotype-specific transcript (HST) per haplotype and copy of its gene, spliced from that haplotype's own sequence.

- A haplotype without a copy of a gene gets no transcript of it; one with two copies gets two, numbered along the walk.
- A gene whose alleles the graph keeps on nodes of their own, as in hyperpolymorphic exons, is placed by its sequence instead of its nodes.
- A gene a deletion cuts below a copy is reported as a partial copy, apart from a haplotype that lacks the gene; its transcripts with every exon left are still written.
- Exons inside a repeat unit of the gene follow the haplotype's own number of units.
- A transcript whose exon the haplotype lacks is reported and not written, so other transcripts of the gene still are.
- Where a deletion joins two copies of a segment inside a gene into one, an exon of either copy is found on it.
- An exon boundary that falls inside an indel is moved to the nearest base the haplotype shares with the reference, and flagged; one inside a substitution of the reference's length keeps its offset, and is noted.
- A haplotype's own annotation (`--hap-annotations`, e.g. HPRC's CAT or Liftoff) adds its transcripts, cut straight from the walk: isoforms and copies the reference annotation lacks. Each HST carries its source: `projected`, `annotated`, or both. Where the two disagree, a clean projection wins and the annotation fills what the projection fails on; an annotated copy lying on a paralog's copy is flagged.

## How it works

1. The GTF's exons are placed on the reference path through the `contig:start-end` in its name. Each exon boundary becomes a reference node and an offset inside it, since boundaries rarely fall between nodes. A transcript reaching past either end of the path keeps its exons on the path, cut at the path's end, so a gene at the edge of a locus still gets transcripts; one with no exon on the path is left out.

2. Each gene's copies are found with one transcript, its marker: the one tagged Ensembl_canonical, else the one with the most exonic bases. Its coding exons are used where it has them, because paralogs are best told apart by what they encode.

3. On each haplotype walk, read in both orientations, the nodes that occur once in the gene's stretch of the reference walk are chained in order. A chain pays for every base of indel between its links, so it stays within one copy instead of jumping to the next; across a tandem stretch (step 8) it pays nothing, as haplotypes carry any number of units there. The gaps between chained nodes are then aligned node by node. Each chain holding a tenth of the best chain's bases is one pass of the walk through the gene; the rest are stray matches of single nodes. A pass that shares at least `--min-exon-share` of the marker's exonic bases with the reference is a candidate copy, and its stretch of the walk is the part matching the marker's span.

4. Candidates are accepted best first, by that share. A candidate loses to an accepted pass of its own gene that it overlaps. It also loses to an accepted pass of a gene apart on the reference (their markers do not overlap) that takes the same stretch of the walk, more than half of the longer of the two. Two genes folded onto the same nodes, as the copies of a duplicated gene are, thus compete for each stretch, and the one the stretch resembles more takes it.

5. A short candidate lying inside a longer accepted pass of a gene apart on the reference is dropped when that pass matches the stretch at least as well as the candidate does, and the rest compete again. Such a stretch is part of the other gene's copy that resembles this gene: a pseudogene fragment, or a gene nested in the other copy of a duplicated gene. A candidate that matches better is kept, as a copy left by a deletion between two paralogs lies inside whatever gene surrounds them.

6. A gene left without a copy on a walk, and not lost to another gene, is looked for by sequence. Node sharing fails where exons lie mostly on allele-specific nodes: every allele is then a node of its own, and a full-length copy shares too few of the marker's bases.
   - Where: the window of the walk that the gene's best chain below the share points to (the one matching most of the marker's bases, then most of the gene's stretch: where every exon is an allele of its own, the introns and flanks alone point to the gene), placed from the chain's anchors nearest the stretch's ends, with a tenth of the stretch and 100 bases more on each side. A window over three times the stretch is given up: the anchors then disagree on where the gene lies.
   - What is aligned: the gene's span of the reference with 100 bases on each side (an indel across an end of the gene is then a gap inside the alignment), within a band of a quarter of its bases, where a copy lies. Failing that (a deletion, or a span over about 45 kb), the marker's span with the same 100 bases on each side, at any distance (up to about 20 kb), which is enough to tell a partial copy. A copy found that way is noted `placed_by_sequence:marker_span`, and transcripts with an exon reaching past the interval aligned (the marker's span and those 100 bases, cut at the reference path's ends) are not written for it (`e1:not_aligned`).
   - Which bases count as aligned: those where the alignment holds, that is where, over the 15 bases on each side, the matches are 0.7 or more of the matches, mismatches and gaps (a gap of any length counting once): unit-cost alignment otherwise sets a stretch the walk lacks against whatever follows it, with many short gaps. It also scatters an indel's bases among others that happen to match, so gaps with under 12 aligned bases between them are set again as one indel: the shorter side against the longer's left or right end, whichever matches more bases.
   - The measure: among the marker's exonic bases, those aligned are the coverage, and those equal among them the identity. Alleles of one gene, hyperpolymorphic ones too, stay within a few percent of each other over their exons, and paralogs mostly do not, so the identity must reach `--min-identity`. Coverage of `--min-exon-share` or more is then a copy, its transcripts noted `placed_by_sequence`; coverage from 0.2 up to the share a partial copy (`partial_copy`), reported and not written; less, or a lower identity, no copy.
   - Boundaries: an exon boundary at an aligned base keeps its place; one at a base the walk lacks moves to the nearer end of that deletion, flagged as in step 7.
   - Competition: placements are accepted best first, whole copies before partial ones. One over the stretch of an accepted copy of a gene apart on the reference is that gene's; so is one whose chain lies inside a longer accepted pass of such a gene that matches the stretch at least as well (step 5).

   Genes with a tandem stretch are placed by nodes only, as their units are.

   A gene still without a copy, and not lost to another gene, can hold whole transcripts: a deletion that takes the marker's exons leaves those that do not use them. Each transcript is projected through the gene's best chain as through a pass (step 7), and one with every exon in place is written. The walk then has one copy of the gene, below the share: its rows are noted `below_share` with the chain's `best_share`, and the transcripts lacking an exon are `missing_exon`. An exon counts there only when the share of its bases on nodes the reference walks nowhere outside the gene's span (not the nodes around it, which may be a neighbour's) reaches `--min-exon-share`, the copy's test taken exon by exon: below it the exon is missing, or noted `e2:shared_nodes` where other nodes bring it to the share. An exon of a tandem stretch (step 8) is placed by sequence, in a window the chain sets, so the same test is taken on the walk's nodes under it (`tandem:shared_nodes` where it fails). A chain below the share can lie on a paralog, and an exon on nodes the paralog shares is not told from the paralog's own, so a gene deleted beside a paralog stays without a copy; and a stretch diverged from the gene is still judged by its sequence alone.

7. Within each copy, an exon boundary is carried to the walk through the node holding it, at the same offset. When the walk does not pass that node (the boundary lies inside a variant), the walk's unmatched stretch there is compared with the reference's. Of the same length, it is a substitution (clustered SNPs the graph folds into one node): every base keeps its place, so the boundary keeps its offset and is noted `in_variant`; the splice sites are checked as elsewhere. Otherwise the boundary moves to the nearest position the walk shares with the reference, and a tie keeps the variant inside the exon. An exon sharing no base with the walk is missing.

   One stretch of a walk can be two stretches of the reference. Where a segment lies twice inside a gene and the graph folds its copies onto the same nodes, a deletion from one copy into the next leaves one copy with alleles of both. The pass follows one reference copy, so the exons of the other would be missing though the walk holds them. A second chain of the gene's nodes, on reference steps the pass leaves out and inside the pass on the walk, is taken as that other copy: the walk passes nodes that only it has. The steps between that chain's anchors are aligned node by node, and an exon the pass lacks, lying between them, is carried through that alignment (`e1:merged_copy`) when it then lies in its place among the transcript's exons on the walk. The same bases can so be an exon of two transcripts, one from each reference copy. A walk's copy without such nodes is the copy the pass follows alone, and the other copy's exons are missing.

8. Exons that repeat with a tandem unit inside the gene start at the same node and offset in each unit the reference walks, so a transcript's exons starting where another one starts mark its tandem stretch, when the units between those starts are of comparable length (the longest at most twice the shortest) and each holds its exons whole; otherwise the equal starts are a node walked twice for another reason, and the exons are left to their nodes. The stretch is placed unit by unit, by sequence. On the reference, the window between the exons around the stretch is cut into a head, its units (each starting where the walk passes the stretch's first exon start, with or without exons of its own) and a tail. On a haplotype, the head and the tail are aligned at the two ends of its window. The units tile the rest: each next unit starts where the walk best holds a reference unit's first bases, or ends where it best holds a unit's last bases, so an indel near one end of a unit leaves the boundary found from the other. Each unit then takes the reference unit sharing most 15-mers with it (ties: the unit of the same rank when the counts agree, else the one nearest in relative position), and that unit's annotated exons, perhaps none, are lifted through their alignment. The alignment counts only in runs of matched bases, with each indel kept in one piece, so an exon with an insertion holds it. The reference walk gives back its annotation; a transcript it does not is reported. The marker's tandem stretch does not count towards its share, as haplotypes carry the unit a different number of times.

9. The exons are cut from the walk with the haplotype's alleles and joined in transcript order; a minus-strand transcript is reverse-complemented. The two splice-site bases at each end of every intron are compared with the reference's; an intron the haplotype shortens below those four bases is lost, and its exons are joined.

10. A haplotype's own annotation is in its assembly's coordinates, and the haplotype's path names its contig and interval (`sample#hap#contig:start-end`, the start 0-based as pggb writes it). Its transcripts with an exon on the path are cut from the walk at those coordinates, strand-aware (a model with coding rows and no exon rows takes its coding rows as exons; exon rows naming a transcript the annotation lacks are skipped, with a warning counting them); one reaching past the path is cut at its ends (`clipped_at_edge`). Each gene record of the annotation is one copy.

11. The annotated transcripts are merged with the projected ones, copy by copy. An annotated transcript's gene is found by its gene id (the source gene's, where the annotation names one; an Ensembl id without its version), and by its name only where the GTF does not know that id, so a gene the GTF has elsewhere is not taken for a projected gene of the same name. Each annotated copy is matched to the projected copy of its gene it overlaps most on the walk. An annotated transcript with the same id and sequence as the projected one of that copy is that record, of both sources (`projected+annotated`). Any other is written beside the projected ones, under these rules:

    - A clean projected record (`ok`, `splice_site_changed`, `boundary_moved`, `clipped_at_edge`) wins, and the annotated one differing from it is `discordant`: the projection carries the reference's exons base by base, while an annotation lifted by alignment can shift a boundary or drop an exon.
    - Where the projection fails on a transcript (a missing exon, a partial copy, or a copy it lacks), the annotation fills it (`annotation_fill`), when it has the reference transcript's exon count. With another count it is another model, not the missing transcript, and is `discordant`.
    - An annotated copy the projection lacks, lying on another gene's projected copy (over half its exonic bases in that copy's stretch, the two genes apart on the reference), is a `paralog_lift`: the annotation put the gene on its paralog. It is no copy. A filling transcript lying on another gene's copy is one too. A gene the GTF has off the path is checked against every gene's copy, as its place on the path is not known.
    - Where the projected record is itself in doubt (`intron_lost`, `not_as_annotated`), both stay.
    - A gene the GTF lacks entirely (by id and by name), such as a novel gene or one of another release, is never judged a paralog lift: where it belongs on the reference cannot be known, and letting its sequence into a recruit index at its own locus costs little. It is written, counted and indexed as annotated, as is a novel isoform of a known gene.

    `discordant` and `paralog_lift` records are written but left out of a recruit index. With `--prefer-annotation` the annotation wins instead: an annotated transcript differing from a projected one leaves that record `discordant`, and no annotated record is `discordant` or a `paralog_lift` (a lift is still noted in `detail`). A gene's copies are the projected ones, followed along the walk by the annotated copies the projection lacks that hold a kept record.

## Required inputs

- `-i, --gfa <graph.gfa>`: a locus GFA (plain or gzipped) with one P or W line per haplotype
- `--gtf <path>`: the annotation in reference coordinates (Ensembl or GENCODE, plain or gzipped); its `exon` rows give the exons, its `CDS` rows the coding parts
- `-r, --reference-path <name>`, with `--gtf`: the reference path; a unique part of the name is enough. Its name must carry `contig:start-end`, as for `bubble --gtf`
- or `--hap-annotations <table.tsv>` alone, for annotated transcripts only (below); both together give both sources
- `-o, --out-prefix <prefix>`: output prefix

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `--hap-annotations <tsv>` | haplotypes' own annotations: one row per haplotype, `haplotype` (`sample#hap`, as in the path names) and `path` (its GFF3 or GTF, plain or gzipped); a header row naming `haplotype` is skipped | — |
| `--prefer-annotation` | with `--hap-annotations`: annotated transcripts win where they differ from projected ones, and none is left out of a recruit index | off |
| `--contig-alias <tsv>` | contig names of the annotations matched to the paths' contigs: each row lists names of one contig (a chromAlias table works); without it a name matches as it is, or once a `sample#hap#` prefix is removed | — |
| `--min-exon-share <F>` | the share of the marker's exonic bases a pass must walk, or a sequence placement must cover, to count as a copy; a lower value also counts partial copies. Below it, the share of each exon of a transcript that must lie on the gene's own nodes for the transcript to be written | `0.8` |
| `--min-identity <F>` | the identity of the marker's exonic bases a stretch must reach to be placed by sequence, as a copy or a partial copy | `0.93` |
| `-t, --threads <N>` | worker threads, one haplotype each (`0` = all cores) | `0` |
| `-q, --quiet` | no progress logs | off |

## Outputs

| file | contents |
|------|----------|
| `<prefix>.hst.fa` | one sequence per HST, 5' to 3', named `transcript|gene|haplotype|copy` (`|annotated:<record>` added to an annotated transcript written beside the projected ones, with the annotation's own id, so two records of one source transcript keep apart), with `source=` and `status=` after the name |
| `<prefix>.hst.tsv` | one row per transcript, haplotype and copy, a row with copy 0 where the haplotype has no copy |
| `<prefix>.hst.summary.tsv` | one row per gene: copies per haplotype and rows per status |

`hst.tsv` (`.` where a value does not apply):

| column | meaning |
|--------|---------|
| `transcript_id`, `transcript_name`, `gene_id`, `gene_name`, `biotype`, `strand` | the transcript, from the GTF |
| `haplotype` | the path |
| `copy`, `copies` | this copy, numbered along the walk as written, and the gene's copies on the walk (with haplotypes' annotations, the copies they fill too) |
| `source` | `projected`, `annotated` (from the haplotype's annotation; `orientation` is `.`, and `ref_length`, `ref_exons` are the projected transcript's of that id, if any) or `projected+annotated` (both, the same sequence) |
| `status` | below |
| `record` | the annotation's own id of the transcript (GFF3 `ID`, GTF `transcript_id`; several where annotated transcripts join one projected record), `.` for a projected one |
| `orientation` | `+` where the walk passes the gene in the reference's orientation, `-` where it passes it reversed |
| `path_start`, `path_end` | the HST's first and last exonic base on the walk as written, 0-based, end excluded |
| `length`, `exons` | the HST's length and exons |
| `ref_length`, `ref_exons` | the same on the reference path (the part on the path, for a transcript reaching past it) |
| `same_as_ref` | whether the HST is the reference transcript's sequence |
| `moved_boundaries` | exon boundaries moved off a variant |
| `changed_splice_sites` | intron ends whose two splice-site bases differ from the reference's |
| `exon_starts`, `exon_lengths` | each exon's first base on the walk (as `path_start`) and length, in transcript order |
| `detail` | notes, `;`-joined: `e3.start:-3` (exon 3's 5' boundary moved, the exon 3 reference bases shorter), `e5.end:in_variant` (exon 5's 3' boundary lies in a substitution and kept its offset), `e2:missing`, `e1:walk_end` (the walk ends before it), `e1:merged_copy` (found where the walk's one copy of a segment is the reference's second copy too), `i1.donor:GT>AT` (a changed splice site, transcript strand), `i2:intron_lost`, `clipped:3':200` (annotated exonic bases past the path at the transcript's 3' end), `tandem:14>46` (the exons of the tandem stretch on the reference and on the haplotype; the transcript's exons and introns are then numbered over its own), `tandem_resized:1` (tandem exons with an indel), `tandem_unplaced:1` (tandem exons a unit's alignment leaves without bases), `rival:<gene>:<share>><share>` (the paralog this copy was taken from, and both shares), `placed_by_sequence` with `best_share:<F>` and `identity:<F>` (a copy placed by sequence: its best share by nodes, and the marker's identity), `below_share` with `best_share:<F>` (a copy counted for its whole transcripts, the marker below the share), `e2:shared_nodes` (on such a copy, an exon walked on nodes the reference walks outside the gene too; `tandem:shared_nodes`, a tandem exon placed on such nodes), `e1:not_aligned` (an exon outside the stretch aligned for such a copy), for a partial copy `best_share:<F>`, `coverage:<F>` and `identity:<F>`, for no copy `best_share:<F>` and `taken_by:<gene>` or `inside:<gene>`, and for an annotated transcript `annotation_fill` (it fills what the projection fails on), `on:<gene>` (it lies on that gene's copy), `exons:<n>/<m>` (its exons against the reference transcript's) or, on a projected record it overrules, `annotated:<record>` |

Statuses, most severe first:

- `no_copy`: the walk has no copy of the gene, nor a whole transcript of it;
- `partial_copy`: the walk holds part of the gene, less than a copy: its marker's exons are covered below the share, at `--min-identity` or more, and no transcript is whole; no sequence is written and the copy count stays 0;
- `missing_exon`: the copy lacks an exon of this transcript; no sequence is written;
- `paralog_lift`: an annotated transcript lying on another gene's copy, the gene put on its paralog (with `--hap-annotations`);
- `discordant`: a record the other source contradicts (with `--hap-annotations`; step 11);
- `not_as_annotated`: on the reference walk only, a transcript whose projection there differs from its annotation (a warning names them; the projection of that transcript is not to be trusted on any haplotype);
- `intron_lost`: the haplotype shortens an intron below its two splice sites, so the exons around it are joined;
- `clipped_at_edge`: the transcript reaches past the reference path, and only its part on the path is projected;
- `boundary_moved`: an exon boundary lies inside an indel and was moved;
- `splice_site_changed`: a splice site differs from the reference's;
- `ok`: every boundary carried exactly.

`hst.summary.tsv` has one row per gene, in reference order: `transcripts` projected, of which `clipped` reach past the reference path, and `outside` (transcripts with no exon on the path, not projected), `haplotypes_by_copies` (`0:2,1:120,2:9` is 2 haplotypes with no copy, 120 with one, 9 with two), `hsts` written, `same_as_ref`, `annotated` (HSTs from the haplotypes' annotations), the rows of each status (`discordant` and `paralog_lift` only with `--hap-annotations`), and `copies_differ`: for a gene some haplotype's annotation has (found as above), each haplotype with an annotation whose own copies of the gene differ from the run's, as `sample#hap:run/annotated`.

## Limitations

- Whole transcripts of a gene without a copy are looked for on its best chain only, so on one stretch of the walk, and counted as one copy. A transcript with an exon on nodes a paralog shares, or walked below the share, is not written there.
- A copy is found through the nodes it shares with the reference copy, or by sequence where a chain of shared nodes points to it. A copy no chain reaches is not found, nor is a second copy on a walk where the first was found by nodes.
- Placement by sequence takes the window of the gene's best chain only. Where that chain lies on a paralog, the gene's own copy elsewhere on the walk is not looked at: the placement is refused by its identity, or is the paralog's where the reference annotates it.
- A stretch reaching `--min-identity` over the marker's exons, where a chain points, is placed as a copy. No identity separates every allele from every paralog: genes whose alleles differ more need a lower value, at which a paralog the reference does not annotate can pass. `identity` in `detail` shows how close each placement was.
- An exon is found on a merged copy only between two nodes of its own reference copy that the walk passes in order; past the last of them it is missing, as is every exon of that copy on a walk carrying the other copy's alleles alone.
- A short exon between diverged introns can fall in a stretch where the alignment does not hold; the transcripts using it are then not written for that copy (`missing_exon`).
- A haplotype's tandem unit cut short (a partial copy) still takes the layout of the reference unit it resembles most; exons that fall outside what it holds are left out and counted (`tandem_unplaced`).
- Tandem units are told apart by their sequence: a haplotype unit takes the layout of the reference unit it resembles most, and among identical units the choice follows their order. `tandem:14>46` gives the exon count against the reference's.
- A copy mixing the differences of two near-identical paralogs goes to the one it shares more coding bases with; `rival` in `detail` gives the margin.
- Only the reference annotation is projected, so isoforms and genes the reference lacks are found only through a haplotype's own annotation. A changed splice site is reported, not judged: its transcript is still written.
- An annotated transcript joins a projected one only by the same transcript id (for CAT, its source transcript) on the copy it overlaps; an annotation from another source (RefSeq ids in Liftoff) stays beside the projected transcripts, unjudged.
- A filling transcript is judged by its exon count and where it lies, not by its boundaries.
- A gene the GTF lacks is never a `paralog_lift`, so one an annotation lifted onto a paralog's copy is kept as a copy of its own and indexed. Using a GTF of the whole genome, rather than of the loci, lets lifts of genes from elsewhere be judged.
- With `--prefer-annotation`, a truncated annotated model replaces the projected transcript it differs from; one shorter than recruit's k-mers leaves that transcript out of the index altogether.

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
