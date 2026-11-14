# Module `transcripts`

CLI: `panvar transcripts`

## What it does

Projects a reference annotation onto every haplotype walk of a locus graph. Each transcript becomes one haplotype-specific transcript (HST) per haplotype and copy of its gene, spliced from that haplotype's own sequence.

- A haplotype without a copy of a gene gets no transcript of it; one with two copies gets two, numbered along the walk.
- Exons inside a repeat unit of the gene follow the haplotype's own number of units.
- A transcript whose exon the haplotype lacks is reported and not written, so other transcripts of the gene still are.
- An exon boundary that falls inside a variant is moved to the nearest base the haplotype shares with the reference, and flagged.
- A haplotype's own annotation (`--hap-annotations`, e.g. HPRC's CAT or Liftoff) adds its transcripts, cut straight from the walk: isoforms and copies the reference annotation lacks. Each HST carries its source: `projected`, `annotated`, or both.

## How it works

1. The GTF's exons are placed on the reference path through the `contig:start-end` in its name. Each exon boundary becomes a reference node and an offset inside it, since boundaries rarely fall between nodes. A transcript reaching past either end of the path keeps its exons on the path, cut at the path's end, so a gene at the edge of a locus still gets transcripts; one with no exon on the path is left out.

2. Each gene's copies are found with one transcript, its marker: the one tagged Ensembl_canonical, else the one with the most exonic bases. Its coding exons are used where it has them, because paralogs are best told apart by what they encode.

3. On each haplotype walk, read in both orientations, the nodes that occur once in the gene's stretch of the reference walk are chained in order. A chain pays for every base of indel between its links, so it stays within one copy instead of jumping to the next; across a tandem stretch (step 7) it pays nothing, as haplotypes carry any number of units there. The gaps between chained nodes are then aligned node by node. Each chain holding a tenth of the best chain's bases is one pass of the walk through the gene; the rest are stray matches of single nodes. A pass that shares at least `--min-exon-share` of the marker's exonic bases with the reference is a candidate copy, and its stretch of the walk is the part matching the marker's span.

4. Candidates are accepted best first, by that share. A candidate loses to an accepted pass of its own gene that it overlaps. It also loses to an accepted pass of a gene apart on the reference (their markers do not overlap) that takes the same stretch of the walk, more than half of the longer of the two. Two genes folded onto the same nodes, as the copies of a duplicated gene are, thus compete for each stretch, and the one the stretch resembles more takes it.

5. A short candidate lying inside a longer accepted pass of a gene apart on the reference is dropped when that pass matches the stretch at least as well as the candidate does, and the rest compete again. Such a stretch is part of the other gene's copy that resembles this gene: a pseudogene fragment, or a gene nested in the other copy of a duplicated gene. A candidate that matches better is kept, as a copy left by a deletion between two paralogs lies inside whatever gene surrounds them.

6. Within each copy, an exon boundary is carried to the walk through the node holding it, at the same offset. When the walk does not pass that node (the boundary lies inside a variant), the boundary moves to the nearest position the walk shares with the reference, and a tie keeps the variant inside the exon. An exon sharing no base with the walk is missing.

7. Exons that repeat with a tandem unit inside the gene start at the same node and offset in each unit the reference walks, so a transcript's exons starting where another one starts mark its tandem stretch. The stretch is cut into units where its first exon starts, and each reference unit's layout is the exons it holds, in order. On the haplotype, the same start cuts the walk into units between the exons around the stretch. A haplotype with the reference's number of units takes the reference's layouts in order, so the reference walk gives back its annotation; otherwise each unit takes the layout that places most of its exons exactly, then that of the reference unit sharing most of its nodes, then the nearest in relative position. Each exon of a layout is placed where the unit holds its first and last nodes at about its length apart; a unit missing one of the two (a variant at the exon's first or last base) gets the exon placed from the other at the reference's length, and a unit missing the start node is placed with the layout again where the previous one ends. The marker's tandem stretch does not count towards its share, as haplotypes carry the unit a different number of times.

8. The exons are cut from the walk with the haplotype's alleles and joined in transcript order; a minus-strand transcript is reverse-complemented. The two splice-site bases at each end of every intron are compared with the reference's; an intron the haplotype shortens below those four bases is lost, and its exons are joined.

9. A haplotype's own annotation is in its assembly's coordinates, and the haplotype's path names its contig and interval (`sample#hap#contig:start-end`, the start 0-based as pggb writes it). Its transcripts with an exon on the path are cut from the walk at those coordinates, strand-aware; one reaching past the path is cut at its ends (`clipped_at_edge`). Its copies are the annotation's: each gene record is one copy, numbered along the walk. An annotated transcript with the same id, copy and sequence as a projected one is that record, of both sources (`projected+annotated`); any other is written beside the projected ones, so where the two disagree both stay.

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
| `--contig-alias <tsv>` | contig names of the annotations matched to the paths' contigs: each row lists names of one contig (a chromAlias table works); without it a name matches as it is, or once a `sample#hap#` prefix is removed | — |
| `--min-exon-share <F>` | the share of the marker's exonic bases a pass must walk to count as a copy; a lower value also counts partial copies | `0.8` |
| `-t, --threads <N>` | worker threads, one haplotype each (`0` = all cores) | `0` |
| `-q, --quiet` | no progress logs | off |

## Outputs

| file | contents |
|------|----------|
| `<prefix>.hst.fa` | one sequence per HST, 5' to 3', named `transcript|gene|haplotype|copy` (`|annotated` added to an annotated transcript written beside a projected one), with `source=` and `status=` after the name |
| `<prefix>.hst.tsv` | one row per transcript, haplotype and copy, a row with copy 0 where the haplotype has no copy |
| `<prefix>.hst.summary.tsv` | one row per gene: copies per haplotype and rows per status |

`hst.tsv` (`.` where a value does not apply):

| column | meaning |
|--------|---------|
| `transcript_id`, `transcript_name`, `gene_id`, `gene_name`, `biotype`, `strand` | the transcript, from the GTF |
| `haplotype` | the path |
| `copy`, `copies` | this copy, numbered along the walk as written, and the gene's copies on the walk |
| `source` | `projected`, `annotated` (from the haplotype's annotation; `orientation` is `.`, and `ref_length`, `ref_exons` are the projected transcript's of that id, if any) or `projected+annotated` (both, the same sequence) |
| `status` | below |
| `orientation` | `+` where the walk passes the gene in the reference's orientation, `-` where it passes it reversed |
| `path_start`, `path_end` | the HST's first and last exonic base on the walk as written, 0-based, end excluded |
| `length`, `exons` | the HST's length and exons |
| `ref_length`, `ref_exons` | the same on the reference path (the part on the path, for a transcript reaching past it) |
| `same_as_ref` | whether the HST is the reference transcript's sequence |
| `moved_boundaries` | exon boundaries moved off a variant |
| `changed_splice_sites` | intron ends whose two splice-site bases differ from the reference's |
| `exon_starts`, `exon_lengths` | each exon's first base on the walk (as `path_start`) and length, in transcript order |
| `detail` | notes, `;`-joined: `e3.start:-3` (exon 3's 5' boundary moved, the exon 3 reference bases shorter), `e2:missing`, `e1:walk_end` (the walk ends before it), `i1.donor:GT>AT` (a changed splice site, transcript strand), `i2:intron_lost`, `clipped:3':200` (annotated exonic bases past the path at the transcript's 3' end), `tandem:14>46` (the exons of the tandem stretch on the reference and on the haplotype; the transcript's exons and introns are then numbered over its own), `tandem_inferred:1` (tandem exons placed from one end), `rival:<gene>:<share>><share>` (the paralog this copy was taken from, and both shares), and for no copy `best_share:<F>` and `taken_by:<gene>` or `inside:<gene>` |

Statuses, most severe first:

- `no_copy`: the walk has no copy of the gene;
- `missing_exon`: the copy lacks an exon of this transcript; no sequence is written;
- `intron_lost`: the haplotype shortens an intron below its two splice sites, so the exons around it are joined;
- `clipped_at_edge`: the transcript reaches past the reference path, and only its part on the path is projected;
- `boundary_moved`: an exon boundary lies inside a variant and was moved;
- `splice_site_changed`: a splice site differs from the reference's;
- `ok`: every boundary carried exactly.

`hst.summary.tsv` has one row per gene, in reference order: `transcripts` projected, of which `clipped` reach past the reference path, and `outside` (transcripts with no exon on the path, not projected), `haplotypes_by_copies` (`0:2,1:120,2:9` is 2 haplotypes with no copy, 120 with one, 9 with two), `hsts` written, `same_as_ref`, `annotated` (HSTs from the haplotypes' annotations), the rows of each status, and `copies_differ`: each haplotype with an annotation whose copies of the gene differ from the projection's, as `sample#hap:projected/annotated`.

## Limitations

- A copy is found through the nodes it shares with the reference copy. A copy the graph keeps on nodes of its own, neither folded onto the reference nor aligned to it, is not found.
- Tandem exons are found by their first and last nodes, so a unit whose exon differs from every reference unit at both ends is not found; `tandem:14>46` gives the count against the reference's.
- A copy mixing the differences of two near-identical paralogs goes to the one it shares more coding bases with; `rival` in `detail` gives the margin.
- Only the reference annotation is projected, so isoforms and genes the reference lacks are found only through a haplotype's own annotation. A changed splice site is reported, not judged: its transcript is still written.
- An annotated transcript joins a projected one only by the same transcript id and copy number (for CAT, its source transcript); an annotation from another source (RefSeq ids in Liftoff) stays beside the projected transcripts.

## Example

```bash
# every transcript on every haplotype
panvar transcripts -i locus.gfa.gz -r GRCh38 --gtf genes.gtf.gz -o locus
# count partial copies too
panvar transcripts -i locus.gfa.gz -r GRCh38 --gtf genes.gtf.gz -o locus --min-exon-share 0.5 -t 2
# with the haplotypes' own annotations beside the projected ones
panvar transcripts -i locus.gfa.gz -r GRCh38 --gtf genes.gtf.gz --hap-annotations hap_annotations.tsv \
  --contig-alias chrom_alias.tsv -o locus
# the haplotypes' own annotations alone
panvar transcripts -i locus.gfa.gz --hap-annotations hap_annotations.tsv -o locus
```

`scripts/transcripts_examples.sh` projects the loci GTF onto the paralog graphs in `tests/real_data` and prints each gene's copies per haplotype, in `tests/examples/transcripts`.
