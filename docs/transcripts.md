# Module `transcripts`

CLI: `panvar transcripts`

## What it does

Projects a reference annotation onto every haplotype walk of a locus graph. Each transcript becomes one haplotype-specific transcript (HST) per haplotype and copy of its gene, spliced from that haplotype's own sequence.

- A haplotype without a copy of a gene gets no transcript of it; one with two copies gets two, numbered along the walk.
- A transcript whose exon the haplotype lacks is reported and not written, so other transcripts of the gene still are.
- An exon boundary that falls inside a variant is moved to the nearest base the haplotype shares with the reference, and flagged.
- Each HST carries its source, `projected`, so transcripts from other annotations can join them later.

## How it works

1. The GTF's exons are placed on the reference path through the `contig:start-end` in its name. Each exon boundary becomes a reference node and an offset inside it, since boundaries rarely fall between nodes. A transcript reaching past the path is not projected.

2. Each gene's copies are found with one transcript, its marker: the one tagged Ensembl_canonical, else the one with the most exonic bases. Its coding exons are used where it has them, because paralogs are best told apart by what they encode.

3. On each haplotype walk, read in both orientations, the nodes that occur once in the gene's stretch of the reference walk are chained in order. A chain pays for every base of indel between its links, so it stays within one copy instead of jumping to the next. The gaps between chained nodes are then aligned node by node. Each chain holding a tenth of the best chain's bases is one pass of the walk through the gene; the rest are stray matches of single nodes. A pass that shares at least `--min-exon-share` of the marker's exonic bases with the reference is a candidate copy, and its stretch of the walk is the part matching the marker's span.

4. Candidates are accepted best first, by that share. A candidate loses to an accepted pass of its own gene that it overlaps. It also loses to an accepted pass of a gene apart on the reference (their markers do not overlap) that takes the same stretch of the walk, more than half of the longer of the two. Two genes folded onto the same nodes, as the copies of a duplicated gene are, thus compete for each stretch, and the one the stretch resembles more takes it.

5. A short candidate lying inside a longer accepted pass of a gene apart on the reference is dropped when that pass matches the stretch at least as well as the candidate does, and the rest compete again. Such a stretch is part of the other gene's copy that resembles this gene: a pseudogene fragment, or a gene nested in the other copy of a duplicated gene. A candidate that matches better is kept, as a copy left by a deletion between two paralogs lies inside whatever gene surrounds them.

6. Within each copy, an exon boundary is carried to the walk through the node holding it, at the same offset. When the walk does not pass that node (the boundary lies inside a variant), the boundary moves to the nearest position the walk shares with the reference, and a tie keeps the variant inside the exon. An exon sharing no base with the walk is missing.

7. The exons are cut from the walk with the haplotype's alleles and joined in transcript order; a minus-strand transcript is reverse-complemented. The two splice-site bases at each end of every intron are compared with the reference's.

## Required inputs

- `-i, --gfa <graph.gfa>`: a locus GFA (plain or gzipped) with one P or W line per haplotype
- `-r, --reference-path <name>`: the reference path; a unique part of the name is enough. Its name must carry `contig:start-end`, as for `bubble --gtf`
- `--gtf <path>`: the annotation in reference coordinates (Ensembl or GENCODE, plain or gzipped); its `exon` rows give the exons, its `CDS` rows the coding parts
- `-o, --out-prefix <prefix>`: output prefix

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `--min-exon-share <F>` | the share of the marker's exonic bases a pass must walk to count as a copy; a lower value also counts partial copies | `0.8` |
| `-t, --threads <N>` | worker threads, one haplotype each (`0` = all cores) | `0` |
| `-q, --quiet` | no progress logs | off |

## Outputs

| file | contents |
|------|----------|
| `<prefix>.hst.fa` | one sequence per HST, 5' to 3', named `transcript|gene|haplotype|copy`, with `source=` and `status=` after the name |
| `<prefix>.hst.tsv` | one row per transcript, haplotype and copy, a row with copy 0 where the haplotype has no copy |
| `<prefix>.hst.summary.tsv` | one row per gene: copies per haplotype and rows per status |

`hst.tsv` (`.` where a value does not apply):

| column | meaning |
|--------|---------|
| `transcript_id`, `transcript_name`, `gene_id`, `gene_name`, `biotype`, `strand` | the transcript, from the GTF |
| `haplotype` | the path |
| `copy`, `copies` | this copy, numbered along the walk as written, and the gene's copies on the walk |
| `source` | `projected` |
| `status` | below |
| `orientation` | `+` where the walk passes the gene in the reference's orientation, `-` where it passes it reversed |
| `path_start`, `path_end` | the HST's first and last exonic base on the walk as written, 0-based, end excluded |
| `length`, `exons` | the HST's length and exons |
| `ref_length`, `ref_exons` | the same on the reference |
| `same_as_ref` | whether the HST is the reference transcript's sequence |
| `moved_boundaries` | exon boundaries moved off a variant |
| `changed_splice_sites` | intron ends whose two splice-site bases differ from the reference's |
| `exon_starts`, `exon_lengths` | each exon's first base on the walk (as `path_start`) and length, in transcript order |
| `detail` | notes, `;`-joined: `e3.start:-3` (exon 3's 5' boundary moved, the exon 3 reference bases shorter), `e2:missing`, `e1:walk_end` (the walk ends before it), `i1.donor:GT>AT` (a changed splice site, transcript strand), `rival:<gene>:<share>><share>` (the paralog this copy was taken from, and both shares), and for no copy `best_share:<F>` and `taken_by:<gene>` or `inside:<gene>` |

Statuses, most severe first:

- `no_copy`: the walk has no copy of the gene;
- `missing_exon`: the copy lacks an exon of this transcript; no sequence is written;
- `boundary_moved`: an exon boundary lies inside a variant and was moved;
- `splice_site_changed`: a splice site differs from the reference's;
- `ok`: every boundary carried exactly.

`hst.summary.tsv` has one row per gene, in reference order: `transcripts` projected and `outside` (reaching past the reference path), `haplotypes_by_copies` (`0:2,1:120,2:9` is 2 haplotypes with no copy, 120 with one, 9 with two), `hsts` written, `same_as_ref`, and the rows of each status.

## Limitations

- A copy is found through the nodes it shares with the reference copy. A copy the graph keeps on nodes of its own, neither folded onto the reference nor aligned to it, is not found.
- Copies of a repeat inside a gene (tandem exons) are not expanded: the HST keeps the reference's exons.
- A copy mixing the differences of two near-identical paralogs goes to the one it shares more coding bases with; `rival` in `detail` gives the margin.
- Only the reference annotation is projected, so isoforms and genes the reference lacks are not found. A changed splice site is reported, not judged: its transcript is still written.

## Example

```bash
# every transcript on every haplotype
panvar transcripts -i locus.gfa.gz -r GRCh38 --gtf genes.gtf.gz -o locus
# count partial copies too
panvar transcripts -i locus.gfa.gz -r GRCh38 --gtf genes.gtf.gz -o locus --min-exon-share 0.5 -t 2
```

`scripts/transcripts_examples.sh` projects the loci GTF onto the paralog graphs in `tests/real_data` and prints each gene's copies per haplotype, in `tests/examples/transcripts`.
