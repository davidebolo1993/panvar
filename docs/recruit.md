# Module `recruit`

CLI: `panvar recruit`

## What it does

Sorts a sample's reads to locus graphs without aligning them. Each locus graph is indexed once against a [background](background.md): every k-mer of its haplotypes is labelled by how often the rest of the genome holds it. Each read, or pair, that touches a locus is then called for that locus:

- `certain`: it carries enough sequence found nowhere else in the background (private k-mers) to belong to the locus, and nothing that points to a homologous copy elsewhere;
- `ambiguous`: it has locus evidence but not enough, or conflicting, with the reason;
- `unrecruited`: it touched the locus's screen but reads as a homolog, or carries no locus evidence.

The output is a table with one row per fragment and locus, the recruited reads as FASTQ if asked, and a summary. Short reads (pairs or single), HiFi and ONT reads are handled alike, from FASTQ, SAM, BAM or CRAM in any sort order.

## How it works

Indexing and scoring are panvar's own; SAM, BAM and CRAM are read with [htslib](https://github.com/samtools/htslib), and the background is queried through `background`'s store or its [SSHash](https://github.com/jermp/sshash) dictionary.

The index, built once per locus graph and kept beside it (`<locus>.plx`, with a report `<locus>.index.tsv`):

1. Every canonical 31-mer of the graph's haplotypes is collected, with the nodes it starts on and its largest copy number within one haplotype. A k-mer once in every haplotype is a backbone k-mer: its depth in a sample is the depth of one copy.
2. The locus is located in every background genome: the reference path's interval (`-r`; its name carries `contig:start-end`) is found through the background's anchors, and the k-mers of each copy found are counted and subtracted from the background's counts. A genome whose copy is broken (a copy reaching a contig end, several copies, or none) is also scanned for stretches made almost entirely of the locus's k-mers (95% of every 5,000 bases), and so is the contig of a copy whose end had to be estimated from the reference (no flank anchor where expected, as past an array longer than the reference's): pieces of the locus the anchors cannot place, such as a tandem array cut by a contig break or running past the estimated end, are subtracted too. A stretch less than 5,000 bases from a contig end reaches that end when at least half of the k-mers in between are the locus's: no window fits there, and assemblies hold more errors near their contig ends. When 0.5% of a locus's k-mers are tandem (3 or more copies within a haplotype) yet left with a few copies outside the located ones, the array has pieces next to copies that look whole, and every genome is scanned for that locus (in parallel, one genome per thread). A paralog 1% diverged keeps about 70% of the locus's k-mers, so it is not taken for a piece.
3. Each k-mer is labelled by what is left (its count outside the locus): private (none), shared (up to `--informative-max` per background haplotype, on average) or repeat (more). Private k-mers of low complexity, or that become shared when the located copies are trimmed at their uncertain ends (unless those ends lie in such a stretch), are labelled shared.
4. Decoys: in the decoy genomes (`--decoy-genomes`), clusters of shared k-mer hits outside the locus's own copy (20 hits, each within 1,500 bp of the next) are homologous windows; their k-mers that are not in the locus and are shared in the background become decoy k-mers, the alleles of the homologs.
5. The screen is the closed syncmers ([Edgar 2021](https://doi.org/10.7717/peerj.10805); 11-mer minimum at either end) among the private and shared k-mers, about one in ten.
6. Calibration: reads are simulated for three classes (2×150 pairs at 0.2% error, HiFi-like reads at 0.1%, ONT-like reads at 3%) from the haplotypes of up to eight samples held out one at a time (the k-mers only that sample carries read as new while its reads are scored), from the decoy windows, and from 5 Mb of the coordinate reference on each side of the locus. They are screened and scored as below, and the index keeps, per class, how many are certain at each cut from one to eight private sites.

Reads, per sample:

7. Every read's closed syncmers are looked up in the screens of all loci; a fragment (a pair, or a single read) is scored for each locus whose screen it hits. In a coordinate-sorted file, a read waits for a mate within 100 kb; a read whose mate lies farther is kept only if it hits a screen, and its mate is fetched at the end, through the index or by a second pass. Mates are fetched in their order on the genome, nearby ones by one query: a CRAM decodes against the reference, and fetching in any other order reloads a reference contig at almost every mate.
8. Each k-mer of a read is labelled private, shared, repeat, decoy, or neither (a sequencing error or a new allele).
9. The read's locus segment is its highest-scoring stretch: locus k-mers count +1, decoy k-mers −2, and a run of k-mers in neither table is free for its first 31 k-mers (one differing base makes that many) and costs after them. When the read leaves the locus at its ends (beyond the segment, few locus k-mers), the segment alone is judged; otherwise the whole read is.
10. Evidence is counted in independent sites: a private (or decoy) k-mer counts only when it starts 31 or more bases after the last one counted, so all k-mers one differing base touches count once. A pair's mates are pooled without counting a site twice.
11. A fragment is certain with at least the cut's private sites (the calibrated cut, below), decoy sites at most 5% of them, at least 0.05 sites per 31 callable bases, and no more runs of unknown k-mers than sequencing errors and about one new allele per kilobase explain (not `divergent`). A pair is judged on its pooled evidence when a mate lies wholly in the locus; failing that, one mate certain on its own makes the pair certain (`rescued`), unless the other mate leans to a homolog away from a locus edge. A fragment that runs out of the locus is `partial`.
12. Otherwise a fragment with locus evidence is ambiguous, with the first reason that holds: `conflict` (private and decoy sites), `uncallable` (mostly unreadable k-mers), `divergent`, `sparse` (enough sites, too thinly spread), `single_site`, `shared_only`; one with only decoy sites is unrecruited as `homolog`, one with nothing as `no_evidence`.
13. The cut: from the first reads, the run's read class (short below 1 kb mean length; long reads HiFi at an error rate up to 1%, else ONT) picks a calibration, and each locus takes the lowest cut of at least two private sites whose simulated off-target reads were certain at most at `--false-certain-rate`. `--min-sites` sets one cut for every locus instead.

## Required inputs

- `--locus <path>`: a locus graph (GFA; PanSN path names) or its index (`.plx`); repeatable
- to index a graph: `--background <store>` and `-r, --reference <name>` (the graph path whose name carries the coordinates of the background's coordinate reference)
- reads: `-1`/`-2` FASTQ pairs, `-i` any one file (FASTQ, SAM, BAM, CRAM), or `--samples <sheet.tsv>` for several samples
- `-o, --out <path>`: the table (or a folder, with `--samples`); not needed with `--index-only`

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `--index-only` | build or check the indexes, read nothing | off |
| `--save-index <dir>` | where indexes are kept and reused | beside each graph |
| `--coordinate-reference <name>` | the background's coordinate reference the graph's coordinates refer to | its only one |
| `--informative-max <N>` | copies per background haplotype outside the locus up to which a k-mer is shared, not repeat | `8` |
| `--decoy-genomes <A,B>` | genomes whose homologous windows give decoys | the coordinate reference |
| `--no-calibration` | build indexes without calibration (the cut is then `--min-sites`, 2) | off |
| `--interleaved` | `-i` FASTQ files hold mates as consecutive records | off |
| `--cram-reference <fa>` | the FASTA (with `.fai`) a CRAM was compressed against | the one its header names |
| `--dry-run` | count each input's records by kind and the time to read them | off |
| `--false-certain-rate <r>` | the off-target certain rate each locus's cut is chosen for | `1e-6` |
| `--min-sites <N>` | one private-site cut for every locus | calibrated |
| `--error-rate <e>` | the reads' error rate per base | from base qualities |
| `--preset <p>` | the error rate of `illumina`, `hifi`, `ont` or `clr` reads | |
| `--reads-out <dir>` | the recruited reads as gzip FASTQ, per locus and fate | off |
| `--certain-only` | with `--reads-out`, certain reads only | off |
| `--summary <path>` | the run's summary | beside the table |
| `--dump-index <plx>=<tsv>` | write every k-mer of an index as text (columns below), and stop | |
| `-t, --threads <N>` | scoring threads; a reader and a quarter as many decompression threads run beside them | Slurm's allocation, else all cores |
| `-q, --quiet` | no progress logs | off |

A sample sheet holds tab-separated lines of sample, file and optionally a mate file (paths relative to the sheet); a sample may take several lines, read in turn. With a sheet, `-o` is a folder holding `<sample>.tsv` per sample.

## Outputs

| file | contents |
|------|----------|
| `<out>` | one row per screened fragment and locus (columns below) |
| `<out stem>.summary.tsv` | the run: error rate and its source, fragments, screened, mean read length; each input's records by kind; per locus, fragments by fate and reason, the backbone depth, the fragments certain for another locus too, and the private-site cut with the calibration behind it |
| `<locus>.plx` | the index |
| `<locus>.index.tsv` | the index's report: the locus (k-mers by label, backbone, pieces, decoys, screen), each copy found per genome with the k-mers subtracted, the pieces found off them, the decoy windows, ceilings (the share of read-length windows of each haplotype with 2 private sites) and the calibration (per read class and cut: held-out fragments certain, off-target certain of those screened, and the 95% bound) |
| `<reads-out>/<locus>.<fate>_1.fq.gz`, `_2.fq.gz` | with `--reads-out`: recruited pairs, as sequenced, with qualities (`<locus>.<fate>.fq.gz` for single reads) |

An index dumped with `--dump-index` has one row per k-mer of the locus: `kmer` (its canonical sequence), `class` (`P`, `S` or `R`), `outside` (its background count outside the locus), `multiplicity` (its largest copy number within one haplotype), `backbone` (1 if once in every haplotype), `screen` (1 if a screen syncmer) and `nodes` (the graph nodes it starts on). Decoy k-mers follow as class `D`, with `.` in the fields only locus k-mers have.

Table columns:

| column | meaning |
|--------|---------|
| `fragment` | the read or pair name |
| `locus` | the locus |
| `fate` | `certain`, `ambiguous` or `unrecruited` |
| `reason` | `.`, or why: `partial`, `rescued` for certain; the ambiguous and unrecruited reasons above |
| `mates` | each mate's place: `in` (wholly in the locus), `part` (a segment covering part of it), `out` |
| `seg1`, `seg2` | each mate's locus segment in read bases (`start-end`), `.` when it has none |
| `sites_p`, `sites_d` | independent private and decoy sites in the segments |
| `decoy_out` | decoy sites outside the segments |
| `runs_n` | runs of k-mers in neither table, in the segments |
| `seg_kmers` | callable k-mers in the segments |
| `n_p`, `n_s`, `n_r`, `n_d`, `n_n` | k-mers per label over whole reads |
| `kmers`, `windows`, `bases` | k-mers read, k-mer windows (callable or not), bases |
| `multi_p` | private k-mers multi-copy within a haplotype |
| `nodes` | the graph nodes of the private k-mers in the segments, `.` if none |
| `other_loci` | the other loci the fragment was scored for |
| `certain_loci` | the loci it is certain for |

## Limitations

- k is 31, as in the background.
- The calibration simulates about 10,000 off-target fragments per locus and class: no false certain fragment among them bounds the rate at a few per 10,000, not one per million, and overlapping fragments share sites. The summary says whether the rate was shown or only none beyond it seen. Rates per million come from whole-genome negatives (see the dev scripts).
- A homologous copy carrying the locus's alleles at two distinguishing bases within one read (gene conversion, a paralogous variant at both) reads as certain; a larger background, where more of those alleles are shared, reduces it.
- Short pairs across a locus edge are less often certain than inside it: the located copies' ends are uncertain, so the k-mers near them are labelled shared.
- Long spliced transcripts cross many junctions and read as divergent.
- A locus missing from a background genome cannot be subtracted there; its k-mers keep that genome's counts, unless the scan for pieces finds it.
- Genomes given to the background as graphs are not scanned for pieces.
- A homologous copy elsewhere close enough to carry three or more of the reference's unique anchors in order is located as a copy of the locus and subtracted, even one the coordinate reference lacks. Past that, the scan for pieces subtracts only stretches tied to the locus (on a located copy's contig, or reaching a contig end): density alone cannot tell a close paralog from a piece, so one with other sequence on both sides, on a contig without a copy, stays outside and its k-mers shared.

## Example

```bash
# the indexes, once per locus graph (kept beside each graph)
panvar recruit --locus locus_a.gfa --locus locus_b.gfa --background bg.bgk -r grch38 --index-only -t 16
# a sample's reads, against both loci
panvar recruit --locus locus_a.plx --locus locus_b.plx -i sample.cram --cram-reference GRCh38.fa \
  -o sample.tsv --reads-out sample_reads -t 16
# several samples from a sheet
panvar recruit --locus locus_a.plx --locus locus_b.plx --samples samples.tsv -o out -t 16
```

`scripts/recruit_examples.sh` indexes the synthetic locus of `scripts/background_examples.sh` against its store and recruits pairs drawn from the locus, from its diverged paralog elsewhere in the reference, and from unrelated sequence, in `tests/examples/recruit`: the locus pairs are certain, the paralog's are not.

A 30× short-read genome takes a few minutes on 8 threads against several loci; nearly every read is rejected by the screen. The simulation benchmark, the real-read checks and the comparisons are in `dev/` ([dev/README.md](../dev/README.md)).
