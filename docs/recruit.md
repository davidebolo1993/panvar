# Module `recruit`

CLI: `panvar recruit`

## What it does

Sorts a sample's reads to locus graphs without aligning them. Each locus graph is indexed once against a [background](background.md): every k-mer of its haplotypes is labelled by how often the rest of the genome holds it. Each read, or pair, that touches a locus gets one score for that locus, from the sequence it carries that is found nowhere else in the background (private k-mers) against the sequence that points to a homologous copy elsewhere (decoy k-mers), and is called:

- `certain`: its score reaches the locus's certain threshold, the lowest at which calibrated off-target reads reach it at most at `--certain-rate`, with at least two private sites;
- `likely`: its score reaches the likely threshold, set the same way for `--likely-rate`;
- `ambiguous`: it has locus evidence but too little, or conflicting, with the reason;
- `unrecruited`: it touched the locus's screen but reads as a homolog, or carries no locus evidence.

The output is a table with one row per fragment and locus, the recruited reads as FASTQ if asked, and a summary. Short reads (pairs or single), HiFi and ONT reads are handled alike, from FASTQ, SAM, BAM or CRAM in any sort order.

## How it works

Indexing and scoring are panvar's own; SAM, BAM and CRAM are read with [htslib](https://github.com/samtools/htslib), and the background is queried through `background`'s store or its [SSHash](https://github.com/jermp/sshash) dictionary.

The index, built once per locus graph and kept beside it (`<locus>.plx`, with a report `<locus>.index.tsv`):

1. Every canonical 31-mer of the graph's haplotypes is collected, with the nodes it starts on and its largest copy number within one haplotype. A k-mer once in every haplotype is a backbone k-mer: its depth in a sample is the depth of one copy.
2. The locus is located in every background genome: the reference path's interval (`-r`; its name carries `contig:start-end`) is found through the background's anchors, and the k-mers of each copy found are counted and subtracted from the background's counts. A copy is certain between its outermost inside anchors; past them its ends are estimated and may run into flank. A locus with fewer than three anchors inside (a short one) takes a copy only between flank anchors on both sides, and with none inside, only the interval projected from those flank anchors is certain. Where a locus's haplotypes differ from the reference, its anchors hold only in part of each copy, so the copy is walked out from them over every 1,000 bases in which at least half of the k-mer starts are the locus's (repeats left out): the copy's own sequence, new alleles included, carries the walk, and a flank stops it. The report gives each copy's anchored stretch and its walked-out ends. A genome whose copy is broken (a copy reaching a contig end, several copies, or none) is also scanned for stretches made almost entirely of the locus's k-mers (95% of every 5,000 bases), and so is the contig of a copy whose end had to be estimated from the reference (no flank anchor where expected, as past an array longer than the reference's): pieces of the locus the anchors cannot place, such as a tandem array cut by a contig break or running past the estimated end, are subtracted too. A stretch less than 5,000 bases from a contig end reaches that end when at least half of the k-mers in between are the locus's: no window fits there, and assemblies hold more errors near their contig ends. When 0.5% of a locus's k-mers are tandem (3 or more copies within a haplotype) yet left with a few copies outside the located ones, the array has pieces next to copies that look whole, and every genome is scanned for that locus (in parallel, one genome per thread). A dense stretch is a piece only when it is tied to the locus: next to a located copy (less than 5,000 bases between them, or as little to another such piece), reaching a contig end, or covering a whole contig shorter than 5,000 bases (with at least half of its k-mers the locus's). One with other sequence on both sides, apart from every copy, may be a close paralog the coordinate reference lacks, and stays outside, also on a copy's contig: an assembly contig can hold a whole chromosome arm, close paralogs and all. A genome indexed for reading by interval (`.fai`, and `.gzi` for bgzip) is read only where a piece can lie, around each located copy and at each contig's ends, widened while a piece reaches the edge of what was read; any other genome is read in full.
3. Each k-mer is labelled by what is left (its count outside the locus): private (none), shared (up to `--informative-max` per background haplotype, on average) or repeat (more). Private k-mers of low complexity, or that become shared when the located copies are cut at their walked-out ends (unless the rest lies in such a stretch), are labelled shared.
4. Decoys: in the decoy genomes (`--decoy-genomes`), clusters of shared k-mer hits outside the locus's own copy (20 hits, each within 1,500 bp of the next) are homologous windows; their k-mers that are not in the locus and are shared in the background become decoy k-mers, the alleles of the homologs.
5. The screen is the closed syncmers ([Edgar 2021](https://doi.org/10.7717/peerj.10805); 11-mer minimum at either end) among the private and shared k-mers, about one in ten.
6. Calibration, for three read classes (2×150 pairs at 0.2% error, HiFi-like reads at 0.1%, ONT-like reads at 3%), screened and scored as below:
   - locus reads, simulated from the haplotypes of up to eight samples held out one at a time (the k-mers only that sample carries read as new while its reads are scored);
   - off-target reads, simulated from the decoy windows and from 5 Mb of the coordinate reference on each side of the locus, read as the index labels them: they gain private sites only by sequencing errors. Long reads are laid over a decoy window at every overlap with it, from one k-mer at an end of the read to the whole read, one start drawn in each stride so that each overlap comes up as often as sequencing gives it: a read that touches a homolog only at its end is judged on a short stretch, where two sites score high, and a threshold must have seen such reads. (The first rules' counts keep reads that end at least 1 kb into a window, simulated apart.) Pairs overlap a window by a read length or more, at a stride of tens of bases, which leaves no overlap out: a window reaches 1,500 bases past the homolog on each side;
   - off-target fragments of the calibration genomes, each held out of the background in turn: its homologous windows with 15 kb of context and 1 Mb on each side of its copies of the locus, read with every locus k-mer that only this genome holds outside the locus labelled private, and every decoy k-mer that only it holds labelled unknown. That is what a new sample's paralog looks like: alleles the locus shares with no genome of the background, and alleles of its own no decoy records. These are not simulated read by read: the genome's sequence is labelled once, fragments of each class are laid over it as the simulated ones are (about 100,000 pairs, or up to 4,000 long reads, over the windows and as many over the flanks; long reads at every overlap with a window), and a fragment holding no private k-mer, which cannot be called, is only counted (as screened with the chance that one of its screen k-mers is read without an error). The others are scored from the labels, with the class's errors drawn along them and every k-mer over an error read as unknown. An error that makes a private k-mer by chance is not looked for here: the simulated off-target reads measure those, and the two add up (step 14).

   The calibration genomes are taken among candidates: the decoy genomes other than the coordinate reference, or those `--calibration-genomes` lists. Every candidate is read once, all loci at once, and each locus takes up to `--calibration-max` (8) of them: first the one whose windows and flanks hold the most distinct homolog k-mers (locus k-mers that are not private, and decoy k-mers), then each time the one adding the most k-mers the taken ones lack, while one adds any (ties: the earlier listed). Many candidates cost little more than few: only their windows and flanks are read. A candidate that is a decoy genome has its own homologous windows; one that is not is not scanned for them: its windows are where the anchors place the coordinate reference's decoy windows, so a homolog the coordinate reference lacks is not seen in it. With no candidate there is no calibration genome, and the index report says so.

   A stretch of a calibration genome that reads as a copy of the locus is not off-target: where 95% of 5,000 k-mer starts are the locus's own k-mers (private once the genome is held out, or shared; repeats do not count), the density the pieces of step 2 take, the fragments over it are left out and the report counts them. The index cannot place such a stretch (no anchors in order, not next to a copy), so its k-mers stay shared; but a sample's reads from it are reads of the locus's sequence, and counting them as false would raise the thresholds for the whole locus.

   From these the index keeps, per class, the site rates of the score (step 11): private and decoy sites per 31 k-mers of judged sequence in the locus reads (`r_in`, `q_in`) and in off-target ones (`r_out`, `q_out`), each with half a site and one 31-k-mer stretch added so that no rate is 0. The off-target rates are each calibration genome's own, to which the private sites sequencing errors make (the simulated off-target reads') are added, and of those the median over the genomes: a genome that holds much of the locus elsewhere does not set them. Without calibration genomes they are the simulated off-target reads'. It keeps histograms of the deciding scores, in bins of 1/16, of the locus reads, the simulated off-target reads and the calibration genomes' fragments; the thresholds are chosen from them when reads are scored, so the rates are a run's choice. It keeps too, for `--rules v1`, how many reads are certain at each private-site cut from one to eight.

Reads, per sample:

7. Every read's closed syncmers are looked up in the screens of all loci; a fragment (a pair, or a single read) is scored for each locus whose screen it hits. In a coordinate-sorted file, a read waits for a mate within 100 kb; a read whose mate lies farther is kept only if it hits a screen, and its mate is fetched at the end, through the index or by a second pass. Mates are fetched in their order on the genome. A CRAM fetches mates within 200 kb of each other by one query, as it decodes against the reference and any other order reloads a reference contig at almost every mate; a BAM queries each mate's position alone, as a query decodes every record between its ends.
8. Each k-mer of a read is labelled private, shared, repeat, decoy, or neither (a sequencing error or a new allele).
9. The read's locus segment is its highest-scoring stretch: locus k-mers count +1, decoy k-mers −2, and a k-mer in neither table is free within 31 k-mers of a locus or decoy k-mer its run touches, on either side (one differing base makes that many), and costs beyond that; free k-mers at either end of the stretch are kept, so a read and its reverse complement get the same segment. When the read leaves the locus at its ends (beyond the segment, few locus k-mers, and a homologous flank on one side at most), the segment alone is judged; otherwise the whole read is. A side is a homologous flank when its decoy sites would count against the segment, two or more and at least the decoy share of the segment's private sites, so a stray decoy k-mer far along a long read does not make the whole read judged.
10. Evidence is counted in independent sites: a private (or decoy) k-mer counts only when it starts 31 or more bases after the last one counted, so all k-mers one differing base touches count once. A pair's mates are pooled without counting a site twice.
11. The score of a judged stretch with n callable k-mers, p private sites and d decoy sites is

    S = p · ln(r_in / r_out) − (n / 31) · (r_in − r_out) + d · ln(q_in / q_out) − (n / 31) · (q_in − q_out)

    the log-likelihood ratio of the counts seen, for a locus read against an off-target one, each a Poisson process of sites along the stretch. It holds length, density and decoys in one number: a long stretch from the locus is expected to carry many private sites, so a few on a long stretch score low; decoy sites weigh against private ones. Sites cluster on paralog variants, so S ranks fragments rather than giving a probability: the thresholds come from calibration, which absorbs that. A read without a segment is judged whole, its decoy sites counted, so its k-mers count in n.
12. A stretch may call a fragment when it has a private site, is callable (not mostly unreadable k-mers) and not `divergent`: no more runs of unknown k-mers than sequencing errors and about one new allele per kilobase explain. The errors expected are the read's own, from its base qualities, where they exceed the run's rate: HiFi reads range from 99% to 99.99% accurate, and the least accurate are not judged by the average. A pair's stretches are the pooled one (both mates' segments, a site the mates share counted once, n the sum of both) and each mate alone. The pooled stretch may call the pair when a mate lies wholly in the locus: two mates that both leave it would pool islands. A mate alone may call it unless the other mate leans to a homolog (more decoy sites than private ones) away from a locus edge; the mates lie an insert apart, so the other one is in the locus or its flank. The fragment's deciding score is the best of the stretches that may call it.
13. A fragment is certain when the pooled stretch, or one mate alone (`rescued`), reaches the certain threshold with at least two private sites: one differing base cannot be told from one sequencing error or one new allele, whatever it scores. It is likely when a stretch reaches the likely threshold. A certain fragment that runs out of the locus is `partial`. A fragment with no private site and two or more decoy sites is unrecruited as `homolog`. Any other fragment that is not certain keeps the first reason that holds: `conflict` (private and decoy sites), `uncallable`, `divergent`, `mates_out` (the pooled stretch would be certain, but both mates leave the locus), `single_site`, `low_score` (two or more private sites, below the threshold), `shared_only`; a likely fragment has one, an ambiguous one too, and one with nothing is unrecruited as `no_evidence`.
14. The thresholds: the run's read class (`--preset`, or from the first reads: short below 1 kb mean length; long reads HiFi at an error rate up to 1%, else ONT) picks a calibration, and each locus takes, per tier, the lowest score bin whose off-target fragments reach it at most at the tier's rate. The off-target share at a bin is the sum of two: the simulated off-target reads reaching it, of those screened (private sites made by sequencing errors), and the calibration genomes' fragments reaching it, of those screened (private sites a new sample holds). They are different causes of a false read, so they add; pooling the fragments would dilute each by the other's count.

    A rate is shown when the 95% bound of that share is within it. When the fragments are too few for that even with none beyond a threshold (the bound is 3 over the number screened, for each kind), the threshold is extrapolated: the scores are taken to go on beyond those seen as an exponential tail, and the threshold is the lowest at which such a tail lies above every score seen. The tail's scale is fitted from what makes a higher score: one more private site adds ln(r_in / r_out), and is as rare as a private site is among the off-target fragments (the share of them that can be called, at its 95% bound), so the scale is ln(r_in / r_out) / ln(1 / that share). The extrapolated threshold is never below the observed one. The summary gives, per tier, the threshold, the held-out recall, the off-target fragments reaching it and their share, the 95% bound, and whether the rate was shown or the threshold extrapolated (with the scale, and the threshold the fragments seen give).
15. The first rules stay for comparison (`--rules v1`), and decide an index without calibration or a run with `--min-sites`: a fragment is certain with at least the cut's private sites, decoy sites at most 5% of them (25% in reads of 1 kb or more), callable and not divergent, pooled and rescued as above. The cut is calibrated per locus, the lowest of at least two private sites (three for ONT) whose simulated off-target reads were certain at most at `--certain-rate`; without calibration it is 2 (3 for ONT); `--min-sites` sets it for every locus.

## Required inputs

- `--locus <path>`: a locus graph (GFA; PanSN path names), its haplotypes (FASTA, one record per path, named as the paths would be), or its index (`.plx`); repeatable. A FASTA locus gives the index its graph would, without the nodes: they are the only thing the index takes from the graph beyond its paths' sequences, so a FASTA locus's `nodes` are `.`.
- to index a graph: `--background <store>` (a store, or a packaged background: see `background --package`) and `-r, --reference <name>` (the graph path whose name carries the coordinates of the background's coordinate reference)
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
| `--calibration-genomes <A,B>` | the candidates for the genomes calibration holds out one at a time (FASTA; not the coordinate reference) | the decoy genomes other than the coordinate reference |
| `--calibration-max <N>` | candidates each locus holds out, at most, taken by the homolog k-mers each adds | `8` |
| `--no-calibration` | build indexes without calibration (fragments are then decided by the private-site cut, `--min-sites`, 2; 3 for ONT reads) | off |
| `--interleaved` | `-i` FASTQ files hold mates as consecutive records | off |
| `--cram-reference <fa>` | the FASTA (with `.fai`) a CRAM was compressed against | the one its header names |
| `--dry-run` | count each input's records by kind and the time to read them | off |
| `--certain-rate <r>` | the off-target rate each locus's certain threshold is chosen for | `1e-4` |
| `--likely-rate <r>` | the same for likely; at least `--certain-rate` | `1e-2` |
| `--min-sites <N>` | decide by one private-site cut for every locus instead of the score | |
| `--error-rate <e>` | the reads' error rate per base, in [0, 0.5); a read whose qualities say more is held to its own | from base qualities |
| `--decoy-share-long <r>` | decoy sites per private site beyond a segment that make a side a homologous flank in reads of 1 kb or more (0.05 in shorter ones), in [0, 1] | 0.25 |
| `--preset <p>` | the read type, `illumina`, `hifi` or `ont`: its calibration class and error rate (0.002, 0.001, 0.03); `--error-rate` overrides the rate. Corrected long reads (HiFi, error-corrected ONT) are `hifi`, any other long reads `ont` | from the reads |
| `--reads-out <dir>` | the recruited reads as gzip FASTQ, per locus and fate (certain, likely, ambiguous) | off |
| `--certain-only` | with `--reads-out`, certain reads only | off |
| `--summary <path>` | the run's summary | beside the table |
| `--dump-index <plx>=<tsv>` | write every k-mer of an index as text (columns below), and stop | |
| `-t, --threads <N>` | scoring threads; beside them run a reader, and for SAM, BAM or CRAM a quarter as many decompression threads and a quarter as many threads preparing records, so the run can use about 1.5 × N cores while it reads (ask a scheduler for N and expect that) | the CPUs available |
| `-q, --quiet` | no progress logs | off |

A sample sheet holds tab-separated lines of sample, file and optionally a mate file (paths relative to the sheet); a sample may take several lines, read in turn. With a sheet, `-o` is a folder holding `<sample>.tsv` per sample.

## Outputs

| file | contents |
|------|----------|
| `<out>` | one row per screened fragment and locus (columns below) |
| `<out stem>.summary.tsv` | the run: error rate and its source, fragments, screened, mean read length; each input's records by kind; per locus, fragments by fate and reason, the backbone depth, the fragments certain for another locus too, the score's rates, and per tier the threshold, the held-out recall, the off-target fragments reaching it of those screened, their share and its 95% bound, and whether the threshold was extrapolated (or, decided by the first rules, the private-site cut with the calibration behind it) |
| `<locus>.plx` | the index |
| `<locus>.index.tsv` | the index's report: the locus (k-mers by label, backbone, pieces, decoys, screen, calibration genomes), each copy found per genome with the k-mers subtracted, the pieces found off them, the decoy windows, ceilings (the share of read-length windows of each haplotype with 2 private sites, counted as the scorer counts them: k bases apart), the first rules' calibration (per read class and cut: held-out fragments certain, off-target certain of those screened, and the 95% bound), each calibration genome (windows, bases read, homolog k-mers, k-mers relabelled private, decoys unlabelled, bases that read as a copy of the locus) and the windows placed in those that are not decoy genomes, the score's rates per read class, each calibration genome's off-target fragments (screened, with a private site, left out over a copy of the locus) and its own rates, and the thresholds, recall and off-target share at rates from 1e-1 to 1e-5, each observed or extrapolated |
| `<reads-out>/<locus>.<fate>_1.fq.gz`, `_2.fq.gz` | with `--reads-out`: recruited pairs (certain, likely, ambiguous), as sequenced, with qualities (`<locus>.<fate>.fq.gz` for single reads); with `--samples`, in `<reads-out>/<sample>/`. An earlier run's files for the same loci are removed first. Each file is a chain of gzip members, compressed by all threads (any gzip reader reads it as one) |

An index dumped with `--dump-index` has one row per k-mer of the locus: `kmer` (its canonical sequence), `class` (`P`, `S` or `R`), `outside` (its background count outside the locus), `multiplicity` (its largest copy number within one haplotype), `backbone` (1 if once in every haplotype), `screen` (1 if a screen syncmer) and `nodes` (the graph nodes it starts on). Decoy k-mers follow as class `D`, with `.` in the fields only locus k-mers have.

Table columns:

| column | meaning |
|--------|---------|
| `fragment` | the read or pair name |
| `locus` | the locus |
| `fate` | `certain`, `likely`, `ambiguous` or `unrecruited` |
| `reason` | `.`, or why: `partial`, `rescued` for certain; why not certain for the others (above) |
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
| `certain_loci` | how many loci it is certain for |
| `score` | the deciding score (the pooled one when no stretch may call the fragment), `.` when the index has no score calibration for the run's reads |
| `tier` | the tier the deciding score reaches, `certain`, `likely` or `.`; with `--rules v1` or `--min-sites`, where the score would put the fragment |

## Limitations

- k is 31, as in the background.
- The simulated off-target reads are about 100,000 short pairs, and up to 4,000 long reads, from the decoy windows and as many from the flanks, fewer when the source is short (two 5 Mb flanks give about 1,300 HiFi or 2,000 ONT reads); fewer still hit the screen. No off-target fragment above a threshold bounds its share at 3 over the number screened (a 95% bound): showing 1e-4 takes 30,000 screened of each kind. Below that the certain threshold is extrapolated (step 14): a model of how the scores go on, not a measurement, and it costs some recall where the fragments seen would have allowed a lower threshold. Overlapping fragments share sites, and the bounds treat them as independent. Rates per million come from whole-genome negatives (see the dev scripts).
- The thresholds are per locus and read class, from all calibration genomes' fragments together: a paralog that only one calibration genome holds, with stretches identical to the locus but too broken or too short to read as a copy of it (step 6), gives fragments that score as the locus's, and the thresholds rise above them for the whole locus. The median keeps the rates from it, not the thresholds. With few genomes in the background, most paralogs are held by one genome only, and certain recall can fall far.
- The constants the score rests on, kept as they are: a side of a read is a homologous flank with two or more decoy sites and at least 5% of the segment's private sites (25% in reads of 1 kb or more); a fragment with no private site and two or more decoy sites is a homolog; a calibration genome is read 1 Mb on each side of its copies; each site rate has half a site and one 31-k-mer stretch added; scores are binned at 1/16. None is fitted to data.
- A homologous copy carrying the locus's alleles at two distinguishing bases within one read (gene conversion, a paralogous variant at both) reads as certain; a larger background, where more of those alleles are shared, reduces it.
- Short pairs across a locus edge are less often certain than inside it: the located copies' ends are uncertain, so the k-mers near them are labelled shared.
- Long spliced transcripts cross many junctions and read as divergent.
- A locus missing from a background genome cannot be subtracted there; its k-mers keep that genome's counts, unless the scan for pieces finds it.
- Genomes given to the background as graphs are not scanned for pieces.
- A homologous copy elsewhere close enough to carry three or more of the reference's unique anchors in order is located as a copy of the locus and subtracted, even one the coordinate reference lacks. Past that, the scan for pieces subtracts only stretches tied to the locus (next to a located copy, or reaching a contig end): density alone cannot tell a close paralog from a piece, so one with other sequence on both sides, apart from every copy, stays outside and its k-mers shared.

## Example

```bash
# the indexes, once per locus graph, kept beside each graph
panvar recruit --locus locus_a.gfa --locus locus_b.gfa --background bg.bgk -r GRCh38 --index-only
# one sample
panvar recruit --locus locus_a.plx --locus locus_b.plx -1 sample_1.fq.gz -2 sample_2.fq.gz -o sample.tsv
# a sample sheet, one table per sample in out/
panvar recruit --locus locus_a.plx --locus locus_b.plx --samples samples.tsv -o out
```

`scripts/recruit_examples.sh` indexes the synthetic locus of `scripts/background_examples.sh` against its store and recruits pairs drawn from the locus, from its diverged paralog elsewhere in the reference, and from unrelated sequence, in `tests/examples/recruit`: the locus pairs are certain, the paralog's are not.
