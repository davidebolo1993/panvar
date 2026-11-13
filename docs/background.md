# Module `background`

CLI: `panvar background`

## What it does

Counts every 31-mer of a set of genomes once, so that `recruit` can tell, for each k-mer of a locus, how often the rest of the genome holds it. The set can be one genome (the coordinate reference), many assemblies, or graphs of many haplotypes. The store also holds anchors: k-mers unique in the coordinate reference, placed in every genome, so any interval of the reference can be found in each genome without aligning.

- Counts are exact: every k-mer with how often it occurs, summed over strands.
- A store grows: genomes can be added later as segments, each counted in its own job; a query sums the segments.
- Once final, a store can be compacted into an [SSHash](https://github.com/jermp/sshash) dictionary, smaller and faster to query; `recruit` uses it when it is current.

## How it works

1. A k-mer and its reverse complement count as one (the canonical k-mer). A base other than A, C, G or T ends the k-mers on either side of it.
2. Consecutive k-mers that share a minimizer (the smallest canonical 11-mer by hash) form a super-k-mer, as in KMC 2 ([Deorowicz et al. 2015](https://doi.org/10.1093/bioinformatics/btv022)), and the minimizer picks one of `--partitions` partitions, so both strands of a k-mer reach the same partition. Readers cut the genomes into chunks; workers write each super-k-mer, 2 bits a base, to the spill file of its partition group under `--tmp-dir`. Identical super-k-mers are stored once, with their multiplicity.
3. Each spill file is counted on its own thread within its share of `--memory`; a partition too large for it is counted in sorted runs on disk and merged. A counted partition is written as its sorted k-mers and their counts. The k-mers are Elias-Fano coded ([Vigna 2013](https://doi.org/10.1145/2433396.2433409)): each 62-bit k-mer is split into high bits, stored as gaps in unary in a bit vector, and low bits, stored at a fixed width, which takes about 2 + log2(universe / n) bits per k-mer and still allows a binary search without decoding. The counts are bit-packed: each takes the same number of bits, the width that fits most of them (often 1 or 2), and the few larger ones are listed apart.
4. A graph genome (GFA) is counted from its nodes and the junctions its walks use, each weighted by the number of walks through it, so its haplotypes are counted without spelling them. Each PanSN haplotype (`sample#haplotype`) counts as one haplotype of the background.
5. Anchors: closed syncmers ([Edgar 2021](https://doi.org/10.7717/peerj.10805)) that occur once in the coordinate reference are sampled about every kilobase. While a genome is counted, where it holds each anchor is recorded. A reference interval is found in a genome by chaining the anchors it holds inside and around the interval.
6. A run writes a segment: a complete store of its genomes, in its own folder. `--add` counts more genomes into a new segment of the store; `--segment-only` counts them into a folder in a separate job, and `--attach` commits such folders. A store changes only when its manifest is rewritten, so a reader never sees half a segment, and a lock keeps two jobs from committing at once.
7. `--merge` streams the segments' sorted partitions into one segment, summing the counts. `--range I/N` splits the partitions among N jobs; the last to finish commits. Merged-away segments stay on disk until `--clean`, since running jobs may still read them.
8. `--compact` spells the merged store's k-mers as a set of strings, each k-mer once, in strings where consecutive k-mers overlap by 30 bases (extended greedily over the store), so a string of n k-mers holds n + 30 bases, 2 bits each, instead of 31n. SSHash ([Pibiri 2022](https://doi.org/10.1093/bioinformatics/btac245)) then indexes the strings by their minimizers: a minimal perfect hash of the minimizers points to the strings (or string pieces) holding each, and a k-mer is found by hashing its minimizer and checking those few positions; its count is a weight stored along the strings. The dictionary is written beside the store as `<store>.sshash`.

## Required inputs

- `-o, --out <dir>`: the store (a folder)
- the genomes, as `--genome NAME=PATH` (repeatable) or `--genomes <tsv>` (lines of name and path): FASTA, plain, gzip or bgzip (with `.fai` and `.gzi` they are read by interval where needed), or a GFA of haplotypes
- `--coordinate-reference NAME`: the genome whose coordinates the locus graphs carry; `recruit` needs it

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `--genome NAME=PATH` | a genome to count; repeatable | |
| `--genomes <tsv>` | genomes as lines of name and path | |
| `--coordinate-reference NAME` | place this genome's anchors in every genome; on an existing store, adds them | |
| `--add` | count the genomes into a new segment of an existing store | off |
| `--segment-only <dir>` | count the genomes into `<dir>` without committing (a parallel job) | |
| `--attach <dir>...` | commit segments built with `--segment-only` | |
| `--merge` | merge the store's segments into one | off |
| `--range I/N` | with `--merge`: only the I-th of N partition ranges (0-based) | all |
| `--compact` | write the merged store as an SSHash dictionary | off |
| `--clean` | remove what no committed segment or running job uses | off |
| `--locate <region>` | print where each genome holds a region of the coordinate reference | |
| `--dump <path>` | every k-mer and its count, segments summed | |
| `--verify <path>` | check that another background (store or `.sshash`) holds the same counts | |
| `--partitions <N>` | partitions of a new store, a power of two | `4096` |
| `--memory <size>` | memory bound for counting, such as `16G` | `2G` |
| `--tmp-dir <dir>` | spill files | `$TMPDIR`, else inside the store |
| `-t, --threads <N>` | worker threads | the CPUs available |
| `-q, --quiet` | no progress logs | off |

## Outputs

| file | contents |
|------|----------|
| `<store>/manifest.json` | the committed segments, their genomes (path, size, checksum) and the coordinate references |
| `<store>/seg-N/` | a segment: its partitions (`parts-0.bin`, or one file per merge range), `segment.tsv` (its genomes and counts) and its anchor places per coordinate reference |
| `<store>/anchors/<reference>/` | the coordinate reference's anchors |
| `<store>.sshash` | with `--compact`: the dictionary `recruit` reads instead of the store while its fingerprint matches |
| `--dump <path>` | one line per k-mer: the k-mer and its count, tab-separated, no header |
| standard output, with `--locate` | one row per copy found in each genome (`genome`, `copy`, `contig`, `start`, `end`, `strand`, the anchors inside and around it, the chains joined); a genome holding no copy gets one row with copy 0 |

`recruit` refuses a store whose genome files changed since they were counted (their size is recorded): the counts would no longer describe the files it reads to locate a locus.

## Limitations

- k is 31 and the minimizer 11; both are fixed in the store format.
- A compacted dictionary cannot grow: add genomes to the store and compact again.
- Counting needs temporary disk of about a byte per genome base for the spill files, and several times that for sorted runs when `--memory` is small; put `--tmp-dir` on a local disk.
- Anchors find a locus in a genome only where the genome holds enough of the interval's anchors in order; a copy broken across contigs, or inverted, is found in pieces.
- A graph genome without PanSN path names counts as one haplotype, which lowers the shared cut of `recruit` for it.

## Example

`scripts/background_examples.sh` counts a synthetic genome set two ways, all at once and grown step by step (a genome added, another counted in its own job and attached, merged, compacted), checks that both hold the same counts, and locates the locus in every genome, in `tests/examples/background`.
