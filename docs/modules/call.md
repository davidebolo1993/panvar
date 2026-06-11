# Call Module (graph-native SV calling)

CLI entrypoint:

- `panvar call`

## What it does

`call` types structural variants directly on the pangenome graph. For each bubble it compares every
haplotype's source→sink walk to a **designated reference path's** walk and reads the differences off a
node-level alignment — no per-cluster minimap2/dotplot machinery. It then fights call fragmentation in
two ways and writes a tidy multi-sample VCF:

1. **Within a haplotype**, fragmented same-type events that sit close together are **coalesced** into one
   (`--merge-distance-bp`).
2. **Across haplotypes**, equivalent events are **merged** by a length-weighted node-set Jaccard
   (`--merge-jaccard`) into a single multi-sample record.

Input is expected to be a **panphorte-normalized GFA**, so a tandem duplication is a single repeat-unit
(`REP`) node traversed N times; copy number then falls straight out of the walk.

### Event types

- **DEL** — reference-only nodes (deleted from the haplotype).
- **INS** — haplotype-only nodes (inserted vs reference). With `--classify-ins`, minimap2 refines a
  subtype `INS_SUBTYPE=NOVEL|DUP` (does the inserted sequence map back to the local reference?). The
  primary `SVTYPE` stays `INS`.
- **INV** — a haplotype run that is the reverse-complement node-walk of a reference run.
- **DUP** — a `REP` node (self-loop) traversed a different number of times than the reference. `REF_CN`
  is the reference copy number; per-sample `CN` is reported in `FORMAT`.

## Required inputs

- `--gfa <graph.gfa>` / `-i` (ideally panphorte-normalized; re-`vg snarls` + `panvar bubble` it first)
- one of:
  - `--bubble-prefix-in <module1-prefix>` (auto uses `<module1-prefix>.bubbles.csv`)
  - `--bubbles-csv-in <module1.bubbles.csv>`
- `--reference-path <name>` — a path present in the GFA, used as the diff baseline
- `-o, --out-prefix <prefix>`

## Key options

- `--min-sv-bp <N>` — minimum event size to report (default `50`)
- `--merge-distance-bp <N>` — coalesce nearby same-type events within a bubble (default `100`)
- `--merge-jaccard <X>` — cross-haplotype node-set Jaccard threshold to merge events (default `0.80`)
- `--classify-ins` — refine INS subtype NOVEL/DUP via minimap2 (`--minimap-preset`, `--minimap-best-n`,
  `--ins-dup-min-identity`)
- `--bubble-id <N>` — restrict to one bubble (repeatable)
- `--no-per-bubble-vcf` — only write the concatenated region VCF
- `--quiet`

## Outputs

- `<prefix>.bubble_<id>.vcf` — one multi-sample VCF per bubble (unless `--no-per-bubble-vcf`)
- `<prefix>.region.vcf` — all bubble records concatenated

Both are VCF 4.2. Samples are the **haplotypes** (every P/W path; haploid). For each record:

- `FORMAT=GT`: `1` = carrier, `0` = traverses the bubble but reference-like, `.` = does not traverse it.
- `FORMAT=CN`: per-sample copy number (DUP records; `.` otherwise).
- `INFO`: `SVTYPE, SVLEN, END, BUBBLE_ID, START_NODE, END_NODE, EVENT_NODES` (the variant node set),
  `NMERGED` (carrier count), `INS_SUBTYPE` (refined INS only), `REF_CN` (DUP), and the event sequence
  (`INSSEQ`/`DELSEQ`/`INVSEQ`, omitted when very long).
- `CHROM`/`POS` come from the reference path's genomic label (e.g. `...chr6:31891045-...`) when
  parseable, else a graph-relative offset; the graph **node** ids in `START_NODE`/`END_NODE` are the
  exact anchors regardless of coordinate availability.

## Algorithm

For each bubble:

1. Take the reference walk and each haplotype's canonical source→sink walk
   (`canonical_bubble_path_steps`). Group identical walks into **distinct alleles** (call once per
   allele, expand genotypes by membership).
2. Mark copy-number nodes (any node traversed > 1× in the reference or an allele) and fold their
   consecutive self-repeats into one alignment token (so a `REP` self-loop becomes a single anchor and
   surfaces as a CN delta, not a spurious INS/DEL).
3. Align the reference and haplotype token walks (diagonal only on equal tokens, else gaps) and read off
   DEL / INS / INV / DUP. INV is detected when a haplotype gap-block is the reverse-complement node-walk
   of the opposing reference gap-block.
4. Coalesce nearby same-type events within the haplotype (`--merge-distance-bp`).
5. Merge events across haplotypes: same type and length-weighted node-set Jaccard `≥ --merge-jaccard`
   (DUP merges on shared `REP` node identity, with per-sample CN). Emit one VCF record per merged group.

## Example

```bash
# Normalize tandems, then re-snarl + re-bubble the normalized graph
panvar panphorte -i graph.flp.gfa --bubble-prefix-in out/bubble -o out/panphorte --min-similarity 0.90
vg snarls ... > out/panphorte.normalized.snarls.jsonl
panvar bubble -i out/panphorte.normalized.gfa -o out/bubble2 --snarls-in out/panphorte.normalized.snarls.jsonl

# Call structural variants
panvar call \
  -i out/panphorte.normalized.gfa \
  --bubble-prefix-in out/bubble2 \
  --reference-path "GRCh38#0#chr6:31891045-32123783" \
  -o out/call \
  --merge-distance-bp 100 --merge-jaccard 0.80 --classify-ins
```
