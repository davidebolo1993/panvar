# Call Module (graph-native SV calling)

CLI entrypoint:

- `panvar call`

## What it does

`call` types structural variants directly on the pangenome graph. For each bubble it compares every
haplotype's source→sink walk to a **designated reference path's** walk and reads the differences off a
node-level alignment. It then fights call fragmentation in
two ways and writes a tidy multi-sample VCF:

1. **Within a haplotype**, fragmented same-type events that sit close together are **coalesced** into one
   (`--merge-distance-bp`).
2. **Across haplotypes**, equivalent events are **merged** when they share an anchor window and either
   their node sets overlap (length-weighted Jaccard ≥ `--merge-jaccard`) **or** their sequences are
   similar (identity ≥ `--merge-seq-identity`). The sequence key consolidates events that are the same
   biologically but thread different graph nodes (e.g. a microsatellite tangle).
3. **Sub-threshold rescue**: events are merged *before* the size filter (down to `--rescue-min-bp`), and
   a joint re-scan checks every haplotype at each called locus — so a 49 bp deletion in one haplotype is
   rescued (`GT=1`) by a 51 bp call in another instead of being dropped and genotyped `0`. A merged
   record is reported only if its representative (largest member) reaches `--min-sv-bp`.

Input is expected to be a **panphorte-normalized GFA**, so a tandem duplication is a single repeat-unit
(`REP`) node traversed N times; copy number then falls straight out of the walk.

### Event types

- **DEL** — reference-only nodes (deleted from the haplotype).
- **INS** — haplotype-only nodes (inserted vs reference). With `--classify-ins`, minimap2 refines a
  subtype `INS_SUBTYPE=NOVEL|DUP` (does the inserted sequence map back to the local reference?). The
  primary `SVTYPE` stays `INS`.
- **INV** — a haplotype run that is the reverse-complement node-walk of a reference run.
- **DUP** — a copy-number gain. Two sources:
  - a `REP` node (self-loop) traversed a different number of times than the reference (panphorte's
    collapsed tandem arrays);
  - with `--cn-from-multiplicity`, a **folded duplication panphorte left intact** (no self-loop):
    a bubble where a haplotype's **peak node-traversal multiplicity** exceeds the reference's peak — the
    extra paralog copy is folded onto shared nodes, so it is traversed once more than the reference does.
    `SVLEN` is the duplicated content (Σ node_len × excess traversals); the peak (not per-node excess) is
    what isolates real gene dosage from cluster background. This is how the GSTM1 gene duplication, buried
    in a ~77 kb segmental-duplication cluster with no adjacent repeat, is called.

  Either way `REF_CN` is the reference copy number and per-sample `CN` is reported in `FORMAT`.

## Required inputs

- `--gfa <graph.gfa>` / `-i` (ideally panphorte-normalized; re-`vg snarls` + `panvar bubble` it first)
- one of:
  - `--bubble-prefix-in <module1-prefix>` (auto uses `<module1-prefix>.bubbles.csv`)
  - `--bubbles-csv-in <module1.bubbles.csv>`
- `--reference-path <name>` — the path used as the diff baseline. Accepts either the full path
  name or a **case-insensitive substring** (e.g. `grch38` → `grch38#1#chr6:...`). An exact name always
  wins; otherwise the substring must match exactly one path, else `call` errors — "not found" or
  "ambiguous" with the candidate list (e.g. `grch38` when both `GRCh38_0` and `grch38_1` exist; pass the
  unambiguous `GRCh38_0`)
- `-o, --out-prefix <prefix>`

## Key options

- `--min-sv-bp <N>` — minimum size of a reported (merged) event (default `50`)
- `--merge-distance-bp <N>` — coalesce nearby same-type events within a bubble (default `100`)
- `--merge-jaccard <X>` — cross-haplotype node-set Jaccard threshold to merge events (default `0.80`)
- `--merge-seq-identity <X>` — cross-haplotype event-sequence identity to merge (default `0.80`)
- `--merge-size-ratio <X>` — length-ratio floor for the sequence merge (default `0` = use
  `--merge-seq-identity`). The sequence merge only compares two events whose shorter/longer length ratio
  clears this floor. **Lower it** to merge same-locus, same-motif events of *different sizes* — e.g. the
  several different-length deletions of one STR/microsatellite, or insertion size-classes — into a single
  record (the record then carries `INFO=SVLEN_RANGE`). The default keeps distinct sizes separate.
- `--min-haplotypes <N>` — drop records carried by fewer than N haplotypes (default `1` = off)
- `--rescue-min-bp <N>` — floor for sub-threshold events kept for merge/rescue (default `min-sv-bp/2`)
- `--classify-ins` — refine INS subtype NOVEL/DUP via minimap2 (`--minimap-preset`, `--minimap-best-n`,
  `--ins-dup-min-identity`)
- `--cn-from-multiplicity` — emit `DUP` from peak node multiplicity for folded bubbles with no self-loop
  (e.g. GSTM1) that panphorte left intact (see Event types and the section below)
- `--bubble-id <N>` — restrict to one bubble (repeatable)
- `--no-per-bubble-vcf` — only write the concatenated region VCF
- `--no-variant-paths` — skip the `<prefix>.variant_paths.tsv` provenance sidecar
- `--quiet`

### Copy-number gains/losses and segmental duplications

- A **tandem array** (panphorte-normalized into a `REP` self-loop) is a `DUP` record with per-sample
  `CN`. A **copy loss** at such a locus shows up as a sample whose `CN` differs from `REF_CN`.
- A **non-tandem extra copy** (a segmental duplication that panphorte did not collapse) can be read two
  ways. With `--cn-from-multiplicity` it becomes a true **`DUP`** record: the copy is folded onto shared
  nodes, so the carrier's peak node multiplicity exceeds the reference's and copy number falls out of the
  walk directly — baseline-corrected, no re-alignment. On the bundled GSTM1 locus (a ~77 kb cluster, no
  adjacent repeat to seed) this emits one `DUP` (`SVLEN≈18.5 kb`, `REF_CN=3`) carried by exactly the two
  copy-2 samples (HG01346, NA19240), `CN=4`, with zero false positives across the other 463 haplotypes.
  Without the flag the same event surfaces as an **INS**; with `--classify-ins` it is labelled
  `INS_SUBTYPE=DUP` because the inserted sequence maps back to the local reference. When
  `--cn-from-multiplicity` emits the DUP, the walk-diff's redundant view of the same extra copy (a
  "duplication insertion" carried by exactly the DUP's carriers, of comparable size) is **dropped**, so
  the event is reported once — as the DUP — not double-counted as both a DUP and an INS.
- Note copy number is **reference-relative**: presence/absence reads as `INS` against a reference that
  lacks the copy and as `DEL` against one that has it. For an absolute per-haplotype count, read CN from
  panphorte's `copies.tsv` / the `DUP` record on a relaxed-similarity normalization (see the panphorte
  docs), not from a single diff against one reference.

## Outputs

- `<prefix>.region.vcf` — all bubble records, **coordinate-sorted** (POS, then END, then ID) with
  **unique IDs**, so it is directly `bgzip` + `tabix -p vcf` / `bcftools index`-able.
- `<prefix>.bubble_<id>.vcf` — one multi-sample VCF per bubble, also sorted (unless `--no-per-bubble-vcf`)
- `<prefix>.variant_paths.tsv` — per-variant path provenance (unless `--no-variant-paths`): one row per
  (variant, carrier haplotype) — `variant_id, bubble_id, svtype, sample, gt, sub_walk` — where `sub_walk`
  is that carrier's realized walk through the event (between the flanking reference nodes) as a GFA-style
  `>node`/`<node` string. Joins 1:1 to the VCF `ID`; carrier rows reconcile with `NMERGED`. This is the
  interpretable bridge for the `describe` module and manual inspection.
- `<prefix>.node_track.tsv` — per-bubble x-axis for plotting (unless `--no-variant-paths`): one row per
  inside node — `bubble_id, order, node_id, length_bp, genomic_pos, in_reference, is_cn` — giving **every**
  inside node of the bubble (the same set and order as the `inspect` `node_counts.tsv` columns), ordered by
  **numeric node id**, which is the reference order because the graph is odgi-sorted along the reference.
  Carries each node's bp length plus reference/CN flags. Feeds the SV colormap below.

## Plotting the calls

`scripts/plot_sv_map.R` draws a per-bubble structural-variant colormap in the style of the node/edge
coverage heatmaps: rows = haplotypes, x = the bubble's nodes ordered by the reference and **length-scaled**,
each cell colored by the SV the haplotype carries there relative to the reference (DEL / INV / INS-NOVEL /
INS-DUP / DUP shaded by copy number), grey where it traverses without a call, white where it does not
traverse the node. It joins this module's VCF + `node_track.tsv` with an `inspect` per-bubble
`node_counts.tsv` (the traversal matrix):

```bash
panvar inspect -i graph.gfa --bubble-prefix-in out/bubble2 --bubble-id 4 -o out/inspect/b4
scripts/plot_sv_map.R \
  --vcf out/call.region.vcf \
  --node-counts out/inspect/b4.bubble_4.node_counts.tsv \
  --node-track out/call.node_track.tsv \
  --bubble-id 4 --out out/sv_bubble_4
```

The VCFs are VCF 4.2. Samples are the **haplotypes** (every P/W path; haploid). For each record:

- `FORMAT=GT`: `1` = carrier, `0` = traverses the bubble but reference-like, `.` = does not traverse it.
- `FORMAT=CN`: per-sample copy number (DUP records; `.` otherwise).
- `INFO`: `SVTYPE, SVLEN, END, BUBBLE_ID, START_NODE, END_NODE, EVENT_NODES`, `NMERGED` (carrier count),
  `SVLEN_RANGE` (min,max member size when a merged record spans differing sizes), `EVENTID` (links a
  co-located DEL+INS substitution), `INS_SUBTYPE` (refined INS only), `REF_CN` (DUP), and the event
  sequence (`INSSEQ`/`DELSEQ`/`INVSEQ`, omitted when very long). `EVENT_NODES` is **deduplicated and
  ordered** (reference nodes by genomic position, haplotype-only nodes after them in walk order), so
  `START_NODE`→`END_NODE` reads as a coherent progression; the full per-carrier node walk lives in
  `variant_paths.tsv`. A substitution (same locus deletes a reference block and inserts a haplotype block)
  yields two records (DEL + INS) sharing one `EVENTID`.
- `CHROM`/`POS` come from the reference path's genomic label (e.g. `...chr6:31891045-...`) when
  parseable, else a graph-relative offset.

## Algorithm

For each bubble:

1. Take the reference walk and each haplotype's canonical source→sink walk
   (`canonical_bubble_path_steps`). Group identical walks into **distinct alleles** (call once per
   allele, expand genotypes by membership).

   **Empty-interior alleles (single-node indels).** The shared bubble-path machinery
   (`find_best_bubble_path_interval`) locates a path's crossing of a bubble by requiring it to traverse
   at least one **interior** node (`inside_count >= 1`) — sensible for `inspect`/`panphorte`, which act
   on interior content. But a *pure* deletion's allele goes source→sink directly (no interior node), and
   a *pure* insertion's reference allele does too — so that allele scored `inside_count = 0` and was read
   as "does not traverse the bubble." The consequence: a single-node deletion lost its deleting allele
   entirely (no DEL emitted), and a single-node insertion made the reference itself look absent (the whole
   bubble was skipped). Complex bubbles hid this, because their alternate alleles still traverse other
   interior nodes — so it only bit clean two-anchor / one-variable-node sites (common in STR/microsatellite
   regions). `call` adds a local fallback (`bubble_steps`): when the interval finder returns nothing but
   the path has the bubble's source immediately adjacent to its sink, it uses the empty-interior allele
   `[source, sink]`. This is confined to the caller; the shared interval finder is unchanged.
2. Mark copy-number nodes (any node traversed > 1× in the reference or an allele) and fold their
   consecutive self-repeats into one alignment token (so a `REP` self-loop becomes a single anchor and
   surfaces as a CN delta, not a spurious INS/DEL).
3. Align the reference and haplotype token walks and read off DEL / INS / INV. To scale to large bubbles
   (e.g. a ~2500-node gene-presence locus), the two walks are first split at shared **anchor** nodes
   (tokens appearing exactly once in both walks, chained in a common monotonic order); the quadratic DP
   then runs only inside the small segments between anchors. This removes the whole-bubble size cap, and
   because the breakpoints fall on the same anchors for every haplotype, the same large event (a gene
   deletion/insertion) is identical across haplotypes and merges cleanly. INV is detected when a
   haplotype gap-block is the reverse-complement node-walk of the opposing reference gap-block.
4. Coalesce nearby same-type events within the haplotype (`--merge-distance-bp`).
5. Merge events across haplotypes (keeping all events down to `--rescue-min-bp`): same type, anchors
   within a window, and node-set Jaccard `≥ --merge-jaccard` **or** sequence identity
   `≥ --merge-seq-identity`. The position window **scales with event size** (`--merge-distance-bp` plus
   the smaller event's length), so a large insertion/deletion whose breakpoint floats by kilobases
   across haplotypes still merges — the sequence/Jaccard gate prevents over-merging distinct events.
   The largest member represents the record. DUP merges on shared `REP` node identity, with per-sample CN.
6. Joint re-scan every haplotype at each merged locus to set `GT` (rescuing sub-threshold carriers),
   keep records whose representative reaches `--min-sv-bp` and carrier count reaches `--min-haplotypes`,
   and emit one VCF record per group.

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
