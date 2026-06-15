# Call Module (Module 3 — graph-native SV calling)

CLI entrypoint:

- `panvar call`

## What it does

`call` types structural variants on the pangenome graph. For each bubble it compares every
haplotype's source→sink walk to a **designated reference path's** walk and reads the differences off a
node-level alignment. It then fights call fragmentation in
two ways and writes a tidy multi-sample VCF:

1. **Within a haplotype**, fragmented same-type events that sit close together are **coalesced** into one
   (`--merge-distance-bp`). The gap is measured in **either** reference space **or** the haplotype's own
   sequence space — so two insertions that are far apart on the reference because a deletion sits between
   them, yet contiguous in the sample, still coalesce.
2. **Across haplotypes**, equivalent events are **merged** by **transitive connected-component**
   clustering: events are nodes, an edge joins two that share an anchor window and either overlap in node
   set (length-weighted Jaccard ≥ `--merge-jaccard`) **or** in sequence (identity ≥
   `--merge-seq-identity`), and each connected component becomes one record. Transitivity matters — if A
   matches B and B matches C but A does not match C, all three still collapse into one record instead of
   fragmenting across haplotypes (the old greedy first-fit could split them). The sequence key
   consolidates events that are the same biologically but thread different graph nodes (e.g. a
   microsatellite tangle); the largest member represents the record.
3. **Graph-level force-call (sub-threshold rescue)**: events are clustered *before* the size filter (down
   to `--rescue-min-bp`), then every non-carrier haplotype is interrogated at each called locus against
   its **own walk-diff** — so a 49 bp deletion in one haplotype is rescued (`GT=1`) by a 51 bp call in
   another instead of being dropped and genotyped `0`. Because the representative is fixed and carriers
   only accrue, this is a single pass (no re-iteration) and adds no re-alignment. A merged record is
   reported only if its representative (largest member) reaches `--min-sv-bp`.

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

- `--gfa <graph.gfa>` / `-i` (ideally panphorte-normalized; use `panphorte --reference-path`'s
  `.normalized.sorted.gfa` so it is already sorted + re-snarled — no `vg`/`odgi` needed)
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
- `--merge-distance-bp <N>` — coalesce nearby same-type events within a bubble (default `100`). The gap
  is checked in **both** reference and haplotype sequence space (whichever is closer wins), and it also
  sets the base width of the cross-haplotype merge window.
- `--merge-jaccard <X>` — cross-haplotype node-set Jaccard threshold to merge events (default `0.80`)
- `--merge-seq-identity <X>` — cross-haplotype event-sequence identity to merge (default `0.80`)
- `--merge-size-ratio <X>` — length-ratio floor for the sequence merge (default `0` = use
  `--merge-seq-identity`). The sequence merge only compares two events whose shorter/longer length ratio
  clears this floor. **Lower it** to merge same-locus, same-motif events of *different sizes* — e.g. the
  several different-length deletions of one STR/microsatellite, or insertion size-classes — into a single
  record (the record then carries `INFO=SVLEN_RANGE`). The default keeps distinct sizes separate.
- `--min-haplotypes <N>` — drop records carried by fewer than N haplotypes (default `1` = off)
- `--min-maf <X>` — drop records with allele frequency `AF = AC/AN` below X, where `AN` is the number of
  haplotypes that **traverse the bubble** (so `.`-genotyped haplotypes are excluded from the denominator).
  Complements the count-based `--min-haplotypes`. Default `0` = off.
- `--rescue-min-bp <N>` — floor for sub-threshold events kept for merge/rescue (default `min-sv-bp/2`)
- `--classify-ins` — refine INS subtype NOVEL/DUP via minimap2 (`--minimap-preset`, `--minimap-best-n`,
  `--ins-dup-min-identity`)
- `--cn-from-multiplicity` — emit `DUP` from peak node multiplicity for folded bubbles with no self-loop
  (e.g. GSTM1) that panphorte left intact (see Event types and the section below)
- `--cn-from-coverage` — emit **total-module copy number** on folded paralog clusters where the reference
  itself traverses the module ≥2× (e.g. CYP2D6/2D7); see the copy-number section below
- `--multiallelic-loci` — collapse a bounded locus (e.g. an STR/VNTR) into ONE multiallelic record
  (`REF` + `ALT1,ALT2,…` explicit sequences, per-sample `GT` indexing the allele); `--multiallelic-max-bp`
  (default 5000) bounds it so large SVs keep their typed per-event records
- `--bubble-id <N>` — restrict to one bubble (repeatable)
- `--no-per-bubble-vcf` — only write the concatenated region VCF
- `--no-variant-paths` — skip the `<prefix>.variant_paths.tsv` provenance sidecar
- `--quiet` — disable the per-bubble progress bar (printed to stderr; the run summary stays on stdout)

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
- A **folded paralog cluster** where the paralogs are too similar to separate (e.g. CYP2D6/2D7, ~95%
  identical) collapses into one tangled bubble in which the *reference itself* traverses the module ≥2×.
  Here peak-multiplicity fails (a gene deletion need not touch the single peak node, so deletions stay at
  `REF_CN`). `--cn-from-coverage` instead reads **total module copy number** from coverage:
  `copies ≈ (bp the haplotype spells through the bubble) / (module unit bp)`, with the unit derived
  graph-only as `ref_spelled_bp / ref_fold`. Because it uses all the traversed sequence it recovers
  **deletions** (fewer bp → fewer copies) as well as gains. On CYP2D6 it separates the 9 deletions
  (`CN=2`) from normal (`CN=3`) from gains/trip (`CN≥4`) — monotonic with the pangene total. It is
  reference-relative and reports the **total module** copy number, not a per-gene (2D6-vs-2D7) count;
  per-subtype resolution is a future refinement. When it fires for a bubble it replaces that bubble's
  (folded, unreliable) walk-diff with the single CN record.

  **The three CN detectors are complementary, applied per bubble in this precedence:** (1) panphorte
  **self-loop DUP** (clean adjacent tandems → exact loop count); else (2) `--cn-from-coverage` total-module
  CN on folded clusters (`ref_fold ≥ 2`); else (3) `--cn-from-multiplicity` peak DUP (a single extra copy
  folded onto a high-multiplicity node, e.g. GSTM1). Note: in a bubble panphorte *partially* collapsed (a
  self-loop exists), the peak/coverage paths are skipped, so any residual un-collapsed copies surface as
  `INS`/`INS_SUBTYPE=DUP` rather than a DUP.

  **Passing both `--cn-from-coverage` and `--cn-from-multiplicity` is safe and gives the widest recall** —
  they do not conflict, they *compose* per bubble. Coverage fires only where the reference itself folds
  (`ref_fold ≥ 2`) and, where it does, takes precedence and replaces that bubble's walk-diff; peak
  multiplicity then handles the remaining folded bubbles where the reference does **not** fold. So coverage
  catches the tangled paralog clusters (CYP2D6) and peak multiplicity catches the single-extra-copy folds
  (GSTM1), with no double counting.

  **Roadmap — per-paralog resolution.** Both folded-cluster detectors report **total module** copy number,
  not per-paralog (2D6 vs 2D7, C4A vs C4B), because pggb folds the conserved backbone onto shared nodes.
  The intended route to per-copy DEL/INS/CN is an **optional pangene/segdup BED** that *unfolds* the cluster
  — labelling which nodes belong to which paralog so the cluster can be split and each copy typed
  independently. That is deliberately not built yet; a dedicated dispersed-segdup module was considered and
  set aside in favor of the BED-driven unfold.
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
  Carries each node's bp length plus reference/CN flags.
- `<prefix>.variant_nodes.tsv` — per-variant node set (unless `--no-variant-paths`): one row per variant —
  `variant_id, bubble_id, svtype, node_ids` (the deduplicated `EVENT_NODES`). The handoff for the
  `describe` module: restrict k-mer markers to the nodes that participate in called variation.

The VCFs are VCF 4.2. Samples are the **haplotypes** (every P/W path; haploid). For each record:

- `FORMAT=GT`: `1` = carrier, `0` = traverses the bubble but reference-like, `.` = does not traverse it.
- `FORMAT=CN`: per-sample copy number (DUP records; `.` otherwise).
- `INFO`: `SVTYPE, SVLEN, END, BUBBLE_ID, START_NODE, END_NODE, EVENT_NODES`, `NMERGED` (carrier count),
  population fields **`NS`** (haplotypes with data = traversing the bubble), **`AN`** (= NS), **`AC`**
  (carrier count) and **`AF`** (= AC/AN, the allele frequency `--min-maf` filters on),
  `SVLEN_RANGE` (min,max member size when a merged record spans differing sizes), `EVENTID` (links a
  co-located DEL+INS substitution), `INS_SUBTYPE` (refined INS only), `REF_CN` (DUP), `NALLELES`
  (multiallelic records from `--multiallelic-loci`), and the event
  sequence (`INSSEQ`/`DELSEQ`/`INVSEQ`, omitted when very long). `EVENT_NODES` is **deduplicated and
  ordered** (reference nodes by genomic position, haplotype-only nodes after them in walk order), so
  `START_NODE`→`END_NODE` reads as a coherent progression; the full per-carrier node walk lives in
  `variant_paths.tsv`. A substitution (same locus deletes a reference block and inserts a haplotype block)
  yields two records (DEL + INS) sharing one `EVENTID`.
- `CHROM`/`POS` come from the reference path's genomic label (e.g. `...chr6:31891045-...`) when
  parseable, else a graph-relative offset.

## Plotting

`scripts/plot_sv_map.R` draws a **node-level structural-variant map** in the same `odgi viz` style as the
`inspect` node/edge-coverage heatmaps: rows = haplotypes, **columns = the bubble's internal nodes in
reference-sorted order** (the same node set and order as the `inspect` `node_counts.tsv` columns), each
cell colored by the **called event** the haplotype carries on that node. Because each event is painted
only on its own nodes, an event occupies just its node columns (in genomic order) — a deletion reads as a
contiguous colored gap, an insertion as a block at its anchor — and the figure lines up 1:1 with the
coverage heatmap of the same bubble.

Cell colors: white = node not traversed; grey = traversed, reference-like; **DEL** red, **INS-NOVEL**
green, **INS-DUP** purple, **INV** orange, **multiallelic** teal (shaded by allele index). A **DUP**
carrier's whole traversed module is shaded **blue by per-node multiplicity** (the local copy count read
from the coverage substrate) — so a copy-number gene like GSTM1 or CYP2D6 shows its dosage as blue
intensity across the module, with DEL/INS painted on top. (DUP is encoded with a single representative
node in the VCF, so this whole-module shading is what makes the copy number visible.)

It reads the per-bubble `inspect` `node_counts.tsv` (substrate + node order) and the `call` VCF (events,
carriers, `GT:CN`). The bubble id is inferred from the `node_counts.tsv` filename (or pass `--bubble-id N`).
Optional: `--node-lengths <...node_lengths.tsv>` scales column widths by node bp; `--clusters <...clusters.tsv>`
orders and color-bars the rows by walk cluster (to see whether clusters track the calls);
`--reference-path NAME` pins that haplotype (substring match) as the top row; `--max-nodes`/`--max-paths`
cap the figure.

```bash
panvar inspect -i graph.gfa --bubble-prefix-in out/bubble --bubble-id 4 --cluster -o out/inspect/b4
Rscript scripts/plot_sv_map.R \
  --node-counts out/inspect/b4.bubble_4.node_counts.tsv \
  --vcf out/call.region.vcf \
  --node-lengths out/inspect/b4.bubble_4.node_lengths.tsv \
  --clusters out/inspect/b4.bubble_4.clusters.tsv \
  --reference-path grch38 --out out/sv_map_b4
```

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
4. Coalesce nearby same-type events within the haplotype, merging when the gap is within
   `--merge-distance-bp` in **reference space OR this haplotype's own sequence space** (the latter
   catches events split on the reference by an intervening deletion but contiguous in the sample).
5. Merge events across haplotypes by **transitive single-linkage clustering** (union-find connected
   components; keeping all events down to `--rescue-min-bp`). Two events are linked when they share a
   type, sit within a position window, and have node-set Jaccard `≥ --merge-jaccard` **or** sequence
   identity `≥ --merge-seq-identity`; each component is one record. The position window **scales with
   event size** (`--merge-distance-bp` plus the smaller event's length), so a large insertion/deletion
   whose breakpoint floats by kilobases across haplotypes still merges — the sequence/Jaccard gate
   prevents over-merging distinct events, and connected components make the merge order-independent and
   transitive (A–B–C collapse even when A and C do not directly match). Edge-building is a
   position-sorted windowed sweep, so it stays near-linear. The largest member represents the record.
   DUP merges on shared `REP` node identity, with per-sample CN (handled separately from this clustering).
6. **Graph-level force-call**: interrogate every non-member haplotype at each merged locus against its
   own walk-diff (including its sub-threshold events) to set `GT` and rescue sub-threshold carriers — a
   single pass, since the representative is fixed and carriers only accrue. Keep records whose
   representative reaches `--min-sv-bp` and carrier count reaches `--min-haplotypes`, and emit one VCF
   record per group.

## Example

No external tools (`odgi`/`vg`) are needed: `bubble` and `panphorte` sort and snarl internally. Tandem
normalization with `panphorte --reference-path` repositions the appended REP nodes along the reference
(so `numeric node id == reference order` again) and re-snarls, emitting a call-ready graph + bubbles.

```bash
# 1. sites on the raw graph (internal sort + cactus snarls)
panvar bubble -i graph.gfa --reference-path grch38 -o out/bubble

# 2. normalize tandems + internally re-sort + re-snarl (writes normalized.sorted.gfa + bubbles.csv)
panvar panphorte -i out/bubble.sorted.gfa --bubble-prefix-in out/bubble \
  --reference-path grch38 -o out/panphorte --min-similarity 0.90

# 3. call structural variants on the panphorte sorted output
panvar call \
  -i out/panphorte.normalized.sorted.gfa \
  --bubble-prefix-in out/panphorte \
  --reference-path grch38 \
  -o out/call \
  --merge-distance-bp 100 --merge-jaccard 0.80 --classify-ins
```
