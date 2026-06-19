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
   fragmenting across haplotypes. The sequence key
   consolidates events that are the same biologically but thread different graph nodes (e.g. a
   microsatellite tangle); the largest member represents the record.
3. **Graph-level force-call (sub-threshold rescue)**: events are clustered *before* the size filter (down
   to `--rescue-min-bp`), then every non-carrier haplotype is interrogated at each called locus against
   its **own walk-diff** — so a 49 bp deletion in one haplotype is rescued (`GT=1`) by a 51 bp call in
   another instead of being dropped and genotyped `0`. Because the representative is fixed and carriers
   only accrue, this is a single pass. A merged record is
   reported only if its representative (largest member) reaches `--min-sv-bp`.

Which graph to call on depends on the locus topology (see the copy-number section):

- **Tandem-repeat regions** (e.g. LPA KIV-2): call on the **panphorte-normalized GFA**, so the variable
  tandem is a single repeat-unit (`REP`) node traversed N times and copy number falls straight out of the
  walk.
- **PGGB-collapsed paralog clusters / CNV loci** (e.g. C4, CYP2D6, GSTM1): call on the **`bubble` sorted
  GFA** (the unfolded graph). Here the copies live as **node multiplicity** on shared nodes;
  `panphorte` folding would hide that signal, so it is *not* used as the call substrate for these — its
  collapse never feeds their variant calling.

### Event types

- **DEL** — reference-only nodes (deleted from the haplotype).
- **INS** — haplotype-only nodes (inserted vs reference). With `--classify-ins`, minimap2 refines a
  subtype `INS_SUBTYPE=NOVEL|DUP` (does the inserted sequence map back to the local reference?). The
  primary `SVTYPE` stays `INS`.
- **INV** — a haplotype run that is the reverse-complement node-walk of a reference run.
- **DUP** — a copy-number gain or loss. Three sources (tried per bubble in a fixed precedence; see the
  copy-number section below):
  - a `REP` node (self-loop) traversed a different number of times than the reference (panphorte's
    collapsed tandem arrays);
  - with `--cn-from-coverage`, a **PGGB-collapsed paralog cluster** the reference itself traverses ≥2× —
    copy number is the total sequence each haplotype spells over the **full bubble walk** divided by one
    copy's bp, so it recovers losses as well as gains, and it counts copies that are collapsed onto shared
    nodes (interleaved, not a contiguous tandem block);
  - with `--cn-from-multiplicity`, a **folded duplication panphorte left intact** (no self-loop) where the
    reference does *not* fold: a haplotype's **peak node-traversal multiplicity** exceeds the reference's
    peak. `SVLEN` is the duplicated content (Σ node_len × excess traversals); the peak (not per-node
    excess) isolates real dosage from cluster background.

  In every case `REF_CN` is the reference copy number and per-sample `CN` is reported in `FORMAT`.

## Required inputs

- `--gfa <graph.gfa>` / `-i` (ideally panphorte-normalized; use `panphorte --reference-path`'s
  `.normalized.sorted.gfa`)
- one of:
  - `--bubble-prefix-in <panphorte-prefix>` (auto uses `<panphorte-prefix>.bubbles.csv`)
  - `--bubbles-csv-in <panphorte-prefix.bubbles.csv>`
- `--reference-path <name>` — the path used as the diff baseline. Accepts either the full path
  name or a **case-insensitive substring** of it (a short unique tag, e.g. a sample or contig name). An
  exact name always wins; otherwise the substring must match exactly one path, else `call` errors —
  "not found" or "ambiguous" with the candidate list
- `-o, --out-prefix <prefix>`

## Key options


```bash
panvar call -i <graph.gfa> (-b <prefix> | -c <bubbles.csv>) -r <name> -o <prefix> [options]
```

Running `panvar call` with no arguments prints this help. 

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
- `--min-maf <X>` — drop records whose carrier frequency `AF = AC/AN` is below X, where `AN` is the number
  of haplotypes that **traverse the bubble** (`.`-genotyped haplotypes are excluded from the denominator).
  `AF` is the **ALT (carrier) frequency and is not folded**, so this drops near-absent variants but keeps
  near-fixed ones (`AF≈0.99` passes a `0.05` cut). Complements the count-based `--min-haplotypes`.
  Default `0` = off.
- `--rescue-min-bp <N>` — floor for sub-threshold events kept for merge/rescue (default `min-sv-bp/2`)
- `--classify-ins` — refine INS subtype NOVEL/DUP via minimap2 (`--minimap-preset`, `--minimap-best-n`,
  `--ins-dup-min-identity`)
- `--cn-from-multiplicity` — emit `DUP` from peak node multiplicity for folded bubbles with no self-loop that panphorte left intact (see Event types and the section below)
- `--cn-from-coverage` — emit **total-module copy number** on folded paralog clusters where the reference itself traverses the module ≥2×; see the copy-number section below
- `--multiallelic-loci` — collapse a bounded locus (e.g. an STR/VNTR) into ONE multiallelic record `REF` + `ALT1,ALT2,…` explicit sequences, per-sample `GT` indexing the allele); `--multiallelic-max-bp`
  (default 5000) bounds it so large SVs keep their typed per-event records
- `--bubble-id <N>` — restrict to one bubble (repeatable)
- `--no-per-bubble-vcf` — only write the concatenated region VCF
- `--no-variant-paths` — skip the `<prefix>.variant_paths.tsv` provenance sidecar
- `--quiet` — disable the per-bubble progress bar (printed to stderr; the run summary stays on stdout)

### Copy-number gains/losses and segmental duplications

Copy number is read off the graph in three ways. Which one applies depends on how the locus is
represented after `panphorte`, so the detectors are tried **per bubble in a fixed precedence** and never
double-count.

**1. Clean tandem array → self-loop `DUP` (always on).** When `panphorte` collapsed an adjacent tandem
into a `REP` self-loop, copy number is the exact loop count: a `DUP` record with `REF_CN` and per-sample
`CN`. A **copy loss** is simply a sample whose `CN` is below `REF_CN`; a **gain** is one above it.
`panphorte` folds **single copies** of the unit too (within a confirmed array), so a one-copy haplotype
traverses the `REP` node once and reads `CN = 1` — not `0`. The only haplotype that reads `CN = 0` is one
`panphorte` could not fold at all (its unit is too divergent to align to the consensus at
`--min-similarity`); that is an honest "unresolved", not a deletion call.

**2. PGGB-collapsed paralog cluster → total-module CN from coverage (`--cn-from-coverage`).** When PGGB
collapses identical paralog copies (e.g. the C4 long-long / short-short RCCX modules, or CYP2D6/2D7/2D8P)
onto **shared nodes**, the copies are carried as **node multiplicity** — a 2-copy haplotype re-traverses
those nodes twice — not as a contiguous tandem block in the spelled sequence, and the **reference itself**
traverses the module two or more times. Copy number is read from how much sequence each haplotype spells
over the **full bubble walk**, normalised to one copy:

```text
copies      ≈ (bp the haplotype spells over the FULL bubble walk) / (one-copy bp)
one-copy bp  = (bp the reference spells over the full walk) / (times the reference folds over it)
```

The **full walk** is the widest source→sink span with **all repeats included** — this is the crux. The
ordinary (minimal-span) bubble walk covers each distinct inside node once and so collapses the repeated
copies onto a single traversal, flattening every haplotype to the same bp; counting over the full span
preserves multiplicity == copy number. Because it uses *all* the traversed sequence it recovers
**losses** (fewer bp → fewer copies) as well as **gains**, monotonically, and works whether or not the
collapse left a self-loop node. It reports the **total module** copy number, not a per-paralog count. When
it fires for a bubble it is the **authority** for that bubble — the self-loop and walk-diff paths are
skipped. The per-sample `CN` is written for **every** haplotype that traverses the module (its absolute
module count), but `GT=1` (and therefore `AC`/`AF`) marks only the **carriers** — haplotypes whose count
differs from `REF_CN` (a gain or a loss) — so a copy-invariant module is not emitted as a variant.
Validated: C4 total CN 131/131 exact; CYP2D6 concordant against the collapsed D6+D7 truth (the residual
misses carry an extra unannotated CYP2D8P / 2D7-hybrid module the gene BED does not count).

**3. Single folded extra copy → peak-multiplicity `DUP` (`--cn-from-multiplicity`).** When an extra copy
was folded onto shared nodes but the *reference* does **not** fold (1 or no copy), the carrier's peak
node-traversal multiplicity exceeds the reference's, and copy number falls straight out of the walk — no
re-alignment. `SVLEN` is the duplicated content (Σ node_len × excess traversals). Keying on the **peak**
multiplicity, rather than any per-node excess, separates true gene dosage from cluster background (per-node
excesses only reflect which paralog is present). When this emits a `DUP`, the walk-diff's redundant view
of the same copy — a "duplication insertion" carried by exactly the DUP's carriers, of comparable size —
is **dropped**, so the event is reported once.

**Precedence and composition.** Per bubble: (1) **coverage CN** if `--cn-from-coverage` and the reference
folds over the full walk (≥2×) — it is the authority and the other two are skipped; else (2) self-loop
`DUP` if a `REP` self-loop exists; else (3) peak-multiplicity `DUP` (`--cn-from-multiplicity`). Passing
both CN flags is safe and gives the widest recall — they *compose* rather than conflict: coverage handles
the bubbles where the reference folds, peak multiplicity handles folded duplications the reference does
not share, with no overlap.

**When no CN flag fires**, an extra copy that `panphorte` left intact surfaces through the ordinary
walk-diff as an **`INS`**; with `--classify-ins` it is labelled `INS_SUBTYPE=DUP` when the inserted
sequence maps back to the local reference.

**Absolute vs reference-relative.** The DEL/INS walk-diff is reference-relative (presence/absence reads as
`INS` against a reference that lacks the copy, `DEL` against one that has it). The two folded-cluster
detectors instead report an **absolute** per-haplotype copy number: `--cn-from-coverage` divides the
haplotype's full-walk bp by one copy's bp, so the reported `CN` is the haplotype's own module count, not a
difference against the reference. The reference choice (e.g. GRCh38) sets the unit-bp denominator but does
not change a haplotype's count — picking a different reference yields the same per-haplotype `CN`.

> **Per-paralog resolution.** Both folded-cluster detectors report total-module copy number
> because PGGB folds the conserved backbone onto shared nodes. They cannot say *which* paralog a copy is
> (e.g. CYP2D6 vs an extra CYP2D8P): the residual CYP2D6 mismatches are haplotypes carrying an unannotated
> CYP2D8P / 2D7-hybrid module that the total-module count includes but the gene BED does not. The intended
> route to per-copy typing is an optional pangene BED that *unfolds* the cluster — labelling which nodes
> belong to which paralog. In progress.

### Multiallelic loci (`--multiallelic-loci`)

By default every event at a bubble is its own VCF record (one DEL, one INS, …). At a small, bounded locus
that varies mainly by **which sequence** a haplotype carries — e.g. a short tandem repeat with several
length alleles — that scatters one site across many records. `--multiallelic-loci` instead collapses such
a bubble into a **single record** with explicit sequences: `REF` plus `ALT1,ALT2,…`, one per distinct
interior spelling, with per-sample `GT` indexing the allele each haplotype carries (`NALLELES` counts
them). This is **not** copy-number specific — it represents ordinary sequence alleles (insertions,
deletions, substitutions); copy-number `DUP` records keep their `REF_CN` / `CN` form. It fires only when
the locus is small enough (`--multiallelic-max-bp`, default 5000) and shows real variation (the largest
allele differs from `REF` by at least `--min-sv-bp`); otherwise the bubble falls back to per-event
records, so large SVs keep their typed representation.

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
  (carrier count) and **`AF`** (= AC/AN). `AF` is the **carrier (ALT) frequency** — the numerator counts
  carriers, the denominator counts every haplotype that traverses the bubble (so the reference and other
  reference-like haplotypes are in `AN` but not `AC`). It is **not folded**, so a near-fixed variant reads
  as `AF≈0.99`, not `0.01`; `--min-maf` drops records below this value.
  Also: `SVLEN_RANGE` (min,max member size when a merged record spans differing sizes) and the
  cross-haplotype merge evidence **`MERGE_JACCARD`** (strongest node-set overlap that joined a member),
  **`MERGE_SEQID`** (strongest sequence identity, when overlap did not decide it) and
  **`MERGE_SIZE_RATIO`** (smallest/largest merged member size) — all present only on records that actually
  merged ≥2 events; `EVENTID` (links a co-located DEL+INS substitution), `INS_SUBTYPE` (refined INS only),
  `REF_CN` and `RU_LEN` (DUP: reference copy number, and repeat-unit length in bp for one copy),
  `NALLELES` (multiallelic records from `--multiallelic-loci`), and the event sequence
  (`INSSEQ`/`DELSEQ`/`INVSEQ`, omitted when very long). `EVENT_NODES` is **deduplicated and ordered**
  (reference nodes by genomic position, haplotype-only nodes after them in walk order), so
  `START_NODE`→`END_NODE` reads as a coherent progression; the full per-carrier node walk lives in
  `variant_paths.tsv`. A substitution (same locus deletes a reference block and inserts a haplotype block)
  yields two records (DEL + INS) sharing one `EVENTID`.
- `CHROM`/`POS` come from the reference path's genomic label (e.g. `...chr6:31891045-...`) when
  parseable, else a graph-relative offset.

## Plotting

### Whole-VCF variant map (`scripts/plot_vcf_map.R`)

The headline figure is a single **oncoprint-style map of the entire VCF**: rows = haplotypes,
**columns = the called variants grouped by bubble** (one facet per bubble, each column labelled by its
variant ID), each cell colored by the **called event** the haplotype carries at that variant. It reads
like the VCF itself rather than like the graph, so it is the at-a-glance answer to "what did `call`
find, and who carries it?".

Cell colors: grey = reference-like (haplotype does not carry the variant); **DEL** red, **INS-NOVEL**
green, **INS-DUP** purple, **INV** orange, **multiallelic** teal (shaded by allele index). A **DUP** is
shaded **blue by the haplotype's absolute copy number** (`FORMAT:CN`) for **every** haplotype — loss
(light) / reference (mid) / gain (dark) — so the DUP column reads as a copy-number gradient. Rows are
sorted oncoprint-style so haplotypes with the same event pattern group together (event-free haplotypes
are dropped); `--reference-path NAME` pins that haplotype on top. It needs only the `call` VCF.

```bash
Rscript scripts/plot_vcf_map.R \
  --vcf results/real_data/c4/call/call.region.vcf \
  --reference-path grch38 \
  --title "c4 variant map" \
  --out results/real_data/c4/plots/c4_vcf_map
```

Optional row controls mirror the coverage heatmaps: `--clusters <...clusters.tsv>` keeps only the cluster
**representative** rows, `--cluster-by <...clusters.tsv>` groups/orders rows by walk cluster (a thin
separator between clusters), and `--max-paths N` caps the row count. Two layout flags: `--flip` transposes
the map (variants on Y, haplotypes on X, legend at the bottom); `--scale` draws each variant's rectangle
**proportional to its size** along the variant axis — `|SVLEN|` for DEL/INS/INV/multiallelic, and the
repeat-unit length `RU_LEN` (one copy) for a DUP, so a high-copy DUP is sized by its unit rather than
ballooning with copy number (`--scale-transform raw|sqrt|log1p`, default `sqrt`). `RU_LEN` is a new DUP
INFO field carrying the exact repeat-unit bp `call` computes (e.g. LPA KIV-2 = 5547 bp).

### Node-level SV map (`scripts/plot_sv_map.R`)

For the **graph-structure** view of a single bubble, `scripts/plot_sv_map.R` draws a **node-level
structural-variant map** in the same `odgi viz` style as the
`inspect` node/edge-coverage heatmaps: rows = haplotypes, **columns = the bubble's internal nodes in
reference-sorted order** (the same node set and order as the `inspect` `node_counts.tsv` columns), each
cell colored by the **called event** the haplotype carries on that node. Because each event is painted
only on its own nodes, an event occupies just its node columns (in genomic order) — a deletion reads as a
contiguous colored gap, an insertion as a block at its anchor — and the figure lines up 1:1 with the
coverage heatmap of the same bubble.

Cell colors: white = node not traversed; grey = traversed, reference-like; **DEL** red, **INS-NOVEL**
green, **INS-DUP** purple, **INV** orange, **multiallelic** teal (shaded by allele index). A **DUP** is
shaded **blue by the haplotype's absolute copy number** (the VCF `FORMAT:CN`), for **every** haplotype
that traverses the locus — so losses (`CN < REF_CN`, light), reference (`= REF_CN`, mid) and gains
(`> REF_CN`, dark) all read as a sequential shade — with DEL/INS painted on top. Two encodings are
handled automatically: a **self-loop / REP DUP** paints the CN on the repeat-unit node (a clean
copy-number strip — e.g. LPA KIV-2 at the `REP` node), while a **whole-module coverage DUP** (event node
is the bubble boundary, e.g. C4 / CYP2D6 called on the `bubble` graph) paints the CN across the whole
traversed module, so each haplotype's copy number reads as a horizontal band.

It reads the per-bubble `inspect` `node_counts.tsv` (substrate + node order) and the `call` VCF (events,
carriers, `GT:CN`). The bubble id is inferred from the `node_counts.tsv` filename (or pass `--bubble-id N`).
Optional: `--node-lengths <...node_lengths.tsv>` scales column widths by node bp; the row-clustering flags
mirror the coverage heatmaps — `--clusters <...clusters.tsv>` keeps only the cluster **representative**
rows, while `--cluster-by <...clusters.tsv>` keeps all rows but **groups/orders them by walk cluster**
(representative first, a thin separator between clusters, no color sidebar); `--reference-path NAME` pins
that haplotype (substring match) as the top row; `--max-nodes`/`--max-paths` cap the figure. Note the
`inspect` clusters group haplotypes by **walk similarity**, which is not the same as the called **copy
number** — the latter is shown directly by the DUP blue-intensity shading regardless of row order.

```bash
./build/panvar inspect -i results/real_data/lpa/panphorte/panphorte.normalized.sorted.gfa \
  --bubble-prefix-in results/real_data/lpa/panphorte/panphorte --bubble-id 7 \
  --cluster \
  --cluster-similarity 0.97 \
  -o results/real_data/lpa/inspect_panphorte/inspect
Rscript scripts/plot_sv_map.R \
  --node-counts results/real_data/lpa/inspect_panphorte/inspect.bubble_7.node_counts.tsv \
  --vcf results/real_data/lpa/call/call.region.vcf \
  --node-lengths results/real_data/lpa/inspect_panphorte/inspect.bubble_7.node_lengths.tsv \
  --clusters results/real_data/lpa/inspect_panphorte/inspect.bubble_7.clusters.tsv \
  --cluster-by results/real_data/lpa/inspect_panphorte/inspect.bubble_7.clusters.tsv \
  --reference-path grch38_1 \
  --out results/real_data/lpa/call/call.bubble_7
```

## Algorithm

Each step below is traced on a tiny worked dataset in
[algorithm_example.md](../algorithm_example.md). For each bubble:

1. **Alleles.** Take the reference walk and each haplotype's source→sink walk; group identical walks into
   distinct **alleles** and call once per allele (genotypes expand by membership). A *pure* single-node
   indel — an allele that goes straight from source to sink with no interior node — is handled explicitly,
   so a lone deletion or insertion at a two-anchor site is not mistaken for "does not traverse the bubble".
2. **Copy-number nodes.** Mark any node traversed more than once (in the reference or an allele) and fold
   its consecutive self-repeats into one alignment token, so a `REP` self-loop becomes a single anchor and
   surfaces as a copy-number change, not a spurious INS/DEL.
3. **Diff against the reference.** Split the reference and haplotype walks at shared **anchor** nodes
   (each appearing once in both walks, in a common order) and align only the segments between anchors.
   This bounds cost on large bubbles and gives the same breakpoints for every haplotype, so identical
   events merge cleanly. Each gap block becomes a **DEL** (reference-only nodes), **INS** (haplotype-only
   nodes), **INV** (a haplotype run that is the reverse-complement node-walk of the reference run), or a
   **substitution** (a co-located DEL + INS sharing one `EVENTID`).
4. **Coalesce within a haplotype.** Merge consecutive same-type events whose gap is within
   `--merge-distance-bp`, measured in **either** reference space **or** the haplotype's own sequence space
   — the latter catches events split on the reference by an intervening deletion but contiguous in the
   sample.
5. **Merge across haplotypes.** Cluster equivalent events by **transitive single-linkage** (connected
   components), keeping events down to `--rescue-min-bp`. Two events link when they share a type, fall in
   a **position window** (`--merge-distance-bp` widened by the smaller event's size, so a breakpoint that
   floats across haplotypes still merges), and either overlap in **node set** (length-weighted Jaccard
   `≥ --merge-jaccard`) **or** in **sequence** (identity `≥ --merge-seq-identity`, comparing only sizes
   whose shorter/longer ratio clears `--merge-size-ratio`). Connected components make the merge transitive
   and order-independent — A–B–C collapse into one record even when A and C never match directly. The
   largest member represents the record, and the evidence that joined the members is reported in
   `MERGE_JACCARD` / `MERGE_SEQID` / `MERGE_SIZE_RATIO`. Copy-number records merge separately, on shared
   `REP` node identity.
6. **Force-call, then filter.** Re-test every non-carrier haplotype at each record against its own diff
   and add it as a carrier when it matches — so a sub-threshold event in one haplotype is genotyped `1`
   when a larger matching event is called elsewhere, instead of `0`. One pass suffices (the representative
   is fixed). Finally keep records whose representative reaches `--min-sv-bp` and whose carrier count
   reaches `--min-haplotypes`.

## Example

```bash
./build/panvar call \
  -i results/real_data/lpa/panphorte/panphorte.normalized.sorted.gfa \
  --bubble-prefix-in results/real_data/lpa/panphorte/panphorte \
  --reference-path grch38_1 \
  -o results/real_data/lpa/call/call \
  --merge-distance-bp 100 --merge-jaccard 0.80 --classify-ins \
  --min-maf 0.05 \
  --cn-from-multiplicity --cn-from-coverage
```
