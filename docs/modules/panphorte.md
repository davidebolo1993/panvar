# Panphorte Module (Module 2 — bubble repeat normalization)

CLI entrypoint:

- `panvar panphorte`

## What it does

`panphorte` reworks tandem-repeat bubbles into a compact, copy-number-explicit form and writes a new GFA that can be re-processed by `bubble`/`inspect`/downstream. In a repeat expansion, a tandem array may show up as a long linear run of nodes whose copies have *identical/very similar sequences* (often as distinct node ids); this hides copy number behind raw length and bloats the bubble. The `panphorte` module detects each tandem array and **collapses it to a single repeat unit node with a self-loop**, so copy number becomes the number of self-loop traversals.

> **Scope.** `panphorte` is for **tandem-repeat** loci (a unit repeated in-line within a haplotype, e.g. LPA KIV-2). For **PGGB-collapsed paralog clusters** (e.g. C4/CYP2D6) the copies are already folded onto shared nodes as node *multiplicity* — they are not a contiguous tandem to normalize — so their copy number is recovered directly in `call --cn-from-coverage` on the `bubble` graph, **not** from `panphorte`. Running `panphorte` on such a region is harmless but its collapsed graph is not the call substrate there, and its `copies.tsv` is not the copy-number source.

### Differences from the original panphorte

This is panvar's own re-implementation of the [panphorte](https://github.com/GenoGra/Panphorte) idea. While the original
tool works in the full base sequence space, this re-implementation
differs deliberately:

| | Original panphorte | This re-implementation |
|---|---|---|
| Repeat detection | scans full base sequence (`O(bp^2)`) | scans the **node-walk**, comparing each step by the sequence it spells (`O(node_count^2)` — a ~1000x reduction on large bubbles) |
| Aligner | external | **none**; a banded aligner is built in for the approximate mode |
| Threading | — | per-haplotype parallel seed scan + copy detection (`--threads`) |
| Modes | single collapse | explicit **exact (sequence-preserving)** vs **approximate (lossy)** via `--min-similarity` |


In exact mode (default, `--min-similarity 1.0`) every rewrite is
**sequence-preserving** — each path spells exactly the same sequence through the normalized graph as through the original; the tool enforces this internally and aborts if any path's sequence would change. Approximate mode (`--min-similarity < 1.0`) is intentionally **lossy** (see below): copies are canonicalized to a representative unit, so haplotype sequences change within collapsed repeat regions and the invariant is not applied there. The scope of this module is to rewrite repeat stretches so that these can be called downstream as DUP. Arrays of identical copies **and** arrays with short interruptions between copies (e.g. `CGG (A) CGG CGG CGG`) are collapsed; the interrupting segment is kept as a nested step between repeat-unit traversals, so the result is still exactly sequence-preserving (`--max-interruption-frac` bounds how much interruption an array may contain).

## Required inputs

- `--gfa <graph.gfa>` / `-i <graph.gfa>` (W- or P-line GFA; the output preserves the input line type)
- one of:
  - `--bubble-prefix-in <module1-prefix>` (auto uses `<module1-prefix>.bubbles.csv`)
  - `--bubbles-csv <module1.bubbles.csv>`

## Key options

```bash
panvar panphorte -i <graph.gfa> (-b <prefix> | -c <bubbles.csv>) -o <prefix> [options]
```

- `-o, --out-prefix <prefix>`: writes the normalized GFA and `<prefix>.panphorte.report.tsv` (see Outputs
  for which GFA, depending on `--reference-path`)
- `-r, --reference-path <name>`: sort + re-snarl the normalized graph along this reference so the output
  is call-ready in one step (see Outputs)
- `--bubble-id <N>`: restrict to a bubble ID; repeatable
- `--min-unit-bp <N>`: minimum repeat-unit span to normalize (default `50`)
- `--min-copies <N>`: minimum tandem copies that must occur in **some** haplotype for a bubble to be treated as an array (default `2`). It defines the array, not the per-haplotype fold threshold: once a bubble is a confirmed array, **every** haplotype with `≥ 1` aligned copy of the unit is folded — so a single-copy haplotype gets a `copies = 1` row and traverses the `REP` node once.
- `--max-interruption-frac <f>`: max fraction of an array's bp that may be interruptions (default `0.25`)
- `--min-similarity <f>`: minimum identity to treat a block as a copy of the repeat unit (default `1.0` = exact, sequence-preserving). Values `< 1.0` enable **approximate, lossy collapse** of near-identical (divergent) copies.
- `--threads <N>`: worker threads for the approximate seed scan and copy detection (`0` = auto).
- `--gtf <path>`: reference-coordinate GTF; **after re-sorting**, project its genes onto the normalized graph's reference nodes and write `<prefix>.bandage_genes.csv` (needs a PanSN `--reference-path`). This is emitted separately from `bubble --gtf` because `panphorte`'s collapse renumbers nodes, so the gene positions differ on the normalized graph. lncRNAs are skipped.
- `--quiet`: disable the progress bar / logs


## Approximate collapse (`--min-similarity < 1.0`)

Exact collapse only merges byte-identical adjacent copies, so it misses **divergent** repeats — copies that differ by SNVs or small indels. With `--min-similarity <f>` `panphorte` does **single-block, lossy** collapse: one representative repeat unit per bubble, whose near-identical copies are found by **banded sequence alignment**. A small worked trace is in [algorithm_example.md](../algorithm_example.md).


1. **Seed one representative unit** per bubble from the exact tandem detector: the unit must come from a clean **adjacent identical pair** somewhere in the cohort, so it recovers the true repeat period. The most-supported unit across all paths wins.
2. **Find copies** per path by aligning the unit (and its reverse complement) into the path's spelled bubble sequence: short k-mer anchors propose copy starts, then a **banded global alignment** decides each copy, its extent, and its orientation. The band width is `(1 − f)·|unit|` (uncapped), so a copy that differs from the unit by a **large internal indel** still aligns when `f` is low enough. Copies need **not be adjacent** within a haplotype; sequence between accepted copies is kept as literal nodes, so only the copies themselves lose detail.
3. **Collapse** the copies to **one REP node** (the representative unit) traversed once per copy — a self-loop when copies are adjacent, edges through flanking literals when they are not. This is **lossy**: within-copy SNVs/small indels are discarded (collapsed copies become identical), appropriate when copy number is the variant of interest, not the within-copy differences. **Per-copy orientation is preserved** (three copies, two forward + one reverse → `REP +,−,+`). Folding applies to **single copies** of the seeded unit as well: a haplotype carrying one copy is routed through the `REP` node once (`copies = 1`), so it is not left un-normalized.

**Limitation:** the unit must occur as `≥ 2` **adjacent** identical copies in at least one haplotype to be seeded. Once seeded, copies in any haplotype may be non-adjacent/divergent. A duplication with **no** adjacent identical pair anywhere — e.g. paralogs separated by other sequence, or a copy embedded in a larger segmental-duplication cluster — is **not** seeded here, so `panphorte` leaves that bubble intact. Such folded duplications are recovered downstream by `call`, which reads copy number from the un-collapsed bubble's node-traversal multiplicity and emits a `DUP` record. See [call](call.md).


## Outputs

The rewritten graph collapses each array into a repeat-unit node with a self-loop (`L U + U + 0M`); the supporting paths route through it `copies` times, and duplicate-copy nodes are dropped when no longer referenced. **What is written depends on `--reference-path`:**

- **Without `--reference-path`** → `<prefix>.normalized.gfa` only. S-lines keep the input's node order minus the dropped copies, and the new REP nodes are **appended at the end**, not at their genomic locus.Downstream `call` treats **numeric node id == reference order**, so this graph must still be re-sorted before calling.

- **With `--reference-path <name>`** → the output is made **call-ready** in one step: the normalized graph is sorted + flipped along the reference (repositioning the appended REP nodes into reference order and renumbering ids) and re-snarled. Only the sorted form is written (the unsorted `normalized.gfa` is not):
  - `<prefix>.normalized.sorted.gfa` — the sorted, call-ready graph
  - `<prefix>.bubbles.csv` — the re-snarled bubbles
  - `<prefix>.bandage_nodes.csv` — node colors for Bandage inspection
  - with `--gtf <path>` (and a PanSN `--reference-path`): `<prefix>.bandage_genes.csv` (`Name,Colour,Gene`) projecting genes onto the **normalized** graph. `panphorte` renumbers nodes when it collapses, so this is distinct from `bubble`'s gene CSV.

- `<prefix>.panphorte.report.tsv`, one row per bubble:
  - `bubble_id`, `normalized` (yes/no), `unit_bp`, `paths_normalized`, `min_copies`, `max_copies`, `nodes_collapsed`
- `<prefix>.panphorte.copies.tsv` (approximate mode), one row per (haplotype, collapsed array) — the provenance for downstream calling: `path_name, sample, bubble_id, copies, unit_bp, orientations, mean_identity, region_bp, from_node, to_node`. Reads as "in this haplotype, original walk nodes `from_node..to_node` (`region_bp` bp) collapsed into the unit looped `copies` times (orientation pattern `orientations`, mean alignment identity `mean_identity`)".
- Run summary on stdout: bubbles seen/normalized, paths rewritten, nodes removed/added, edges added.

## Algorithm overview

For each bubble, for each path crossing it:

1. Take the path's source-to-sink node-walk and tokenize each step by a hash of the sequence it spells (orientation-aware), so identical-sequence copies match even when they are distinct node ids.
2. Find non-overlapping tandem arrays: establish the unit from a clean adjacent pair, then extend bidirectionally collecting further copies, **tolerating short interruptions** between them (kept as nested literal steps). A run is accepted when the unit spans `>= --min-unit-bp`, repeats `>= --min-copies` times, and interruptions stay within `--max-interruption-frac` (copies verified by actual sequence equality; node tokens are precomputed so detection avoids re-hashing sequences).
3. Collapse each array to a repeat-unit node (deduplicated across paths by canonical unit sequence) with a self-loop; reroute the path through the unit `copies` times and wire the flanking edges; drop now-unused duplicate nodes.

## Example

```bash
./build/panvar panphorte \
  -i results/real_data/lpa/bubble/bubble.sorted.gfa \
  --bubble-prefix-in results/real_data/lpa/bubble/bubble \
  --reference-path grch38#1 \
  -o results/real_data/lpa/panphorte/panphorte \
  --min-similarity 0.97
```
