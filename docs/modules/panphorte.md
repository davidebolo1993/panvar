# Panphorte Module (bubble repeat normalization)

CLI entrypoint:

- `panvar panphorte`

## What it does

`panphorte` reworks tandem-repeat bubbles into a compact, copy-number-explicit form and writes a new
GFA that can be re-processed by `bubble`/`inspect`/downstream. In a repeat expansion, a tandem array
shows up as a long linear run of nodes whose copies have *identical sequences* (often as distinct
node ids); this hides copy number behind raw length and bloats the bubble. `panphorte` detects each
tandem array and **collapses it to a single repeat unit node with a self-loop**, so copy number
becomes the number of self-loop traversals.

It is panvar's own efficient, dependency-free take on the "panphorte" idea: it reuses panvar's own
bubbles (no BubbleGun), detects repeats on the **node-walk** (comparing each step by the sequence it
spells) rather than scanning the full base sequence (`O(node_count^2)` instead of `O(bp^2)` — a
~1000x reduction on large bubbles), and uses **no external aligner**.

**Correctness:** in exact mode (default, `--min-similarity 1.0`) every rewrite is
**sequence-preserving** — each path spells exactly the same sequence through the normalized graph as
through the original; the tool enforces this internally and aborts if any path's sequence would
change. Approximate mode (`--min-similarity < 1.0`) is intentionally **lossy** (see below): copies
are canonicalized to a representative unit, so haplotype sequences change within collapsed repeat
regions and the invariant is not applied there.

**Scope:** only tandem repeats (the DUP class) are rewritten — that is the class where linear
expansion obstructs calling. DEL / INS / INV are already encoded in the walks (a haplotype skipping
reference nodes is a deletion, extra nodes an insertion, reversed orientation an inversion) and are
*read* by variant calling, not rewritten here. Arrays of identical copies **and** arrays with short
interruptions between copies (e.g. `CGG (A) CGG CGG CGG`) are collapsed; the interrupting segment is
kept as a nested step between repeat-unit traversals, so the result is still exactly
sequence-preserving (`--max-interruption-frac` bounds how much interruption an array may contain).

## Required inputs

- `--gfa <graph.gfa>` / `-i <graph.gfa>` (W- or P-line GFA; the output preserves the input line type)
- one of:
  - `--bubble-prefix-in <module1-prefix>` (auto uses `<module1-prefix>.bubbles.csv`)
  - `--bubbles-csv <module1.bubbles.csv>`

## Key options

- `-o, --out-prefix <prefix>`: writes `<prefix>.normalized.gfa` and `<prefix>.panphorte.report.tsv`
- `--bubble-id <N>`: restrict to a bubble ID; repeatable
- `--min-unit-bp <N>`: minimum repeat-unit span to normalize (default `50`)
- `--min-copies <N>`: minimum tandem copies to normalize (default `2`)
- `--max-interruption-frac <f>`: max fraction of an array's bp that may be interruptions (default `0.25`)
- `--min-similarity <f>`: minimum identity to treat a block as a copy of the repeat unit (default `1.0`
  = exact, sequence-preserving). Values `< 1.0` enable **approximate, lossy collapse** of
  near-identical (divergent) copies.
- `--quiet`: disable progress logs

## Approximate collapse (`--min-similarity < 1.0`)

Exact collapse only merges byte-identical adjacent copies, which misses divergent repeats (e.g. LPA
KIV-2 copies differ by SNVs; C4 modules differ by small indels and come in long/short forms). With
`--min-similarity <f>` panphorte does **single-block, lossy** collapse: one representative repeat
unit per bubble, whose near-identical copies are found by **banded sequence alignment**.

1. **Seed one representative unit** per bubble from the exact tandem detector: the unit must come from
   a clean **adjacent identical pair** somewhere in the cohort (so we recover the true repeat period —
   the whole ~32 kb C4 module, the ~5.5 kb KIV-2 unit — not a small segment that merely recurs). The
   most-supported unit across all paths wins.
2. **Find copies** per path by aligning the unit (and its reverse complement) into the path's spelled
   bubble sequence: multi-seed 16-mer anchors propose copy starts, then a **banded global alignment**
   (`fit_align`) decides each copy, its extent, and its orientation. The band is `(1 − f)·|unit|`
   (uncapped), so a copy that differs from the unit by a **large internal indel** still aligns when
   `f` is low enough — this is what lets a C4 **short** module (missing the ~6.4 kb HERV-K, ≈ 20 % of
   the long unit) align to the long unit at `--min-similarity 0.70`. Copies need **not be adjacent**
   within a haplotype; sequence between accepted copies is kept as literal nodes (so only the copies
   themselves lose detail).
3. **Collapse** the copies to **one REP node** (the representative unit) traversed once per copy —
   a self-loop when copies are adjacent, edges through flanking literals when they are not. This is
   **lossy**: within-copy SNVs/small indels are discarded (collapsed copies become identical),
   appropriate when copy number is the variant of interest, not the within-copy differences.
   **Per-copy orientation is preserved** (three copies, two forward + one reverse → `REP +,−,+`).

Per-path detection runs on `--threads` workers (default = hardware concurrency); haplotypes are
processed independently, and results are deterministic regardless of thread count (each path writes
its own slot; the graph mutation is serial and in path order). Exact mode (`--min-similarity 1.0`)
is unchanged and sequence-preserving.

On the bundled data (threaded, ~1 min/locus):

- **LPA KIV-2** collapses to one **5547 bp** unit at **mean 19.3 / max 32** occurrences across all
  465 haplotypes — matching `LPA_repeats.tsv`.
- **C4** collapses to the **32738 bp (long) RCCX** unit. At `--min-similarity 0.90` only the long
  modules align (43 haps × 2–3). Lowering to **0.70** pulls the **short** modules into the same REP
  (they align across the HERV-K indel), giving **114 haps** (98 × 2, 16 × 3) — the 98 × 2 tracks the
  `chr6_…plot.bed.gz` count of ~98 haplotypes with two C4A+C4B modules. Note the long/short (and
  C4A/C4B) distinction is **not preserved**: every copy is canonicalized to the long-unit sequence.
  Separating long vs short per copy would require splitting the unit into a core + HERV node (light
  node-splitting from the long-vs-short alignment gap), which is not done here.

**Limitation:** the unit must occur as `>= 2` **adjacent** identical copies in at least one haplotype
to be seeded. Once seeded, copies in any haplotype may be non-adjacent / divergent.

## Outputs

- `<prefix>.normalized.gfa`: the rewritten graph. Each collapsed array becomes a repeat-unit node
  with a self-loop (`L U + U + 0M`), the supporting paths route through it `copies` times, and
  duplicate-copy nodes are dropped when no longer referenced by any path. Re-run `vg snarls` on this
  GFA before re-running `panvar bubble`; `inspect` only needs the new GFA + a bubbles CSV.
- `<prefix>.panphorte.report.tsv`, one row per bubble:
  - `bubble_id`, `normalized` (yes/no), `unit_bp`, `paths_normalized`, `min_copies`, `max_copies`,
    `nodes_collapsed`
- `<prefix>.panphorte.copies.tsv` (approximate mode), one row per (haplotype, collapsed array) —
  the provenance for downstream calling: `path_name, sample, bubble_id, copies, unit_bp, orientations,
  mean_identity, region_bp, from_node, to_node`. Reads as "in this haplotype, original walk nodes
  `from_node..to_node` (`region_bp` bp) collapsed into the unit looped `copies` times (orientation
  pattern `orientations`, mean alignment identity `mean_identity`)".
- Run summary on stdout: bubbles seen / normalized, paths rewritten, nodes removed / added, edges
  added.

## Algorithm overview

For each bubble, for each path crossing it:

1. Take the path's source-to-sink node-walk (reusing the bubble-path machinery) and tokenize each
   step by a hash of the sequence it spells (orientation-aware), so identical-sequence copies match
   even when they are distinct node ids.
2. Find non-overlapping tandem arrays: establish the unit from a clean adjacent pair, then extend
   bidirectionally collecting further copies, **tolerating short interruptions** between them (kept
   as nested literal steps). A run is accepted when the unit spans `>= --min-unit-bp`, repeats
   `>= --min-copies` times, and interruptions stay within `--max-interruption-frac` (copies verified
   by actual sequence equality; node tokens are precomputed so detection avoids re-hashing sequences).
3. Collapse each array to a repeat-unit node (deduplicated across paths by canonical unit sequence)
   with a self-loop; reroute the path through the unit `copies` times and wire the flanking edges;
   drop now-unused duplicate nodes.

## Example

```bash
# 1) bubbles on the (orientation-normalized) graph
panvar bubble -i graph.flp.gfa -o out/bubble --snarls-in graph.flp.snarls.jsonl

# 2) normalize tandem-repeat bubbles
panvar panphorte -i graph.flp.gfa --bubble-prefix-in out/bubble -o out/panphorte

# 3) re-snarl + re-bubble the normalized graph, then continue downstream
vg snarls ... > out/panphorte.normalized.snarls.jsonl   # recompute snarls
panvar bubble -i out/panphorte.normalized.gfa -o out/bubble2 --snarls-in out/panphorte.normalized.snarls.jsonl
```
