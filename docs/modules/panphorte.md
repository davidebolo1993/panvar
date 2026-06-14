# Panphorte Module (bubble repeat normalization)

CLI entrypoint:

- `panvar panphorte`

## What it does

`panphorte` reworks tandem-repeat bubbles into a compact, copy-number-explicit form and writes a new
GFA that can be re-processed by `bubble`/`inspect`/downstream. In a repeat expansion, a tandem array
may show up as a long linear run of nodes whose copies have *identical/very similar sequences* (often as distinct
node ids); this hides copy number behind raw length and bloats the bubble. The `panphorte` module detects each
tandem array and **collapses it to a single repeat unit node with a self-loop**, so copy number
becomes the number of self-loop traversals.

### Differences from the original panphorte

This is panvar's own re-implementation of the [panphorte](https://github.com/GenoGra/Panphorte) idea. While the original
tool works in the full base sequence space, this re-implementation
differs deliberately:

| | Original panphorte | This re-implementation |
|---|---|---|
| Repeat detection | scans full base sequence (`O(bp^2)`) | scans the **node-walk**, comparing each step by the sequence it spells (`O(node_count^2)` — a ~1000x reduction on large bubbles) |
| Aligner | external | **none**; banded `fit_align` is built in for the approximate mode |
| Threading | — | per-haplotype parallel seed scan + copy detection (`--threads`) |
| Modes | single collapse | explicit **exact (sequence-preserving)** vs **approximate (lossy)** via `--min-similarity` |


**Correctness:** in exact mode (default, `--min-similarity 1.0`) every rewrite is
**sequence-preserving** — each path spells exactly the same sequence through the normalized graph as
through the original; the tool enforces this internally and aborts if any path's sequence would
change. Approximate mode (`--min-similarity < 1.0`) is intentionally **lossy** (see below): copies
are canonicalized to a representative unit, so haplotype sequences change within collapsed repeat
regions and the invariant is not applied there.

**Scope:** tandem repeats (DUP events) are rewritten — that is the class where linear
expansion obstructs calling. DEL / INS / INV are already encoded in the walks (a haplotype skipping
reference nodes is a deletion, extra nodes an insertion, reversed orientation an inversion) - and are
*read* by downstream variant calling. Arrays of identical copies **and** arrays with short
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
- `--threads <N>`: worker threads for the approximate seed scan and copy detection (`0` = auto).
- `--quiet`: disable the progress bar / logs


## Approximate collapse (`--min-similarity < 1.0`)

Exact collapse only merges byte-identical adjacent copies, which misses divergent repeats (e.g. LPA
KIV-2 copies differ by SNVs; C4 A/B modules differ by small indels and come in long/short forms). With
`--min-similarity <f>` panphorte does **single-block, lossy** collapse: one representative repeat
unit per bubble, whose near-identical copies are found by **banded sequence alignment**.

1. **Seed one representative unit** per bubble from the exact tandem detector: the unit must come from
   a clean **adjacent identical pair** somewhere in the cohort (so it recovers the true repeat period — the ~5.5 kb KIV-2 unit in LPA, for instance, or the whole ~32 kb in C4). The
   most-supported unit across all paths wins (for instnace the C4 **long** form is
   seeded and the short modules then align into it).
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



**Limitation:** the unit must occur as `>= 2` **adjacent** identical copies in at least one haplotype
to be seeded. Once seeded, copies in any haplotype may be non-adjacent / divergent. A duplication with
no adjacent identical pair anywhere (e.g. paralogs separated by other sequence, or a copy embedded in a
larger segmental-duplication cluster like GSTM1) is **not** seeded here, so the `panphorte` module leaves that bubble
intact. Such folded duplications are recovered downstream by `call`, which reads
copy number from the **peak node-traversal multiplicity** of the un-collapsed bubble relative to the
reference path and emits a `DUP` record. See [call](call.md).


## Outputs

- `<prefix>.normalized.gfa`: the rewritten graph. Each collapsed array becomes a repeat-unit node
  with a self-loop (`L U + U + 0M`), the supporting paths route through it `copies` times, and
  duplicate-copy nodes are dropped when no longer referenced by any path. Re-run `vg snarls` on this
  GFA before re-running `panvar bubble`; `inspect` only needs the new GFA + a bubbles CSV.

  **Node ordering:** S-lines are emitted in the input's node order minus the dropped duplicate-copy
  nodes, so every surviving original node keeps its input position — if the input was `odgi sort`ed
  along the reference, the normalized graph stays sorted. The new REP nodes are the **only** exception:
  they are **appended at the end**, not placed at their genomic locus.

  **Re-sort the normalized graph with `odgi sort` before re-snarling.** Downstream `call` writes
  `node_track.tsv` ordered by **numeric node id** and treats that as reference order ("the graph is
  `odgi sort`ed along the reference"); the appended REP nodes break that until the graph is re-sorted.
  So the intended pipeline is **panphorte → `odgi sort` → `vg snarls` → `bubble` → `call`**:

  ```bash
  odgi build -g out/panphorte.normalized.gfa -o - \
    | odgi sort -i - -o - -p Ygs -P \
    | odgi view -i - -g > out/panphorte.normalized.sorted.gfa
  ```

  Then run `vg snarls` + `panvar bubble` + `panvar call` on `...normalized.sorted.gfa`. (`odgi sort`
  repositions the REP nodes topologically along the reference and renumbers ids in sort order.)
- `<prefix>.panphorte.report.tsv`, one row per bubble:
  - `bubble_id`, `normalized` (yes/no), `unit_bp`, `paths_normalized`, `min_copies`, `max_copies`,
    `nodes_collapsed`
- `<prefix>.panphorte.copies.tsv` (approximate mode), one row per (haplotype, collapsed array) —
  the provenance for downstream calling: `path_name, sample, bubble_id, copies, unit_bp, orientations,
  mean_identity, region_bp, from_node, to_node, n_long, n_short`. Reads as "in this haplotype, original
  walk nodes `from_node..to_node` (`region_bp` bp) collapsed into the unit looped `copies` times
  (orientation pattern `orientations`, mean alignment identity `mean_identity`)". `n_long`/`n_short`
  split those `copies` by per-copy length: a copy spanning `< 0.90 ×` the (long) unit is `short`: `n_long + n_short = copies`.
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
# normalize tandem-repeat bubbles
panvar panphorte \
  -i tests/real_data/c4.gfa \
  --bubble-prefix-in tests/results/c4/bubble \
  -o tests/results/c4/panphorte

# re-sort the normalized graph so REP nodes are repositioned along the reference
odgi build -g tests/results/c4/panphorte.normalized.gfa -o - \
  | odgi sort -i - -o - -p Ygs -P \
  | odgi view -i - -g > tests/results/c4/panphorte.normalized.sorted.gfa

# re-snarl (see `bubble` docs) then re-bubble on the SORTED graph
panvar bubble \
  -i tests/results/c4/panphorte.normalized.sorted.gfa \
  --snarls-in tests/results/c4/panphorte.normalized.sorted.snarls.jsonl \
  -o tests/results/c4/bubble_normalized
```
