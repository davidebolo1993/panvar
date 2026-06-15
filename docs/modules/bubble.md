# Bubble Module (Module 1)

Date: 2026-04-20

CLI entrypoint:

- `panvar bubble`

## What it does

`bubble` turns any GFA into `panvar` bubble sites for downstream modules, **with no external tools**
(no `odgi`, no `vg`). By default it:

1. sorts + flips the graph along the reference internally (restores `numeric node id == reference order`)
2. finds **snarls** internally with a vendored cactus / 3-edge-connected decomposition (matches
   `vg snarls`; see [snarls vs superbubbles](#snarls-vs-superbubbles-background))
3. infers bubble-internal nodes from path intervals between each snarl's source/sink
4. computes path support and internal sequence span per candidate
5. applies base site filters (`min-variant-bp`, `min-path-support`)
6. optionally merges nearby bubbles (`--merge-nearby-bp`)
7. writes the **sorted GFA** + module-handoff + visualization outputs

## Required inputs

1. `--gfa <graph.gfa>` (any GFA; it is sorted internally)
2. `--reference-path <name>` — reference path name or unique case-insensitive substring, used to order
   the internal sort/flip and snarl finder. (Not needed if you pass `--snarls-in`, the legacy override.)

## Outputs

- `*.bubbles.csv`: refined bubble/site table used by module 2/3
- `*.sorted.gfa`: the internally sorted (+flipped) graph — **use this for downstream `panphorte`/`call`**
- `*.bandage_nodes.csv`: node colors for Bandage
- optional `--snarl-debug-tsv <path>`: candidate-level diagnostics
- optional `--emit-snarls-jsonl <path>`: the internal snarls in vg-style JSONL (for inspection)

Output directories are auto-created when missing.

## Input graph (from pggb)

A graph is typically produced with `pggb`; no manual `odgi sort`/`odgi flip`/`vg snarls` step is needed
anymore — `bubble` does the equivalent internally:

```bash
pggb -i <haplotypes.fa> -o <pggb.outdir>
odgi view -i <pggb.outdir>/*smooth.final.og -g > <graph.gfa>   # just GFA conversion
panvar bubble -i <graph.gfa> --reference-path <reference.id> -o out/bubble
```

**Legacy override.** To reproduce an exact external `vg snarls` run, pass `--snarls-in <snarls.jsonl>`
(from `vg snarls -A integrated <graph.gfa> | vg view -R -j -`). In that mode the graph is used **as-is**
(no internal sort), so the JSONL node ids match. The two modes are mutually exclusive — never combine an
external snarls file with internal sorting, since sorting renumbers nodes.

No manual re-sort is needed after `panphorte` either: `panphorte --reference-path …` re-sorts and
re-snarls its normalized output itself (see [panphorte](panphorte.md)).

## Snarls vs. superbubbles (background)

`bubble` consumes **snarls**, not superbubbles. The difference matters for what variation can be
represented:

- A **superbubble** (what tools like BubbleGun enumerate) is a single-source / single-sink subgraph that
  is **directed and acyclic**: every internal node is reachable from the source and reaches the sink, and
  the only way in or out is through the two boundaries. By construction a superbubble **cannot** contain a
  cycle or an inversion.
- A **snarl** (Paten et al.; what `vg snarls` emits) is the more general structure: a pair of boundary
  nodes whose removal separates an internal subgraph from the rest, defined on the **bidirected** graph.
  Snarls therefore also capture **inversions, tandem cycles, and tangles** that a superbubble can't, and
  they **nest** (the snarl tree; the acyclic subclass — "ultrabubbles" — is essentially the superbubble).
  `-A integrated` reports the nested decomposition; `bubble` takes the **top-level** snarls.

panvar uses snarls because real pangenome variation (the inversions and tandem expansions this toolkit
targets) lives precisely in the cyclic/inverted sites superbubbles omit.

**Internal finder.** panvar reproduces `vg snarls` internally by vendoring vg's cactus / 3-edge-connected
decomposition (`integrated_snarls.cpp`), so the default top-level snarl set matches vg (validated
bubble-for-bubble on the bundled c4/lpa/gstm1 graphs). Pass `--superbubbles` to instead emit only the
acyclic superbubble subset — useful to *see* the difference: on a locus with a tandem cycle or inversion,
the default snarl mode reports the site while `--superbubbles` omits it.

### Why a snarl need not be crossed by every path

A snarl is a property of the graph **topology**, computed independently of the paths. It is bounded by two
specific nodes (`source`, `sink`), and a haplotype only "supports" that snarl if its walk visits **both**
boundaries. In a pangenome many haplotypes don't:

- a haplotype may take a different route that **bypasses a boundary** (e.g. a large deletion that removes
  the boundary region, or an alternative local structure);
- an assembly may be **fragmented/partial** and simply not span the locus;
- with nested snarls, a path can cross the parent yet take a branch that **skips a child** snarl.

So `path_support` is the count of paths whose walk actually crosses `source → sink` (found by
`find_best_bubble_path_interval`), and it is normally **less than the total number of P/W paths** — that
is expected, not an error. Downstream, `call` marks a sample that doesn't cross a bubble with genotype
`.` (vs `0` for "crosses but reference-like" and `1` for a carrier).

## Algorithm overview

For each top-level snarl candidate:

1. get `(source, sink)` boundary from the internal cactus finder (or `--snarls-in` JSONL)
2. find path intervals crossing source->sink (or sink->source, then canonicalize)
3. collect internal nodes seen between boundaries
4. measure:
   - number of crossing paths
   - internal bp support across crossing paths
5. apply base filters (`min-path-support`, `min-variant-bp`)
6. optionally merge nearby surviving bubbles
7. assign bubble IDs

## Filter behavior

- `--min-variant-bp <N>`
  - keep only candidates where at least one supporting path has internal sequence span `>= N`
  - default `50`
  - `0` disables this filter

- `--min-path-support <N>`
  - keep only candidates with at least `N` crossing paths
  - default `0` (disabled)

## Optional Nearby Merge

- `--merge-nearby-bp <N>` (default `0`, disabled)
- merge is applied after `--min-variant-bp` and `--min-path-support` filters
- if enabled, consecutive surviving bubbles are merged when the shortest-path distance
  from previous `sink` to next `source` is `<= N` bp
- distance is computed from node lengths across graph connectivity
Example:

If bubble A ends at node `2527` and bubble B starts at node `2527`, and node `2527` has length `10 bp`,
then with `--merge-nearby-bp 20` those bubbles are merged.

## Key options

- `--reference-path <name>` — reference for the internal sort/flip + snarl finder (required unless
  `--snarls-in`)
- `--superbubbles` — emit only acyclic superbubbles instead of all snarls
- `--no-flip` — skip reorienting nodes to the reference forward strand (still sorts)
- `--sorted-gfa-out <path>` — where to write the sorted GFA (default `<prefix>.sorted.gfa`)
- `--emit-snarls-jsonl <path>` — also write the internal snarls as vg-style JSONL
- `--snarls-in <path>` — legacy override: use an external `vg snarls` JSONL (graph used as-is, no sort)
- `--min-variant-bp <N>`
- `--min-path-support <N>`
- `--merge-nearby-bp <N>`
- `--snarl-debug-tsv <path>`

## Bubble CSV Columns

Current schema:

- `bubble_id`
- `source`
- `sink`
- `inside_node_count`
- `total_node_count`
- `path_support`
- `min_inside_bp`
- `max_inside_bp`
- `inside_nodes`


Two extra per-bubble metrics are also computed — they drive the
`--min-variant-bp` filter and are surfaced in the optional debug TSV:

- `long_path_support`: count of supporting paths with inside span `>= --min-variant-bp`
- `inversion_signal`: true when at least one internal node is observed in both orientations across
  supporting paths (such a bubble is kept even when no path reaches `--min-variant-bp`)

## Debug TSV columns

When `--snarl-debug-tsv` is enabled, each snarl candidate becomes one row:

- `candidate_id`: candidate index in scan order
- `source`, `sink`: boundary nodes
- `inside_node_count`: number of inferred internal nodes
- `n_paths`: number of crossing paths
- `min_inside_bp`: smallest internal bp span among crossing paths
- `long_path_support`, `inversion_signal`: the two metrics above
- `accepted`: `1` if the candidate is kept after all filters (and merging), else `0`

Final bubbles produced by `--merge-nearby-bp` have no original candidate row, so one extra
`accepted=1` row is appended for each.

## Example

```bash
# Internal pipeline (no vg/odgi): sort + flip + cactus snarls.
./build/panvar bubble \
  -i tests/real_data/c4.gfa \
  -o tests/results/c4/bubble \
  --reference-path grch38 \
  --merge-nearby-bp 20
# downstream uses tests/results/c4/bubble.sorted.gfa

# Legacy override with an external vg snarls file (graph used as-is, not re-sorted):
./build/panvar bubble -i tests/real_data/c4.gfa -o tests/results/c4/bubble \
  --snarls-in tests/real_data/c4.snarls.jsonl
```

## Bandage Colors

Bandage node colors provide context and retained calls:

- blue: nodes in non-SNP candidate bubbles (pre-filter context)
- red: nodes in retained output bubbles (`*.bubbles.csv`)

