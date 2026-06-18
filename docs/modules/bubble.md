# Bubble Module (Module 1)

Date: 2026-04-20

CLI entrypoint:

- `panvar bubble`

## What it does

`bubble` turns a pangenome graph (GFA) into `panvar` bubble sites for downstream modules. By default it:

1. sorts + flips the graph along the reference internally
2. finds **snarls** internally with a vendored cactus / 3-edge-connected decomposition (matches closely
   `vg snarls`; see [Backround](#snarls-vs-superbubbles-background))
3. infers bubble-internal nodes from path intervals between each snarl's source/sink
4. computes path support and internal sequence span per candidate
5. applies base site filters (`--min-variant-bp`, `--min-path-support`)
6. optionally merges nearby bubbles (`--merge-nearby-bp`)
7. writes the **sorted GFA** + **CSV of bubbles** + **visualization outputs**

## Required inputs

1. `--gfa <graph.gfa>` (any GFA)
2. `--reference-path <name>` — reference path name or unique case-insensitive substring, used to order the internal sort/flip and snarl finder (not needed if you pass `--snarls-in`).

## Outputs

- `*.bubbles.csv`: refined bubble/site table used by downstream modules (`inspect`, `panphorte`)
- `*.sorted.gfa`: the reference-sorted/-flipped graph to be used as input for downstream `panphorte`
- `*.bandage_nodes.csv`: node colors for Bandage visualization
- optional `--snarl-debug-tsv <path>`: candidate-level diagnostics
- optional `--emit-snarls-jsonl <path>`: the internal snarls in vg-style JSONL

Output directories are auto-created when missing.

## Input graph (from pggb)

A graph is typically produced with `pggb`:

```bash
pggb -i <haplotypes.fa> -o <pggb.outdir>
panvar bubble -i <pggb.outdir>/*smooth.final.gfa --reference-path <reference.id> -o <panvar.outdir>/bubble
```

When run with `--snarls-in <snarls.jsonl>` (from `vg snarls -A integrated <graph.gfa> | vg view -R -j -`), the graph is used **as-is**. 

## Background: Snarls vs. superbubbles

`bubble` consumes **snarls** by default, not superbubbles.

A **superbubble** is a single-source / single-sink subgraph that is **directed and acyclic**: every
internal node is reachable from the source and reaches the sink, and the only way in or out is through
the two boundaries. By construction it **cannot** contain a cycle or an inversion. A **snarl** (*i.e.* what `vg snarls` emits) is the more general structure: a pair of boundary nodes whose removal
separates an internal subgraph from the rest, defined on the **bidirected** graph. Snarls therefore
capture the **inversions, tandem cycles, and tangles** that a superbubble cannot, and they
**nest**; their acyclic subclass — "ultrabubbles" — is essentially the superbubble. `panvar` uses snarls
because the pangenome variation this toolkit targets (including inversions and tandem expansions) often
lives in the cyclic or inverted sites that superbubbles omit.

Internally, `panvar` reproduces `vg snarls` by mirroring vg's cactus / 3-edge-connected decomposition, so
the default top-level snarl set matches vg; passing `--superbubbles` instead emits only the acyclic
superbubble subset. This acyclic subset is the same superbubble concept used by tools like 
[BubbleGun](https://doi.org/10.1093/bioinformatics/btac448), whose detector finds superbubbles
directly, via the Onodera algorithm (BubbleGun additionally *sub-labels* its
superbubbles into "simple bubbles", "insertions", and "super" — so its narrow "super" excludes
simple/insertion sites, whereas panvar's `--superbubbles` keeps **all** of them).

A snarl is a property of the graph **topology**, so it is computed independently of the paths: it is
bounded by two nodes (`source`, `sink`), and a haplotype only "supports" it when its walk visits **both**
boundaries. In a pangenome some haplotypes may not, because:

- a haplotype may take a route that **bypasses a boundary** (for example a large deletion that removes
  the boundary region, or an alternative local structure);
- an assembly may be **fragmented or partial** and simply not span the locus;
- with nested snarls, a path can cross the parent yet take a branch that **skips a child** snarl.


**References.**

- Onodera, Sadakane, Shibuya. *Detecting Superbubbles in Assembly Graphs.* WABI 2013.
  (the superbubble concept)
- Paten, Eizenga, Rosen, Novak, Garrison, Hickey. *Superbubbles, Ultrabubbles, and Cacti.*
  J. Comput. Biol. 25(7):649–663, 2018. <https://doi.org/10.1089/cmb.2017.0251>
  (snarls / ultrabubbles / the cactus decomposition behind `vg snarls` and panvar's finder)
- Dabbaghie, Ebler, Marschall. *BubbleGun: enumerating bubbles and superbubbles in genome graphs.*
  Bioinformatics 38(17):4217–4219, 2022. <https://doi.org/10.1093/bioinformatics/btac448>


## Algorithm overview

For each top-level snarl candidate:

1. get `(source, sink)` boundary from the internal cactus finder (or `--snarls-in` JSONL)
2. find path intervals crossing source->sink (or sink->source, then canonicalize)
3. collect internal nodes seen between boundaries
4. measure:
   - number of crossing paths
   - internal bp support across crossing paths
5. apply base filters (`--min-path-support`, `--min-variant-bp`)
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


## Key options

```text
panvar bubble -i <graph.gfa> -r <name> [-o <prefix>] [--superbubbles] [options]
```

- `-i, --gfa <path>` — input GFA (required)
- `-r, --reference-path <name>` — reference for the internal sort/flip + snarl finder (required unless
  `--snarls-in`)
- `-o, --out-prefix <prefix>` — output prefix (default `bubble_calls`)
- `-s, --superbubbles` — emit only acyclic superbubbles instead of all snarls
- `--no-flip` — skip reorienting nodes to the reference forward strand (still sorts)
- `--sorted-gfa-out <path>` — where to write the sorted GFA (default `<prefix>.sorted.gfa`)
- `--emit-snarls-jsonl <path>` — also write the internal snarls as vg-style JSONL
- `--snarls-in <path>` — legacy override: use an external `vg snarls` JSONL (graph used as-is, no sort)
- `--min-variant-bp <N>`
- `--min-path-support <N>`
- `--merge-nearby-bp <N>`
- `--snarl-debug-tsv <path>`
- `-q, --quiet` — disable the progress bar

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
- `inversion_signal`: true when at least one internal node is observed in both orientations across supporting paths (such a bubble is kept even when no path reaches `--min-variant-bp`)

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
./build/panvar bubble \
  -i tests/real_data/lpa.gfa \
  -o tests/results/lpa/bubble/bubble \
  --reference-path grch38_1

#or override with an external vg snarls file
./build/panvar bubble -i tests/real_data/lpa.gfa -o tests/results/lpa/bubble/bubble --snarls-in tests/real_data/c4.snarls.jsonl
```

## Bandage Colors

Bandage node colors provide context and retained calls:

- blue: nodes in non-SNP candidate bubbles (pre-filter context)
- red: nodes in retained output bubbles (`*.bubbles.csv`)

