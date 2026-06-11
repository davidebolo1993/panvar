# Bubble Module (Module 1)

Date: 2026-04-20

CLI entrypoint:

- `panvar bubble`

## What it does

`bubble` refines precomputed snarl boundaries into `panvar` bubble sites for downstream modules.

In practice:

1. load top-level snarl boundaries from `--snarls-in` JSONL
2. infer bubble-internal nodes from path intervals between source/sink
3. compute path support and internal sequence span per candidate
4. apply base site filters (`min-variant-bp`, `min-path-support`)
5. optionally merge nearby bubbles (`--merge-nearby-bp`)
6. write module handoff and visualization outputs

## Required inputs

1. `--gfa <graph.gfa>`
2. `--snarls-in <snarls.jsonl>`

## Outputs

- `*.bubbles.csv`: refined bubble/site table used by module 2/3
- `*.bandage_nodes.csv`: node colors for Bandage
- optional `--snarl-debug-tsv <path>`: candidate-level diagnostics

Output directories are auto-created when missing.

## Preparing snarl JSONL (offline)

`snarls.jsonl` is generated before running `panvar bubble`:

```bash
#graph building
pggb -i <haplotypes.fa> -o <pggb.outdir>
#manipulating
odgi paths -i <pggb.outdir>/*smooth.final.og -L | grep <reference.id> > <pggb.outdir>/ref.path.txt
odgi sort -i <pggb.outdir>/*smooth.final.og -Y -H <pggb.outdir>/ref.path.txt -o - | odgi flip -i - --ref-flips <pggb.outdir>/ref.path.txt -o - | odgi view -i - -g | sed 's/_inv$//g'>  <graph.gfa>
#snarls calling
vg snarls -A integrated <graph.gfa>  | vg view -R -j - >  <graph>.snarls.jsonl
```

## Algorithm overview

For each top-level snarl candidate:

1. read `(source, sink)` boundary from JSONL
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

Compatibility:

- readers still accept legacy columns (`site_mode`, `type`, `nesting_level`, `parent_id`, `long_path_support`, `inversion_signal`) when present; fields not used by current logic are ignored.

Internal metrics (computed, not emitted in simplified CSV):

- `long_path_support`: count of supporting paths with inside span `>= --min-variant-bp`
- `inversion_signal`: true when at least one internal node is observed in both orientations across supporting paths

## Debug TSV columns

When `--snarl-debug-tsv` is enabled, each candidate row includes:

- `candidate_id`: candidate index in scan order
- `source`, `sink`: boundary nodes
- `inside_node_count`: number of inferred internal nodes
- `n_paths`: number of crossing paths
- `min_inside_bp`: smallest internal bp span among crossing paths
- `accepted`: `1` if kept after filters, else `0`
- `reason`: acceptance/rejection reason (`accept:final`, `reject:min-variant-bp`, etc.)

Common `reason` values:

- `accept:snarl-jsonl`: valid parsed snarl candidate before final filters
- `accept:final`: candidate survives all enabled filters and merge handling
- `accept:final-merged`: synthetic accepted row for a final merged bubble endpoint that did not exist as an original single snarl candidate
- `reject:min-path-support`
- `reject:min-variant-bp`
- `reject:merged-nearby`: candidate survived base filters but was fused into a nearby merged bubble
- `reject:filtered-out`: fallback when removed but no specific reason was recorded
- parse/import rejects such as `reject:parse-start-end`, `reject:endpoint-not-in-graph`, `reject:empty-inside`, `reject:duplicate-endpoints-smaller-inside`

## CLI summary terms

- `simple|super|insertion`: topology labels computed on retained bubbles for summary only

## Example

```bash
./build/panvar bubble \
  -i tests/real_data/c4.gfa \
  -o tests/results/c4/bubble \
  --snarls-in tests/real_data/c4.snarls.jsonl \
  --merge-nearby-bp 20
```

## Bandage Colors

Bandage node colors provide context and retained calls:

- blue: nodes in non-SNP candidate bubbles (pre-filter context)
- red: nodes in retained output bubbles (`*.bubbles.csv`)

For clustering behavior details and usage modes, see:

- `docs/modules/allele.md`

## Module handoff

Use bubble output as input to module 2/3:

- `panvar allele --bubble-prefix-in <prefix>`
- `panvar call --bubble-prefix-in <prefix>`

Equivalent explicit-file form (module 3 only):

- `panvar call --bubbles-csv-in <prefix>.bubbles.csv`
