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
4. apply site filters (`min-variant-bp`, `min-path-support`, nesting cap)
5. write module handoff and visualization outputs

It does not execute `vg` at runtime.

## Required inputs

1. `--gfa <graph.gfa>`
2. `--snarls-in <snarls.jsonl>`

## Outputs

- `*.bubbles.csv`: refined bubble/site table used by module 2/3
- `*.bandage_nodes.csv`: node colors for Bandage
- optional `--snarl-debug-tsv <path>`: candidate-level diagnostics

## Preparing snarl JSONL (offline)

`snarls.jsonl` is generated before running `panvar bubble`:

```bash
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
   - nesting metadata from JSONL hierarchy
5. apply filters
6. keep accepted sites and assign bubble IDs / nesting levels

## Filter behavior

- `--min-variant-bp <N>`
  - keep only candidates where at least one supporting path has internal sequence span `>= N`
  - default `50`
  - `0` disables this filter

- `--min-path-support <N>`
  - keep only candidates with at least `N` crossing paths
  - default `0` (disabled)

- `--max-nesting-level <N>`
  - keep only candidates up to nesting level `N`
  - `0` means no cap

## Key options

- `--min-variant-bp <N>`
- `--min-path-support <N>`
- `--max-nesting-level <N>`
- `--snarl-debug-tsv <path>`

## Debug TSV columns

When `--snarl-debug-tsv` is enabled, each candidate row includes:

- `candidate_id`: candidate index in scan order
- `source`, `sink`: boundary nodes
- `inside_node_count`: number of inferred internal nodes
- `n_paths`: number of crossing paths
- `min_inside_bp`: smallest internal bp span among crossing paths
- `nested`: `1` if this top-level snarl has child snarls in JSONL, else `0`
- `accepted`: `1` if kept after filters, else `0`
- `reason`: acceptance/rejection reason (`accept:final`, `reject:min-variant-bp`, etc.)

## Example

```bash
./build/panvar bubble \
  -i tests/real_data/c4.gfa \
  -o tests/results/c4/bubble \
  --snarls-in tests/real_data/c4.snarls.jsonl
```

For clustering behavior details and usage modes, see:

- `docs/modules/allele.md`

## Module handoff

Use bubble output as input to module 2/3:

- `panvar allele --bubbles-csv-in <prefix>.bubbles.csv`
- `panvar call --bubbles-csv-in <prefix>.bubbles.csv`
