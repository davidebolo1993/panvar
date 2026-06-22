# Bubble Module (Module 1)

CLI: `panvar bubble`

## What it does

Turns a pangenome graph (GFA, typically from `pggb`) into `panvar` **bubble sites** for the downstream
modules. It sorts/flips the graph along the reference, finds [snarls](../algorithms/bubble.md#terms)
internally (a vendored cactus / 3-edge-connected decomposition matching `vg snarls`), infers each bubble's
internal nodes from path intervals, scores path support + internal span, applies size/support filters, and
writes the sorted GFA + a bubble CSV + Bandage visualization files.

Mechanism, filters, and a worked trace: **[algorithms/bubble.md](../algorithms/bubble.md)**.

## Required inputs

- `-i, --gfa <graph.gfa>` — any GFA.
- `-r, --reference-path <name>` — reference path name or unique case-insensitive substring; orders the
  internal sort/flip and snarl finder. Not needed with `--snarls-in`.

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `-o, --out-prefix <p>` | output prefix | `bubble_calls` |
| `-s, --superbubbles` | emit only acyclic [superbubbles](../algorithms/bubble.md#terms) instead of all snarls | off (all snarls) |
| `--min-variant-bp <N>` | keep a bubble only if some path's internal span ≥ N (`0` = off) | `50` |
| `--min-path-support <N>` | keep only bubbles crossed by ≥ N paths | `0` (off) |
| `--merge-nearby-bp <N>` | merge consecutive bubbles ≤ N bp apart (after filters) | `0` (off) |
| `--no-flip` | sort but don't reorient nodes to the reference strand | off |
| `--snarls-in <path>` | use an external `vg snarls` JSONL (graph used as-is, no sort) | — |
| `--emit-snarls-jsonl <path>` | also write the internal snarls as `vg`-style JSONL | — |
| `--snarl-debug-tsv <path>` | per-candidate diagnostics ([columns](../algorithms/bubble.md#debug-tsv---snarl-debug-tsv)) | — |
| `--gtf <path>` | project a reference-coordinate GTF's genes onto reference nodes → `<prefix>.bandage_genes.csv` (needs a [PanSN](https://github.com/pangenome/PanSN-spec) `--reference-path`) | — |
| `-q, --quiet` | disable the progress bar | off |

## Outputs

| file | contents |
|------|----------|
| `<prefix>.bubbles.csv` | the bubble table consumed by `inspect`/`panphorte`/`call` (cols below) |
| `<prefix>.sorted.gfa` | reference-sorted/flipped graph — the input for downstream `panphorte` |
| `<prefix>.bandage_nodes.csv` | Bandage node colors (blue = candidate context, red = retained bubbles) |
| `<prefix>.bandage_genes.csv` | (with `--gtf`) `Name,Colour,Gene` per bubble |

`bubbles.csv` columns: `bubble_id, source, sink, inside_node_count, total_node_count, path_support,
min_inside_bp, max_inside_bp, inside_nodes`.

## Example

Matches `scripts/genes/lpa.sh` (run via the per-gene drivers in `scripts/genes/`):

```bash
./build/panvar bubble \
  -i tests/real_data/lpa.gfa.gz \
  -o results/real_data/lpa/bubble/bubble \
  --reference-path GRCh38 --gtf tests/real_data/Homo_sapiens.GRCh38.116.gtf.gz
```

Algorithm & worked example: see [algorithms/bubble.md](../algorithms/bubble.md). References:
[references.md](../references.md#bubble).
