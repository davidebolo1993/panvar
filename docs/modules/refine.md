# Module `refine`

CLI: `panvar refine`

## What it does

`refine` fixes graph artifacts (e.g. spurious `INS` + `DEL` events) at the graph level. It runs after `panphorte`, on the normalized graph, and for each selected bubble re-aligns the actual per-haplotype interior sequences with POA ([abPOA](https://github.com/yangao07/abPOA), affine-gap partial-order alignment) — collapsing artifacts into a single clean event. While this is safe in bubbles with no duplication signal:

- in a bubble that carries an unfolded copy-number signal - the reference or any haplotype revisits a non-REP interior node ≥2× (e.g. a paralog collapse or a private duplication `panphorte` left unfolded) — `refine` does nothing - i.e. the bubble is skipped entirely, because POA would linearize those copies and destroy the CN signal `call` reconstructs from them;
- in a bubble with a folded tandem (a `REP` node from `panphorte`), `refine` splits each haplotype's interior at every REP block, copies the `REP×n` run (per-haplotype copy count and orientation is preserved), and POA-aligns only the residual flanks around it.
  
The exact distinct-walk collapse means abPOA aligns only the handful of distinct interior sequences, so the binding cost is mainly interior length, not haplotype count.

Algorithm and worked trace: [algorithms/refine.md](../algorithms/refine.md).


## Required inputs

- `-i, --gfa <graph.gfa>` — a `panphorte` normalized and sorted GFA (`<prefix>.normalized.sorted.gfa`).
- `-b, --bubble-prefix-in <prefix>` (or `--bubbles-csv-in <path>`) — the input bubbles (`<prefix>.bubbles.csv`).
- `-r, --reference-path <name>` — reference path name or unique case-insensitive substring (orders the re-sort/re-snarl).
- `-o, --out-prefix <prefix>`.

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `--gtf <path>` | reference-coordinate GTF; also write `<prefix>.bandage_genes.csv` | — |
| `--bubble-id <id[,id...]>` | refine only these bubble ids (targeted mode) instead of auto over the whole locus | auto (all) |
| `--max-poa-bp <N>` | skip a residual segment whose median interior exceeds this | `5000` |
| `--max-walks <N>` | skip a residual segment with more than this many distinct walks | `500` |
| `--min-bubbles <N>` | only rebuild regions fusing ≥ N bubbles (`2` = cross-bubble over-splits only, leaves single clean bubbles untouched) | `1` |
| `--no-flip` | do not reorient nodes to the reference forward strand | — |
| `-q, --quiet` | disable progress logs | — |

## Outputs

| file | contents |
|------|----------|
| `<prefix>.normalized.sorted.gfa` | the refined, reference-sorted/flipped graph |
| `<prefix>.bubbles.csv` | bubble decomposition on the refined graph (the sorted/flipped bubbles) |
| `<prefix>.bandage_nodes.csv` | Bandage node colouring |
| `<prefix>.bandage_genes.csv` | Bandage gene track (only with `--gtf` and a PanSN reference) |

## Example

See the [LPA walkthrough](../walkthrough.md) for this module in a full end-to-end run.
