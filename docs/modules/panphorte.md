# Panphorte Module (Module 2 — bubble repeat normalization)

CLI: `panvar panphorte`

## What it does

Rewrites **tandem-repeat** bubbles into a compact, copy-number-explicit form: each detected tandem array is
collapsed to a single repeat-unit (`REP`) node with a self-loop, so copy number becomes the number of
self-loop traversals (which `call --cn-from-multiplicity` reads directly). Writes a new GFA for `inspect`/
`call`.

Mechanism, exact-vs-approximate collapse, and a worked trace:
**[algorithms/panphorte.md](../algorithms/panphorte.md)**.

> **Scope.** For **tandem** loci (a unit repeated in-line, e.g. LPA KIV-2). **Paralog clusters** (C4,
> CYP2D6) are folded by PGGB onto shared nodes — not a contiguous tandem — so their copy number is recovered
> by `call --cn-from-coverage` on the `bubble` graph, *not* from panphorte (see the
> [CN-topology table](../README.md#copy-number-one-method-per-locus-topology)). Running panphorte there is
> harmless but its graph is not the call substrate.

## Required inputs

- `-i, --gfa <graph.gfa>` — the `bubble` `*.sorted.gfa` (W- or P-line; output preserves the line type).
- one of `-b, --bubble-prefix-in <bubble-prefix>` (auto-uses `<prefix>.bubbles.csv`) or `-c, --bubbles-csv <path>`.

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `-o, --out-prefix <p>` | output prefix (GFA + report; which GFA depends on `--reference-path`) | — |
| `-r, --reference-path <name>` | sort + re-snarl the normalized graph along this reference → call-ready in one step | — |
| `--min-similarity <f>` | identity to treat a block as a copy of the unit; `1.0` = exact (sequence-preserving), `< 1.0` = approximate/lossy collapse of divergent copies | `1.0` |
| `--min-unit-bp <N>` | minimum repeat-unit span to normalize | `50` |
| `--min-copies <N>` | tandem copies needed (in some haplotype) to treat a bubble as an array; once an array, every haplotype with ≥1 copy folds | `2` |
| `--max-interruption-frac <f>` | max fraction of an array's bp that may be interruptions | `0.25` |
| `--threads <N>` | workers for the approximate seed scan / copy detection (`0` = auto) | `0` |
| `--gtf <path>` | after re-sorting, project genes → `<prefix>.bandage_genes.csv` (needs PanSN `--reference-path`); separate from `bubble --gtf` because collapse renumbers nodes | — |
| `-q, --quiet` | disable progress/logs | off |

## Outputs

| file | contents |
|------|----------|
| `<prefix>.normalized.sorted.gfa` | sorted, call-ready normalized graph (**with** `--reference-path`) |
| `<prefix>.normalized.gfa` | unsorted normalized graph (**without** `--reference-path`; must be re-sorted before `call`) |
| `<prefix>.bubbles.csv` | re-snarled bubbles (with `--reference-path`) |
| `<prefix>.panphorte.report.tsv` | per bubble: `bubble_id, normalized, unit_bp, paths_normalized, min_copies, max_copies, nodes_collapsed` |
| `<prefix>.panphorte.copies.tsv` | (approximate) per (haplotype, array): `path_name, sample, bubble_id, copies, unit_bp, orientations, mean_identity, region_bp, from_node, to_node` — the CN provenance for `call` |
| `<prefix>.bandage_nodes.csv` | Bandage node colors |

## Example

Matches `scripts/genes/lpa.sh`:

```bash
./build/panvar panphorte \
  -i results/real_data/lpa/bubble/bubble.sorted.gfa \
  --bubble-prefix-in results/real_data/lpa/bubble/bubble \
  --reference-path GRCh38 -o results/real_data/lpa/panphorte/panphorte \
  --min-similarity 0.97
```

Algorithm & worked example: see [algorithms/panphorte.md](../algorithms/panphorte.md). References:
[references.md](../references.md#panphorte).
