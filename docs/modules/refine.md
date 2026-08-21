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
| `--max-poa-bp <N>` | skip a residual segment whose longest sequence exceeds N bp, measured over the DISTINCT sequences abPOA receives rather than over every carrier | `5000` |
| `--max-poa-work <N>` | estimated abPOA budget in DP cells (longest × total bases) over those distinct sequences; bounds total work independently of the length bound | `0` (no bound) |
| `--max-walks <N>` | skip a residual segment carrying more than N distinct sequences (identical carriers collapse first, so this counts sequences, not walks) | `500` |
| `--min-bubbles <N>` | only rebuild regions fusing at least N bubbles | `1` |
| `--resnarl-min-variant-bp <N>` | interior-span filter applied to the re-snarled bubbles CSV | `50` |
| `--partial-path-policy <p>` | `skip` or `retain` a region only some paths fully traverse. Retaining keeps their old nodes AND the old edges between them, so refined and unrefined topology coexist; `retain` is experimental | `skip` |
| `--no-flip` | do not reorient nodes to the reference forward strand | off |

## Outputs

| file | what it is |
|------|------------|
| `<prefix>.normalized.sorted.gfa` | the refined, reference-sorted/flipped graph |
| `<prefix>.bubbles.csv` | bubble decomposition on the refined graph (the sorted/flipped bubbles) |
| `<prefix>.refine.report.tsv` | one row per region: the decision taken and why (see below) |
| `<prefix>.bandage_nodes.csv` | Bandage node colouring |
| `<prefix>.bandage_genes.csv` | Bandage gene track (only with `--gtf` and a PanSN reference) |

`refine.report.tsv` columns:

| column | meaning |
|--------|---------|
| `region` | the region's number in this run |
| `n_bubbles` | how many input bubbles were fused into it |
| `source`, `sink` | its outer anchor nodes |
| `decision` | `rebuilt` or `skipped` |
| `reason` | which guard fired, and the sizes behind it |

## Limitations

- A region carrying an unfolded copy-number signal is skipped entirely, because re-aligning it would linearize the copies and remove the signal `call` reconstructs from them. Refinement therefore does not reach the sites where duplications are typed.
- A region only some haplotypes traverse fully is skipped by default. Rewriting the rest while retaining those haplotypes' nodes would also retain the edges between them, so refined and unrefined topology would coexist at the same site — which the losslessness check cannot detect, since every haplotype still spells the same bases. `--partial-path-policy retain` rebuilds anyway and is experimental.
- The resource guards decide per region, so a locus can be partly refined. The report names which guard fired and the sizes behind it, so a skip can be acted on without re-running.
- Bubble ids are reassigned when the refined graph is re-snarled, so an id in the output does not refer to the same site as that id in the input.

## Example

See the [walkthrough](../walkthrough.md) for this module in a full end-to-end run.
