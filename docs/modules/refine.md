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
| `--max-poa-bp <N>` | skip a residual segment whose longest residual sequence, measured over the DISTINCT sequences abPOA receives rather than over every carrier.normalized.sorted.gfa` | the refined, reference-sorted/flipped graph |
| `<prefix>.bubbles.csv` | bubble decomposition on the refined graph (the sorted/flipped bubbles) |
| `<prefix>.bandage_nodes.csv` | Bandage node colouring |
| `<prefix>.bandage_genes.csv` | Bandage gene track (only with `--gtf` and a PanSN reference) |

## Example

See the [LPA walkthrough](../walkthrough.md) for this module in a full end-to-end run.

## Guards, policy and the decision report

`--max-poa-bp` bounds the LONGEST residual sequence in a segment; `--max-poa-work` optionally bounds the
estimated DP-cell budget (longest x total bases). Both are measured over the DISTINCT sequences abPOA is
handed, so replicating an identical haplotype cannot change a decision without changing the POA input.
`--max-walks` counts those distinct sequences, not walks.

A region that some paths traverse only partly cannot have those paths rewritten. Retaining their nodes
also retains the old edges between them, so pre-refinement topology survives beside the refined one --
which sequence losslessness cannot detect, because every path still spells the same bases.
`--partial-path-policy` defaults to `skip`; `retain` rebuilds anyway and is experimental. Measured on C4
and LPA the two policies are identical (3 rebuilt / 2 skipped, and 7 / 5), because no region at either
locus has a partial traverser.

`--resnarl-min-variant-bp` sets the interior-span filter for the re-snarled call-ready CSV, which
otherwise silently reapplies `bubble`'s own 50 bp default to an input built with a different one.

`<prefix>.refine.report.tsv` records one row per region: region number, bubble count, anchors, decision
and reason. Reasons name the guard that fired and the sizes behind it (carriers, distinct sequences,
longest, distinct total bp), so a skip can be acted on without re-running.

Losslessness is an acceptance criterion, not an assumption: every path's spelled sequence is compared
against its input, and every consecutive step pair must be joined by a link that exists in the
orientation walked, both before anything is written. The whole output family is staged and committed
only once that passes.
