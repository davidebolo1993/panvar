# Making figures for talks

Self-contained, slide-ready figures that explain what each panvar module does, built only from outputs
the pipeline already writes. Two helper scripts cover the graph schematics; the population-scale view
reuses the `call` heatmap.

- `scripts/plot_graph_schematic.R` — tube-map schematics (`--mode panphorte`, `--mode sv`).
- `scripts/plot_sv_map.R` — the population node-level event heatmap (see [call](modules/call.md)).

All need R + `ggplot2` (`conda install -c conda-forge r-ggplot2`). Each writes a `.png` and a `.pdf`.

## bubble — which nodes are the sites

`bubble` already writes `<prefix>.bandage_nodes.csv`. Open the original GFA in **Bandage**, load that CSV
(`File → Load CSV`), and color by it to show the picked bubble nodes on the graph. No extra tooling.

## panphorte — tandem normalization (one slide per event)

`--mode panphorte` reads `<prefix>.panphorte.copies.tsv` (written by **approximate**-mode panphorte,
`--min-similarity < 1.0`) and draws, per normalized bubble: **before** = the repeat unit drawn `copies`
times (orientation per copy; short copies narrower), **after** = a single `REP` box with a self-loop
annotated `×copies`. It picks the reference plus the min/max-copy haplotypes as the 3-4 examples.

```bash
panvar bubble    -i graph.gfa -o out/bubble --snarls-in graph.snarls.jsonl
panvar panphorte -i graph.gfa --bubble-prefix-in out/bubble -o out/panphorte --min-similarity 0.70
Rscript scripts/plot_graph_schematic.R --mode panphorte \
  --copies out/panphorte.panphorte.copies.tsv --bubble-id <N> --out out/slide_panphorte_<N>
```

The four reference scenarios, each its own slide:

| locus | what to show | how |
|---|---|---|
| **LPA** (KIV-2) | a clean tandem VNTR collapses to one REP unit | `--min-similarity 0.90` |
| **C4** | the long/short modules only fold together at low similarity (the short copy is missing the ~6.4 kb HERV-K) | `--min-similarity 0.70` (try `0.90` to show it *not* folding) |
| **GSTM1** | a dispersed segdup with no adjacent identical pair is **not** seeded → stays intact for `call` | the bubble is absent from `copies.tsv`; show it via the `--mode sv` DUP panel instead |
| **CYP2D6/2D7** | a folded paralog cluster the reference itself traverses ≥2× → stays intact, recovered by `call --cn-from-coverage` | as GSTM1 — use the `call` side |

## call — how each SV type looks in the graph (one panel per type)

`--mode sv` reads the `call` VCF + `variant_paths.tsv` + `node_track.tsv`, picks one representative record
per SV type in a bubble, and draws the reference (top) vs one carrier (bottom), zoomed to the event:
**DEL** = a red deleted region on the reference, carrier shorter; **INS** = a green inserted block on the
carrier; **INV** = an orange reversed block; **DUP** = the unit repeated `×CN` on the carrier. Segment
width is `sqrt(bp)` so kb-scale events stay readable.

```bash
panvar call -i graph.sorted.gfa --bubble-prefix-in out/bubble2 \
  --reference-path grch38 -o out/call --cn-from-multiplicity --cn-from-coverage
Rscript scripts/plot_graph_schematic.R --mode sv \
  --vcf out/call.region.vcf --variant-paths out/call.variant_paths.tsv \
  --node-track out/call.node_track.tsv --bubble-id <N> --out out/slide_sv_<N>
# or a single record:  --variant-id bubble4_DEL_2019
```

## call — population view (who carries what)

`scripts/plot_sv_map.R` is the companion population-scale figure: a node-level event heatmap
(haplotypes × bubble nodes, colored by called event, CN-shaded for duplications) in the same `odgi viz`
style as the `inspect` coverage heatmaps. Good for the "this variant segregates across the cohort" slide.

```bash
panvar inspect -i graph.gfa --bubble-prefix-in out/bubble --bubble-id <N> --cluster -o out/inspect/b<N>
Rscript scripts/plot_sv_map.R \
  --node-counts out/inspect/b<N>.bubble_<N>.node_counts.tsv \
  --vcf out/call.region.vcf --clusters out/inspect/b<N>.bubble_<N>.clusters.tsv \
  --reference-path grch38 --out out/slide_map_<N>
```
