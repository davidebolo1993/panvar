# Module `rebuild`

CLI: `panvar rebuild`

## What it does

Re-induces a locus graph (GFA — Graphical Fragment Assembly — from [pggb](https://github.com/pangenome/pggb)) before bubble decomposition, for the minority of graphs whose construction left them too fragmented to decompose sanely. It:
- runs a cheap degree gate over the input GFA to decide whether the graph is pathological, and passes healthy graphs through untouched;
- if pathological, rebuilds the locus by progressive graph generation: haplotypes are added most-complete-first, each aligned to the graph built so far and augmented with the parts that do not match;
- recovers a walk for every haplotype and writes a plain GFA ready for `bubble`.

pggb induces its graph with [seqwish](https://github.com/ekg/seqwish), whose transitive closure over all-pairs alignments can collapse sequences that recur across the locus — repeat copies, or the same short segment reused in different places — onto shared, very high-degree nodes, complicating downstream variant calling. Rebuilding avoids this because progressive re-alignment never closes over every pairwise match: it threads each haplotype colinearly, so recurring sequence branches locally instead of merging globally. In `rebuild`, graph generation is delegated to [minigraph](https://github.com/lh3/minigraph), while the module decides when to intervene, the haplotype order that generation consumes, and the emitted GFA.

Algorithm and worked trace: [algorithms/rebuild.md](../algorithms/rebuild.md).


## Required inputs

- `-i, --gfa <graph.gfa>` — the locus graph, with one `P` line per haplotype.
- `-o, --out <path>` — output GFA.


## Key options

| flag | what it does | default |
|------|--------------|---------|
| `--hub-degree <N>` | node degree (distinct link neighbours) at which a node counts as a hub | `50` |
| `--min-hubs <N>` | number of hubs at or above which the graph is called pathological | `10` |
| `--force` | rebuild even when the gate says healthy; for testing and for small inputs the gate is not calibrated on | off |
| `--kmer <N>` | k for the k-mer richness metric that orders haplotypes: distinct k-mers first, total k-mers breaking ties | `21` |
| `--min-var <N>` | minimum variant length augmented into the graph; smaller differences are left at the backbone allele | `50` |
| `--tmp-dir <path>` | parent directory for the per-haplotype FASTA scratch; a dedicated subfolder is created under it and removed on exit | beside `--out` |
| `-t, --threads <N>` | worker threads (`0` = auto) | `0` |
| `-q, --quiet` | disable progress/logs | off |

## Outputs

| file | contents |
|------|----------|
| `<out.gfa>` | the rebuilt graph: `S` lines, `L` lines carrying orientation, one `P` line per recovered haplotype |


## Notes

The rebuilt graph is structural, not base-exact. Differences shorter than `--min-var` are not augmented, so a haplotype's `P` line spells the backbone allele at those positions rather than its own sequence. Expect per-haplotype identity below `1.0`, with the shortfall concentrated in substitutions and short indels; lowering `--min-var` does not fully remove it, because the underlying generator does not augment single-base differences at any practical setting. This matters for anything that assumes a path spells its haplotype exactly. 
