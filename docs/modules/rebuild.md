# Module `rebuild`

CLI: `panvar rebuild`

## What it does

Re-induces a locus graph (GFA — Graphical Fragment Assembly) before bubble decomposition, for the minority of graphs whose construction left them too fragmented to decompose sanely. It:
- runs a degree gate over the input graph to decide whether it is pathological, and passes healthy graphs through untouched;
- if pathological, rebuilds the locus by progressive graph generation: haplotypes are added most-complete-first, each aligned to the graph built so far and augmented with the parts that do not match;
- recovers a walk for every haplotype, checks each against the original sequence, and writes a plain GFA ready for `bubble`.

Graphs induced by transitive closure over all-pairs alignments can collapse sequence that recurs across a locus — repeat copies, or a short segment reused in different places — onto shared, very high-degree nodes, which makes downstream variant calling hard. Progressive re-alignment instead threads each haplotype colinearly, tending to keep recurrence in local branches rather than merging it globally. Graph generation is delegated to [minigraph](https://github.com/lh3/minigraph); this module decides when to intervene, in what order haplotypes are added, whether the result is acceptable, and what is emitted.

Algorithm and worked trace: [algorithms/rebuild.md](../algorithms/rebuild.md).

## Required inputs

- `-i, --gfa <graph.gfa>` — the locus graph, with one `P` or `W` walk per haplotype; accepted rebuilds are emitted with `P` lines
- `-o, --out <path>` — output GFA

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `--hub-degree <N>` | node degree at which a node counts as a hub, measured per node side | `50` |
| `--min-hubs <N>` | number of hubs at or above which the graph is called pathological | `10` |
| `--force` | rebuild even when the gate says healthy; for testing and for inputs the gate is not calibrated on | off |
| `--kmer <N>` | k for the richness metric that orders haplotypes: distinct k-mers first, total k-mers breaking ties | `21` |
| `--min-var <N>` | minimum variant length augmented into the graph; smaller differences stay at the backbone allele | `50` |
| `--min-align-len <N>` | minimum alignment length that may contribute events; `0` scales it from the haplotype lengths | `0` (auto) |
| `--min-recovered-identity <X>` | required identity of each recovered walk against its original haplotype | `0.98` |
| `--min-matched-cover <X>` | required fraction of each haplotype covered by matching bases | `0.95` |
| `-r, --reference-path <name>` | this path must be recovered within the same bounds; exact match wins, else a unique substring | — |
| `--allow-loss` | accept a rebuild that fails the contract, recording what it violated | off |
| `--audit <path>` | per-path audit TSV | `<out>.rebuild_audit.tsv` |
| `--tmp-dir <path>` | parent directory for scratch; a dedicated subfolder is created under it and removed on exit | beside `--out` |
| `-t, --threads <N>` | worker threads (`0` = auto) | `0` |
| `-q, --quiet` | disable progress logs | off |

## Outputs

| file | contents |
|------|----------|
| `<out>` | the rebuilt graph: `S` lines, `L` lines carrying orientation, one `P` line per recovered haplotype. If the rebuild is not accepted this is the input graph, unchanged |
| `<out>.rebuild_audit.tsv` | one row per haplotype, plus the run's verdict |

Audit columns:

| column | meaning |
|--------|---------|
| `path` | the haplotype's name in the input graph |
| `original_bp` | its length in the input |
| `recovered_steps` | how many nodes its recovered walk visits |
| `envelope_cover` | fraction of the haplotype between the first and last aligned base, which does not see internal gaps |
| `matched_cover` | fraction of the haplotype covered by bases that actually match, which does |
| `chain_identity` | identity within the aligned region |
| `walk_identity` | identity of the walk re-spelled from the rebuilt graph against the original haplotype |
| `status` | `ok`, or why this haplotype failed: `not_recovered`, `low_cover`, `low_identity`, `identity_unavailable` |

The trailing lines record the run's disposition (`#verdict`), the reason if it was refused (`#reason`), and the thresholds that applied.

## Limitations

- The rebuilt graph is structural, not base-exact. Differences shorter than `--min-var` are not augmented, so a haplotype's walk spells the backbone allele at those positions rather than its own sequence. Expect per-haplotype identity below 1, concentrated in substitutions and short indels.
- Lowering `--min-var` recovers part of that at the cost of a denser graph, but not all of it: the generator does not augment single-base differences at any practical setting.
- Acceptance checks fidelity, not improvement. A rebuild that faithfully reproduces every haplotype but does not reduce tangling is accepted, and the structural change is reported separately so it can be judged.
- Sequence a haplotype does not share with the graph is dropped rather than added as private nodes, so a haplotype that cannot be recovered within the bounds causes the whole rebuild to be refused rather than partially represented.

## Example

`rebuild` runs before `bubble` and only on graphs the gate flags, so it is not part of the main walkthrough run. See the [rebuild section of the walkthrough](../walkthrough.md#rebuild--re-inducing-a-tangled-locus).
