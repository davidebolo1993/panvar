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
| `--min-align-len <N>` | minimum alignment length that may contribute events; `auto` scales it to half the seed haplotype | `auto` |
| `--tmp-dir <path>` | parent directory for the per-haplotype FASTA scratch; a dedicated subfolder is created under it and removed on exit | beside `--out` |
| `-t, --threads <N>` | worker threads (`0` = auto) | `0` |
| `-q, --quiet` | disable progress/logs | off |

## Outputs

| file | contents |
|------|----------|
| `<out.gfa>` | the rebuilt graph: `S` lines, `L` lines carrying orientation, one `P` line per recovered haplotype |


## Notes

The rebuilt graph is structural, not base-exact. Differences shorter than `--min-var` are not augmented, so a haplotype's `P` line spells the backbone allele at those positions rather than its own sequence. Expect per-haplotype identity below `1.0`, with the shortfall concentrated in substitutions and short indels; lowering `--min-var` does not fully remove it, because the underlying generator does not augment single-base differences at any practical setting. This matters for anything that assumes a path spells its haplotype exactly.

Lowering `--min-var` does recover a useful part of it. On a locus of a few tens of kb, moving from `50` to `10` roughly doubles the number of haplotypes that still spell a distinct sequence, at the cost of a denser graph; below `10` the curve flattens. The default stays at `50` because that is the scale the downstream calling is tuned for.

The generator drops an alignment entirely when it is shorter than `--min-align-len`, before `--min-var` is ever consulted, and its native default assumes chromosome-scale input. On a locus graph that threshold can reject every alignment, in which case nothing is augmented and the result collapses to the seed alone. `auto` avoids this by scaling the threshold to the locus, and as a backstop a rebuild that recovers no variation at all is discarded and the input graph is passed through unchanged.

## Example

`rebuild` runs before `bubble` and only on graphs the gate flags, so it is not part of the main LPA run. The [rebuild section of the walkthrough](../walkthrough.md#rebuild--re-inducing-a-tangled-locus) shows it on MYOM2 (chr8) instead — a tangled locus before and after re-induction.

## What the coverage numbers mean

Three different questions, which one number used to answer:

| reported | definition | what it detects |
|----------|------------|-----------------|
| `envelope cover` | `(qe − qs) / haplotype length` | where the alignment begins and ends — the chain's outer span |
| `matched cover` | matching bases / haplotype length | how much of the haplotype is **genuinely** aligned; an internal gap shows up here where the envelope hides it |
| `chain identity` | matches / alignment-block length | how well it matches where it does align |
| `recovered-walk identity` | edit identity of the walk re-spelled from the **rebuilt** graph against the original haplotype | what a caller actually gets back — an error anywhere in mapping, chain choice or emission |

The distinction is not academic. On the C4 locus (131 haplotypes, forced rebuild) the envelope reads **0.9999** while matched cover is **0.9910** and the worst recovered walk is **0.9740** — so a `98–99%` acceptance threshold read off the envelope would pass a haplotype that comes back 97.4% correct. `recovered-walk identity` is the quantity such a threshold has to be stated in, and the run reports how many haplotypes fall below 0.99.

## The acceptance contract

`rebuild` drives minigraph, which augments variation **above `--min-var`** — sub-threshold differences (SNPs, short indels) are collapsed by construction, so a recovered walk is never byte-identical to the haplotype it came from. Losslessness is therefore the wrong contract. What is checked instead is structural, plus a threshold:

- every input path comes back as a walk;
- every consecutive pair of steps on that walk is joined by a link that exists in the emitted graph;
- matched cover ≥ `--min-matched-cover` (default `0.95`) and recovered-walk identity ≥ `--min-recovered-identity` (default `0.98`);
- if `-r/--reference-path` is given, that path must be present and meet the same bounds.

A recovered-walk identity that could not be **computed** counts as a failure, not a pass. Absence of evidence is not evidence of a good recovery, and treating it as one would let exactly the unverifiable cases through the check meant to catch them.

**If any of these fails, the rebuilt graph is discarded and the original is written unchanged**, with the reason stated. A graph that silently drops or truncates a haplotype is worse than no rebuild, because everything downstream would then agree with it. `--allow-loss` accepts anyway and records what was violated.

**Acceptance proves fidelity, not untangling.** It says every haplotype came back within the bounds; it says nothing about whether the graph got simpler, which is the reason for rebuilding. The run therefore reports hubs, maximum handle degree and self-loops before and after — on the same per-end measure — and labels the outcome `untangled` or `NOT untangled`. A rebuild can be perfectly faithful and still not have helped.

The reference clause is about **recovery, not seeding**. Which haplotype seeds the graph stays richness-driven — that is deliberate — but a reference that did not come back invalidates every reference-relative bubble and call downstream, so it is checked.

### Audit sidecar

`<out>.rebuild_audit.tsv`, one row per path: `original_bp`, `recovered_steps`, `envelope_cover`, `matched_cover`, `chain_identity`, `walk_identity`, and a `status` of `ok` / `not_recovered` / `low_cover` / `low_identity` / `identity_unavailable`. Trailing `#`-prefixed rows carry the global verdict, the rejection reason and the bounds that applied, so the file explains itself without the log. It is staged and renamed like the graph, and a failure to write it is an error rather than a silent omission. The verdict can then be read rather than trusted — on C4 at the default bounds, 130 paths are `ok` and one is `low_identity` at 0.9740, which is the single haplotype that rejects the run.

| flag | what it does | default |
|------|--------------|---------|
| `--min-recovered-identity <X>` | per-path recovered-walk identity bound | `0.98` |
| `--min-matched-cover <X>` | per-path matching bases / haplotype length | `0.95` |
| `-r, --reference-path <name>` | this path must be recovered within the bounds | none |
| `--allow-loss` | accept a rebuild that fails the contract, recording why | off |
| `--audit <path>` | where to write the sidecar | `<out>.rebuild_audit.tsv` |

## Reproducibility

Three properties a rebuild has to have, none of which is visible in a single run's output:

- **Deterministic.** Haplotypes are ordered by richness, then by path **name**, then by original index. `std::sort` is not stable, so equal richness previously produced an implementation-defined order — and ties are the norm in a panel of near-duplicate haplotypes, not an edge case. Two runs could seed from different haplotypes and build different graphs.
- **Thread-invariant.** One thread and eight produce byte-identical output, audit included.
- **Strand-independent.** k-mer richness is counted over **canonical** k-mers (a k-mer and its reverse complement are one), so reverse-complementing the whole graph does not change the seed. Without that, richness followed whichever strand the GFA happened to store.

k-mer counting also only counts windows free of ambiguity, and `total` counts exactly the windows `distinct` does — previously `total` included windows containing `N` while `distinct` skipped them, so the redundancy figure `total − distinct` was wrong by however much ambiguity a haplotype carried. `k > 31` now behaves identically to `k ≤ 31`; it used to keep every window verbatim, ambiguity and strand included.

## Chain selection

When a haplotype produces several graph chains, one has to win — concatenating them would stitch an incoherent walk. The choice is by **matching bases**, then identity within the aligned block, then mapping quality, with a deterministic tie-break. Ranking by longest query *span* let a chain that reaches further through a large internal gap beat one that genuinely aligns more.

## The pathology gate

Degree is measured **per handle**. A bidirected node has two ends, and pooling them conflates two different shapes: 25 neighbours on each side is a clean two-sided branch, while 50 on one side is the tangle the gate exists to find — pooled, both read 50. A node's degree is the larger of its two handle degrees, and self-loops (how a tandem array appears after folding, and not pathology) are counted separately. On the six reference loci this lowers `maxdeg` — c4 11 → 9, LPA 20 → 17, ANKRD36C 21 → 19 — without changing any verdict, since all sit far below the default threshold of 50.
