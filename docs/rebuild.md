# Module `rebuild`

CLI: `panvar rebuild`

## What it does

Rebuilds a tangled locus graph so its sites can be found. Graphs made from all-pairs alignments ([pggb](https://github.com/pangenome/pggb)) can pile sequence that recurs across a locus onto a few nodes with hundreds of neighbours. `rebuild` builds the locus again one haplotype at a time, so recurring sequence stays in local branches, and folds duplicated genes so `call` can count their copies.

- A graph that is not tangled is copied through unchanged.
- The rebuilt graph is written only if every haplotype can still be read back from it; otherwise the input is written unchanged.
- The graph is structural: differences shorter than `--min-var` (50 bp) are not added, so a path spells there the allele of the node it walks, the backbone's or another haplotype's.

## How it works

Graph generation and mapping are [minigraph](https://github.com/lh3/minigraph)'s; realignments use [minimap2](https://github.com/lh3/minimap2) and [edlib](https://github.com/Martinsos/edlib).

1. A hub is a node with `--hub-degree` (50) or more neighbours on one side. A graph with `--min-hubs` (10) hubs or more is rebuilt; `--force` rebuilds any graph.
2. The haplotype with the most distinct k-mers or the reference (`-r`) goes first as the backbone; the others follow by distinct k-mers.
3. minigraph adds the haplotypes in that order: each is aligned to the graph so far and its differences of `--min-var` or more are added.
4. Each haplotype is mapped to the new graph (minigraph) and its walk read out. An assembly gap long enough to break the alignment in two is bridged with a node holding the haplotype's own N bases, so the walk spells the gap and not reference sequence. A large inversion breaks it in three (forward, reversed, forward); the three are joined into one walk with turns at the breakpoints, kept only if that walk spells the haplotype within the bounds.
5. Insertions and deletions of `--min-var` or more that minigraph left out, found by aligning a haplotype to its walk (minimap2), are written into the walk. Near a haplotype's end too few bases lie beyond an indel to align, so there the bases past it, `--min-var` or more, must match the walk's end; the rest is the indel. Sequence past the graph's ends is kept as tip nodes, nodes are cut where haplotypes end, and unused nodes are dropped.
6. minigraph keeps a duplication's copies apart, so a copy change would look like a long insertion. Duplicated stretches (repeated at 95% identity or more, 3 kb apart or more: `--fold-identity`, `--fold-min-copy-bp`) are found by aligning sequences to themselves (`--mode folded`, the default); each haplotype's copies are placed, and all copies are made to walk one shared unit (a `unit graph` from the copies). The unit is cut wherever a copy differs from it, so copies share every node but their own differences. Tandem repeats with a unit under 3 kb are left to `panphorte`.
7. Every haplotype's path, as written (after folding), must cover its sequence and spell it closely, within the bounds below. With `-r`, the reference must come back exactly.

## Required inputs

- `-i, --gfa <graph.gfa>`: the locus graph, one `P` or `W` walk per haplotype (plain or gzipped, not a pipe)
- `-o, --out <path>`: output `.gfa` or `.gfa.gz`

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `-r, --reference-path <name>` | reference path: seeds the graph and must come back exactly. A unique part of the name is enough | richest haplotype |
| `--mode <folded\|minigraph>` | fold copy-number modules, or write minigraph's graph as it is | `folded` |
| `--fold-identity <X>` | identity of a repeat that makes a module, 0.80 to 1; a copy is held to at least 0.90, or to this value if lower | `0.95` |
| `--fold-min-copy-bp <N>` | shortest copy that makes a module, 1000 or more; under 3000 it warns, as shorter tandem units are `panphorte`'s | `3000` |
| `--min-var <N>` | shortest difference added (minigraph `-L`) | `50` |
| `--min-align-len <N>` | shortest alignment that may add variation (minigraph `-l`); `0` scales it to the locus | `0` |
| `--base-align` | base-level alignment while generating (minigraph `-c`): finds more, up to twice as slow | off |
| `--clip-to-reference` | with `-r`: drop sequence running past the reference's ends instead of keeping it as a tip node | off |
| `--kmer <N>` | k-mer size for the order | `21` |
| `--hub-degree <N>`, `--min-hubs <N>` | the gate | `50`, `10` |
| `--force` | rebuild even a graph the gate passes | off |
| `--min-matched-cover <X>` | share of each haplotype's bases that must align and match | `0.95` |
| `--min-recovered-identity <X>` | identity each path must reach against its haplotype | `0.98` |
| `--allow-loss` | write the rebuild even when a check fails, recording why (never with an inexact reference) | off |
| `--tmp-dir <path>` | where scratch files go | beside `--out` |
| `-t, --threads <N>` | worker threads (`0` = all cores) | `0` |
| `-q, --quiet` | no progress logs | off |

A haplotype passes when its matched cover and identity reach the two thresholds, its path neither runs `--min-var` or more past its ends nor stops `--min-var` or more short of them, and its length is within max(1000, 4×`--min-var`, 0.5%) bp of the haplotype's.

## Outputs

| file | contents |
|------|----------|
| `<out>` | the rebuilt graph (the input unchanged when not accepted). Every segment has an `mn:Z:` tag: minigraph's segment name, or what `rebuild` added (tips, gaps, repaired insertions, the fold's unit and pieces) |
| `<prefix>.rebuild.report.tsv` | one row per haplotype, then `#` lines with the verdict, the thresholds (`#max_end_overrun_bp`, `#max_end_short_bp`, ...), the counts of each repair, and one `#module` line per copy-number module (where, unit length, copy range, folded or why not) |

Main report columns:

| column | meaning |
|--------|---------|
| `path` | the haplotype |
| `matched_cover` | matching bases / haplotype length, from minigraph's mapping; decides acceptance |
| `walk_identity` | the walk spelled and aligned end to end with the haplotype: 1 − edits / the longer one; decides acceptance |
| `length_diff_bp` | walk length minus haplotype length |
| `start_overrun_bp`, `end_overrun_bp` | bases spelled past the haplotype's ends |
| `start_short_bp`, `end_short_bp` | ACGT bases of the haplotype past the alignment at its start or end that the path does not hold (non-ACGT bases excluded) |
| `after_fold` | for a folded path: `exact`, `within` (inside the bounds) or `differs`; `.` if not folded |
| `status` | `ok`, `ok_after_fold`, `fails_after_fold`, or why it failed (`low_cover`, `low_identity`, `length_mismatch`, `end_overrun`, `end_short`, `not_recovered`, `reference_inexact`, ...) |

The other columns (`original_bp`, `recovered_steps`, `envelope_cover`, `chain_identity`, `gap_bp`, `clipped_bp`, `inversion_joined`) give the detail behind those.

### Visualization

To compare graphs, such as the input and the rebuilt one, `scripts/compare_graphs.py` writes per-graph and per-path tables (nodes, hubs, sequence, how each path changed) and `scripts/plot_compare_graphs.R` draws them, one panel per measure; both list their options under `--help`:

```
python3 scripts/compare_graphs.py -i <input.gfa> -i <rebuilt.gfa> -l input,rebuilt -o <compare prefix>
Rscript scripts/plot_compare_graphs.R --prefix <compare prefix>
```

The panels:

| panel | what it shows |
|-------|---------------|
| nodes, edges | the graph's size |
| sequence (Mb) | bases stored in all nodes |
| 1 bp nodes (%) | share of nodes one base long: SNPs and short indels cut the graph this finely |
| hubs, max degree | nodes with `--hub-degree` (50) or more neighbours on one side, as the gate counts them, and the most neighbours any node has |
| nodes/kb | nodes per kb of the longest path |
| graph bp / median haplotype bp | about 1 when the haplotypes share their sequence, higher when repeats are left unmerged, lower when a repeat is collapsed |
| self-loops | nodes linked to themselves: a repeat folded onto one node (a unit graph of several nodes has none) |
| paths as in first (%) | paths spelling the same sequence as in the first graph: 100% after a lossless step, low after a rebuild |
| revisits per path | steps minus distinct nodes: passes through a folded copy or a repeat node again |
| steps per kb | nodes a path walks per kb of its sequence |
| length vs first (bp) | each path's length minus its length in the first graph |


## Limitations

- Not base-exact: differences under `--min-var` are not kept, and a path takes the allele of the node it walks, which can be another haplotype's. Inside a folded module a path spells its own sequence only where both ends of its stretch place uniquely in its haplotype (`exact_stretches` on the `#module` line); otherwise it keeps the rebuilt walk's spelling. Lowering `--min-var` keeps more, with a denser graph.
- One haplotype that cannot be read back within the bounds refuses the whole rebuild (unless `--allow-loss`). A walk follows the haplotype's single best mapping chain, so a large difference the generator left out near an end can leave that end off the walk.
- Only arrays whose copies are `--fold-identity` identical or more and `--fold-min-copy-bp` or more apart are folded (95%, 3 kb), and only with 2 kb of unique reference sequence within 50 kb of each end; inverted copies are left as they are. Of two arrays a few kb apart, the second may be refused (`overlaps a module already folded`) or not found. A copy under `--fold-identity` a few kb from an array can be counted as one of its copies, and the module then spans the sequence between; when that copy lies before the array, the array may not be folded at all.
- Generation is serial: about a minute for a 200 kb locus, up to about 25 minutes for the largest panels.

## Example

`rebuild` runs before `bubble`, on graphs the gate flags.

`scripts/rebuild_examples.sh` rebuilds a tangled locus and four paralog loci from `tests/real_data` twice, seeded by the reference (`-r`) and without it, compares both rebuilt graphs with the input and plots the comparison, all in `tests/examples/rebuild`. The paralogs are rebuilt with `--force`, since their graphs pass the gate.
