# Design: `panphorte --refine`, realigning the residue around folded arrays

**Branch:** `panphorte-refine-design`, from `main` `e91b73c`.

**Status:** implemented with the changes in §6, which take precedence over §2–§5 where they differ. The user
documentation is `docs/modules/panphorte.md` and `docs/algorithms/panphorte.md` §5.

## 1. Problem

After panphorte folds a tandem array into a REP node, the bases outside the folded copies keep the input graph's
representation:
- partial unit copies;
- flank variation;
- the fragments at copy boundaries.

**At chr3:46654531–46940288 (bubble 5)**, the same 30 bp half copy `GTAGATATATATCTCTCTAGAGATATATAT` is spelled through three
different node chains, inherited unchanged from the input graph.
- Of the 156 nodes in the site, 69 are 1 bp.
- The separate `refine` step does not repair this. It skips the site (`REP skeleton mismatch`), because haplotypes
  differ in their REP arrangement:
  - 114 carry no array;
  - 345 carry one block;
  - 7 carry two blocks.

**What is wanted:** an opt-in step, inside panphorte, that realigns that residue at each site it folded. It must:
- include the haplotypes with no array copy;
- keep the result only when it is simpler;
- run only when nothing else forbids realignment.

## 2. Behaviour

`panvar panphorte … --refine` (off by default).
- Without the flag, panphorte's output is **byte-identical** to today's.
- With it, every bubble panphorte normalized (`normalized=yes`) is a candidate site.
- Bubbles panphorte did not fold are left to the standalone `refine`, which is unchanged.

### 2.1 When a site is refused (the existing refine guards, applied unchanged)

A site is left exactly as folded, with a recorded reason, when any of these holds:

| reason | why |
|---|---|
| `unfolded_dup` | a path revisits a non-REP interior node two or more times. **This is the duplication guard:** realigning would linearize copies whose count is the signal `call` reads |
| `partial_paths` | a path touches the interior without spanning both anchors |
| `ambiguous_anchor_traversal` | a path has two equally good traversals between the anchors |
| `reference_not_traversing` | the reference does not cross the site |
| `segment_too_long` / `poa_work` / `too_diverse` | the existing `--max-poa-bp`, `--max-poa-work` and `--max-walks` limits, measured over the distinct residual rows |
| `not_a_clean_chain` | the site's anchors do not reduce to one source and one sink |

**REP nodes, their self-loops, their orientation and every path's REP copy count are never modified.** Only the
non-REP residue is rewritten.

**The one refine guard dropped** is `REP skeleton mismatch`. Removing it is the point of this step (§2.2).

### 2.2 One joint alignment per site, with REP cut points

This replaces refine's "same skeleton, segment by segment" alignment, so that mixed REP arrangements and
**non-carriers** align together.

1. **Flatten each traversing path's interior** (oriented to the reference, as refine does) into:
   - its residual sequence: the interior with every REP run removed;
   - its **cut offsets**: the positions in that residual where a REP run sat, with the run copied verbatim.

   Examples:
   - a non-carrier: its whole interior, with no cut;
   - a one-block carrier: pre + post, with one cut;
   - a two-block carrier: pre + mid + post, with two cuts.
2. **The distinct rows** are keyed by (residual sequence, cut offsets). Identical haplotypes therefore share a row by
   construction, and get identical walks.
3. **abPOA** aligns the distinct residual sequences: all of them together, carriers and non-carriers. The
   non-carriers' bases share columns with the carriers' flanks wherever they align.
4. **Compaction** as in refine (a run of columns traversed by exactly the same rows becomes one node), with one addition:
   **a forced node boundary at every row's cut column.** No rebuilt node can span a REP insertion point, so each
   carrier's walk can be spliced back as `nodes(pre) · REP run · nodes(mid) · … · REP run · nodes(post)`.
5. **The splice** works as in refine: the old non-REP interior nodes are dropped (if no path still uses them), the
   rebuilt nodes are added, every traversing path's interior is rewritten, and missing adjacencies are added as `0M`
   links.

### 2.3 "Keep only if simpler" (per site)

The site is measured before and after, over its non-REP interior:
- **N:** distinct non-REP nodes;
- **S:** the sum, over traversing paths, of their non-REP interior steps;
- **P:** nodes used by exactly one path (reported, not decisive).

**The rebuilt site is accepted** iff `N_after ≤ N_before` **and** `S_after ≤ S_before`, with at least one strict.
Otherwise the folded site is kept unchanged and recorded as `kept_original:not_simpler`. Simpler never means "the
aligner ran".

**Example:** a plain realignment of this repeat without folding went from 160 to 374 nodes. Under this rule it would be
rejected.

### 2.4 Acceptance (whole run, before anything is written)

These are checked after every accepted site has been spliced. A failure is a hard error, and nothing is written.

| check | requirement |
|---|---|
| lossless against the **folded** graph | every path spells exactly what it spelled after folding. Approximate folding is already lossy inside copies by design; this step must add no change |
| REP invariance | every path's ordered list of (REP node, orientation, copy count) is unchanged |
| links | every consecutive step pair has a link in the orientation walked |
| identical haplotypes | paths whose folded walks spelled identical sequence over a refined site walk identical node routes there |
| the existing "represented twice" guard | re-run over the refined graph |

### 2.5 Report

The panphorte report (`<prefix>.panphorte.report.tsv`) gains, per bubble:
- `refine_status`: `not_requested`, `not_folded`, `refined`, `kept_original:not_simpler` or `skipped:<reason>`;
- `refine_nodes_before`, `refine_nodes_after`, `refine_steps_before`, `refine_steps_after`, `refine_single_path_nodes_before`,
  `refine_single_path_nodes_after`.

The run summary gains counts of refined, kept and skipped sites.

### 2.6 Options

| option | values | default | effect |
|---|---|---|---|
| `--refine` | flag | off | enable the step |
| `--refine-max-poa-bp` | N | 5000 | as refine's `--max-poa-bp` |
| `--refine-max-walks` | N | 500 | as refine's `--max-walks` (distinct residual rows) |
| `--refine-max-poa-work` | N (0: none) | 0 | as refine's `--max-poa-work` |

## 3. Implementation plan

1. **Extract refine's region machinery** into a shared internal module (`src/refine_core.cpp`, internal header):
   interior extraction, the guards, the distinct rows, `poa_msa`, `build_interior_graph`, the splice and the
   losslessness check.
   - The standalone `refine` keeps its exact current behaviour, including the skeleton-equality requirement.
   - **Its outputs must be byte-identical before and after the extraction:** `refine_stats` passes, and the saved chr3
     refine outputs are reproduced byte for byte.
2. **Add the joint mode** (§2.2): flattening with cut offsets, distinct rows keyed by (sequence, cuts), and forced
   boundaries in `build_interior_graph`, which is exercised only by the new mode.
3. **Add the panphorte step:** after folding, the obsolete-route removal and the "represented twice" check; before
   sort/flip and re-snarl. It runs for each folded bubble (anchors = the bubble's source and sink), followed by the
   simpler rule and the acceptance checks. The report columns and options follow.

## 4. Verification

**Fixtures** (in `tests/panphorte_stats.sh`):

| fixture | expectation |
|---|---|
| a folded site with flank variation split into 1 bp nodes, carriers and non-carriers | `refined`; fewer nodes and steps; identical haplotypes on identical routes; lossless; REP counts unchanged |
| mixed arrangements: no array, one block, two blocks | `refined`; every cut respected (REP runs spliced at their original residual offsets); lossless |
| an unfolded duplication in the site | `skipped:unfolded_dup`; the site unchanged byte for byte |
| a repeat whose realignment produces more nodes | `kept_original:not_simpler`; the site unchanged |
| no `--refine` | output byte-identical to the current panphorte |

**Locus check (chr3, the colleague's pipeline + `--refine`):**
- lossless and REP invariance across all 466 paths;
- the site `refined` or `kept_original` with before/after figures;
- whether the three half-copy chains merge;
- the standalone `refine` afterwards still behaves.

**Regression:**
- the full enabled suite (including `panphorte_stats` and `refine_stats`);
- the standalone `refine` byte-identical;
- `panphorte` without `--refine` byte-identical on the test loci and on chr3.

## 5. Risks and limits

**The cut boundaries can add nodes.** Forced boundaries at each carrier's cut can split nodes that would otherwise merge.
If the cut offsets vary a lot between carriers, the simpler rule will keep the original.

**POA in tandem repeats** places indels by heuristic (no left-normalization). The joint alignment may still stagger
identical insertions. The simpler rule bounds the harm, but merging the three chains at chr3 is a hope to be measured,
not a guarantee.

**Cost:** one abPOA per folded site over the distinct residual rows (on the order of 50 at chr3 bubble 5; up to about 2.3 kb for the one haplotype carrying a ~1.5 kb insertion, below
the 5,000 bp limit),
within the existing limits.

**Integration:**
- **Where it lives:** `main` (panphorte and refine live there).
- **The genotype candidate** (M15-I) never touches `src/panphorte.cpp`, `src/refine.cpp` or `tests/panphorte_stats.sh`,
  so the planned merge stays independent.
- **The shared module:** a new file. The genotype code does not use refine.

## 6. As implemented (2026-10-22)

### 6.1 Changes from the proposal

| proposal | implemented | why |
|---|---|---|
| one joint alignment of every crossing haplotype, non-carriers included, with forced node boundaries at each REP cut (§2.2) | **per side**: one alignment of the carriers' flanks before the array, one of those after it; non-carriers are left as folded | measured at chr3 bubble 5 (§6.2): the joint alignment was rejected by the simpler rule (12,090 → 50,452 steps) |
| REP runs as cut points, one residual segment per gap between runs | the **array block**: a carrier's first to last REP step, verbatim, fragments between copies included | at lower `--min-similarity`, folded arrays carry interned fragments between copies (`REP GT REP GT REP`), which the per-run duplication guard read as an unfolded duplication |
| the duplication guard over the whole interior | over the **flanks** only | the same reason; nothing inside the block is realigned, so nothing in it can be linearized |
| a shared module extracted from `refine` | a separate file, `src/panphorte_refine.cpp`, with its own copies of the interior, compaction and splice helpers | the standalone `refine` is untouched, so it is byte-identical by construction |
| new columns in `panphorte.report.tsv` | a separate `<prefix>.panphorte.refine.tsv` | the report is written inside the folding loop; a separate file keeps the output without `--refine` byte-identical |
| — | `--refine` needs `--reference-path` | sites are aligned in reference orientation and the reference must cross each one |
| — | the old route check is per rewritten carrier; links walked only by old routes are dropped | non-carriers keep old flank nodes, so the old route could otherwise stay walkable through them |

### 6.2 Evidence (chr3:46654531–46940288, the colleague's pipeline, `--min-similarity 0.97`)

**Why non-carriers are left as folded.** The 114 haplotypes without a folded copy still carry their array as literal sequence:
- 6–10 copies at 87.5–96.4% identity to the folded unit, so none folds at 0.97;
- their interiors are 470–824 bp, against 218–296 bp of carrier flanks.

POA of tandem repeats of varying length gives no consistent alignment, as the measurements before → after realignment show:

| what is realigned | nodes | steps |
|---|---|---|
| non-carriers alone | 115 → 126 | 6,863 → 8,715 |
| carriers, flanks glued into one sequence | 52 → 39 | 5,227 → 6,948 |
| carriers, per side | 52 → 30 | 5,227 → 4,275 |
| everyone jointly (§2.2) | 155 → 572 | 12,090 → 50,452 |

**Result at bubble 5:**
- `refined`: 153 → 140 nodes and 12,054 → 11,115 steps outside the blocks;
- 1 bp nodes in the site: 69 → 53;
- the 30 bp half copy in HG00126#1, HG00128#2 and HG00140#1 is now one shared 5-node chain, where it was three different chains;
- spelled sequence and REP copies identical on 466 of 466 paths;
- 45 distinct walks for 45 distinct sequences;
- the run is deterministic, and the standalone `refine` downstream makes the same decisions as before.

**Lower thresholds do not help.**
- At 0.95 folding is refused (a rewritten route would still be walkable).
- At 0.93 / 0.90 / 0.87, 462 of 466 haplotypes fold, but 3 / 3 / 7 distinct haplotype sequences are erased. The newly folded arrays start at a different offset (left flank 63 vs 78 bp) and keep 0–3 divergent copies in their right flank (114–295 bp). `--refine` then rejects the site as not simpler.
- The copies are a mix of variants: 52–58 bp, differing by AT-run indels; five units cover 89% at 0.97. Folding them would need several REP nodes per site, which is not designed here.

**Other loci.** Without `--refine`, the output is byte-identical to `main` (`e91b73c`) at chr3 and at chr17, ACOT, ANKRD36C, C4, CYP2D6, GSTM1 and LPA (similarity 1.0 and 0.97). With `--refine`:
- nothing folds at chr17, ACOT, C4, CYP2D6 or GSTM1;
- ANKRD36C bubble 6 and LPA bubble 7 (0.97) are skipped as `unfolded_dup`;
- MYOM2 was not run (its graph is to be rebuilt).

All locus results above were re-measured with the final code.
