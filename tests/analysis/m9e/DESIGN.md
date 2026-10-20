# M9-E: why the three improving neighbours were missed (bounded reachability audit)

**Cases:** M9-D's E1 (LOO d1, sv:5), E3 (LOO d3, sv:5) and E6 (LOO d6, sv:5). The improving alternatives are frozen as
M9-D defined them: context-A truth replacements whose total exceeds the reported winner (from M9-D's
`construct/<case>.candidates.tsv` and `analyse/contrasts.tsv`). No new truth-derived candidates are selected.

**Not done:** caller change, new search, tuning, confidence remapping, clustering, promotion, full caller rerun, merge,
push or cleanup.

**Results:** `tests/results/m9e_reachability_20261020/`.

## Sources and code facts used

**Saved traces from the M9-C runs:**
- `repair_rounds.tsv`: each trunk round's move, its `local_pair_rank`, and the final scan's
  generated/scored/duplicate/suppressed/self/prior-rejected counts;
- `tie_branches.tsv`: branch terminals, rounds and stop;
- `mosaic_candidates.tsv`: every scored label pair, with `scored_in` (`union`, `roundN`, `branchK_roundR`);
- `locus_genotype.tsv`: the reported winner.

**Neighbourhood** (`genotype_repair.cpp`, `CoordinateRepairSearch`):
- **The local pairs are ranked once per block,** by proposal emission descending, then `(allele_x, allele_y)` ascending,
  and the top `repair_local_pairs` (32, from the run's status) are kept.
- **A neighbour sets one block to one local pair** on both haplotypes, in either assignment.
- **A neighbour is then skipped** if it is self, already visited, or a duplicate within the round; it is **rejected** if
  its sequence prior is −∞; otherwise it is **scored**.
- **The list does not depend on the refinement context.** The full score does.

**Markers** (`marker_model.cpp`, `build_marker_panel`):
- markers are k = 31 syncmers counted **inside each block allele's own sequence**;
- markers shared by several blocks are dropped;
- so an allele shorter than 31 bp carries no proposal marker, and a junction-spanning k-mer is never a proposal marker.
- There is no existing per-marker dump. Marker availability is therefore read from allele lengths (the M9-D model dump)
  and from the emission tie structure. Individual markers are **not** enumerated.

## Units

| unit | timeout | purpose |
|---|---|---|
| `fixture` | 10 s | context-aware inconsistency classification on six synthetic rows (A against the reported winner; B carrier versus its called insertion) |
| `ranks.LOO.d1`, `ranks.LOO.d3`, `ranks.LOO.d6` | 60 s each | `genotype_probe score-set --rank-blocks sv:5,ctx:2,ctx:4,ctx:6,ctx:8,ctx:12,ctx:16,ctx:18,sv:13,sv:17 --emissions-only 1`: the run's model, curation, reads and library, with `repair_local_pairs` from the run's status |
| `analyse` | 30 s | see below |

**What `ranks` runs.** The probe calls the caller's own functions: `build_marker_panel(model, 31, 15)`,
`count_marker_observations`, `estimate_marker_depth` with the library's declared background (none is declared, so it is
estimated in both), `proposal_marker_parameters`, `build_proposal_emissions`, and the identical sort.
- k = 31 and s = 15 equal the runs' settings (`kmer_size`, `syncmer_s`, preset defaults).
- **The probe does not print its depth fit.** Equivalence is therefore gated behaviourally: every recorded trunk move's
  block, pair and `local_pair_rank` must be reproduced exactly by the probe's ranking. A failure stops interpretation.
- The blocks ranked are sv:5 and every block of a recorded move.

**What `analyse` reports:**
1. **Terminal contexts.**
   - Does the reported winner equal the trunk's last moved state (E1, E3)? Does it equal a converged tie-branch terminal
     (E6, with the trunk stated as well)?
   - The final neighbourhood scan of the trunk: its counts from the last round.
   - For each branch, the number of candidates scored in its terminal round.
2. **The rank gate.**
3. **Reachability of each frozen improving pair at sv:5:**
   - rank, emission, cutoff emission (rank 32) and gap;
   - whether it is in the list;
   - ties at its emission and at the cutoff;
   - the tie rule.
4. **Marker availability:**
   - called and improving allele lengths;
   - internal 31-mer count;
   - the number of sv:5 alleles shorter than 31 bp;
   - the number of distinct pair emissions;
   - the top tie class;
   - the composition of the 32-pair list.
5. **The generation trace,** per frozen improving alternative: excluded by the shortlist (never generated), or in the list
   and therefore generated at each scan.
6. **Sequence-level scored-pool membership.** M9-D's inventory has no within-block duplicates and cross-block
   equivalences only within sv:5..ctx:6. So only pool candidates equal outside sv:5..ctx:6 can spell the same pair; they
   are compared on the sv:5 + ctx:6 bytes.
7. **The context-aware reclassification of all M9-D contrasts.** An inconsistency is claimed only for a pool member
   scoring above the reported winner. A context-B carrier beating its called insertion is not, by itself, an
   inconsistency.

**Where it is recorded:** the correction is appended to M9-D's `RESULT.md`. Its tables are untouched.

**Budget:** 300 s, including 2 s of disclosed pre-freeze inspection. One process at a time.

**If existing records cannot resolve a stage, the result says so.** No generation history is inferred from aggregate
counters.
