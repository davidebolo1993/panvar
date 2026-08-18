# Deferred module review follow-ups

This is the project-level ledger for issues deliberately left until every module has completed its
review pass. It is not a live scratchpad: add or revise a module only when its pass is considered closed
and the project owner explicitly asks to **record** it.

Items here are non-blocking unless marked otherwise. Completed work should be removed rather than kept
as history; Git already provides that history.

## Associate

### Missing capability

- **Rare-variant aggregate association.** The implemented scope is common, single-feature association.
  Add and calibrate gene/bubble-level burden, collapsing, and SKAT-style tests. Firth and SPA improve a
  rare feature tested on its own; they do not provide rare-variant aggregation.

### Known statistical limits

- **Rare binary far-tail calibration.** Even after SPA, the measured type-I error is about 1.7 times
  nominal at `p < 0.001` for the rarest features. These results remain exploratory. Revisit with an exact
  or better-calibrated small-count method when rare-variant support is designed.
- **Rare binary mid-distribution distortion.** Under the covariate-preserving parametric null, the rarest
  features have a non-uniform middle of the p-value distribution even though error at 0.05 and 0.01 is
  near nominal. Preserve this case in future calibration work rather than relying only on phenotype
  permutation.
- **Feature-tier `Meff` is biological, not statistical.** Counting distinct bubbles is a useful grouping,
  but not a formal effective-number-of-tests estimator. Raw Bonferroni and BH remain the defensible
  corrections unless a phenotype-blind genotype-correlation estimator is implemented for this tier.

### Validation and cleanup debt

- **Make the GEMMA comparison assert absolute agreement.** `tests/gwas/validate_gemma.sh` currently
  reports correlations of beta and `-log10(p)` and can skip successfully when GEMMA or inputs are absent.
  Add per-feature absolute/relative tolerances for beta, standard error, and p-value for both linear and
  LMM models, with a pinned runnable environment or committed reference output. Correlation alone can
  hide a systematic p-value difference; until this exists, do not describe the LMM as independently
  validated to numerical tolerances.
- **Remove stale SPA wording in the implementation.** The warning comment near the association summary
  still says SPA is not implemented even though it is; keep the warning about residual anti-conservatism,
  but make the source comment match the actual method.

## Bubble

### Retained-site disjointness

- **DONE, with the policy reversed on measurement.** `bubble` now runs a deterministic conflict pass
  before assigning IDs, so the emitted interiors are pairwise disjoint and Panphorte's preflight passes.
  ANKRD36C previously failed the pipeline outright at `panphorte: bubbles 1 and 2 both claim interior
  node 7467`.

  The recorded policy was to retain the enclosing snarl and drop the smaller ones. The comparison this
  entry asked for was run, and it does not support that:

  | retained | sites | `call` records | reconstruction gap closed | worse than baseline |
  |---|---|---|---|---|
  | enclosing | 1 | 21 | 3.98% | 279 |
  | **smaller** | **10** | **21** | **89.54%** | **245** |

  Same record count either way; one large site simply describes a haplotype far worse than ten
  well-localised ones. The allele VCF is lossless under both (0 bp residual), so this affects only the
  interpreted output. `smaller` is now the behaviour. Only ANKRD36C is affected — the other five loci
  emit no overlapping sites and are byte-identical.

  Still owed: a synthetic nested-snarl fixture, and an ANKRD36C regression asserting pairwise-disjoint
  interiors plus a successful Panphorte preflight.

### Bounded traversal limitations

- **Remove or expose the 64-endpoint traversal cap.** The repeated-boundary interval search examines only
  the first 64 later endpoint occurrences for each start. A focused 70-repeat path therefore reports an
  interior span of 64 rather than 70. Replace the quadratic search with an exact bounded-time selection,
  or at minimum make truncation explicit in diagnostics; high-copy tandem loci must not be silently
  clipped.
- **Expose graph-interior traversal truncation.** Graph-derived interior discovery abandons its handle
  traversal after `2^20` visited handles and silently falls back to the path-derived interior. Report this
  condition and either reject the candidate or label the result incomplete, so the stated graph-derived
  interior contract does not quietly become panel-derived on a very large or leaky pair.

### Legacy-mode contract

- **Decide whether reference-free `--snarls-in` is diagnostic-only.** When a reference is supplied, the
  external door now validates it and resolves the same exact/unique-case-insensitive alias as internal
  mode. Without one, it still exits successfully after warning: imported pairs have no reliable reference
  order, boundary orientations default, reference-coordinate merging is skipped, and no sorted graph is
  emitted. Either require a reference for downstream-ready output, preserve and consume the oriented vg
  handles, or state prominently that this form is diagnostic-only; the module page currently says a
  reference is simply “not needed with `--snarls-in`”.

### Transactional and documentation cleanup

- **Reject colliding output destinations.** Input/output aliasing is checked, but two outputs may still
  name the same file. For example, identical `--bubbles-csv` and `--bandage-csv` destinations succeed and
  leave only the Bandage CSV. Preflight pairwise-distinct final paths before staging them.
- **Strengthen the multi-file commit contract.** The output family is staged, but destinations are
  installed sequentially and a later commit failure can leave a partially updated family. Add rollback
  or a manifest/disposition protocol if these files need true family-level atomicity.
- **Bring the public description up to the final contract.** Update the module summary and primary CSV
  table for graph-derived interiors, boundary orientations and allele-support columns; add
  `--min-alt-support` to the key-option table; document the `--snarls-in` limitation chosen above; and
  include `--emit-snarls-jsonl` in both the public output list and the command's `wrote:` summary.

## Call

### Missing cross-module capability

- **Consume Panphorte's REP provenance to aggregate phase-rotated REP nodes.** When one tandem array is
  folded at several rotations, Panphorte emits a REP node per phase and `call` counts each separately,
  so a site's total copy number is split across records. The aggregation should join on
  `output_rep_node` (not `created_rep_node`, which is pre-sort), determine the current bubble from node
  membership rather than trusting `input_bubble_id` (re-snarling reassigns ids), group by current site
  **plus** `canonical_motif` (never motif alone, or the same motif at two sites would merge), and sum
  `traversal_count × copy_quantum` across the phase nodes for every traverser including reference-like
  and zero-copy samples. Emit the grouped node set so the aggregation is auditable.

  **This fixes no current call.** Across all six reference loci the provenance table holds one row
  (LPA) and no site carries more than one REP node, so there is nothing to aggregate; any test would be
  synthetic. Recorded as a capability gap rather than a correctness defect, and worth building when a
  locus actually produces phase-rotated REPs. Tests to write with it: two rotations at one site
  aggregate; the same motif at two sites does not; missing or duplicate provenance rows fail loudly;
  `copy_quantum` affects the dosage; and CNBP plus the allele VCF are unchanged by aggregation.

### Accepted limitations

- **`PEAK` has no external truth validation, and cannot get one from the current data.** Every live
  `PEAK` record describes an 11–92 bp micro-module while the available truth counts whole gene copies,
  so the comparison is definitionally invalid — it produces numbers that look like accuracy and are
  not. The route is tested on a synthetic fixture with a known micro-module copy number and keeps
  `CN_CONFIDENCE=HEURISTIC`.
- **Repeated module boundaries remain an explicit span choice.** `CN`, `CNBP`,
  `CN_MODULE_REF_BP`, and the record's `POS`/`END` now use one occurrence-aware resolver. When a source
  or sink recurs, it deliberately selects the widest first-source-to-last-sink span and reports the
  number of affected traversers in `CN_SPAN_AMBIGUOUS`; that span can include content between separate
  visits and should not be read as uniquely resolved. Add a focused forward/reverse fixture that pins
  the warning count and proves all five quantities describe the same chosen span. The ordinary
  repeated-event fixtures exercise occurrence-aware calls, but currently do not produce a CN record
  with a repeated boundary.
- **Merged-record representative sequence.** A merged record still hands every carrier the
  representative's sequence; `MERGE_DIAMETER` now measures how far that reaches (0.0000 on a
  132-member LPA record) but does not fix it. Per-carrier sequence needs the allele VCF.

## Benchmark

### Corrected measurements — what the old numbers meant

- **Every QV figure recorded for this project before the ledger is the `graph` column**, which is an
  optimistic upper bound: a block is substituted when it shares a node with ANY call at the bubble, and
  no genotype is read at all. It is kept, unchanged and relabelled. The `called` reconstruction beside
  it is the strict form.
- **`over_threshold_bp` read 0 at all six reference loci and could not do otherwise.** The residual was
  split by contiguous alignment-run length, and edlib's co-optimal edit path distributes a gap over
  chance matches: a clean 60 bp deletion of non-repetitive sequence returns as fourteen runs of 1–10
  bases. Any claim resting on that split, and on the carrier truth flag built from it, is withdrawn.
- **"Recall is 100% everywhere; FN = 0 at all six loci" is withdrawn.** It followed from the truth flag
  almost never firing. On the walk-derived ledger, ACOT has 145 genotyping false negatives and CYP2D6
  has 3. One ACOT case is a 56,889 bp divergent block at bubble 7 whose haplotype is genotyped `0` at
  every one of the nine records there.
- **Carrier precision was understated, substantially.** ANKRD36C 3.9% → 69.2%, CYP2D6 55.7% → 100%,
  LPA 77.4% → 94.7%, ACOT 76.1% → 98.6%.
- **No above-threshold truth event is left untouched by some record** at any of the six loci: 0 `missed`.
  That is a discovery *ceiling*, not a proof of correct discovery — `called` means a record shares at
  least one node with the block, not that it covers or represents it. The firm half of the result is
  the negative: had a record been absent entirely, the event would have shown as `missed`.
- **Where the loss lives is representation.** The genotype residual is partitioned exactly, over one
  common observation set, into five consecutive steps. Across the six loci: `discovery_or_attribution`
  is **0 bp everywhere** (the retained calls restore exactly as much as a perfect in-scope caller
  would), `carrier_missed` is 0 at four loci and at most 15.6 kb, `false_positive_damage` is small and
  sometimes negative, and `representation` runs 43 kb to 8.57 Mb. Out-of-scope sub-threshold variation
  is 13–78 kb and is reported separately, since `call` was never asked to emit it.

  Three accounting faults were found in the first version of that partition and are fixed: comparative
  totals were taken over different populations (a partial VCF join put `carrier` ABOVE `called` and
  produced a loss of minus 20 bp); false-positive carrier damage was charged to representation; and
  below-threshold variation was charged to discovery. Each has a fixture.

### Accepted limitations

- **The allele VCF has no per-call attribution, and none is attempted.** `call --allele-vcf` writes one
  row per BUBBLE (`bubbleN_ALLELES`) spelling every allele; `variant_nodes.tsv` is one row per CALL.
  They share no id, so the `carrier` level and the per-call `loss_bp` terms do not exist in that mode
  and are reported `NA` / `not_applicable`. A bubble/allele-level attribution model would restore them
  — substitute a block when the haplotype's GT at that bubble names a non-reference allele — but it
  answers a coarser question than the region-VCF version, and giving the same column two meanings
  depending on the input file is the kind of silent semantic drift this campaign keeps paying for.
  Build it as separately named columns if it is wanted.
- **`genotype` is not bounded by `carrier`.** `graph ≥ called ≥ carrier` are nested ceilings; the
  genotype level applies every record a haplotype carries rather than only attributed eligible blocks,
  so it can beat the walk-based ceiling locally. The last two `loss_bp` terms are signed for this
  reason, and the sweep asserts only the nested part.
- **The reconstruction substitutes the whole divergent block, not the record's own interval.** A
  single overlapping node authorises the entire maximal block between shared anchors, so `called` and
  `carrier` do not prove the record delimits that sequence. Restricting the substitution to the
  record's `START_NODE`..`END_NODE` span would tighten both.
- **`called` means a record's node lies in the block, not that the record covers or reproduces it.**
  That is stronger than the previous bubble-wide node union but still weaker than comparing the
  record's reconstructed effect against the block's sequence. The `called` reconstruction inherits the
  same limit: it splices in the haplotype's TRUE block, so it measures what the retained records would
  achieve if each reproduced its block exactly. Closing it needs a per-record effect comparison — apply
  the record's own REF/ALT (or CNBP) within the block and align that against the block's true sequence.
- **No `filtered_other` class.** A size-eligible truth event removed by `--min-haplotypes`, an AF
  filter, a tangle guard or a resource limit is indistinguishable from one never found. Separating them
  needs a decision audit emitted by `call` — every raw truth-side event with its decision and reason —
  which is a `call` capability, not benchmark's.
- **`false_positive_damage` is bubble-grained.** An emitted call with no compatible truth event is
  attributed per (haplotype, bubble): if the haplotype has any eligible truth event at that bubble the
  observation counts as representation, even when a second record there is spurious. Record-level
  attribution needs a per-record effect comparison, the same thing the `called` limitation above needs.
- **DUP reconstruction is heuristic and labelled so.** Both `--dup-model` branches tile an inferred
  reference span, so the length is right and the sequence only approximately; the count is reported as
  `gt_records/heuristic`. **Measured cost at CYP2D6: genotype `variation_recovered` is 86.1% without
  `--cn` and 39.3% with it, at the same threshold on the same graph and bubbles.** The copy-number
  records are therefore where that locus's reconstruction loss lives, not the DEL/INS ones. Reconstruct
  a DUP from its allele sequence rather than by tiling, or accept and state the ceiling.
- **The threshold sweep is not automated.** `--min-sv-bp` reclassifies the ledger and gates the `called`
  reconstruction, but does not re-run discovery. A real experiment runs `call` at each threshold and
  `benchmark` on that run's output at the same threshold, over a fixed bubble set, and reports where
  monotonicity fails — those failures are where merging and rescue interact.
- **The allele VCF is a serialization ceiling.** It reaches 0 residual because it spells every allele
  out, which demonstrates lossless serialization, not call sensitivity.

## Describe

### Corrected measurements

- **A pure deletion produced no features at all.** All three step-building sites required an interior
  bubble node, so a path taking the direct source→sink edge was not counted as a traverser; the node
  that discriminates the deletion was then dropped as non-discriminative *because only its carriers had
  been observed*. Real effect, like-for-like: cyp2d6 +245 path observations (+38%), lpa +496 (+13%),
  c4 +26, and c4's 1340 "NA" cells were dropped traversals rather than missing data. Any dosage matrix
  built before this is missing its deletion carriers.
- **Both association gates re-run and pass.** The LPA KIV-2 signal is unmoved (bubble 7, beta −0.0193,
  p_bonf 6.1e-06) on a 13% larger substrate, and a matched null (same flags, only the code differing)
  gives type-I 0.05052→0.05082 at 0.05 and 0.00081→0.00081 at 0.001, lambda 0.9094→0.9116, with the
  worst per-feature KS p and the maxT threshold bit-identical.
- **The ~23% per-feature non-uniformity in that null is pre-existing**, present identically before and
  after, and belongs to `associate`'s recorded limits rather than to this pass.

### Accepted limitations

- **`--variant-flank-bp` has two granularities on purpose.** Bases for k-mers, whole nodes for graph
  dosage, because a node dosage is a property of the whole node. It therefore admits more nodes than
  bases. Splitting it into two flags would be clearer; one flag with a documented asymmetry was chosen
  over adding CLI surface.
- **Node and edge features share one id namespace.** They are told apart by the sidecar's `encoding`
  column, not by prefixed ids. A collision needs a node literally named like an oriented-edge key,
  while renaming feature ids would break every downstream join against `call`'s `variant_nodes.tsv`.
- **Pooling is locus-wide and stays that way.** Measured for the double-counting hazard: c4 has exactly
  one node claimed by two bubbles, cyp2d6 none, and zero emitted features localise to a shared node at
  either locus. The features spanning several bubbles (9 at c4, 31 at cyp2d6) are one k-mer sequence at
  genuinely different nodes. Revisit only if a locus produces features on a shared boundary.
- **Pooled carrier maps are held in memory across all bubbles.** `--max-wide-features` bounds the dense
  matrices, not this. Streaming bubble-qualified rows or an external sort/merge is the fix if a panel
  ever exceeds memory; no current locus does.
- **The missingness rules are latent at the reference loci.** Requiring every contributing bubble and
  every assigned haplotype to be observable is correct, but graph features never span more than one
  bubble at these loci and no incomplete observation occurs, so nothing measured today changes.

## Rebuild

### Deferred capability

- **R3: patch unmapped sequence into private nodes.** Rebuild currently rejects and rolls back a result
  that cannot recover every haplotype within the fidelity contract. That is safe. Measure rejection rates
  on real loci first; implement private-node patching only if losses occur often enough that it would
  rescue useful rebuilds without creating misleading graph structure.

### Validation and test debt

- **Preserve decision precision in the audit.** `fmt_exact()` uses ten significant digits, while per-path
  cover and identity values still use four decimals. Use `max_digits10` (or an explicitly justified
  precision) for both thresholds and the values compared with them so a decision can be reproduced from
  the audit alone.
- **Add direct unit fixtures for hard-to-reach branches.** Extract/test the chain comparator, the
  `identity_unavailable` rejection, and per-handle degree at the hub threshold. The current integration
  suite exercises the surrounding workflow but does not pin these branches independently.
- **Make the remaining integration fixtures exercise their claims.** The unchanged-structure fixture
  currently permits the degenerate early-return path, so it need not reach the `unchanged` label. Add a
  non-degenerate faithful rebuild with unchanged metrics. Also add an ambiguity-containing `k > 31`
  fixture; C4 proves strand canonicalization at `k=41` but contains no ambiguous bases.
- **Round-trip the accepted output through downstream consumers.** Parse an accepted rebuilt GFA again
  and run at least bubble detection on it in a focused fixture, so in-memory validation and emitted-file
  semantics cannot drift apart.

### Policy and transactional hardening

- **Decide whether “faithful but not untangled” should be accepted in production.** The current contract
  deliberately validates fidelity and only reports structural improvement. If real data shows accepted
  rebuilds that do not improve hubs/handle degree, consider making non-improvement a rollback condition
  rather than a warning.
- **Treat graph plus audit as one logical transaction.** Each file is staged individually, but the audit
  is committed before the graph. A graph commit failure can therefore leave an audit describing an
  output that was not installed. Either commit the graph first and clearly mark audit failure, or add a
  small manifest/transaction protocol that makes the pair's final disposition unambiguous.

## Refine

### Biological validation debt

- **Add the defining before/after `call` fixture.** The structural contracts are now well covered, but
  there is still no end-to-end example in which a reproducible graph-builder artifact calls as a split
  `INS`+`DEL` before refinement and as the intended clean event afterwards. Build this on a realistic
  pggb-style artifact and assert the exact records on both sides; that is the test that demonstrates the
  module's biological purpose rather than only its losslessness and graph validity.
- **Measure the stricter POA guards on every real locus.** `--max-poa-bp` now uses the longest allele and
  `--max-poa-work` independently bounds total distinct-sequence work, which is the honest resource
  contract. Record how many regions each guard skips on the six reference loci so safer accounting is
  not mistaken for unchanged refinement coverage.

### Documentation and portability cleanup

- **Repair the module page's options/output tables.** `docs/modules/refine.md` currently has output rows
  spliced into the `--max-poa-bp` option row, leaving both sections malformed. Also tighten the oriented
  traversal sentence in the algorithm page and document the final `--max-poa-work` and partial-path
  decision/report semantics in one coherent table.
- **Close and check staged streams before committing.** The report stream is flushed but remains open
  when `StagedOutputs::commit()` renames the family. This works on Unix but is not portable to platforms
  that refuse renaming open files, and `flush()` alone does not turn a late write failure into an error.
- **Declare `PANVAR_SLOW_TESTS` as a CMake option.** It currently works as an undeclared cache variable.
  Add an `option(...)` with help text and a default, and make clear that the default build must be
  reconfigured before the rebuild test can be selected explicitly.

## Inspect

### Clustering accuracy and scale

- **Use exact Jaccard when neither walk sketch is truncated.** The union-bottom-k estimator fixes the
  old length-ratio bias, but it deliberately subsamples even when both complete shingle multisets fit in
  the 512-element sketches. The worked example therefore estimates `J=0.5` where the exact value already
  available in memory is `0.667`; near a clustering threshold this can split or join a pair needlessly.
  Retain the exact shingle cardinality/truncation state and use exact multiset Jaccard for complete pairs,
  reserving the estimator for pairs where at least one sketch was actually truncated.
- **Bound or accelerate all-pairs clustering.** With `U` distinct walks, Inspect stores a `U x U` double
  distance matrix and compares every pair. This is reasonable for the reviewed panels but can become the
  dominant memory and runtime cost on a large cohort with many distinct haplotypes. Add a documented
  guard/diagnostic or an LSH/banding candidate stage before claiming cohort-scale clustering.
- **Make the representative independent of GFA path order.** Identical walks retain the first path name
  encountered, and that name becomes `representative_path`; reordering otherwise identical `P`/`W`
  records can therefore change the clusters TSV. Choose the representative by a stable rule such as the
  lexicographically smallest member after the structural medoid has been selected.

### Validation and transactional cleanup

- **Check more than one surviving crossing when matching the CSV to the graph.** Rejecting absent nodes
  and zero crossing paths closes the dangerous cases, but a stale graph with the same node IDs and only
  some of the original paths still succeeds. Compare the number of emitted crossings with the CSV's
  `path_support` (and, where the transformation promises to preserve it, the allele-support summary), or
  report an explicit mismatch warning.
- **Pin the two final second-pass branches directly.** The derived-output collision preflight and the
  present-nodes-but-no-crossing rejection are implemented and work, but `inspect_stats.sh` has no focused
  assertions for either case. Add them so the exact regressions fixed by the final pass cannot return.
- **Finish and commit the output family defensively.** Check every TSV stream after its final write/close;
  an open succeeding does not prove that later writes reached storage. Inspect also shares the sequential
  `StagedOutputs::commit()` limitation recorded for Bubble, so a late rename/copy failure can install only
  part of the per-bubble family unless the shared helper gains rollback or a manifest protocol.

## Panphorte

### Deferred capability and sensitivity

- **Recover safely when reused topology preserves a replaced route.** The per-edit acceptance check now
  rejects a normalization when the original local route remains walkable because all of its nodes and
  links are still needed elsewhere. That is the correct safe result. If this occurs often on real loci,
  clone the reused context before rewriting so the normalized site can be separated without breaking the
  paths that still need the original topology; until then, retain the hard refusal.
- **Make exact folding invariant to graph node segmentation.** Default exact mode still seeds tandem
  structure from repeated node-step runs. Two byte-identical copies split into different node boundaries
  can therefore be missed even though folding them would be lossless; only approximate mode reaches the
  base-level period fallback. Add an exact base-sequence fallback that proves byte identity before
  rewriting, with fixtures whose identical copies have deliberately different step partitions.
- **Reduce or quantify approximate-mode seed false negatives.** Copy alignment is attempted only after
  finding an exact shared 16-mer. A copy above `--min-similarity` can remain invisible when distributed
  substitutions disrupt every 16-mer. Consider multiple shorter/spaced seeds, minimizer candidates, or a
  bounded seedless fallback, and calibrate sensitivity and runtime before changing the default detector.

### Validation and cleanup debt

- **Register a real folding regression.** The default real-data smoke test uses C4, where its Panphorte
  rows do not normalize anything. Add an opt-in LPA regression that pins the KIV-2 site, all 466 paths,
  sequence/topology acceptance, delivered REP provenance, and representative copy numbers. Once Bubble
  emits a disjoint ANKRD36C set, add that locus as the nested-site regression as well.
- **Reconcile the partial-boundary guard with its reachable behavior.** Fragment emission currently makes
  `copies_declined_partial_boundary`, `paths_with_partial_boundary`, `partial_boundary`, and
  `--allow-partial-boundary` structurally unreachable on the reviewed detector paths. Either add a valid
  fixture proving the defensive branch can occur, or remove the dead option/status columns and update the
  stale header comments. Also rename the surviving-route diagnostic from “rewritten path(s)” to
  “replaced span(s)”, since the counter is now per edit.
- **Strengthen family-level commit semantics.** Panphorte closes and checks its staged streams, but the
  shared `StagedOutputs::commit()` still installs outputs sequentially. A late rename/copy failure can
  leave only part of the normalized GFA, report, provenance, bubbles, Bandage and optional GTF family;
  address this through the shared rollback or manifest protocol rather than module-specific handling.
