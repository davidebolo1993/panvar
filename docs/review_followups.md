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

- **Resolve overlapping retained snarls before assigning final bubble IDs.** On ANKRD36C, `bubble`
  emits one snarl whose interior contains the interiors of all ten other retained snarls; Panphorte then
  sees the same 5,616 bp repeat array at two scales and correctly refuses the CSV rather than rewriting
  overlapping spans. Apply a deterministic conflict-resolution pass after filtering/merging. The chosen
  project policy is to retain the enclosing/larger snarl and drop the overlapping smaller snarls; report
  the dropped candidates and their conflict so the loss of finer site resolution is visible. Add a
  synthetic nested-snarl fixture plus an ANKRD36C regression asserting pairwise-disjoint emitted
  interiors and successful Panphorte preflight. Because the larger site can combine otherwise distinct
  subevents, also compare its downstream `call` records with the current smaller-site calls before
  treating this policy as biologically neutral.

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
