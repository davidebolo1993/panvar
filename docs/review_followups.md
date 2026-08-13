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

### Bounded and legacy-mode limitations

- **Remove or expose the 64-endpoint traversal cap.** The repeated-boundary interval search examines only
  the first 64 later endpoint occurrences for each start. A focused 70-repeat path therefore reports an
  interior span of 64 rather than 70. Replace the quadratic search with an exact bounded-time selection,
  or at minimum make truncation explicit in diagnostics; high-copy tandem loci must not be silently
  clipped.
- **Expose graph-interior traversal truncation.** Graph-derived interior discovery abandons its handle
  traversal after `2^20` visited handles and silently falls back to the path-derived interior. Report this
  condition and either reject the candidate or label the result incomplete, so the stated graph-derived
  interior contract does not quietly become panel-derived on a very large or leaky pair.
- **Define the supported `--snarls-in` contract.** Without a resolvable reference, imported boundary pairs
  are emitted with no reliable reference order and default handle orientations; reference-coordinate
  merging is also skipped. Either require a reference (and resolve it by the same exact/case-insensitive
  rules as internal mode), preserve and consume the oriented vg handles, or clearly restrict this path to
  diagnostic use rather than downstream-ready Bubble output.

### Transactional and documentation cleanup

- **Reject colliding output destinations.** Input/output aliasing is checked, but two outputs may still
  name the same file. For example, identical `--bubbles-csv` and `--bandage-csv` destinations succeed and
  leave only the Bandage CSV. Preflight pairwise-distinct final paths before staging them.
- **Strengthen the multi-file commit contract.** The output family is staged, but destinations are
  installed sequentially and a later commit failure can leave a partially updated family. Add rollback
  or a manifest/disposition protocol if these files need true family-level atomicity.
- **Bring the public description up to the final contract.** Update the module summary and primary CSV
  table for graph-derived interiors, boundary orientations and allele-support columns; document the
  `--snarls-in` limitation chosen above, and include `--emit-snarls-jsonl` in the reported output list.

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
