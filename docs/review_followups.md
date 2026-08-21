# Deferred module review follow-ups

This is the release ledger after the non-genotype review pass. Completed work is removed rather than
kept as history; corrected historical claims belong in
[`reports/module-review-summary.md`](reports/module-review-summary.md).

Tags have deliberately different meanings:

- **[RELEASE]** — close this, or explicitly narrow the release contract, before tagging a release.
- **[LIMIT]** — accepted behaviour. It is safe to ship only while the public documentation says it.
- **[TEST]** — behaviour is implemented and exercised on the reviewed loci, but the focused regression
  named here is still missing.
- **[LATER]** — capability or hardening work that is safe to defer.

`genotype` is excluded from the default build and this ledger. Its open work is in
[`reports/genotype-round2-verification.md`](reports/genotype-round2-verification.md); development builds
restore it with `-DPANVAR_ENABLE_EXPERIMENTAL_GENOTYPE=ON`.

## Release-wide items

### Output transaction correctness

- **[RELEASE] Repair the shared rollback state machine.** `cli::StagedOutputs::commit()` moves an old
  destination to a reserve, calls `commit_staged()`, and only then appends that destination to
  `installed`. If installation fails after the set-aside, the current reserve is absent from
  `installed`, so rollback never restores it. The existing fault injection fails *before* the
  set-aside and therefore cannot expose this path. Register the current destination and reserve as
  soon as the set-aside succeeds, track whether the new file was installed, and inject a failure after
  set-aside/before installation. Assert that old files are restored, new-only files are absent, and no
  `*-prev`/`*-tmp` files remain.

  `commit_staged()` also falls back from a sibling rename to copying directly over the final path and
  then removes the staged file even if the copy reported an error. That fallback is not atomic and can
  damage the destination it is supposed to protect. Because staging files are siblings of their
  destinations, either refuse an unexpected cross-filesystem rename or copy to another sibling and
  atomically rename that complete copy.

- **[RELEASE] Give `associate` the same output contract as the other commands.** It writes
  `.assoc.tsv` and `.summary.tsv` directly, has no output/input collision preflight, and does not
  explicitly close/check either stream. A late failure can leave a mixed family, and an output can
  alias a genotype, phenotype, covariate, annotation, or kinship input. Use `StagedOutputs`,
  `reject_output_collisions`, and checked close before the shared commit; cover replacement and failure
  rollback in `associate_stats.sh`.

- **[RELEASE] Finish Describe's independent directory transaction.** Its restore loop restores old
  owned entries but does not remove a newly installed entry that had no predecessor if a later install
  fails; restore errors are ignored. Track every installed destination, remove this run's new entries
  in reverse order, restore every reserve, and report any reserve that could not be restored. Add fault
  injection. Also reject an input whose canonical path is an owned output such as
  `bubble_<digits>` inside `--out-dir`, or the commit can consume it and then replace it.

### Packaging and platform contract

- **[RELEASE] Decide whether the install is CLI-only or a supported C++ package.** The installed CLI
  was smoke-tested successfully. The exported library package is not consumable as shipped:
  `panvarTargets.cmake` refers to non-exported `abpoa` and `minigraph` targets and to `ZLIB::ZLIB`
  without a package config that discovers dependencies. Experimental genotype headers are also
  installed in the default build although their implementations are omitted. The small, honest release
  fix is to install only the executable. If a C++ API is intended, add `panvarConfig.cmake`, dependency
  discovery/export, an intentional public-header set, and a clean external consumer compile/link test.

- **[RELEASE] Run a clean Linux/ELF build and test job, or state that this release is macOS-only.** The
  clean AppleClang release build and all 10 enabled tests pass; an AddressSanitizer/UndefinedBehaviorSanitizer
  build also passes those tests. Linux remains important because minigraph uses a separate ELF symbol
  localization path. CI should configure from an empty tree, build with genotype off, run the default
  suite, run the slow rebuild test at least in a scheduled/release job, install, and execute the
  installed binary.

### Small release cleanup

- **[LATER] Remove the clean-build warning in benchmark.** `tot_carrier_aln` is accumulated but never
  read (`src/benchmark_command.cpp`); either report the denominator it was intended to preserve or
  delete it. This is not a behavioural defect.

## Associate

### Accepted statistical and scale limits

- **[LIMIT] Rare binary far-tail calibration.** Even with SPA, measured type-I error for the rarest
  features is about 1.7 times nominal at `p < 0.001`. Treat those single-feature results as
  exploratory.
- **[LIMIT] Rare binary mid-distribution distortion.** The covariate-preserving parametric null remains
  non-uniform in the middle for the rarest features even though error near 0.05 and 0.01 is close to
  nominal. Future calibration must retain this null rather than relying only on permutation.
- **[LIMIT] Feature-tier `Meff` is biological, not statistical.** Counting bubbles is a grouping, not a
  formal effective-number-of-tests estimate. Raw Bonferroni and BH are the defensible corrections at
  this tier.
- **[LIMIT] Variant-tier LD/Meff is region-scale.** LD clumping computes all feature pairs and the Li-Ji
  estimator materializes an `n_tests × n_tests` correlation matrix followed by a full symmetric
  eigendecomposition: quadratic memory/work to form it and cubic eigensolve time, with no cap. Add a
  projected-memory/work guard and an auditable fallback (raw Bonferroni or block/clump Meff) before
  describing this as cohort/genome-scale.

### Deferred capability and validation

- **[LATER] Add rare-variant aggregate tests.** Burden, collapsing, and SKAT-style gene/bubble tests are
  not implemented. Firth and SPA improve a rare feature tested alone; they are not rare-variant
  aggregation.
- **[LATER] Validate LMM numerics against GEMMA with absolute tolerances.** The current comparison uses
  correlations and can skip when dependencies are absent. Pin a runnable reference and assert
  beta/SE/p-value absolute or relative agreement. Until then, do not claim independent numerical
  validation of the LMM.

## Bubble

### Accepted limits

- **[LIMIT] Size filters measure interior span, not allele/reference divergence.** A 1 bp substitution
  inside a 1 kb allele has a 1 kb span and can pass a 50 bp threshold.
- **[LIMIT] `--snarls-in` without `--reference-path` is diagnostic-only.** Imported boundaries are
  unordered, orientations are defaults, reference-coordinate merging is skipped, and no sorted graph
  is written. It must not feed coordinate-bearing `call` output.
- **[LIMIT] Only top-level snarls are retained.** Descending into children was measured and reverted:
  it did not help C4, exploded ACOT, and did not finish at LPA. Keeping the larger snarl is acceptable,
  provided downstream consumers continue to require disjoint selected interiors.

### Validation debt

- **[TEST] Add the real ANKRD36C disjointness regression.** The synthetic nested-snarl test is good, but
  an opt-in real fixture should assert pairwise-disjoint interiors and successful Panphorte preflight.

## Call

### Accepted limits

- **[LIMIT] `PEAK` remains a heuristic micro-module CN route.** Current truth counts whole genes and
  cannot validate the 11–92 bp modules on which PEAK fires. Keep `CN_CONFIDENCE=HEURISTIC`; the
  synthetic known-copy fixture is its present evidence.
- **[LIMIT] Repeated module boundaries require a policy choice.** CN, CNBP,
  `CN_MODULE_REF_BP`, POS, and END share the occurrence-aware resolver, which deliberately chooses the
  widest first-source-to-last-sink span and reports `CN_SPAN_AMBIGUOUS`. That may include sequence
  between separate visits and is not a uniquely resolved biological interval.
- **[LIMIT] A merged record uses one representative sequence for every carrier.** `MERGE_DIAMETER`
  exposes the spread but does not provide per-carrier sequence; use the allele VCF when exact carrier
  alleles matter.

### Deferred capability and validation

- **[LATER] Aggregate phase-rotated Panphorte REP nodes.** Join provenance on `output_rep_node`, map it
  to the current re-snarled site, group by current site plus `canonical_motif`, and sum
  `traversal_count × copy_quantum`. Test two phases at one site, the same motif at two sites, malformed
  provenance, zero/reference-like traversers, and invariance of CNBP/allele serialization. No reviewed
  locus currently has more than one REP node at a site, so this is a capability gap rather than a
  current miscall.
- **[TEST] Add a repeated-boundary MODULE_BP/CNBP fixture.** Exercise forward and reverse traversal,
  pin `CN_SPAN_AMBIGUOUS`, and prove CN, CNBP, `CN_MODULE_REF_BP`, POS, and END all describe the same
  chosen occurrence span. Existing repeated-node event tests do not create a CN record with a repeated
  source/sink.

## Benchmark

### Accepted interpretation limits

- **[LIMIT] The allele VCF result is a serialization ceiling.** Its 0 bp residual at all six loci proves
  that every bubble allele can be serialized and re-implanted, not that region-level discovery or
  genotyping is perfect.
- **[LIMIT] Allele-VCF mode has no per-call attribution.** It has one row per bubble while
  `variant_nodes.tsv` has one row per call; `carrier` and per-call loss terms therefore remain
  `NA/not_applicable`, not zero.
- **[LIMIT] `graph ≥ called ≥ carrier` are nested truth-assisted ceilings; `genotype` is not bounded by
  `carrier`.** Genotype applies every record assigned to a haplotype and can locally outperform or
  damage the truth-block ceiling.
- **[LIMIT] `called`/`carrier` substitute the haplotype's whole true divergent block.** A call-node
  overlap authorizes the block; it does not prove that the record delimits or reproduces it. A stricter
  score must apply each record's own REF/ALT/CNBP effect inside its actual interval and align that to
  truth.
- **[LIMIT] No `filtered_other` truth class exists.** An eligible event removed by AF/support/tangle or
  resource policy is indistinguishable from one never discovered until `call` emits a decision audit.
- **[LIMIT] False-positive damage is attributed at haplotype/bubble level.** A second spurious record in
  a bubble containing one valid truth event is not isolated without per-record effect matching.
- **[LIMIT] DUP genotype reconstruction is heuristic.** It tiles an inferred reference span. The
  measured CYP2D6 result (86.1% without CN versus 39.3% with CN) shows this is a material representation
  limit, not cosmetic reporting.

The threshold sweep is now the correct experiment: it reruns both `call` and `benchmark`, holds truth
events fixed, and verifies that lowering the threshold cannot increase `below_threshold`. Do not claim
that genotype reconstruction itself is monotone; at C4 additional overlapping/merged records make it
worse at lower thresholds, and the validator reports that as a call interaction.

## Describe

### Accepted limits

- **[LIMIT] `--variant-flank-bp` has base granularity for k-mers and whole-node granularity for graph
  dosage.** The asymmetry is intentional but must remain explicit.
- **[LIMIT] Node and edge features share one ID namespace.** A node ID containing the edge separator is
  refused; it is not repaired or annotated after collision.
- **[LIMIT] Pooling is locus-wide.** It has not produced a shared-boundary emitted feature on C4 or
  CYP2D6, but revisit if a larger locus does.
- **[LIMIT] Pooled carrier maps remain resident across all bubbles.** The wide-matrix cap does not bound
  this memory. Stream qualified rows or external-sort them if larger panels approach memory limits.
- **[LIMIT] Multi-bubble missingness is not exercised by current real loci.** The rule is implemented,
  but present features do not span several bubbles with incomplete observations.

### Validation debt

- **[TEST] Strengthen the deletion-positive-marker fixture.** It currently proves only that the k-mer
  matrix is non-empty. Pin the observed ALT-exclusive junction syncmers: carried by every bypass
  deletion haplotype and by no reference/full/truncated path. Also report ALT-exclusive junction-marker
  counts for real deletions; repetitive junctions cannot be guaranteed unique.

## Rebuild

### Accepted/deferred behaviour

- **[LATER] Patch unrecovered sequence into private nodes only if measurements justify it.** Current
  behaviour rejects and rolls back any unverifiable/lossy recovery, which is safe. Implement R3 only if
  real rejection rates show useful graphs are being lost.
- **[LATER] Decide whether faithful-but-not-untangled should remain accepted.** Fidelity is the acceptance
  gate and structural improvement is reported separately. Promote non-improvement to a rollback policy
  only if it occurs materially on real inputs.

### Validation debt

- **[TEST] Unit-test hard-to-reach decisions.** Extract/pin the chain comparator,
  `identity_unavailable`, and per-handle degree exactly at the hub threshold.
- **[TEST] Make integration fixtures non-degenerate.** Require a rebuild that really runs and reports
  unchanged metrics, and add an ambiguity-containing `k > 31` seed fixture. Current N handling and
  `k=41` strand canonicalization are tested separately.
- **[TEST] Round-trip an accepted rebuilt file through downstream parsing and bubble detection.** This
  is the remaining check that emitted GFA semantics cannot drift from in-memory acceptance.

Audit precision, graph-first/audit-second disposition, checked Refine streams, and declaration of
`PANVAR_SLOW_TESTS` are closed and intentionally absent from this list.

## Refine

### Biological validation debt

- **[TEST] Add the defining before/after-call fixture.** Construct a reproducible pggb-style artifact
  that calls as split INS+DEL before refinement and as the intended clean event afterwards; assert both
  exact record sets and unchanged haplotype spellings.
- **[TEST] Measure POA guard coverage on all six loci.** Record regions skipped separately by
  `--max-poa-bp`, `--max-poa-work`, and `--max-walks`, so safer accounting is not mistaken for unchanged
  biological coverage.

## Inspect

### Accepted scale limit

- **[LIMIT] Clustering is dense all-pairs.** It warns above 2,000 distinct walks and refuses above
  25,000; the reviewed maximum is 119. LSH/banding candidates are the route to cohort-scale clustering.

### Validation and hardening debt

- **[LATER] Compare surviving graph support with the CSV.** Present nodes and at least one crossing are
  required, but a stale graph with the same IDs and only some original paths can still pass. Compare or
  warn on emitted crossings versus CSV `path_support`.
- **[TEST] Pin derived-output collision and present-nodes/no-crossing directly.** Both branches exist but
  lack focused assertions in `inspect_stats.sh`.
- **[LATER] Explicitly close and check every TSV stream.** The writer's local streams are destroyed before
  family commit, but destructor close cannot report a late disk error.

## Panphorte

### Accepted/deferred behaviour

- **[LATER] Clone reused context only if hard refusals become common.** A normalization is safely refused
  when another site needs topology whose old local route would otherwise survive. Cloning is the
  recovery, but no reviewed locus currently justifies its complexity.
- **[LATER] Make exact folding invariant to node segmentation.** Exact mode still seeds repeated
  node-step runs. Byte-identical copies split at different node boundaries can be missed; an exact
  base-sequence fallback must prove spelling identity before rewriting.
- **[LATER] Quantify approximate seed false negatives.** Approximate candidates require an exact shared
  16-mer. Calibrate shorter/spaced/minimizer seeds or a bounded seedless fallback before changing the
  default detector.

### Validation and interface debt

- **[TEST] Add opt-in real folding regressions.** LPA should pin the KIV-2 site, 466 path spellings,
  topology acceptance, REP provenance, and representative CNs. ANKRD36C should pin the nested-site
  preflight once paired with Bubble's real disjointness test.
- **[TEST] Resolve the unreachable partial-boundary interface.** Either construct a valid detector path
  that reaches `copies_declined_partial_boundary`/`--allow-partial-boundary`, or remove the dead option,
  status columns, and stale comments. Rename the surviving-route diagnostic from “rewritten paths” to
  “replaced spans” if it still counts edits.
- **[TEST] Add a Panphorte-family rollback fixture after the shared transaction repair.** The present
  injected rollback test runs through Bubble only.

## Verification snapshot (2026-08-21)

- Fresh macOS AppleClang release configure/build with genotype disabled: pass.
- Default tests from that fresh tree: 10/10 enabled pass; slow rebuild test disabled by policy.
- macOS AddressSanitizer + UndefinedBehaviorSanitizer build: 10/10 enabled tests pass. Vendored C
  libraries are not fully instrumented by that CMake flag.
- Installed CLI smoke (`panvar --help` outside the build tree): pass.
- Slow rebuild suite: reported passing in the closing pass, but not rerun in this audit.
- Linux/ELF and installed C++ consumer: not validated; the latter is currently structurally incomplete
  as described above.
