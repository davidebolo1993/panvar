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

All five release items from the 2026-08-21 audit are closed. Recorded here only so the next reader
knows what was checked rather than assumed:

- Shared rollback registers a destination as soon as the SET-ASIDE succeeds, so the window between
  moving the old file aside and installing the new one is owned. `commit_staged()` no longer copies
  onto the destination; it copies to a second sibling and renames that complete copy into place.
  `PANVAR_TEST_FAIL_COMMIT_AFTER_SETASIDE` reaches the window that the previous fixture could not.
- `associate` uses `reject_output_collisions`, `StagedOutputs` and checked closes.
- Describe's directory transaction removes this run's new entries on rollback and reports a reserve it
  could not restore; an input under an owned output name is refused.
- The install is CLI-only (`bin/panvar`). No C++ package is advertised.
- CI configures from an empty tree on Linux and macOS, asserts genotype is absent from the default
  build, installs and runs the installed binary; a second job builds with genotype ON; the slow suite
  runs on a schedule and on demand.
- Benchmark's clean-build warning is gone: `tot_carrier_aln` was the denominator that turns
  `carrier_recon delta` into an identity, and is now reported as `carrier_recon aln_len`.

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
  eigendecomposition: quadratic memory to form it and cubic time to solve it. A guard is now in place —
  a matrix-size note above 5,000 tested variants and, above 20,000, a fall back to LD-clumping Meff or
  raw Bonferroni with `meff_method` recording which was used. Exceeding the cap is not an error because
  Meff only refines the correction. This is still not a genome-scale estimator.

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

- **[TEST] DONE** — `tests/real_regressions.sh` asserts pairwise-disjoint ANKRD36C interiors across 10
  sites and a successful Panphorte preflight on that set. Opt-in under `PANVAR_SLOW_TESTS`.

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
- **[TEST] DONE, and it corrected a documented invariant.** `call_stats.sh` now builds a reference
  visiting both module boundaries twice. `END − POS = CN_MODULE_REF_BP` was stated as the placement
  check for module records; it holds only while every boundary is visited once. POS/END are the widest
  span, CN_MODULE_REF_BP and CNBP sum node length × multiplicity over the module nodes, and with a
  repeated boundary the two differ by exactly the boundary visits between the outermost occurrences —
  measured 700 against 600. `CN_SPAN_AMBIGUOUS` is the flag that says the identity does not apply.
  `docs/algorithms/call.md` is corrected.

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

- **[TEST] DONE** — the fixture pins the ALT-exclusive junction syncmers (carried by all three bypass
  haplotypes, by no path keeping the interior) and reports how many there are: three here. Counting
  rather than thresholding, because a repetitive junction cannot be guaranteed unique.

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
- **[TEST] DONE** — `real_regressions.sh` rebuilds C4, re-reads the accepted GFA with `bubble`, and
  asserts it decomposes with pairwise-disjoint interiors.

Audit precision, graph-first/audit-second disposition, checked Refine streams, and declaration of
`PANVAR_SLOW_TESTS` are closed and intentionally absent from this list.

## Refine

### Biological validation debt

- **[TEST] Add the defining before/after-call fixture.** Construct a reproducible pggb-style artifact
  that calls as split INS+DEL before refinement and as the intended clean event afterwards; assert both
  exact record sets and unchanged haplotype spellings.
- **[TEST] DONE, reported rather than thresholded** — `real_regressions.sh` prints refine's decisions
  per locus (lpa 7 rebuilt/5 skipped, c4 3/2, gstm1 8/1, cyp2d6 5/4, acot 7/2, ankrd36c 5/5) and fails
  only if no locus rebuilds anything. The per-guard breakdown is still owed; today the report gives the
  decision, not which bound produced it.

## Inspect

### Accepted scale limit

- **[LIMIT] Clustering is dense all-pairs.** It warns above 2,000 distinct walks and refuses above
  25,000; the reviewed maximum is 119. LSH/banding candidates are the route to cohort-scale clustering.

### Validation and hardening debt

- **[LATER] DONE** — emitted crossings are compared against the CSV's `path_support` and a mismatch is
  named. Not an error, since a legitimate path-dropping transform produces one.
- **[TEST] DONE** — both are pinned, plus a stale-panel case: a graph carrying fewer paths than the CSV
  records now warns (2 crossings against `path_support=3`) and a matching panel stays silent.
- **[LATER] DONE** — all four inspect TSV streams are closed and checked before the family commit.

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
- **[TEST] DONE** — `panphorte_stats.sh` pins both failure windows on the panphorte family, with the
  graph and report byte-checked and no residue.

## Verification snapshot (2026-08-21)

- Fresh macOS AppleClang release configure/build with genotype disabled: pass.
- Default tests from that fresh tree: 10/10 enabled pass; slow rebuild test disabled by policy.
- macOS AddressSanitizer + UndefinedBehaviorSanitizer build: 10/10 enabled tests pass. Vendored C
  libraries are not fully instrumented by that CMake flag.
- Installed CLI smoke (`panvar --help` outside the build tree): pass.
- Slow rebuild suite: reported passing in the closing pass, but not rerun in this audit.
- Linux/ELF and installed C++ consumer: not validated; the latter is currently structurally incomplete
  as described above.
