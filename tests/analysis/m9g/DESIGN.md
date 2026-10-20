# M9-G: default-off caller integration of junction-aware proposals (JP1.1)

**Scope:** integrate M9-F's proposal mechanism into the existing adaptive search, then validate its mechanics and its
development-set behaviour. No additional ACOT evaluation families are used.

**Not done:** a wider shortlist, prior change, joint search, confidence remapping, new filter, default promotion, merge,
push or cleanup.

**Results:** `tests/results/m9g_junction_integration_20261020/`.

## Implementation contract

**Setting.** Profile key `repair_junction_proposals off|experimental_v1`, default **off**.
- It is recorded in the settings only when set.
- Status rows `repair_junction_proposals`, `junction_*` appear only when the experiment is on.

**Helper:** `include/panvar/junction_proposals.hpp` and `src/junction_proposals.cpp`, about 300 lines, tested
independently by `tests/junction_proposals_test.cpp`.

**Search integration** (`genotype_repair.*`): a single accessor, `CoordinateRepairSearch::local_for(b, first, second)`.
- **Its two consumers:** the round neighbourhood, and the homologue-swap invariance check, which asks for the same
  unordered incumbent.
- **Otherwise unchanged:** neighbours, scoring, the objective and the move rule.
- **Caller wiring** (`genotype_command.cpp`) is about 40 lines:
  - the key and its refusal;
  - one sample-scoped read table;
  - report accumulation over the trunk, legs and tie branches;
  - the status rows.

**Activation.** Only at blocks whose existing allele-pair proposal emissions are **all exactly equal** and that have
**more pairs than `repair_local_pairs`** (32).
- With 32 or fewer pairs, every pair is already generated.
- Nothing is hard-coded: no block id, locus, donor or allele.

**Rule JP1.1** (`junction_proposals.hpp`). It is M9-F's JP1, with one declared change to the residual tie-break.
- **Windows and markers:** the incumbent-context windows L + allele + R, with K − 1 = 30 bases each side, and all
  canonical 31-mers. Canonical means the minimum of the forward and reverse 2-bit codes, the caller's convention.
- **Expected count:** e(m) counts copies elsewhere in both incumbent haplotypes. Read counts y(m) cover both mates of
  every read. λ = rate · copies + background, using the run's own `MarkerDepthEstimate` (`rate_per_copy`,
  `background`), unchanged.
- **Score:** the Poisson kernel, as an exact delta summed in ascending inventory order. This is a deterministic,
  single-threaded reduction.
- **Pair ranking:** max over the two assignments.
- **Shortlist:** 32 pairs, with both assignments generated as before.
- **Ties:** |ΔS| ≤ 1e-9, chained.
- **The declared change — residual tie-break:** inside a tie class the order is by **content**, the sorted pair of
  allele sequences compared lexicographically, instead of JP1's allele indices. It is permutation-invariant.
  - Membership changes against frozen JP1 are measured and disclosed (`gate.parity`). The tolerance is not retuned.

**The incumbent.** It is always retained independently of its proposal rank. The repair scores the incumbent every round
as its current state, so a low proposal rank for the incumbent's own pair never removes it. A neighbour equal to the
incumbent is self-skipped, as before.

**Projection semantics.** Unchanged. JP1.1 only chooses which single-block neighbours are generated; the pooled decision,
its sequence groups and its projections are untouched. Identical sequences under different labels score and tie
identically (unit gate F8).

**Empirical rate check (M9-F) — integration treatment.** Not integrated. It was diagnostic evidence that the fitted
syncmer rate applies to all K-mers on these data; it is not a calibrated production test. The run's fitted rate and
background are used unchanged, and that is stated in the status (`junction_rate_source`).

## Cache ownership and invalidation

**Read table.** `JunctionReadKmers` holds the sorted canonical-K-mer counts of every read of the run.
- It is built once per run from the run's own fragments, only when the experiment is on.
- It is shared read-only by that run's searches through `RepairSearchOptions::junction_reads`, a `shared_ptr` to a
  `const` table.
- Its size is bounded by the read data and reported (`junction_read_table_bytes`).
- It is never persisted. The prepared or persisted locus index is unchanged, and nothing crosses samples: each run builds
  its own table.

**Proposal memo.** `JunctionProposalCache` is owned by each `CoordinateRepairSearch`: the trunk, every leg and every tie
branch has its own.
- **Key:** the **complete dependency** of an invocation, which is the block plus the incumbent's **unordered pair of whole
  label paths**. That pair fixes both flanks and every elsewhere multiplicity in both haplotypes; the reads, rate and
  background are fixed per search.
- **Consequences:** a changed elsewhere multiplicity under identical flanks gets a different key; so does a
  neighbouring-block move; a homologue swap gets the same key.
- **Bound:** 16 entries, oldest evicted.
- **Why this is safe:** resumed legs and branches cannot reuse incompatible state, because each search owns its memo and
  the key is complete.

**Not supported:** diagnostic tie branching (it enumerates the index-ordered lists) and the refinement replay. Both are
refused before any output, as is `repair_rounds 0`.

## Gates (all before real caller runs; any failure refuses every call)

| gate | checks |
|---|---|
| `build` | panvar, genotype_probe and the three test binaries |
| `gate.unit` (`junction_proposals_test`) | F1–F9 (short / empty allele, both-boundary K-mers, elsewhere repeats, reverse-complement, homologue swap, permutation, identical sequences, uninformative ties) against a whole-haplotype oracle; F4 same flanks with changed elsewhere multiplicity refreshes the ranking, and the cache key misses; C1 a neighbouring-block change refreshes, and the key changes; T1 a tie crossing rank 32 selects the same **sequence-content set** under 5 allele-index permutations and 4 concurrent invocations; C2 eligibility; C3 bounded eviction |
| `gate.sequence_report`, `gate.e2e_sequence_first` | existing v2 fixtures, relinked |
| `gate.refusal` | the CLI refuses experimental_v1 with `--repair-rounds 0`, with no output written |
| `parity.<cell>` + `gate.parity` | the probe `junction-proposals` (the caller's helper) on all eight M9-F contexts: every pair score within 1e-6 of the M9-F Python reference; the top-32 set equals the JP1.1 reference (Python scores and tie classes, content order); membership changes against frozen JP1 are reported |

**After the off arm,** `gate.offcheck` must pass before any on run. The off arm (M9-C's own profile, unchanged) must
reproduce the archived M9-C outputs:
- the sequence-first files, locus, block, genotype, score, candidate, equivalence, dosage and settings tables must be
  byte-identical after normalising only each arm's own output prefix;
- the rounds and tie-branch tables must match on their scientific columns (declared allowlists; timing, memory and
  cache-work columns excluded);
- the status must match on its scientific keys;
- there must be no junction rows.

## Paired development runs (frozen matrix)

**Cells:** LOO d1, d3, d6 and LZO d1, d2, the original M9-C inputs.

**Arms:**
- **off:** M9-C's own profile (validated combined route, sequence_first_v2, r = 0.01, joint search off);
- **on:** the same profile plus `repair_junction_proposals experimental_v1`.

**Common to both:** `--dosage-mode audit`, four threads, the 6 GiB and 900 s caller guards, and a 930-s process-group
watchdog. One caller process at a time.

**Order:** the off arm for all five cells, then `gate.offcheck`, then the on arm, then `freeze`. The predictions of both
arms are frozen before any evaluation; then 10 evaluations, then `analyse`.

**Analysis** (`m9g.py analyse`):
- prediction class, sequence status, ties and provisional or refused outcomes;
- whether the sequence pair is the same;
- every block whose v2 genotype changed (not just sv:5), with v2 exactness against the truth-closest alleles and the
  legacy distances;
- locus distance and excess over floor;
- winner fragment / exposure / prior / total;
- trunk moves and tie branches;
- junction eligibility, invocations, cache hits, cutoff-tie crossings and seconds;
- the read table's bytes and seconds;
- wall time and peak memory;
- any deterioration in the LZO controls.

M9-F's offline one-step outcomes are **not** assumed to equal the complete iterative caller outcomes.

## Budget

**Total: 2,400 s,** including 45 s of disclosed pre-freeze development (one compile at 34 s plus 1.5 s, and one unit-test
run).

**Admission:** spent + the unit's timeout + the declared timeouts of every pending freeze/eval/analyse unit (for builds,
gates and calls) ≤ 2,400. Refused units are recorded NOT RUN.
