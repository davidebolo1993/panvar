# Mosaic genotyper: state audit and plan for phases 4-8

Companion to `genotype_mosaic.md`, which remains the authoritative contract. This document records
what is actually built, what its own gates still owe, and the work order to reach the Phase 8
go/no-go.

## 1. Verified state

Commit `d842deb` plus an uncommitted working tree. Build: Release, `-DPANVAR_BUILD_TESTING=ON`.
Both genotype test binaries pass.

| Phase | Component | Code | Gates |
|---|---|---|---|
| 1 model foundation | `genotype_model`, `ls_hmm` | committed `fefd4f9` | pass |
| 2 immutable panel/block | `genotype_model` | committed `715e8ed` | pass |
| 3 strict spelling | `mosaic_spelling` | committed `d842deb` | pass on hand fixture; **graph checkpoint owed** |
| 4 marker proposal | `marker_model` (462 lines) | uncommitted | pass on hand fixture; **gate 3 owed** |
| 5 deterministic top-K | `mosaic_beam` (393 lines) | uncommitted | pass against exhaustive oracle |
| 6 fragment likelihood | — | **absent** | — |
| 7 command integration | — | **absent** | — |

`panvar genotype` is not registered in `main.cpp`, as the contract requires at this stage.

The k-way merge over already-sorted allele lists is present in `build_marker_panel`
(`src/marker_model.cpp:198`), replacing the second N log N sort.

`tests/genotype_mosaic_proposal.cpp` is substantive, not decorative: it carries an independent
exhaustive enumerator over ordered template-pair histories and compares the beam against it at
r=0, r=0.25/0.6, and r=1, asserting best score, retained log-sum, emission and transition
components, and history counts separately. Non-vacuity is asserted explicitly
(`summed_histories_are_live`), and chunk-size variation stands in for scheduling nondeterminism.

## 2. The C4 block structure, measured

Panvar is an SV caller, so `panvar bubble`'s default `--min-variant-bp 50` is the correct block
alphabet: it keeps bubbles whose interior span is SV-sized. On `tests/real_data/c4.gfa.gz` that is
**5 SV bubbles** out of 1005 snarls in the complete decomposition.

The other 1000 snarls are all polymorphic, and every one of them lies in one of the six spans
*between* the SV bubbles — none is nested inside one:

```text
gap  0: 101 snarls    gap  6: 313 snarls
gap  2: 132 snarls    gap  8: 100 snarls
gap  4:  10 snarls    gap 10: 344 snarls
```

So at C4 no inter-block span is invariant, and `invariant_segments` — one shared string per gap —
cannot respell 131 templates byte-for-byte. Phase 3's round-trip gate forces the adapter to **tile**
the locus instead: the six context spans become blocks in their own right, alternating with the five
SV blocks. That yields **11 blocks**, which is the figure the contract's performance paragraph
assumes.

Measured per-block census after canonicalising orientation:

| Block | Kind | Alleles | Min bp | Max bp |
|---|---|---|---|---|
| 0 | context | 68 | 12291 | 12313 |
| 1 | SV | 68 | 277 | 7947 |
| 2 | context | 100 | 20784 | 20808 |
| 3 | SV | 16 | 203 | 506 |
| 4 | context | 10 | 274 | 324 |
| 5 | SV | 8 | 88 | 205 |
| 6 | context | 118 | 47756 | 47808 |
| 7 | SV | 119 | 26667 | 98508 |
| 8 | context | 62 | 13857 | 13862 |
| 9 | SV | 5 | 387 | 456 |
| 10 | context | 126 | 63407 | 63478 |

700 alleles over 11 blocks. The contract's "C4 contains blocks with 118--126 alleles" matches blocks
6, 7 and 10 exactly, which confirms the tiled decomposition is the one it was written against.

Two further structural facts the adapter must handle, both already anticipated by the contract:

- **0 partial frames.** All 131 panel paths traverse every one of the ten SV boundary nodes exactly
  once. No template is refused.
- **60 of 131 templates are antiparallel** to the reference chain. `mosaic_spelling`'s reverse-frame
  handling is load-bearing on real data, not a defensive branch.

130 of the 131 templates spell a distinct whole-locus sequence; the locus is 194,089--265,891 bp.

## 3. Defect found and fixed

`run_beam` stored each retained layer by moving `best`, which after `flush()`'s `best.swap(pool)`
owns the expansion pool's buffer — capacity `expansion_chunk_size`, default 262144 nodes. Layers are
retained for every block, so one chunk-sized allocation stayed alive per block. Storing exactly the
retained nodes instead, measured at C4's shape (131 templates, 11 blocks, r=0.01):

| K | Extensions | RSS before | RSS after | Wall |
|---|---|---|---|---|
| 64 | 11,000,201 | 121 MB | 26 MB | 0.08 s |
| 256 | 43,949,321 | 121 MB | 26 MB | 0.20 s |

Both test binaries remain green after the fix.

Two things worth recording. First, the contract's estimate — "131 squared ordered template states,
K=256, and ten transitions imply roughly 44 million extensions" — is exact: the measured figure is
43,949,321. The proposal stage is not a performance risk, and the stay/switch factorization the
contract permits is not needed for C4.

Second, the waste is proportional to block count, so the fix matters more elsewhere than here. On a
locus with a thousand blocks the same code reached 4.58 GB at K=64, over the Phase 8 memory gate on
its own. At C4's eleven blocks it is 121 MB, which no gate would have caught. It is fixed either way.

## 4. Work order

### Step 1 — graph-facing adapter (`genotype_locus_input`)

The single blocking item. Both Phase 3 and Phase 4 defer an acceptance gate to it, and nothing after
Phase 5 can touch real data without it.

It builds `StructuralLocusInput` from a GFA plus a reference path name, **tiling** the locus: SV
bubbles from `panvar bubble` at its production `--min-variant-bp`, interleaved with context blocks
covering the spans between them. A context span that is genuinely invariant across all templates
collapses to a single-allele block and costs nothing; at C4 all six vary. Block records carry whether
they are SV or context, because only SV blocks are reporting units. It populates the existing types;
it does not rebuild them.

Gates:
1. `verify_complete_panel_round_trip` over all 131 C4 templates, byte-for-byte, with the 60
   antiparallel templates canonicalised.
2. A node-split and a node-merge equivalent C4 graph produce an identical prepared-model
   fingerprint (which this step also introduces — the contract requires it in command output).
3. The census in section 2 is reproduced exactly: 11 blocks, 700 alleles, the per-block allele
   counts above, 0 partial frames, 60 reverse-oriented templates.
4. A template that does not traverse every block reports `Partial` and is refused, proven on a
   synthetic truncated path since C4 has no such case.

### Step 2 — Phase 4 gate 3 on real C4

Run `build_marker_panel` over the real panel and report `scored_alleles == total_alleles` for every
block — 700 alleles, no cap.

This step also measures `dropped_multiblock_markers`. Dropping a marker that varies in two blocks
from both blocks is a correctness-preserving choice the implementation made that the contract does
not describe, and section 2 makes its risk concrete: context block 6 (47.8 kb) and context block 10
(63.5 kb) flank the C4A/C4B modules and are largely near-duplicate sequence, as are SV block 7's
longer alleles. Markers shared between them are dropped from all of them. The drop rate is reported
before the marker model is frozen.

### Step 3 — Phase 6 candidate-neutral fragment likelihood

The largest step, the dominant runtime cost, and the one the previous attempt failed on. Four
commits:

- **3a** `LibraryModel` extension: error, insert, background, outlier and fragment-start intensity
  with provenance, plus validation. Read lengths stay per-fragment observations.
- **3b** Exhaustive substitution-only oracle, test-only. Independently re-derived, not lifted from
  `archive/genotype-factor-2026-09-11:src/genotype_fragments.cpp` — that file is 9514 lines of
  capped, truncated, top-K approximations the contract explicitly bars from being a production base.
  Its recorded traps become permanent regression fixtures: signed binned starts flipping strand;
  `FragmentState` ordering omitting edit counts; a run-wide read length raising a short fragment's
  insert floor; several seeds rediscovering one placement; distinct repeat copies with identical
  statistics.
- **3c** Production scorer: inverted pigeonhole piece index over candidates (not a per-candidate
  position hash), unique-haplotype interning so 2K haploid mosaics of ~230 kb are scored once, exact
  log-multiplicity state compression with separate physical-origin / statistical-class / mass
  counters, analytic exposure, and diploid combination by column pair.
- **3d** Certification: lower/upper score intervals, contender-only refinement, a winner only when
  one lower bound beats every competing upper bound, otherwise an equivalence set or
  `RESCORE_INCOMPLETE`.

All 14 named gates, the four prior mutations that must fail, and the ten runtime counters land with
3d. Compression is accepted only after it reproduces the oracle's mass exactly.

### Step 4 — Phase 7 command integration

`panvar genotype` with the eight declared options and no others. Orchestration only: load, count
markers, propose, spell and intern, recompute the Li-Stephens prior independently of the beam,
score, certify, project. Five outputs, five statuses, model fingerprint in the output. Per-block
output covers SV blocks; context blocks are modelled, not reported as calls.

### Step 5 — Phase 8 one-donor C4 go/no-go

K=64 under 90 s and 4 GB with per-stage wall time, then the K = 8..256 ladder with truth-in-candidate
rank, per-block true edlib distance, score gap, and winner stability. Proposal is measured at 0.20 s
and 26 MB, so effectively the whole budget belongs to Step 3. No cohort work starts before this
passes.

The beam stay/switch factorization is **not** on this path. The contract permits it; C4 does not need
it. Revisit only if a later locus has a block count in the hundreds.

## 5. Reviewer questions

1. **Tiling.** Confirm that promoting the six polymorphic context spans to blocks is the intended
   reading of "ordered blocks and their invariant contexts". At C4 no gap is invariant, so the
   alternative — one shared context string per gap — cannot pass Phase 3's byte-for-byte gate, and
   would inject systematic false bases into every spelled candidate that the fragment likelihood
   would then have to absorb. Related: confirm context blocks are modelled but excluded from
   per-block genotype output, since they are not SV calls.
2. **Long blocks and switch resolution.** SV block 7 reaches 98.5 kb and context block 10 reaches
   63.5 kb. Inside such a block the mosaic cannot recombine at all, so the caller is effectively
   choosing a whole template over those spans. Is that the intended resolution for C4, or should
   long blocks be sub-divided at a declared coordinate? This changes the legal switch set and so is
   a model change needing its own gate.
3. **`dropped_multiblock_markers`.** Is dropping a marker that varies in two blocks from both blocks
   the intended semantics? It is safe against double-counting but potentially lossy exactly across
   the C4A/C4B duplication. Decision requested after Step 2 reports the magnitude.
4. **Oracle provenance.** Confirm that re-deriving the Phase 6 oracle independently, rather than
   porting the archived exhaustive scorer, is the intended reading of "the archived exhaustive
   substitution-only fragment model is the oracle".
