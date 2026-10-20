# M9-G result: default-off junction-aware proposals integrated; off is identical to the archived caller, and on improves locus excess in all three affected development cells with both correct controls unchanged

**Run.** Branch `m9-g-junction-integration`.
- **Commits:** implementation and freeze `262ed6f`; Amendments 1–3 `878e514`, `085a502`, `7b215db`.
- **Results:** `tests/results/m9g_junction_integration_20261020/`.

**Execution: 814.1 of 2,400 s,** including 45 s of disclosed pre-freeze development. No unit was refused.

**Recorded failures** (each my own defect, each fixed by a declared amendment before any further execution):
- **A1:** the runner failed at import; nothing had run.
- **A2:** `gate.refusal` was mis-designed. The v2 schema's own refusal fired first, so the gate was rerun with the legacy
  schema.
- **A3:** `gate.offcheck` treated `dosage_status.tsv` timing rows as scientific. They were the only differences in 84
  checks.

## Implementation (default off)

**Setting:** profile key `repair_junction_proposals off|experimental_v1`.
- **Activation:** only at blocks whose existing pair emissions are all exactly equal and that have more pairs than 32.
- **What it replaces:** that block's index-ordered local list, with the JP1.1 shortlist for the current incumbent.
- **Unchanged:** 32 pairs, both assignments, the incumbent and the objective.

**Helper:** `junction_proposals.{hpp,cpp}`, with a unit test against an independent oracle.

**Search:** one accessor, `CoordinateRepairSearch::local_for`, with two consumers (the neighbourhood and the swap check).

**Caller:** about 40 lines:
- the key and its refusal (`repair_rounds 0`, diagnostic tie branching, the refinement replay);
- one sample-scoped read table;
- accumulation over the trunk, legs and tie branches;
- status rows recorded only when on.

**Caches:**
- **Read table:** one per run, `shared_ptr` to a `const` table, shared read-only by that run's searches, never persisted.
  It is 3.1–3.4 MB and takes 25–46 ms to build.
- **Proposal memo:** one per search, keyed by block plus the incumbent's unordered pair of whole label paths, the
  complete dependency. Bounded at 16.
- The persisted locus index is unchanged.

**JP1.1 = JP1 + a content-based residual tie-break.** Inside a tie class pairs are ordered by their sorted allele
sequences, not by allele indices. **Shortlist membership versus frozen JP1: 0 changes** in all eight M9-F contexts
(`gate.parity`).

**The M9-F empirical rate check is not integrated.** The run's own depth fit is used, and the status says so.

## Gates (all passed before any real caller run)

| gate | result |
|---|---|
| `gate.unit` | 24/24: F1–F9 against the whole-haplotype oracle; same flanks with changed elsewhere multiplicity refreshes the ranking and the cache key misses; a neighbouring-block change refreshes; a tie crossing rank 32 selects the same sequence-content set under 5 permutations and 4 concurrent invocations; eligibility; bounded eviction |
| `gate.sequence_report`, `gate.e2e_sequence_first` | pass (relinked) |
| `gate.refusal.a2` | refused, `repair_junction_proposals experimental_v1 needs the adaptive repair`, with no output written |
| `gate.parity` | C++ against the M9-F Python reference on all eight contexts × 6,786 pairs: max |Δscore| ≤ 1.4e-12; identical tie classes; top 32 identical in order to the JP1.1 reference |
| `gate.offcheck.a3` | **84/84**: the off arm (M9-C's own profile, the new binary) reproduces M9-C's sequence-first, locus, block, score, candidate, equivalence, dosage and settings outputs byte for byte (output prefix normalised; dosage timing rows excluded), the rounds and branch tables on their scientific columns, and the status on its scientific keys; no junction rows |

## Paired development runs (frozen before evaluation; every cell completed both arms)

| cell | off → on prediction | locus excess (floor) | blocks changed | sv:5 (off → on), legacy distance | Δtotal on − off (Δfragment / Δexposure / Δprior) | tie branches off → on | wall s / peak GiB off → on |
|---|---|---|---|---|---|---|---|
| **LOO d1** | RESOLVED → RESOLVED (different pair) | **27,279 → 27,277** (634) | sv:5 | 5/71 → 6/90: 4 → 2 | **+16.68** (+38.86 / 0.00 / −22.18) | 0 → 0 | 43.4 → 46.1 / 1.71 → 2.28 |
| **LOO d3** | RESOLVED → RESOLVED | **47 → 35** (21) | sv:5, ctx:6, ctx:8 | 45/107 → 34/93: 9 → 2 (ctx:6 4 → 3; ctx:8 36 → 32) | **+22.31** (+61.85 / +0.09 / −39.63) | 0 → 0 | 39.7 → 51.9 / 2.07 → 1.53 |
| **LOO d6** | RESOLVED → RESOLVED | **69 → 28** (42) | sv:5 | 0/21 → 18/21: 43 → 2 | **+17.61** (+17.03 / −0.51 / +1.10) | **2 → 0** | 72.8 → 58.1 / 2.00 → 1.88 |
| **LZO d1 (correct)** | RESOLVED, **identical** | 0 → 0 | none | 6/77, exact in both | 0 (identical scores) | 0 → 0 | 28.2 → 27.7 / 1.97 → 2.17 |
| **LZO d2 (correct)** | RESOLVED, **identical** | 0 → 0 | none | 76/81, exact in both | 0 | 0 → 0 | 26.2 → 26.2 / 2.12 → 2.04 |

**Trajectories (on arm):**
- **d1:** the off trunk plus one final sv:5 move to 6/90. The gain equals M9-F's offline one-step result, +16.68.
- **d3:** sv:5 → 34/93, after which ctx:2, ctx:6 and ctx:8 moves follow. The gain (+22.31) exceeds M9-F's offline one-step
  +19.41: the iterative caller continued from the new state.
- **d6:** sv:5 → 18/21 early in the trunk, followed by 7 further moves. It converges without the tie that sent the off arm
  into 2 branches. The final genotypes differ from off only at sv:5 (+17.61, equal to M9-F's offline one-step).
- **Exactness:** no changed sv:5 call is exact. Each ends 2 bp from the truth.

**Junction work (on arm):**
- **Eligible:** exactly 1 block per run (sv:5), found by the general rule.
- **Invocations:** 1–13 per run, and every one was followed by a cache hit (the swap check). Total proposal time is ≤ 0.14
  s per run.
- **Cutoff-tie crossings:** 0–7 invocations per run had a tie class crossing rank 32, resolved by content.

**v2 block confidence at changed blocks.** The sv:5 error-mass bounds of the new, still inexact calls are 0.20 (d1),
0.08 (d3) and 0.43 (d6), against 1e-38, 2e-24 and 9e-7 off. Better competitors now enter the scored pool, so the
pool-conditional masses no longer report near-certainty for inexact calls. Both LZO controls stay highly confident and
exact (5e-16 and 2e-18 error mass).

**Cost and memory.** Wall time changes by −15 to +12 s per cell, the increases coming from extra trunk rounds where
the new state keeps improving. Peak memory varies in both directions within 1.5–2.3 GiB, inside the 6 GiB guard. This is
consistent with the known allocator noise; no systematic increase is shown.

## Reading and recommendation

**What the evidence shows:**
- **Off is safe:** scientifically identical to the archived caller.
- **On, in all three development cells with a known missed improvement:** the caller reaches states the unchanged
  objective prefers (+16.7 to +22.3 nats) and reduces locus excess (−2, −12, −41 bp).
- **Controls:** both correct LZO controls are unchanged.

**Limits:**
- All five cells are inspected development cells, and these three were selected because they failed.
- The improved sv:5 calls remain 2 bp from truth.
- One locus, one eligible block, one coverage.
- Excess changes are development evidence, not an accuracy estimate.

**Recommendation: ready for a separate, predeclared evaluation on additional ACOT families** that have not been used in M9.
That evaluation should be specified before any data is used and should include:
- off/on pairing;
- the predeclared primary criterion: exact block agreement and locus excess at eligible blocks, plus control stability;
- refusal and resource accounting.

**Not started here.** No default change.
