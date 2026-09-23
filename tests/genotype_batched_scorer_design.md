# Batched block-local scorer: design and operation-count model

This is a design, not an implementation, and nothing here is wired into production. Inputs are the NA20509 30× complete search at `0b370e9`, a CPU profile of a 10× round-1 enumeration, and the recorded demand ledgers. The run record is `tests/genotype_within_round_retention_na20509.tsv`.

## What the search costs now (measured, 30×, single donor)

| component | work | seconds |
|---|---|---|
| startup (model, reads, marker fit) | — | 4.1 |
| round-1 seed placement + promotion | 2 full placements | 16.1 |
| per-round scoring, combination, reduce, promotion | 9 rounds | ~19 |
| window enumeration, first time (cold) | 344 distinct windows | 94.4 |
| window enumeration, repeated | 552 re-enumerations | 202.7 |
| **total** | 2,908 window requests → 896 enumerations | **336.9** |

- Cost per enumeration is 0.332 s. Enumeration is 88% of wall time.
- The delta path uses about one core of the eight declared: user time was 349 s against 336.9 s of wall time.
- The repeats come from windows the 384 MiB cache cannot hold:
  - 31 windows (326 MiB) were never admitted: 78.5 s;
  - 98 windows were evicted and then requested again: 124.1 s.

A complete cache is not an option. The 344 windows need 1.23 GiB, and the peak footprint is already 1.83 GiB at 384 MiB. Their representation compresses only 1.24 records per state (`fb19a43`).

## Where one enumeration's time goes (profile)

| function | share of samples |
|---|---|
| `gapped_emission` — mostly the exact feasibility prefilter | 71% |
| sorting candidate starts (`pair<uint,uint>` introsort/partition) | ~12% |
| `collect_form_placements` | ~6% |
| origin building, cache, and the rest | ~11% |

An earlier census at 30× found that 993 M of 1,189 M emission requests end as prefilter refusals, at about 1.6 µs each. Keying refusals costs more than evaluating them, so they cannot be cached.

## The batched block-local scorer as specified — what it would buy

The proposal:

1. group neighbours by distinct (block, local sequence context, end conditions, length shift);
2. enumerate each context once and score every requester from the same exact affected aggregate;
3. reuse the incumbent's unaffected contribution;
4. build positional origins only for the promoted winner.

Measured against that:

- **Distinct contexts ≈ distinct windows.** A window's key is already its local context plus end flags and replacement interval. Within a round, `0b370e9` already enumerates each distinct window once. That is 265–324 distinct per round against 322–325 requests, so grouping has nothing left to take.
- **Sharing flanks across alleles is worth ≤ 4.5%.** Read starts overlapping the replaced allele are a median 95% of a window's placements, and they carry 95.5% of its enumeration time. C4 alleles are long relative to the halo.
- **Reusing the incumbent's unaffected contribution and building the census only for the winner are already in place** (stage 6, `fb19a43`).

Proposed work = distinct contexts × enumeration + candidate assembly. With unbounded memory that is 344 × 0.332 s + about 39 s fixed ≈ **133 s**. Under the real 384 MiB bound it is ≈ **320 s**. **The batched scorer alone cannot reach 60–90 s.** Its operation count is what the current code already executes.

## What could reach 60–90 s (exact levers, not yet measured)

To finish in 60–90 s with today's 896 enumerations, each must cost 0.023–0.057 s, 6–14× less than now. If only the 344 distinct windows are enumerated, 0.06–0.15 s is enough, 2.2–5.4× less. Three exact levers target the measured costs:

1. **An exact pigeonhole pre-screen before the Myers prefilter.** A read split into d disjoint pieces of at least 12 bp, placed within band k, keeps at least d − k pieces intact. So a candidate start needs at least that many exact piece hits on diagonals consistent within ±(band + max gap). This is a necessary condition: it only discards placements the prefilter would refuse. It targets the 71%.
   - **Unmeasured:** the fraction of candidate starts with fewer consistent pieces than the bound.
2. **Deterministic parallel enumeration within a window, across fragments**, on the eight threads already declared, with the ordered reduction the scorer already uses. It keeps one window scratch (145 MiB at 30×) instead of one per thread, so it does not multiply memory the way parallelising across windows would.
   - **Unmeasured:** parallel efficiency.
3. **Bucketed candidate starts instead of the comparison sort** (~12%).

With a 2.5× gain on the emission share, 5× parallel efficiency, and the fixed part halved by the same parallelism, the model gives

`297 s × (0.29 + 0.71 / 2.5) / 5 + ~20 s ≈ 55 s`

(`0.29` is the non-emission share of enumeration time), and 75 s at 3× parallel efficiency. That is plausible, but it rests on two unmeasured multipliers.

**Next step:** measure both multipliers before implementing either:

- a census of candidate starts by count of consistent disjoint pieces, over one 30× round;
- the parallel efficiency of a single window's enumeration at 1, 2, 4 and 8 threads.

Every lever must preserve exact scoring:

- homologue assignment;
- StateKey multiplicity;
- gapped emissions;
- exposure;
- locus-end normalisation;
- the deterministic reduction order.

The adaptive production caller (panel → grouped-prefix mosaic → exact ranking → mandatory neighbour test → accept the best improving repair → repeat until no neighbour improves, failing closed on a tie or refusal) is the consumer. None of this changes that loop.
