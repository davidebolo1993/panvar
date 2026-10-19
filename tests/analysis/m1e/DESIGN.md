# M1-E: residual-sequence sensitivity and allocation robustness

Branch `m1-e-residual-sensitivity` from `1386f31`. This is diagnostic only: no production caller, model, prior, search or
output-contract change, and no truth input to a caller. Only `tests/genotype_probe.cpp` changes: `score-set --contrast`
adds compact fragment-delta summaries computed from the scorer's masses, with no origin recording.

Results: `/Users/davide.bolognini/panvar/tests/results/m1e_residual_sensitivity_20261019`.

## Scope carried forward

M1-D's 264 / 810 nats is the deficit of the selected oracle composition against truth. It is not a proven panel-wide
deficit; a qualification is appended to the M1-D README.

## Restoration set (`construct.py`; frozen before scoring)

These are ORACLE sensitivity interventions, not inferred corrections.

**Starting point.** The M1-D composition:
- P = prefix_u6, which equals both truth prefixes exactly;
- truth1: U1 = period_u13 (1 bp edit distance) and T1 = tail_u44 (23 bp);
- truth2: U2 = period_u35 (3 bp) and T2 = tail_u33 (2 bp).

These distances are edit distances, not counts of variant sites.

**Patch sets.** Each patch replaces one *whole* component with its truth bytes (from the M1-D truth partition by the
panel anchor). The sets are: none, U1, T1, U2, T2, all.

**Arrangements.** Each patched inventory is arranged both ways:
- `a22` (2+2): P+U1+T1 | P+U2+T2;
- `e13` (1+3, as in M1-D): P+T1 | P+U1+U2+T2.

**Controls.** `reported`, with its saved total as expectation, and `truthblock`. That makes 14 declared diplotypes;
`truthblock` equals `a22_all` and is interned once.

**Asserted:**
- inventory SHA-256 identities and the truth components;
- `a22_none` equals the M1-D composition;
- matched inventories and total lengths per patch set;
- each side equals reported left flank + block + reported right flank;
- `a22_all` equals the truth-block reference.

**Recorded:** lengths, SHA-256 and the 62-bp context of every join.

## Scoring

The explicit-sequence scorer (validated in M1-C) is used with unchanged library, band, exposure and options, one batch
per coverage (10x, 30x), and no origin recording.

`--contrast` summaries:
- each patch against none, per arrangement;
- e13 against a22 per patch set;
- the composition and truthblock against reported.

Prior support is reported separately: unsupported sequences get no prior and no total.

## Reporting (`summary.py`, `content.py`)

- Fragment term, exposure and likelihood for each candidate.
- Each single-component gain, the all-component gain, and the interaction (all minus the sum of singles). Additivity is
  not assumed.
- The allocation contrast (e13 − a22) before and after correction, with bit-equality.
- Pooled block content against truth (k31, k151), and the content difference between the matched arrangements (which
  is what their different joins create).

## Budget

600 s cumulative, one process at a time, at most 4 threads, at most 180 s per unit, 6 GiB envelope. The manifest is
frozen and committed before execution, and units run one at a time, stopping on a failure.
