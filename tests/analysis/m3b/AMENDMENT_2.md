# M3-B Amendment 2: intermediate lookup in `rank.py` (reporting bug)

Declared after `eval.a1.*` (417.0 of 600 s spent) and before any further execution.

## The bug

`rank.py` identified the single-block intermediates of the best two-block move by splitting each move component at its
first `:`. Block IDs themselves contain a colon (`ctx:0`, `sv:1`), so the lookup matched the first single move at that
block instead of the same allele change.

**Example.** The 2× joint move `ctx:0:36>4 + sv:1:67>24` was decomposed with `c0001` (ctx:0 36>0) instead of the
single move 36>4. The 5× decomposition was wrong in the same way.

## What is unaffected

- The winner, the incumbent's rank and the count above the incumbent: selection does not use this parsing.
- `eval.a1` outputs for the winner, the incumbent and the top 10.

## Change

- **Fix.** Components are split from the right as `<block id>:<from>><to>:r<rank>`.
- **Units:**
  - `rank.a2` re-runs the truth-free selection. Its `selection.tsv` must be byte-identical to `rank.a1`'s; only
    `decomposition.tsv` may change.
  - `eval.a2.*` re-evaluates so the correct intermediates are included.
- **Records kept.** The `rank.a1` and `eval.a1` outputs stay as recorded. Only the A2 decomposition is used.

**Unit timeouts:** `rank.a2` 30 s, `eval.a2.mosaic_2x` 60 s, `eval.a2.mosaic_5x` 60 s. 417.0 + 30 + 60 + 60 ≤ 600.
