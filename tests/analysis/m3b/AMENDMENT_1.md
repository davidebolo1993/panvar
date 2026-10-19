# M3-B Amendment 1: split scoring for the 5× cells; one control within budget

Declared after `score.mosaic_2x` (139.7 of 600 s spent) and before any further execution.

## Observation

`score.mosaic_2x` scored 1,849 candidates (1,850 haplotypes) in 126.0 s, with a peak RSS of 1.92 GB. The incumbent's
total reproduced.

Scoring time scales with fragments × haplotypes:

| cell | fragments | projected scoring time |
|---|---|---|
| `mosaic_2x` | 2,214 | 126 s (measured) |
| `mosaic_5x` | 5,517 | ≈ 314 s |
| `lzo_5x` | 5,215 | ≈ 297 s |

Each 5× cell is above the 180-s per-unit cap, and both together are above the remaining budget (about 460 s).

## Change

- **Split.** Each scored 5× catalogue is split deterministically (`split_catalogue.py`): catalogue order, with the first
  and second halves of the non-incumbent rows as two parts.
  - Both parts carry the incumbent with its expected total, so each part must reproduce it.
  - A candidate's score depends only on its own haplotypes, the fragments and the library, so the union of the two
    parts is the full catalogue scored once.
  - `rank.py` merges the parts and asserts that the incumbent rows agree exactly.
- **Budget choice.** The budget admits **one** complete 5× control: the **mosaic 5× twin**, the direct counterpart of the
  failing cell, with the same simulated mosaic and a correct incumbent.
- **LZO control.** `lzo_5x` scoring is **not run**. Its frozen catalogue (1,848 candidates) stays on disk, and the gap is
  reported: deterioration of the LZO control is **not assessed** in this batch.
- **Superseded units.** From manifest B, `score.mosaic_5x`, `score.lzo_5x`, `rank` and `eval.*` are not run.
  - The proposal rule, top-4, the 2,048 limit, the scorer and the selection rule are unchanged.

**Units (manifest A1):**

| unit | declared timeout |
|---|---|
| `split.mosaic_5x` | 30 s |
| `score.mosaic_5x.part1` | 180 s |
| `score.mosaic_5x.part2` | 180 s |
| `rank.a1` | 30 s |
| `eval.a1.mosaic_2x` | 60 s |
| `eval.a1.mosaic_5x` | 60 s |
