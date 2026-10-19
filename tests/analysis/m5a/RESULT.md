# M5-A result: the pooled winner is already optimal under the ordinary single-block move

**Run.** Branch `m5-a-pooled-winner-refinement`: freeze `5cec0c6`, M4-B closure `e300b1c`. Results in
`tests/results/m5a_pooled_winner_20261019/`.

**Budget.** 38.8 of 600 s. 13 units, all exit 0 and first-time, with no amendments.

**No caller change.** Everything ran through existing probes from the saved M4-B constant-arm runs.

## 1. M4-B closure (appended to its record)

- Six unchanged prediction pairs; five evaluated accuracy comparisons. ACOT accuracy remains unassessed.
- Hazard matching is a geometry experiment, not calibration. Fitting the copying length stays deferred.
- Distance mode stays experimental and default-off.
- **The transition and geometry state is process-global** (`g_transitions`, `g_geometry`). Per-run ownership is
  required before any concurrent multi-sample use; it was not redesigned here.

## 2. Corrected trunk comparison (M4-B `three_switches_1` 2×)

Both terminal incumbents were scored under both priors, with each run's own terminal reproducing its saved total.
Fragment and exposure are identical across the two specifications, so only the prior changes.

| terminal | fragment | exposure | likelihood | prior (constant) | prior (distance) | total (constant) | total (distance) |
|---|---|---|---|---|---|---|---|
| C_term (constant arm) | −29697.874 | −2129.041 | −31826.915 | −44.485 | −44.562 | **−31871.400** | **−31871.477** |
| D_term (distance arm) | −29733.619 | −2128.221 | −31861.839 | −27.502 | −27.535 | −31889.341 | −31889.375 |

**Class: DIFFERENT_BASIN.**
- C_term is above D_term by 17.94 under the constant prior and by 17.90 under the distance prior. The likelihood gap is
  34.92.
- The objective barely changes for these two candidates: their priors move by less than 0.08 between specifications.
- So the distance arm's trunk stopped at a point its **own** objective ranks 17.9 below the constant arm's terminal. It
  is a different, lower basin, not a re-ranking.
- D_term is a single-block local optimum under the distance prior, since the leg stopped there. That the path to C_term
  passes through an intermediate move the distance prior penalises more (the constant arm's round-3 sv:3 step) is
  plausible, but that neighbour was not scored here.

## 3. Pooled-winner continuation

**Setup.** One ordinary leg per root, constant mode, with the caller's refinement settings: shortlist 32, planned origin
cache, `repair_max_neighbours`, cumulative caps, round bound 32. Every root reproduced its saved total in round 1.

| cell | roots | leg | round-1 best neighbour gap | improving neighbour? | combined | genotype / phase change | truth excess (orig → combined) | leg wall / peak RSS |
|---|---|---|---|---|---|---|---|---|
| C4 LZO HG00171 5× | 1 | 1 round, `no_improvement`, converged | −8.80 (639 neighbours scored) | **no** | UNIQUE = original | none / none | 0 → 0 | 8.1 s / 0.75 GiB |
| CYP2D6 LZO HG02769 5× | 1 | 1 round, `no_improvement` | −8.73 (908 scored) | **no** | UNIQUE = original | none / none | 0 → 0 | 9.0 s / 0.74 GiB |
| C4 `three_switches_1` 2× | 2 (W1, W2) | each 1 round, `no_improvement` | −2.37 each (638 scored) | **no** | GENOTYPE_EQUIVALENT = original (both roots retained) | none / none | 8 → 8 | 4.8 + 4.9 s / 0.52 GiB |

- **Ties and refusals:** none. Moves: 0 everywhere.
- **Score components:** the original and terminal are identical (F, E and prior in `eval.*`).
- **Cost:** each leg is a fresh replay with its own setup and placements. That cost is not an in-caller incremental
  overhead estimate, and it is not charged to the original run's budget.

**Why `three_switches_1` does not move.**
- M3-G's +0.613 one-block substitution inside W1/W2 sets ctx:0 to (40,52). That pair ranks **723rd of 2,346** by
  proposal emission (M3-G), far outside the existing 32-pair shortlist, so the ordinary engine cannot propose it.
- The best move inside the shortlist loses by 2.37.
- The reported winner is therefore already a local optimum of the existing single-block neighbourhood, even though it
  was never refined.
- M3-G already showed the +0.613 move goes to another wrong pair, and the truth pair ranks 90th, so widening the
  shortlist would not fix this cell either.

## 4. Decision (frozen rule): **PARK**

- **No change anywhere.** No root has an improving neighbour, no combined decision changes a genotype or phase, and no
  excess changes. There is no scientific benefit to weigh against the added refinement work, here about 5–9 s and
  0.5–0.75 GiB per root as a replay.
- **The limitation needs restating.** "The final pooled winner may never have been refined" stays true as a fact about
  the pipeline. On these three cells, though, refining it with the existing move set changes nothing. The +0.613
  improvement lies outside the shortlist, so a pooled-winner continuation is not where the M3 errors come from.
- **Not tested:** other cells, a widened shortlist, or joint moves (all out of scope). Three cells do not show the
  winner is always single-block optimal.

The structural/count-model comparison remains the next major scientific task.
