# M5-A: bounded pooled-winner refinement pilot (diagnostic; frozen before any unit)

**Ancestry.** Branch `m5-a-pooled-winner-refinement` from `2aba389` (M4-B), plus the M4-B closure (`e300b1c`).

**No caller change.** No production option and no default change. The binary is the M4-B build. Everything runs through
existing probe interfaces:
- `genotype_probe repair-search`, the ordinary single-block engine `CoordinateRepairSearch`;
- `genotype_probe score-set`.

Orchestration lives in this directory. Sources are the saved M4-B **constant-arm** runs: fresh caller outputs on the
current binary, in constant mode.

## Part X: cross-score of the M4-B `three_switches_1` 2× trunks (truth-free)

**Candidates.**
- **C_term:** the constant arm's terminal incumbent, `40,30,40,4,4,6,84,47,35,1,67 | 52,46,19,7,4,5,52,67,3,1,0`, saved
  total −31871.3998675 under the constant prior.
- **D_term:** the distance arm's terminal, `40,30,57,9,9,4,84,47,35,1,67 | 52,46,19,7,4,5,52,67,3,1,0`, saved total
  −31889.3745955 under the distance prior with the frozen C4 λ_HM.

**Scoring.** score-set scores both under the constant prior (on the constant run) and under the distance prior (on the
distance run, same builder, digest checked).
- Each run's own terminal is declared with its saved total, which must reproduce.
- Fragment and exposure must agree across the two scorings, which recomputes only the priors in effect.

**Classification** (`xscore_table.py`), under the distance prior:
- CHANGED_OBJECTIVE if D_term ≥ C_term − 1e-6;
- DIFFERENT_BASIN if C_term > D_term + 1e-6.

The constant-prior comparison is reported beside it. Cross-prior score differences are never read as worse search.

## Part P: continuation from the reported winner (constant mode; joint search off)

**Cells, controls first.**
1. C4 LZO HG00171 5× (correct control).
2. CYP2D6 LZO HG02769 5× (correct control).
3. C4 `three_switches_1` 2× (the pooled-winner limitation).

**Roots (`ROOTS.tsv`, frozen).** The distinct sequence-pair representatives of each saved winning equivalence set, in
rank order.

| cell | roots |
|---|---|
| HG00171 5× | 1 |
| CYP2D6 5× | 1. Its reported call equals the pooled winner: the dosage audit did not apply, and the `.original` and reported equivalence sets are identical. |
| `three_switches_1` 2× | 2: W1 (rank 1) and W2 (rank 2), one unordered genotype in two arrangements |

All cells have two roots or fewer, so all roots are explored. None is chosen as a representative.

**One leg per root.**
- **Engine:** `repair-search --start <root> --start-score <saved total>`. The engine's round-one guard refuses unless
  the start reproduces its saved total within 1e-6.
- **Settings, the caller's refinement settings from the source run's status:**
  - `local_pairs` 32 (the existing shortlist);
  - context, delta and planned origin cache;
  - `repair_origin_cache_bytes`, `repair_max_neighbours` (670 for C4, 958 for CYP2D6);
  - the cumulative repair caps;
  - `switch_probability` 0.01.
- **Bounds:**
  - round bound **32** (the profile's `repair_max_rounds`), from round 1;
  - an empty visited set: a new leg, not the original trunk;
  - footprint guard = the cell profile's (C4 4 GiB, CYP2D6 6 GiB);
  - time limit = declared timeout − 10 s.
- **Tie rule (existing).** The engine moves only to a unique strictly improving best neighbour. A tie stops the leg with
  its stop reason. Nothing breaks the tie, and no branching follows.
- **Not used:** joint moves, a new shortlist, restarts.

**Combination (`combine.py`, truth-free, frozen; detail in the script header).**
- **Candidates:** every original root, retained as a scored candidate, plus every completed leg's terminal.
- **Top set:** every candidate within 1e-6 of the best total.
- **Outcome:** AMBIGUOUS when the top set spans more than one unordered genotype; GENOTYPE_EQUIVALENT for several
  sequence pairs with one genotype; UNIQUE otherwise.
- **INCOMPLETE** when any leg refused, timed out or hit the round bound. Ambiguity between roots is never discarded.
- The selection digest freezes the predictions.

**Evaluation, after the freeze (`eval.<cell>`).** score-set with `--truth-gfa` scores every candidate (originals with
their saved totals, terminals with their leg totals, all as expectations). It gives F, E, prior, total, and distance,
floor and excess against truth, per candidate and per block.

**Reported per root and cell.**
- Start and final score components.
- Moves, convergence, ties and refusals.
- Whether the original call had an improving neighbour (round 1 moved).
- The predicted sequence pair, block-genotype changes and phase-only changes.
- Truth excess.
- Runtime and peak memory (ledger and leg summary).

This is a fresh replay with its own setup and placements. Its cost is not the eventual incremental in-caller overhead,
and it is not charged to the original run's cumulative budget.

**Recommendation rule (`summary.py`, frozen).**
- **BENEFIT:** the combined genotype changes and excess falls.
- **HARM:** the combined genotype changes and excess rises, or any control's genotype changes.
- **Decision:**
  - INTEGRATE if at least one BENEFIT and no HARM;
  - INVESTIGATE if any HARM, or a genotype change at equal excess;
  - PARK if no genotype changes anywhere.

**Expectation from M3-G (not a selection criterion).** The known +0.613 move at `three_switches_1` 2× ctx:0 goes to
(40,52), which is still wrong.

## Bounds

- 600 s counting every attempt, 4 threads, unchanged resource settings.
- **Declared timeouts:**
  - xscore: 60 s each, table 30 s;
  - legs: HG00171 120 s, CYP2D6 120 s, `three_switches_1` 90 s per root;
  - combine 30 s; evaluations 90 s each; summary 30 s.
- Not done: no ACOT, fitting, simulation, structural model, support-rule change, merge or push.
- Inputs are hashed in INPUTS.tsv.
