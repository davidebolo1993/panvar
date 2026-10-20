# M8-D result: given its declared runtime, HG01975 LOO completes in both arms (~160 s); r = 0.1 calls a different, closer genotype (excess 49 → 28 bp) because the search reached a different state, not because the ranking changed

**Run.** Branch `m8-d-hg01975-loo`, freeze `f65ffec`, with no amendments. Results are in
`tests/results/m8d_hg01975_loo_20261020/`.

**Unchanged:** M8-C, its timed-out runs and its tables.

**Budget:** 429.4 s of the 2,100-s budget for the units, plus a disclosed pre-freeze dry run of the hashcheck (about 3 s).
- No caller guard and no watchdog bound.
- No failures were recorded.

## Inputs (`hashcheck`: all pass)

- **Executables:** panvar `98f339e8…` and genotype_probe `af2ed281…` (full SHA-256 equal to M8-C `prepare.tsv`).
- **Graph:** `acot.gfa` `05a1a82a…`.
- **Reads and library:** HG01975 cov5 reads_1 `47dfcb74…`, reads_2 `e9afc4f4…`, library `d8b4ee7f…`, equal to
  `prepare.tsv`.
- **Curation and the two profiles:** M8-C did not hash these, so they were re-derived with M8-C's frozen prepare rule and
  are byte-identical.
- **Commands:** M8-C's recorded command lines reference exactly these files. Only `--out-prefix` differs, pointing to new
  directories.

## Outcomes (frozen before any truth was read)

| run | outcome | switch_probability (settings) | wall / caller seconds | peak RSS | guard / watchdog |
|---|---|---|---|---|---|
| HG01975 LOO r = 0.01 | **COMPLETE** (strict_call, filter PASS) | 0.01 | 162.9 / 158.5 | 1.86 GiB | neither bound |
| HG01975 LOO r = 0.1 | **COMPLETE** (strict_call, filter PASS) | 0.1 | 161.8 / 157.7 | 2.52 GiB | neither bound |

- **M8-C's 120-s unit bound** stopped both runs at roughly three quarters of the runtime they needed. Tie branching alone
  takes 83–89 s.
- **Determinism:** the M8-D repair trajectories reproduce the M8-C partial traces move for move.

## Evaluation (original evaluator, LOO, completed outputs only)

| | r = 0.01 | r = 0.1 |
|---|---|---|
| unordered sequence pair (hap SHA-256 / bp) | `…` 224,255 + 252,277 | `…` 224,256 + 252,256 (**different**) |
| switches | 4 + 6 (mosaic) | 4 + 6 (mosaic) |
| projection ambiguity | none (winner and every scored candidate) | none |
| availability floor | 22 | 22 |
| **unordered distance** | **71** | **50** |
| **excess over floor** | **49** | **28** |
| exact blocks (SV) | 7/19 (5/9) | **8/19** (5/9) |
| context / SV distance | 49 / 22 | 28 / 22 |
| read-supported blocks (exact; excess) | 4 (2; 19) | 3 (1; 19) |
| phase | PRIOR_ONLY; phased distance NA (truth not spellable under LOO) | same |
| separation gap / candidate-set probability | 0.044 nats / 0.20 | 0.043 nats / 0.16 |

**The calls differ in two blocks, and r = 0.1 is better in both.** Both blocks have PRIOR_ONLY genotype support:

| block | r = 0.01 call → excess | r = 0.1 call → excess | truth-closest alleles |
|---|---|---|---|
| ctx:2 | 74/76 → 20 bp | 74/73 → **0** | 74/73 |
| ctx:4 | 306/242 → 1 | 307/269 → **0** | 306/269 (floor 3) |

The other ten non-exact blocks are identical in both arms, including the read-supported errors at sv:5 and ctx:8.

## Search comparison and decomposition

**Moves** (accepted repair moves, in order):

| run | moves | round-4 delta |
|---|---|---|
| r = 0.01 | ctx:4, ctx:6, ctx:8, **ctx:0** | +29.5 |
| r = 0.1 | ctx:4, ctx:6, ctx:8, **sv:13**, ctx:0 | +32.3 |

- Both arms then hit a tie at ctx:2 and branch twice (`tie_branches_run 2`, two distinct terminal genotypes, complete).
  Each winner is a tie-branch terminal.
- "Search work identical", as M8-C reported, holds only for its completed pairs. Here the candidate pools also differ:
  103,206 vs 104,039 sequence candidates.

**Same contrast, rescored under each run's objective** (`score-set`, same model, curation and reads, the run's own r; each
caller total reproduced):

| scored under | candidate | reached by that run | fragment | exposure | prior | total |
|---|---|---|---|---|---|---|
| r = 0.01 | winner (r = 0.01) | yes | −77,433.192 | −5,947.927 | −98.009 | −83,479.128 |
| r = 0.01 | winner (r = 0.1) | **no** | −77,434.029 | −5,947.677 | −95.437 | **−83,477.143** |
| r = 0.1 | winner (r = 0.01) | yes (tie branch 0, round 8) | −77,433.192 | −5,947.927 | −79.600 | −83,460.718 |
| r = 0.1 | winner (r = 0.1) | yes | −77,434.029 | −5,947.677 | −77.078 | **−83,458.784** |

Each row of the next table is the r = 0.1 winner minus the r = 0.01 winner:

| scored under | fragment | exposure | prior | total |
|---|---|---|---|---|
| r = 0.01 | −0.84 | +0.25 | +2.57 | +1.99 |
| r = 0.1 | −0.84 | +0.25 | +2.52 | +1.93 |

**This is different reachability, not changed ranking:**
- Both objectives prefer the r = 0.1 winner by about 2 nats.
- Under r = 0.01 that winner was never reached: it is absent from the run's scored candidates and from its repair-round
  states.
- The r = 0.1 run did reach the r = 0.01 winner, as a scored state in tie branch 0, round 8, and ranked it lower.
- What r changed is the search path, starting with sv:13 before ctx:0 in round 4. It did not change the ordering of these
  two states.

**What favours the closer call:** the prior (+2.5 nats) and the exposure term (+0.25). The reads slightly favour the
farther call (fragment −0.84). The margins are small: the separation gap is 0.04 nats in both runs.

## Corrected missing-result semantics

`same_call` is **UNKNOWN** whenever either run is incomplete or missing. When both are complete it compares unordered
haplotype-sequence pairs.

**Fixtures (`fixtures.tsv`, 5/5 PASS):**
- both missing → UNKNOWN;
- one missing → UNKNOWN;
- one incomplete → UNKNOWN;
- swapped → True;
- different → False.

M8-C's historical `summary.a2/paired.tsv`, where the timed-out pair reads `same_call = True`, is unchanged. The corrected
table is `runs/summary/paired_corrected.tsv`.

## Closure: the combined experiment, eight completed runs

| donor × scheme | same_call | excess r = 0.01 → r = 0.1 | floor | verdict |
|---|---|---|---|---|
| HG03209 LZO | True | 0 → 0 | 0 | exact in both (non-discriminating) |
| HG03209 LOO | True | 5 → 5 | 9 | equal |
| HG01975 LZO | True | 0 → 0 | 0 | exact in both (non-discriminating) |
| HG01975 LOO | **False** | **49 → 28** | 22 | **r = 0.1 better** |

M8-C's two TIMEOUT rows (120-s unit bound) are kept in `cells.tsv` as historical. M8-D's two completed outcomes supersede
them for the comparison only.

**Decision (frozen rule, predictions differ): report without retuning; no default promotion.**
- **Benefit:** observed in one of two LOO cells (−21 bp excess, +1 exact block), with no deterioration in any cell. Both
  LZO controls stay exact.
- **Mechanism:** the benefit comes from a changed search trajectory, which reaches a state that r = 0.01's own objective
  would also prefer.
- **Fragility:** a near-tie winner (gap 0.04 nats) in PRIOR_ONLY blocks, decided through a ctx:2 tie branch, and slightly
  disfavoured by the reads.
- **Not shown:** that r = 0.1's objective ranks genotypes better, or that the gain generalises. The evidence is one donor
  at one coverage on one locus.
- **Status:** r = 0.1 stays a candidate parameter with one reachability-mediated LOO improvement. It is not promoted.
