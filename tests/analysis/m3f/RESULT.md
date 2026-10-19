# M3-F result: bounded validation beyond the original joint-search recovery

**Run.** Branch `m3-f-validation`. Freeze `7c83812`, Amendments 1–2 (`a1b52af`, `9a53a5c`), post-hoc diagnostic
`a3d2f0b`. Binary sha `a728eb71…` built from `f31b9e6` sources, with no source change since. Results are in
`tests/results/m3f_validation_20261019/`.

**Budget.** 631.5 of 1,200 s, counting every attempt, at 4 threads, one process at a time. The largest caller footprint was
4.40 GB, inside the 6 GiB envelope.

**Scope.** Development validation only: not held out, not a promotion gate.

## Matched results (every attempted cell; refusals kept)

| cell | role | status off → on | excess off → on | blocks / SV exact off → on | joint pass (on) | outcome |
|---|---|---|---|---|---|---|
| LZO HG00171 5× | control (exact) | COMPLETE → COMPLETE | 0 → 0 | 11/5 → 11/5 | complete, no improvement | UNCHANGED |
| `three_switches_1` 5× | control (GE twin) | GE → GE | 0 → 0 | 11/5 → 11/5 | **skipped** (trunk tie-branched) | UNCHANGED |
| `dense_switching_1` 10× | control (GE twin) | GE → GE | 0 → 0 | 11/5 → 11/5 | **skipped** (trunk tie-branched) | UNCHANGED |
| `dense_switching_1` 2× | error | COMPLETE → COMPLETE | 81 → 81 | 5/2 → 5/2 | **skipped** (trunk tie-branched) | UNCHANGED |
| `three_switches_1` 2× | error | GE → GE | 8 → 8 | 10/5 → 10/5 | complete, no improvement | UNCHANGED |
| NA19650 LOO 10× | error (structural) | COMPLETE → **RESOURCE_REFUSAL** | 65,506 → no call | 4/3 → – | **refused inside the pass** (footprint 4.40 GB > 4 GiB guard) | **REFUSAL (harmful)** |

- **Controls:** the gate passed; no control worsened.
- **Predictions:** frozen before truth. Control digest `cd53548b…`; all-cell digest `8b04b7af…`.
- **Selection:** in every completed or skipped pair, the selected unordered pair, the block projections, the winner's
  score and its F/E/P terms are identical between arms.

## Recoveries and harmful changes

- **Recoveries:** none (0 of 3 error cells).
- **Harmful:** one. In NA19650 LOO 10× the on arm turned a (wrong) COMPLETE call into a no-call.
  - The trunk was identical up to the pass: 9 rounds and the same incumbent score, −153310.308864, at about 1.8 GB.
  - The joint pass then pushed the physical footprint past the v6 profile's 4 GiB guard, so the caller refused.
  - The guard was not raised.
- **No harmful move:** neither completed pass moved; there was no score-driven worsening.

## Attribution of the unchanged errors (diagnostic; never fed back)

| cell | frozen class | evidence |
|---|---|---|
| `dense_switching_1` 2× | PASS_SKIPPED_OR_INCOMPLETE | The trunk stopped at an improving-move tie (`rescore_incomplete`) and tie branching resolved it, so by the M3-D contract the pass does not run. 6 wrong blocks; v2 labels 5 of them DATA_SUPPORTED/PASS. |
| `three_switches_1` 2× | OUTSIDE_TWO_BLOCK_NEIGHBOURHOOD | The phased truth differs from the incumbent at 4 blocks on each homologue, and the truth is absent from the 1,848-member catalogue. See the post-hoc note below. |
| NA19650 LOO 10× | PASS_SKIPPED_OR_INCOMPLETE | Refused inside the pass. Off arm: availability floor 113 (truth not spellable), and no scored candidate attains it (closest total-score gap 280.8). sv:7 truth has 4 units; off calls 3+1, an ALLOCATION_ERROR (the M1-B finding). |

**Post hoc** (`posthoc_genotype_level.py`, declared after the freeze; it changes no class). The `three_switches_1` 2×
call is unphased, and its only genotype error is ctx:0: called (4,52), truth (14,56).
- Both homologues must change at one block. That is a cross-homologue move, outside the same-homologue catalogue at
  genotype level too.
- No catalogue member carries the true unordered genotype.
- The needed allele 4→14 ranks 11th in the single-block shortlist, beyond k = 4. Its partner 52→56 ranks 1st.

## Genotype-data support (v2) kept separate from phase

| off-arm blocks, all 6 cells | exact | wrong |
|---|---|---|
| DATA_SUPPORTED / PASS | 53 | **9** |
| PRIOR_DEPENDENT / LOW_QUALITY | 1 | 3 |

- **Wrong blocks labelled PASS:** 9 of the 12 wrong blocks. These are 5 in `dense_switching_1` 2× and 4 in NA19650,
  including sv:7 with an envelope gap of 47.3.
- **Why:** the labels are conditional on the evaluated candidates, as contracted. When the truth was never evaluated,
  DATA_SUPPORTED certifies a wrong genotype. They are not calibrated accuracy.
- **Joint candidates:** in the two completed passes, adding the joint candidates changed no label but shrank 16 envelope
  gaps. For example, `three_switches_1` 2× sv:1 went from 7.00 to 0.75 and ctx:0 from −0.76 to −0.06.
- **Phase:** phase support is unchanged in every pair.

## Costs (on − off; ledger and status)

| cell | wall | joint pass | retained census peak | peak footprint |
|---|---|---|---|---|
| LZO HG00171 5× | 18.5 → 48.1 s (+160%) | 29.8 s | 1.00 GB | 1.32 → **3.55 GB** |
| `three_switches_1` 2× | 17.8 → 30.2 s (+70%) | 12.5 s | 0.33 GB | 0.93 → 2.02 GB |
| NA19650 LOO 10× | 38.7 → 62.5 s (refused) | – | – | 1.80 → **4.40 GB (refused)** |
| skipped passes (3 cells) | ±1 s | ~0 | 0 | ±0.1 GB |

The pass's footprint grows with coverage and catalogue size. At 5× LZO it already reaches 3.55 GB against a 4 GiB guard.

## Recommendation: **restrict** (keep default-off and experimental; do not widen or promote)

- **Scientific reach in this sample.** The pass is applicable in 3 of 6 cells and completes in 2. The pass skips every
  tie-branched trunk (3 of 5 mosaic cells), and it recovers none of the 3 new errors.
- **What the errors need instead.**
  - A cross-homologue change at one block, with one needed allele beyond k = 4.
  - A tie-branched trunk.
  - A LOO allocation the panel cannot spell.
- **What it does well.** Where the pass completed it was safe: no move, and no score-driven harm.
- **Where the harm came from.** Memory, not selection.
- **Next step, if any.** A bounded-memory design for the pass, so it refuses before it starts. Not a wider campaign,
  larger k, or tie-branched trunks; those would be algorithm changes outside this batch.
