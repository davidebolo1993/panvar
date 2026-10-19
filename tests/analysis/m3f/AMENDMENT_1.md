# M3-F Amendment 1 (declared after `accuracy` failed and before any further unit)

**State at declaration.** 613.6 s spent. Recorded so far:
- C complete; `control.gate` PASS.
- E: all six error runs and the `compare` freeze (digest `8b04b7af…`) are done; the six evaluations and `accuracy` have
  been attempted.
- `accuracy` failed (exit 1): a script defect.

**What happened.**
1. `run.loo_na19650_10.on` ended RESOURCE_REFUSAL (`footprint_limit`) in phase `joint search`, with a physical footprint
   of 4,401,073,992 bytes against the v6 profile's 4,294,967,296-byte guard. The trunk up to that point equals the off
   arm's: 9 rounds, the same final incumbent score of −153310.308864.
2. As a result, `eval.loo_na19650_10.on` failed (exit 1, `map::at`): there is no call to evaluate, and it left a header
   plus one partial row in `eval.blocks.tsv`.
3. `accuracy_m3f.py` read that partial row (a field of None) and crashed while writing the support table. The outcome
   logic itself was not reached.

**Changes** (analysis scripts only; no caller, probe, setting, input or prediction changes):
- `accuracy_m3f.py`: a block table is read only when that arm's evaluation summary exists.
- `attribute_m3f.py`:
  - an error cell whose on arm has no evaluable call because it was refused in phase `joint search` is classed
    PASS_SKIPPED_OR_INCOMPLETE;
  - its detail reports the off arm's pool and floor evidence and the off-arm sv:7 structure;
  - "no error in either arm" is made explicit.
- `m3f_run.py`:
  - `structure.loo_na19650_10` passes the off table in both columns when the on evaluation failed, and only the off
    ("historical") column is read;
  - new manifest A1 with one unit, `accuracy.a1` (30 s), run before manifest D. The failed `accuracy` unit stays in the
    ledger.

**Unchanged:** the refusal is a result. The 4 GiB guard is not raised, and the LOO on arm stays REFUSAL_OR_NO_CALL in
every denominator.
