# M11-D amendment A4 (declared before further execution)

**State at declaration:** `report` exited 1 at `float(normalized_excess)`. Spent: 1,941.8 of 3,600 s.
- **Cause:** 34 predicted blocks are **zero-length**: both truth alleles and both called alleles are empty, in CYP2D6 sv:1,
  sv:11 and sv:15, and ACOT sv:3 and sv:11.
- **Consequence:** their normalized excess is undefined (0/0) and was written empty by `outcomes.a2`, as designed. Their
  excess is 0, and their QVfrac is 1 by the L = 0 → QV 60 convention.
- **Fit unaffected:** these blocks are in `fit.a2` with y = 1, as frozen.

**Fix (reporting only):** normalized excess is summarised over defined values, and `zero_length_blocks` is reported as its
own count. The same script was dry-run into the scratchpad first (0.1 s, outputs discarded).

**Disclosed A4 development:** 60 s.

**Re-run unit:** `report.a4` (120 s). The failed `report` unit stays recorded.
