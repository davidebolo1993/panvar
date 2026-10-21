# M12-A amendment A1 (declared before further execution)

**State at declaration:** `build` passed; `gate.fixtures` exited 1 with 48 of 3,976 checks failed. Spent: 153.0 of 900 s.

**What failed, and why:**
- **42 G1 checks.** All are `k=0 boundary oracle=1 brute=1 analytic=0`: templates lying wholly inside a flank (for example,
  starts 0–20 in the 60 bp left flank).
  - These exist on H_0 too, so the oracle and the brute-force count both find exactly one origin.
  - My analytic class "boundary" wrongly expected 0 for every k = 0 boundary fragment.
- **6 G2 checks.** These are the fragment terms of the candidates (0, ·), whose independent score used the same analytic
  counts.
- **Everything else passed,** including every oracle = brute force comparison, every scorer = oracle mass, the exposures,
  G3 (the peak at the true total of 6 copies), G4 and G5.

**This is a defect in the fixture code, not the implementation.**

**Fix (fixture only):**
- a template wholly inside one flank is class `flank` and has 1 origin on every H_k;
- a template that crosses a flank|array boundary stays `boundary` (1 origin for k ≥ 1, 0 for k = 0).

**Units:**
- re-run as `build.a1` (300 s) and then `gate.fixtures.a1` (120 s);
- the case units now need `build.a1`;
- `analyse` reads `gate.fixtures.a1`;
- the failed `gate.fixtures` unit stays recorded.

**No comparison, threshold, case or classification rule changes.**
