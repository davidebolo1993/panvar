# M12-B amendment A2 (declared before any score unit and before any truth comparison)

**State at declaration:** `select` passed and froze its plan. Spent: 990.0 of 1,800 s. No score unit has run, and truth has
not been opened.

**The plan:**
- **Scored:** 18 blocks (CYP2D6 B01–B06, C4 B07–B18).
- **Limitations:** all 11 ACOT blocks (estimates 300–649 s against the 240 s cap).

**The problem:**
- **Scoring time:** from the preflight rates without the 1.5 safety factor, the 18 blocks take about 840 s of actual
  scoring.
- **Admission:** the runner admits a unit while spent + its declared timeout ≤ 1,800, but it reserves nothing for
  `evaluate` (declared 300 s).
- **Consequence:** scoring could leave evaluation inadmissible, so no selection would ever be compared with truth.

**Fix (budget safety only):**
- a `score.*` unit is admitted only if spent + its declared timeout + 300 (evaluate's declared timeout) ≤ 1,800;
- otherwise it is recorded NOT RUN with the budget reason, as a limitation;
- the order (manifest), the estimates, the watchdogs and every rule are unchanged;
- no block is reordered or shortlisted.

**Operator loop:** it continues past a NOT RUN or a failed score unit, because each block is independent and its failure is
that block's recorded result. It stops on any other failure.
