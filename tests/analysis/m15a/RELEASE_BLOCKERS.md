# Release-blocker checklist (short-read release candidate)

| # | blocker | owner / action | status |
|---|---|---|---|
| 1 | Filter procedure approved: no reliable-call class; validity checks A; annotations B; any `REPEATED_UNIT_STRUCTURE` or `DOSAGE_AUDIT_*` formalisation | Davide | **open** |
| 2 | Acceptable per-block excess and normalized excess; catastrophic-error (≥ 1,000 bp) handling | Davide | **open** |
| 3 | `FINAL_EVALUATION_PROTOCOL.md` approved before any reserved input is prepared | Davide / Sol | **open** |
| 4 | Reserved ACOT evaluation run (18 families: read preparation, calls, evaluation) under the approved protocol | next batch | **not started** (outcomes sealed) |
| 5 | RC packaging: named profile against default promotion of the `experimental_*` keys | Davide | **open** |
| 6 | Locus-specific `repair_max_neighbours` (CYP2D6 958, ACOT 1004): per-locus profiles or a general rule; release time-guard default | decision + docs | **open** |
| 7 | ACOT 30× has no reference or ordinary-route equivalence; backend equivalence tested at 5× only | accept as a stated gap, or close | **open** |
| 8 | Full ctest suite on the tagged revision, including the three heavy tests excluded here for budget (`genotype_report`, `genotype_branch_lifetime`, `genotype_block_delta`). `genotype_command_e2e` now takes 175 s, so the stale cost data should be refreshed | before tagging | **open** |
| 9 | Release notes: C4 and CYP2D6 evidence is development-only; simulated reads only; GQ and error mass uncalibrated; dosage audit diagnostic | docs | **open** |
| 10 | Merge, push, tag | Davide | **not done** (out of scope) |

**Closed in M15-A:**
- **Behaviour:** the caller sources are frozen and identical to the checkpoint, and the binary reproduces archived calls
  (15 of 15).
- **Isolation:** M12, M13 and M14 cannot affect predictions (no `src/`/`include/` change).
- **Build:** every target built (23 s).
- **Regression:** the focused regression suite passed. **60 of 60** enabled tests in the frozen selection passed — 58 in
  `ctest`, plus `genotype_command_e2e` and `smoke_real_c4` in `ctest.a1` after amendment A1. 2 tests are disabled in the
  suite, and 3 heavy tests were excluded.
