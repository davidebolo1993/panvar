# Remaining-release checklist (after M15-H)

**Status legend:** done = implemented and tested; open = needs a decision or work before an experimental release.

| item | status |
|---|---|
| no-profile basic command (preset, automatic library, `--dosage-check`, v2 report) | done (M15-D, M15-D.1) |
| declared-library preflight; exit-field agreement | done (M15-D.1) |
| public panel-index build and load | done (M15-E) |
| primary assessment view `assessment_v1` | done (M15-G) |
| help organisation, option inventory, guide | done (M15-H) |
| the `headline_filters_v2` columns: superseded for one experimental release, then removed or moved behind a flag | open (decide timing) |
| block checks for legacy-tied v2 predictions (now NOT_ASSESSED) | open (new scoring and dosage work; separate approval) |
| the caller-only API writes no assessment view | open (document as is, or add) |
| run-guard flags (`--time-limit`, `--memory-limit`) and other aliases (guide §10) | open (separate approval) |
| removals (guide §10: `--dosage-audit`, the `--max-candidates` message, `legacy_i2a`, the diagnostic cache flag spellings) | open (after the release) |
| ACOT 30× route-equivalence gap; CYP2D6 audit unit-library limitation | open (known, documented) |
| release notes: evidence scope (simulated 5×, ACOT reserved, C4/CYP2D6 development only), no calibrated confidence, the two unflagged sv:9 length errors, `--library auto` without accuracy evidence | open (write) |
| merge, tag and publication | open (Davide) |
