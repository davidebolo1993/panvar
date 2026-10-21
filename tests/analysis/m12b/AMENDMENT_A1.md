# M12-B amendment A1 (declared before further execution)

**State at declaration:** all units through `count.B29` passed. Spent: 982.1 of 1,800 s.
- **The failure:** `gate.features` exited 1 in 0.2 s with `KeyError: 'gate.features'`.
- **The cause:** the runner dispatched the unit name as the analysis mode, but the mode is `gate_features`. The gate never
  ran.

**Fix (runner dispatch only):**
- re-run as `gate.features.a1` (120 s; same code path and the same frozen checks);
- the `score.*` units now need `gate.features.a1`;
- the failed unit stays recorded.

**Pre-checked:** the other analysis units' names already equal their modes (`select`, `evaluate`, `preflight_candidates`;
`gate.order` maps to `gate_order`).

**No rule, feature, loss, bound or decision criterion changes.**
