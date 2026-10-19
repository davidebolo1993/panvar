# M3-C Amendment 1: context lifetime in the gate fixture

Declared after `gate` (0.6 of 600 s spent, exit 1) and before any further execution.

**Failure.** `gate` stopped with `another PreparedScoringContext already holds the emission cache`.
- The library owns one emission cache per process, so only one scoring context may be alive at a time.
- The fixture kept the cold reference context alive while creating the incremental one.
- No comparison was made. This is a harness defect.

**Change.** The cold reference is scored inside its own scope and destroyed before the incremental context is created.
- Cases, orders, homologues, the tolerances (1e-9), the byte-identity and same-slot requirements, and the
  zero-full-placement requirement are all unchanged.
- `two-block-lik` already uses one context per process and is unaffected.
- `gate` stays recorded with exit 1.

**Unit:** `gate.a1`, declared timeout 60 s.
