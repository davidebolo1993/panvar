# S2-B: feasibility of pre-placement pruning (measurement only)

**What is here:**
- `DERIVATION.md`: the admissible bound against the implemented likelihood;
- `tests/s2b_preplacement_bound.cpp`: the prototype (fixtures and measure modes);
- `s2b_run.py`: the runner.

**The evaluator.** The diagnostic profile key `diagnostic_union_evidence_dump 1` writes `<prefix>.union_evidence.bin`. It holds the interned union sequences, exact masses, exposures, backgrounds, haploid priors, candidate scores, the top set and the early pair. Only the masses and scores are used, and only as the evaluator; the bound itself reads reads, sequences and the library.

**The two variants:**
- **Pure:** sections 1–3 of the derivation; no read–window DP at all.
- **Lazy:** pure, plus exact emissions at every qualifying anchor of S2-A's early-incumbent haplotypes, reused byte-verified wherever another haplotype has the identical window. That exact store is the lazy variant's own cost. Production's early evaluation verifies only seeded starts.

Results: `tests/results/s2b_preplacement_20261018/README.md`.
