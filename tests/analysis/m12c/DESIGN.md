# M12-C: bounded completion and interpretation of M12-B (frozen before execution)

**M12-B's frozen rule, selections, results and decision (INCONCLUSIVE) are preserved.** No caller integration or redesign.
- **Branch:** `m12-c-dosage-first-completion`, from `m12-b-dosage-first` `e5f8c5c`.
- **Results:** `tests/results/m12c_dosage_first_completion_20261021/`.
- **Budget:** 900 s, including 60 s of disclosed development; 6 GiB; one scoring process at a time.

**Not done:**
- no CYP2D6 rerun and no ACOT scoring;
- no new reads, simulation, caller rerun, validation use, feature change, new sequence, fit, prior change, merge, push or
  cleanup.

## 1. Continuation (the manifest is frozen in `MANIFEST.tsv`)

**Scope:** only M12-B's six budget-lost C4 blocks — B10, B12, B13, B14, B16 and B17 (K1).

**Reused unchanged:**
- M12-B's frozen candidate lists (`runs/select/<block>.candidates.tsv`; SHA-256 recorded): exhaustive within the selected
  state, both assignments, the incumbent with its expected total;
- the inputs; the executable (`genotype_probe` and `m10c_align_eval`, hashes equal to M12-B's manifest);
- the scoring command (`m12b_run.score_set`; 4 threads, unchanged caps) and M12-B's per-unit timeouts.

**Not recomputed:** the count-state selection.

**`bind`** verifies, before any scoring:
- the executable hashes;
- that each candidate list's count equals the selection;
- each block's M12-B plan (admitted, NOT RUN on budget);
- the reads SHA-256 against M12-B's manifest.

It then links M12-B's selection, counts and 12 scored outputs read-only into this package.

**Admission:** a score unit runs only if spent + its timeout + the 300 s evaluation reserve ≤ 900.

## 2. Evaluation (`evaluate`)

**M12-B's own `evaluate()` is reused unchanged** over the augmented set (12 + 6 scored blocks).
- **Gates per new block:** reproduction, bytes, shared candidates, aliases.
- **Reproduction check:** M12-B's 12 rows must come out identical.

**Reported, in addition:**
- **Severity transitions** for every scored block: floor / ≤ 10 / < 1,000 / ≥ 1,000. Worsening while already large is
  flagged separately.
- **The combined coverage ledger:** all 29 original blocks, including the 11 ACOT limitations.
- **The frozen decision rule**, applied transparently to the augmented results. **M12-B's historical decision stays
  INCONCLUSIVE.**

**The K5 record comparison** (existing M12-B records only; no new scoring): the selected pair 104/97 against the
evaluator-best pair 95/10 (and 10/95, and the incumbent). For each:
- per-haplotype pooled-feature states, the total state and the feature-0 total;
- count loss;
- fragment, exposure, prior and total;
- distance and excess.

The truth's per-haplotype states are reported for diagnostics.

**Four quantities are kept distinct:**
- the biological total unit count, where established;
- the complete pooled-feature vector;
- the per-haplotype allocation;
- sequence excess.

## 3. Return

**One recommendation:** pursue the rule, park it, or name the specific remaining obstacle. No new audit or sweep.

**Not claimed:**
- that completing C4 resolves ACOT coverage;
- that the result guarantees promotion.

**Interpretation corrections** to M12-B are appended to its RESULT; nothing sealed is rewritten.
