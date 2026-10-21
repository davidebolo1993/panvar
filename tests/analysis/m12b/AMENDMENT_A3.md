# M12-B amendment A3 (declared before further execution)

**State at declaration:** spent 1,410.7 of 1,800 s.
- **Scored:** 12 blocks (B01–B09, B11, B15, B18).
- **NOT RUN on the A2 budget reserve:** 6 blocks (B10, B12, B13, B14, B16, B17).

**The failure:** `evaluate` exited 1 after writing its per-block lines for B01–B06. Reading B07's feature k-mer list raised
`_csv.Error: field larger than field limit (131072)`.

**Fix (mechanical):** `csv.field_size_limit(sys.maxsize)` in `m12b.py`.
- Re-run as `evaluate.a3` (300 s).
- The failed unit and its partial log stay recorded.
- **No selection, rule, gate or decision criterion changes.** The selections were frozen in `select` and the score units
  before truth was opened.
- **Disclosure:** the partial log shows B01–B06's outcomes, and evaluate.a3 recomputes them identically.
