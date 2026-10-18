# S2-G: complete-reference validation and narrow caller integration (experimental, default off)

- **The route:** profile key `union_fragment_blocks fine_core_4096_experimental`. It runs the combined blocked union route (S2-F) in the caller.
- **Design and contract:** `DESIGN.md`.
- **Gate:** `panvar_genotype_blocked_union_route`.
- **Tools:** runner `s2g_run.py`; caller comparator `s2g_compare.py` (explicit allowances).

**Outcome:**
- **ACOT 30× reference** (ordinary blocked, 1 B origins and state classes) completes. It is bit-identical to the combined route on all 108,939 candidates.
- **Caller 5×** (C4, CYP2D6, ACOT): control and integrated outputs are identical. Footprint peak −3 to −11 %; runtime flat.
- **ACOT 30× integrated caller:** the union completes (2.13 GB at 146.6 s). Refinement then reaches the 4 GiB guard within its first round (no round committed), and the run keeps the union-best checkpoint as PROVISIONAL.

Results: `tests/results/s2g_blocked_integration_20261018/README.md`.
