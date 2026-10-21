# Remaining-release checklist (after M15-C)

| # | item | status |
|---|---|---|
| 1 | Davide's interface review: `docs/guides/genotype_short_read.md`, the help text, and the complete per-locus profiles | **open** |
| 2 | Experimental-release decision (M15-B supports an *experimental* release within its stated scope) | **open** (Davide) |
| 3 | Full enabled ctest on the final tagged revision. The help-corrected binary (`panvar` `5197fd97…`) has run only the CLI contract test and the documented examples; M15-B's 63-of-63 run used `47042a5c…`, whose sources differ only in help strings | **open** (before tagging) |
| 4 | Release notes stating scope and limits: ACOT reserved evidence only; simulated reads with a declared library; C4/CYP2D6 development-only; no calibrated confidence; PASS is not reliability (the two sv:9 errors); `--library auto` and real reads untested | **open** |
| 5 | Packaging: ship the named experimental profiles (`profiles/genotype_experimental_sr_v1_{c4,cyp2d6,acot}.tsv`); defaults unchanged | **done in-tree**, awaiting approval |
| 6 | Locus-specific `repair_max_neighbours`: documented per locus; no general rule | **documented** (decision recorded as "retain") |
| 7 | The ACOT 30× reference-route equivalence gap | **open** (accept as stated, or close) |
| 8 | Merge, push and tag | **not done** (Davide) |
| 9 | Profile layering: repeated `--profile` files override silently, with no merged-profile digest. Documented, not changed | **documented** (possible future interface work) |
| 10 | Panel-sequence recovery uses a repository script (`tests/analysis/s1d/extract_panel_haplotype.py`), not a `panvar` subcommand | **documented** (possible future interface work) |
