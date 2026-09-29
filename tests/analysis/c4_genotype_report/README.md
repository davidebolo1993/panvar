# C4 genotype report: the audited analysis

These scripts build the C4 truth-evaluation report with the generic, locus-independent tool
(`scripts/plot_genotype_results.py`, documented in `docs/modules/genotype_report.md`). They check the report against the
sealed C4 tables. They never genotype, and they only read the sealed packages.

| script | role |
|---|---|
| `build_inputs.py` | Verifies the seals, then joins the sealed sources into `<work-dir>/inputs`: a truth-free package, the run metadata, an evaluator manifest and the mosaic switch plans. |
| `validate.py` | Writes the report to `<work-dir>/report`, then runs the gates: determinism, reconciliation, redraw from `data/`, modality, truth-free mode, and unchanged inputs and seals. It writes `run_record.tsv` and `validation.txt`. |
| `reconcile.py` | Checks the report's `data/` against the 484-job campaign, the 649-cell replay package, the headline audit and the phase evaluation. |
| `c4_common.py` | Roots, seals, re-rooting of recorded paths, and provenance. |

**Roots.** The repository root is this directory's repository, or `--repo-root`. The results root holding the sealed
packages is `--results-root`, by default `<repo>/tests/results`. That directory is gitignored, so a clean worktree
must name it. No path is hard-coded. Absolute paths recorded inside the sealed tables are re-rooted under the results
root when they do not exist as written.

**Run:**

    W=<results root>/c4_genotype_report_<reporting-tool commit>
    python3 tests/analysis/c4_genotype_report/build_inputs.py --results-root <results root> --work-dir "$W"
    python3 tests/analysis/c4_genotype_report/validate.py --results-root <results root> --work-dir "$W"

**What `run_record.tsv` records:**
- the reporting tool's commit, and whether its files match it;
- the genotyper release that produced the outputs (26401b4, the replay of the 11fe681 campaign) and its binary's SHA-256;
- the source of the headline values;
- every sealed package, with its seal file's own SHA-256.

Each report's `manifest.tsv` also hashes every input and output file.
