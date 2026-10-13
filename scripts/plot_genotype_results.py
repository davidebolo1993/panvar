#!/usr/bin/env python3
"""Summary tables, figures and a short report for `panvar genotype` outputs (analysis and plotting only).

Input (--input) is either a directory holding many output prefixes (<prefix>.locus_genotype.tsv, found recursively;
run_id is the prefix path relative to the directory) or a consolidated package directory with runs.tsv and blocks.tsv.

Truth-free QC needs nothing else. Accuracy figures are drawn only from explicitly supplied evaluation tables and are
refused otherwise; truth is never inferred from sample or run names. See docs/modules/genotype.md (Reporting) and --help.

Writes <out-dir>/data/*.tsv, <out-dir>/figures/*.png|pdf, <out-dir>/report.md and <out-dir>/manifest.tsv. The inputs
are only read; an --out-dir inside an input location, or one that is not empty, is refused.

Exit status: 0 on success, 2 on a refusal (malformed or incompatible input), with the reason on stderr.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from genotype_report.common import Refusal  # noqa: E402
from genotype_report.pipeline import run  # noqa: E402


def parse(argv):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, help="a directory of output prefixes, or a package directory (runs.tsv, blocks.tsv)")
    p.add_argument("--out-dir", required=True, help="new or empty output directory")
    p.add_argument("--input-kind", choices=["auto", "directory", "package"], default="auto",
                   help="auto: a package when --input holds runs.tsv and blocks.tsv")
    p.add_argument("--run-metadata", help="TSV: run key (run_id, or dataset and job) plus experiment, coverage and optionally "
                                          "subject, modality, execution_group")
    p.add_argument("--evaluation-blocks", help="TSV of the evaluator's per-block fields (run key, block_id, unordered_distance, "
                                               "availability_floor, excess, truth1_bp, truth2_bp, called1_bp, called2_bp), or a "
                                               "manifest (run key, path) of per-run evaluator block tables")
    p.add_argument("--dosage-content-evaluation", help="TSV: run key, block_id, content_correct[, content_truth_source]")
    p.add_argument("--representative-diagnostics", help="TSV: run key, block_id, original_distance, "
                                                        "representative_phased_distance[, role] (non-authoritative)")
    p.add_argument("--switch-metadata", help="TSV: subject (or run key), switch_class, switch_before_block_hap1, "
                                             "switch_before_block_hap2[, source_a, source_b]")
    p.add_argument("--phase-boundaries", help="phase evaluation boundary table (run key, boundary_index, boundary, status, "
                                              "reason, phase_correct, switch_error)")
    p.add_argument("--phase-calls", help="phase evaluation per-call table (run key, exact_phased_locus, assessed_boundaries, "
                                         "correct_boundaries, switch_errors, all_assessable_correct, any_switch_error)")
    p.add_argument("--phase-projected", help="secondary sequence-projected phase table (run key, boundary_index, "
                                             "projected_status, projected_concordant)")
    p.add_argument("--panel-groups", help="rows of the ordered-block figures: 'NAME=exp1,exp2;NAME2=exp3' (every experiment "
                                          "exactly once); default one row per experiment")
    p.add_argument("--catastrophic-bp", type=float, default=1000.0, help="catastrophic error threshold on excess (default 1000)")
    p.add_argument("--title", default="panvar genotype report")
    return p.parse_args(argv)


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    args = parse(argv)
    try:
        r = run(args, [os.path.basename(sys.argv[0])] + list(argv))
    except Refusal as e:
        print(f"REFUSED: {e}", file=sys.stderr)
        return 2
    print(f"{len(r.figures)} figures, {len(r.tables)} tables, {len(r.skipped)} skipped -> {r.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
