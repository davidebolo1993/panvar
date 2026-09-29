"""Paths, seals and provenance shared by the C4 genotype-report analysis scripts.

No path is hard-coded. The repository root is this file's repository (or --repo-root); the results root holding the
sealed C4 packages is --results-root, by default <repo>/tests/results (gitignored, so a clean worktree must name it).
Absolute paths recorded inside the sealed tables are re-rooted under the results root when they do not exist as
written.
"""
import argparse
import hashlib
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))

# the sealed sources, relative to the results root: (directory, seal file, what it seals)
CAMPAIGN = "genotype_c4_production_11fe681"
FINAL = "genotype_c4_production_11fe681_dosage_repair/final_26401b4"
FRESH_JOBS = "genotype_c4_production_11fe681_dosage_repair/policyC_test/jobs.tsv"
HEADLINE = "c4_headline_filters_26401b4"
PHASE = "c4_phase_eval_5b1f059"
SEALS = [(FINAL, "seal_replay.sha256", "the 649 replay output prefixes (genotyper release 26401b4)"),
         (FINAL, "seal_package.sha256", "the final evaluation package and post-seal evaluator tables"),
         (HEADLINE, "seal_derived_v2.sha256", "the headline_filters_v2 derived tables"),
         (HEADLINE, "seal_eval_v2.sha256", "the headline_filters_v2 truth evaluation"),
         (PHASE, "seal_tables.sha256", "the phase evaluation tables")]
PANEL_GROUPS = "LZO=lzo,fresh_lzo;mosaic LZO=mosaic;LOO=loo,fresh_loo"
TITLE = "C4 short-read genotype report (final_26401b4 replay, headline_filters_v2)"


def parser(description, work=True):
    p = argparse.ArgumentParser(description=description)
    p.add_argument("--repo-root", default=DEFAULT_REPO, help="repository holding scripts/plot_genotype_results.py "
                                                              "(default: this script's repository)")
    p.add_argument("--results-root", help="directory holding the sealed C4 packages (default: <repo>/tests/results)")
    if work:
        p.add_argument("--work-dir", required=True, help="directory for the built inputs and reports")
    return p


def roots(args):
    repo = os.path.abspath(args.repo_root)
    results = os.path.abspath(args.results_root or os.path.join(repo, "tests", "results"))
    for d, _, _ in SEALS:
        if not os.path.isdir(os.path.join(results, d)):
            sys.exit(f"{os.path.join(results, d)} not found: pass --results-root (tests/results is not in a clean worktree)")
    tool = os.path.join(repo, "scripts", "plot_genotype_results.py")
    if not os.path.isfile(tool):
        sys.exit(f"{tool} not found: pass --repo-root")
    return repo, results


def reroot(path, results):
    """A path recorded in a sealed table: as written if it exists, otherwise re-rooted under the results root at its
    first component that names a sealed source directory."""
    if os.path.exists(path):
        return path
    parts = path.replace("\\", "/").split("/")
    for anchor in (CAMPAIGN, FINAL.split("/")[0], HEADLINE, PHASE):
        if anchor in parts:
            return os.path.join(results, *parts[parts.index(anchor):])
    return path


def sha256(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def verify_seals(results):
    ok = True
    for d, seal, _ in SEALS:
        r = subprocess.run(["shasum", "-a", "256", "-c", "--quiet", seal], cwd=os.path.join(results, d), capture_output=True)
        ok &= r.returncode == 0
    return ok


def git(repo, *args):
    r = subprocess.run(["git", "-C", repo] + list(args), capture_output=True, text=True)
    return r.stdout.strip() if r.returncode == 0 else "unknown"


def provenance(repo, results):
    """What a report built here depends on: the reporting tool's commit (and whether its files match it), the genotyper
    release that produced the outputs, and every sealed package with its seal's own digest."""
    tracked = ["scripts/plot_genotype_results.py", "scripts/genotype_report", "tests/analysis/c4_genotype_report"]
    dirty = git(repo, "status", "--porcelain", "--", *tracked)
    rows = [("reporting_tool_commit", git(repo, "rev-parse", "HEAD")),
            ("reporting_tool_files", "clean at that commit" if dirty == "" else "MODIFIED: " + dirty.replace("\n", "; ")),
            ("genotyper_outputs_release", "26401b4 (the final_26401b4 replay of the 11fe681 campaign)"),
            ("genotyper_binary_sha256", sha256(os.path.join(results, FINAL, "bin", "panvar_26401b4"))),
            ("campaign_release", "11fe681 (" + CAMPAIGN + ", 484 jobs)"),
            ("headline_values", "headline_filters_v2 audit derived_v2, which 02f846a's apply_headline reproduces exactly "
                                "(headline_filters_review gate R: 0 mismatches over 7117 blocks and 649 loci)")]
    for d, seal, what in SEALS:
        rows.append((f"sealed:{d}/{seal}", f"{sha256(os.path.join(results, d, seal))} ({what})"))
    return rows


def write_record(path, rows):
    with open(path, "w") as f:
        f.write("key\tvalue\n")
        for k, v in rows:
            f.write(f"{k}\t{v}\n")


def tool_argv(results, inputs, out, meta=None, truth=True):
    """The report command, as a list of arguments to scripts/plot_genotype_results.py."""
    a = ["--input", os.path.join(inputs, "package"), "--out-dir", out, "--run-metadata", meta or os.path.join(inputs, "run_metadata.tsv"),
         "--panel-groups", PANEL_GROUPS, "--title", TITLE]
    if truth:
        a += ["--evaluation-blocks", os.path.join(inputs, "evaluation_manifest.tsv"),
              "--switch-metadata", os.path.join(inputs, "switch_metadata.tsv"),
              "--dosage-content-evaluation", os.path.join(results, HEADLINE, "eval_v2", "dosage_content_blocks.tsv"),
              "--representative-diagnostics", os.path.join(results, FINAL, "package", "layer_D_conditional_representatives.tsv"),
              "--phase-boundaries", os.path.join(results, PHASE, "tables", "boundaries_primary.tsv"),
              "--phase-calls", os.path.join(results, PHASE, "tables", "calls_primary.tsv"),
              "--phase-projected", os.path.join(results, PHASE, "tables", "loo_sequence_projected_boundaries.tsv")]
    return a
