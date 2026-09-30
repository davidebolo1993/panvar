#!/usr/bin/env python3
"""C4-specific adapter for the generic genotype report (scripts/plot_genotype_results.py). It only joins sealed inputs:
it computes no genotype, filter or evaluation value, and it writes only <work-dir>/inputs.

Sources, under --results-root, every seal verified first:
  final_26401b4/replay/{campaign,fresh}  the 649 genotyper output prefixes, read with the tool's own directory loader
  c4_headline_filters_26401b4/derived_v2 the sealed headline_filters_v2 values. The replay outputs predate the headline
                                         columns; 02f846a's apply_headline reproduces these values exactly.
  final_26401b4/package/runs.tsv         experiment, subject, coverage
  evaluator block tables                 per run, as the sealed package chose them (evaluate_integration.py): the
                                         post-seal eval/ table for a dosage-repaired call, else the campaign's own
                                         <out_prefix>.eval.blocks.tsv
and <repo>/tests/genotype_c4_heldout_mosaic_manifest.tsv, the frozen mosaic plans.

Writes <work-dir>/inputs/: package/runs.tsv and package/blocks.tsv (truth-free), run_metadata.tsv,
evaluation_manifest.tsv (paths relative to it), switch_metadata.tsv, inputs_manifest.tsv (every source with SHA-256).
"""
import csv
import os
import sys

import pandas as pd

import c4_common as C


def main():
    args = C.parser(__doc__.split("\n")[0]).parse_args()
    repo, results = C.roots(args)
    sys.path.insert(0, os.path.join(repo, "scripts"))
    from genotype_report import load as L
    from genotype_report.common import read_tsv, write_tsv

    if not C.verify_seals(results):
        sys.exit("a seal does not verify: refusing to build inputs")
    out = os.path.join(os.path.abspath(args.work_dir), "inputs")
    if os.path.exists(out):
        sys.exit(f"{out} exists: choose a fresh --work-dir")
    os.makedirs(os.path.join(out, "package"))
    F, H = os.path.join(results, C.FINAL), os.path.join(results, C.HEADLINE)
    manifest = os.path.join(repo, "tests", "genotype_c4_heldout_mosaic_manifest.tsv")
    jobs = {"campaign": os.path.join(results, C.CAMPAIGN, "jobs.tsv"), "fresh": os.path.join(results, C.FRESH_JOBS)}

    # truth-free runs and blocks from the replay prefixes, with the tool's own loader
    runs, blocks, sources = L.load_directory(os.path.join(F, "replay"))
    db, dl = read_tsv(f"{H}/derived_v2/blocks.tsv"), read_tsv(f"{H}/derived_v2/loci.tsv")
    db["run_id"], dl["run_id"] = db.dataset + "/" + db.job, dl.dataset + "/" + dl.job
    if set(dl.run_id) != set(runs.run_id) or len(dl) != len(runs):
        sys.exit("the headline audit's loci do not match the replay runs")
    loc = dl.set_index("run_id")
    for src, dst in (("GENOTYPE_FILTER", "genotype_filter"), ("GENOTYPE_FILTER_REASON", "genotype_filter_reasons"),
                     ("DOSAGE_FILTER", "dosage_filter"), ("DOSAGE_FILTER_REASON", "dosage_filter_reasons"),
                     ("PHASE_FILTER", "phase_filter"), ("PHASE_FILTER_REASON", "phase_filter_reasons"),
                     ("CALL_SUMMARY", "call_summary")):
        if (runs[dst] != "NA").any():
            sys.exit(f"the replay already carries {dst}")
        runs[dst] = runs.run_id.map(loc[src])
    if (runs.call_class != runs.run_id.map(loc.call_class)).any():
        sys.exit("call classes of the headline audit and the replay disagree")
    called = db[db.call_class != "no_call"].set_index(["run_id", "block_id"])
    key = list(zip(blocks.run_id, blocks.block_id))
    if set(key) != set(called.index) or len(key) != len(called):
        sys.exit("the headline audit's called blocks do not match the replay's reported blocks")
    for src, dst in (("GENOTYPE_FILTER", "genotype_filter"), ("GENOTYPE_FILTER_REASON", "genotype_filter_reasons"),
                     ("DOSAGE_FILTER", "dosage_filter"), ("DOSAGE_FILTER_REASON", "dosage_filter_reasons")):
        blocks[dst] = [called.loc[k, src] for k in key]
    for col in ("zygosity", "genotype_support", "dosage_status", "repair_class"):
        bad = sum(a != b for a, b in zip((called.loc[k, col] for k in key), blocks[col]))
        if bad:
            sys.exit(f"{col} differs between the headline audit and the replay at {bad} blocks")
    write_tsv(runs, os.path.join(out, "package", "runs.tsv"))
    write_tsv(blocks, os.path.join(out, "package", "blocks.tsv"))

    pr = read_tsv(f"{F}/package/runs.tsv")
    pr["run_id"] = pr.dataset + "/" + pr.job
    if set(pr.run_id) != set(runs.run_id):
        sys.exit("the sealed package runs do not match the replay runs")
    meta = pr[["run_id", "experiment", "subject", "coverage"]].copy()
    meta["execution_group"] = "final_26401b4 replay"
    write_tsv(meta, os.path.join(out, "run_metadata.tsv"))

    rows, called_runs = [], set(runs.run_id[runs.call_class != "no_call"])
    for ds, jp in jobs.items():
        with open(jp, newline="") as fh:
            job_rows = list(csv.DictReader(fh, delimiter="\t"))
        for job in job_rows:
            rid = f"{ds}/{job['job']}"
            if rid not in called_runs:
                continue
            lg = read_tsv(os.path.join(F, "replay", ds, job["job"] + ".locus_genotype.tsv")).iloc[0]
            if lg.reported_genotype_source == "caller_with_dosage_repair":
                path = os.path.join(F, "eval", ds, job["job"] + ".blocks.tsv")
            else:
                path = C.reroot(job["out_prefix"], results) + ".eval.blocks.tsv"
            if not os.path.isfile(path):
                sys.exit(f"missing evaluator table {path}")
            rows.append({"run_id": rid, "path": os.path.relpath(path, out)})
    write_tsv(pd.DataFrame(rows), os.path.join(out, "evaluation_manifest.tsv"))

    # mosaic plans, switch positions copied verbatim (a flip at block index i switches between blocks i-1 and i)
    plans = []
    with open(manifest) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[0] == "plan" and p[1] != "name":
                plans.append({"subject": p[1], "switch_class": p[2], "switch_before_block_hap1": p[6] or "none",
                              "switch_before_block_hap2": p[7] or "none", "source_a": p[3], "source_b": p[4]})
    write_tsv(pd.DataFrame(plans), os.path.join(out, "switch_metadata.tsv"))

    extra = [f"{H}/derived_v2/blocks.tsv", f"{H}/derived_v2/loci.tsv", f"{F}/package/runs.tsv", manifest] + list(jobs.values())
    rec = [(os.path.relpath(p, results) if p.startswith(results) else os.path.relpath(p, repo), C.sha256(p))
           for p in sorted(set(sources)) + extra]
    write_tsv(pd.DataFrame(rec, columns=["path", "sha256"]), os.path.join(out, "inputs_manifest.tsv"))
    print(f"runs {len(runs)}, blocks {len(blocks)}, evaluator tables {len(rows)}, plans {len(plans)}, sources {len(rec)}")


if __name__ == "__main__":
    main()
