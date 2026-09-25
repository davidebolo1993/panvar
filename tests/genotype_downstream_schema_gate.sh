#!/usr/bin/env bash
# Downstream schema gate: the campaign classifier and the consolidator on a synthetic five-job campaign.
#
#   tests/genotype_downstream_schema_gate.sh <scratch dir>
#
# A call is decided by the final status, call class and search_conditional -- never by repair_converged,
# which stays false for a tie-branched call. The representation filter annotates calls only: a no-call
# carrying a locus filter is a defect. The package keeps an all-quality and a PASS-only block table, PASS-only
# keeps exactly the blocks whose filter is PASS (fail closed: NOT_EVALUATED and NA are not PASS), the
# high-quality table keeps PASS blocks with READS support, and callability, no-call, LowQual and NOT_EVALUATED
# rates are reported separately, so filtering can never make a call disappear or inflate callability.
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK="${1:?usage: genotype_downstream_schema_gate.sh <scratch dir>}"
rm -rf "$WORK"; mkdir -p "$WORK"
python3 - "$REPO" "$WORK" <<'EOF'
import csv, importlib.util, pathlib, sys
repo, work = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
def load(name):
    spec = importlib.util.spec_from_file_location(name, repo / "scripts" / f"{name}.py")
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module); return module
campaign, consolidator = load("c4_production_campaign"), load("consolidate_c4_production")
failures = []
def same(what, got, want):
    print(("  ok: " if got == want else "  GATE FAILED: ") + f"{what} = {got!r}" + ("" if got == want else f", expected {want!r}"))
    if got != want: failures.append(what)
def write_kv(path, values):
    path.write_text("key\tvalue\n" + "".join(f"{k}\t{v}\n" for k, v in values.items()))
BLOCK_HEADER = ("scheme\tdonor\tmax_mosaic_candidates\tblock_id\tkind\ttruth1_bp\ttruth2_bp\tcalled1_bp\tcalled2_bp\t"
                "unordered_distance\tavailability_floor\texcess\ttruth1_available\ttruth2_available\tgenotype_support\t"
                "called_allele1\tcalled_allele2\ttruth1_closest_allele\ttruth2_closest_allele\ttruth1_closest_distance\t"
                "truth2_closest_distance\tfilter\tquality_filter\tquality_reason\n")
def quality(flt, sup):
    # block_quality_rule_v1, as the caller writes it (the consolidator must read it, never rebuild it)
    reasons = (["representation_conflict"] if flt.startswith("LOWQUAL") else [] if flt == "PASS" else ["representation_not_evaluated"]) + \
              ([] if sup == "READS" else ["no_read_support"])
    return ("PASS", "none") if not reasons else ("LOW_QUALITY", ",".join(reasons))
def blocks(rows):
    return BLOCK_HEADER + "".join(f"loo\tX\t128\t{b}\t{k}\t10\t10\t10\t10\t{d}\t{f}\t{d-f}\ttrue\ttrue\t{sup}\t0\t1\t0\t1\t0\t0\t{flt}\t{quality(flt, sup)[0]}\t{quality(flt, sup)[1]}\n"
                                  for b, k, d, f, sup, flt in rows)
# name: (status keys, files present, eval block rows or None)
jobs = {
    "loo.A.10": ({"status": "COMPLETE", "call_class": "strict_call", "no_call_reason": "none", "phase_support": "PRIOR_ONLY",
                  "search_conditional": "tie_branched", "repair_converged": "false", "repair_stop_reason": "tie_at_best",
                  "locus_filter": "LOWQUAL_REPRESENTATION_CONFLICT", "lowqual_blocks": "sv:7"},
                 [".genotypes.tsv", ".block_genotypes.tsv", ".mosaic_haplotypes.fa", ".repair_rounds.tsv"],
                 [("sv:7", "SV", 32800, 154, "READS"), ("ctx:10", "context", 50, 16, "READS")],
                 {"sv:7": "LOWQUAL_REPRESENTATION_CONFLICT"}),
    "loo.B.10": ({"status": "COMPLETE", "call_class": "strict_call", "no_call_reason": "none", "phase_support": "READS_WITHIN_CANDIDATES",
                  "search_conditional": "none", "repair_converged": "true", "locus_filter": "NOT_EVALUATED",
                  "lowqual_blocks": "none", "not_evaluated_blocks": "ctx:10", "representation_guard_complete": "false"},
                 [".genotypes.tsv", ".block_genotypes.tsv", ".mosaic_haplotypes.fa", ".repair_rounds.tsv"],
                 [("sv:7", "SV", 0, 0, "READS"), ("ctx:10", "context", 9, 5, "READS")], {"ctx:10": "NOT_EVALUATED"}),
    "loo.C.10": ({"status": "CONDITIONAL_GENOTYPE_EQUIVALENT", "call_class": "genotype_equivalent_call", "no_call_reason": "none",
                  "phase_support": "UNRESOLVED", "search_conditional": "tie_branched", "repair_converged": "false",
                  "locus_filter": "NOT_EVALUATED", "lowqual_blocks": "none"},
                 [".genotypes.tsv", ".block_genotypes.tsv", ".equivalence_set.tsv", ".repair_rounds.tsv"],
                 [("sv:7", "SV", 2, 0, "PRIOR_ONLY"), ("ctx:10", "context", 7, 0, "READS")], {"ctx:10": "NA"}),
    "loo.D.10": ({"status": "RESCORE_INCOMPLETE", "call_class": "no_call", "no_call_reason": "divergent_tie",
                  "search_conditional": "tie_branched", "repair_converged": "false", "locus_filter": "NA"},
                 [".repair_rounds.tsv", ".tie_branches.tsv"], None, {}),
    "loo.E.10": ({"status": "REPAIR_INCOMPLETE", "call_class": "no_call", "no_call_reason": "tie_branching_incomplete",
                  "locus_filter": "LOWQUAL_REPRESENTATION_CONFLICT"},
                 [".repair_rounds.tsv"], None, {}),
}
rows = []
for name, (status, files, evalrows, flags) in jobs.items():
    prefix = work / name
    write_kv(pathlib.Path(str(prefix) + ".mosaic_status.tsv"), status)
    for suffix in files: pathlib.Path(str(prefix) + suffix).write_text("x\n")
    if evalrows is not None:
        pathlib.Path(str(prefix) + ".eval.blocks.tsv").write_text(blocks([(b, k, d, f, sup, flags.get(b, "PASS")) for b, k, d, f, sup in evalrows]))
        dist = sum(r[2] for r in evalrows); floor = sum(r[3] for r in evalrows)
        passing = [r for r in evalrows if flags.get(r[0], "PASS") == "PASS"]
        high = [r for r in passing if r[4] == "READS"]
        pdist = sum(r[2] for r in passing); pfloor = sum(r[3] for r in passing)
        write_kv(pathlib.Path(str(prefix) + ".eval.summary.tsv"),
                 {"unordered_distance_total": dist, "availability_floor_total": floor, "excess_over_floor_total": dist - floor,
                  "blocks": 2, "sv_blocks": 1, "blocks_exact": sum(r[2] == 0 for r in evalrows),
                  "sv_blocks_exact": sum(r[2] == 0 for r in evalrows if r[1] == "SV"),
                  "pass_unordered_distance_total": pdist, "pass_excess_over_floor_total": pdist - pfloor,
                  "high_quality_excess_over_floor_total": sum(r[2] - r[3] for r in high)})
    rows.append({"job": name, "campaign": "loo", "subject": name.split(".")[1], "coverage": "10", "out_prefix": str(prefix)})
print("== classifier")
outcomes = {r["job"]: campaign.classify(r)[0] for r in rows}
same("tie-branched COMPLETE with repair_converged false is a strict call", outcomes["loo.A.10"], "strict_call")
same("ordinary COMPLETE", outcomes["loo.B.10"], "strict_call")
same("tie-branched genotype-equivalent call", outcomes["loo.C.10"], "genotype_equivalent_call")
same("divergent_tie is a no-call", outcomes["loo.D.10"], "no_call")
same("a no-call carrying a locus filter is a defect, never a LowQual call", outcomes["loo.E.10"], "defect")
# The package over the four well-formed jobs.
with open(work / "jobs.tsv", "w") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows[:4])
with open(work / "completeness.tsv", "w") as f:
    f.write("job\toutcome\tattempts\twall_seconds\tuser_seconds\tsys_seconds\tpeak_footprint_bytes\tpeak_rss_bytes\n")
    for r in rows[:4]: f.write(f"{r['job']}\t{outcomes[r['job']]}\t1\t1\t1\t0\t1\t1\n")
consolidator.consolidate(work, figure=False)
import pandas as pd
package = work / "package"
runs = pd.read_csv(package / "runs.tsv", sep="\t", dtype=str, keep_default_na=False).set_index("job")
allq = pd.read_csv(package / "blocks_all_quality.tsv", sep="\t", dtype=str)
passq = pd.read_csv(package / "blocks_pass_only.tsv", sep="\t", dtype=str)
highq = pd.read_csv(package / "blocks_high_quality.tsv", sep="\t", dtype=str)
pairs = lambda frame: sorted(map(tuple, frame[["job", "block_id"]].values))
summary = pd.read_csv(package / "summary.tsv", sep="\t").iloc[0]
print("== consolidator")
same("runs keep search_conditional", runs.loc["loo.A.10", "search_conditional"], "tie_branched")
same("runs keep the locus filter of a call", runs.loc["loo.A.10", "locus_filter"], "LOWQUAL_REPRESENTATION_CONFLICT")
same("a no-call has no locus filter", runs.loc["loo.D.10", "locus_filter"], "NA")
same("all-quality blocks: every block of every call", len(allq), 6)
same("PASS-only keeps filter == PASS exactly (LowQual, NOT_EVALUATED and NA out)", pairs(passq),
     [("loo.A.10", "ctx:10"), ("loo.B.10", "sv:7"), ("loo.C.10", "sv:7")])
same("high quality = quality_filter PASS (representation PASS and READS)", pairs(highq), [("loo.A.10", "ctx:10"), ("loo.B.10", "sv:7")])
same("high-quality rows carry quality_filter PASS and reason none", sorted(set(map(tuple, highq[["quality_filter", "quality_reason"]].values))), [("PASS", "none")])
same("a call's unevaluated blocks are named", runs.loc["loo.B.10", "not_evaluated_blocks"], "ctx:10")
same("and the guard is marked incomplete", runs.loc["loo.B.10", "representation_guard_complete"], "false")
same("blocks.tsv is the all-quality table", len(pd.read_csv(package / "blocks.tsv", sep="\t")), 6)
same("callability counts the LowQual call: strict calls", int(summary.strict_calls), 2)
same("genotype-equivalent calls", int(summary.genotype_equivalent_calls), 1)
same("no-calls reported separately", int(summary.no_calls), 1)
same("no-call rate", float(summary.no_call_rate), 0.25)
same("LowQual calls reported separately", int(summary.lowqual_calls), 1)
same("LowQual rate given a call", round(float(summary.lowqual_call_rate_given_call), 6), round(1 / 3, 6))
same("tie-branched calls", int(summary.tie_branched_calls), 2)
same("block filter rate given a call", round(float(summary.block_filter_rate_given_call), 6), round(1 / 6, 6))
same("NOT_EVALUATED calls reported separately", int(summary.not_evaluated_calls), 2)
same("NOT_EVALUATED call rate given a call", round(float(summary.not_evaluated_call_rate_given_call), 6), round(2 / 3, 6))
same("NOT_EVALUATED blocks (NOT_EVALUATED and NA)", int(summary.not_evaluated_blocks), 2)
same("all-quality denominator", int(summary.all_quality_blocks), 6)
same("representation-PASS retained fraction", round(float(summary.representation_pass_fraction), 6), 0.5)
same("high-quality retained fraction", round(float(summary.high_quality_fraction), 6), round(2 / 6, 6))
same("all-quality excess keeps every block", float(summary.excess_total_given_call), float((32800 - 154) + (50 - 16) + 0 + (9 - 5) + 2 + 7))
same("PASS-only excess: PASS blocks only", float(summary.excess_total_pass_only_given_call), float((50 - 16) + 0 + 2))
same("high-quality excess: PASS and READS only", float(summary.excess_total_high_quality_given_call), float((50 - 16) + 0))
print("downstream schema gate: " + ("PASS" if not failures else "FAILED"))
sys.exit(1 if failures else 0)
EOF
