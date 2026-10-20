#!/usr/bin/env python3
"""M11-A family exposure inventory (read-only; existing records only; no genotyping or simulation).

Sources: M8-A families.tsv (per-locus family grouping and donor-level development status), M8-B split.<locus>.tsv (training /
outer-evaluation roles used for prior-parameter scoring), the C4 production campaign jobs.tsv (484 genotyping jobs), the
CYP2D6 pilot reads.tsv (genotyped donors), and the ACOT families named by M8-C/D, M9-C..G and M9-H, plus the development
cells used by the performance and configuration work. Output: one row per locus and family.
  inventory.py <out.tsv>
"""
import collections
import csv
import sys

R = "/Users/davide.bolognini/panvar/tests/results"
FAM = f"{R}/m8a_copying_calibration_20261020/support.full/families.tsv"
SPLIT = f"{R}/m8b_copying_calibration_20261020/runs/inventory/split.{{}}.tsv"
C4JOBS = f"{R}/genotype_c4_production_11fe681/jobs.tsv"
CYPREADS = f"{R}/cyp2d6_shortread_pilot_v1/data/reads.tsv"
ACOT_GENOTYPED = {"HG03209": "M8-C/D", "HG01975": "M8-C/D (timed-out LOO pair completed in M8-D)", "NA19909": "M9-C..G", "NA18747": "M9-C..G",
                  "HG02178": "M9-C..G", "HG00639": "M9-C..G", "NA21106": "M9-C..G", "HG03831": "M9-C..G",
                  "HG02723": "M9-H (L1, Z1)", "HG02074": "M9-H (L2, Z2); M9-I/M10-B..D case", "NA18971": "M9-H (L3)",
                  "HG02083": "M9-H (L4); M9-I/M10-B..D case", "HG00097": "development cell (I-series, P1, S2, M10-A)"}
FAILURE_INSPECTED = {   # donor -> locus-specific failure / outcome inspection beyond aggregate campaign metrics
    ("C4", "NA19650"): "sv:7 allocation failure (M1-B..M6, M10-B..D)", ("C4", "HG03248"): "sv:7 factorized proposal (Batch 3)",
    ("C4", "HG01573"): "census-free branch seeds memory case", ("C4", "HG00096"): "reference donor; development cell; M10 control",
    ("C4", "HG00171"): "P1 sequential-reuse cell",
    ("CYP2D6", "NA19650"): "sv:9 dosage failure, sv:3..sv:5 synonymous case (M1-B, M9-I, M10-B..D)",
    ("CYP2D6", "HG02769"): "development cell (P1, S2, M10-A)",
    ("ACOT", "HG02083"): "M9-H L4 + M9-I/M10 synonymous case", ("ACOT", "HG02074"): "M9-H L2 + M10 few-base case",
    ("ACOT", "NA19909"): "M9-C..G high-confidence sv:5 error attribution", ("ACOT", "NA18747"): "M9-C..G LZO cell",
    ("ACOT", "HG01975"): "M8-C/D copying challenge", ("ACOT", "HG03209"): "M8-C copying challenge",
}


def T(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


def main(out):
    c4_jobs = collections.defaultdict(set)
    for r in T(C4JOBS):
        c4_jobs[r["subject"]].add(f"{r['scheme']}{r['coverage']}")
    cyp_done = {r["sample"] for r in T(CYPREADS)}
    roles = {}
    for loc in ("ACOT", "C4", "CYP2D6"):
        for line in open(SPLIT.format(loc)):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3 and "#" in f[0]:
                roles.setdefault((loc, f[1]), set()).add(f[2])
    rows = []
    for r in T(FAM):
        loc, fam, donors = r["locus"], r["family"], r["donors"].split(",")
        genotyped, inspected = [], []
        for d in donors:
            if loc == "C4" and d in c4_jobs:
                genotyped.append(f"{d}: C4 campaign 11fe681 ({len(c4_jobs[d])} jobs)")
            if loc == "CYP2D6" and d in cyp_done:
                genotyped.append(f"{d}: CYP2D6 pilot")
            if loc == "ACOT" and d in ACOT_GENOTYPED:
                genotyped.append(f"{d}: {ACOT_GENOTYPED[d]}")
            if (loc, d) in FAILURE_INSPECTED:
                inspected.append(f"{d}: {FAILURE_INSPECTED[(loc, d)]}")
        role = ",".join(sorted(roles.get((loc, fam), {"absent"})))
        selection = []
        if role:
            selection.append(f"M8-B prior scoring as {role}")
        if loc == "C4" and genotyped:
            selection.append("C4 headline filters v2 and dosage policy v2 selected on the campaign")
        if loc == "CYP2D6" and genotyped:
            selection.append("CYP2D6 amended profile (958 / 600M) selected on the pilot")
        if loc == "ACOT" and any("M9-C" in g for g in genotyped):
            selection.append("JP1.1 developed on these cells")
        if r["status"] == "REFERENCE_GENOME":
            proposed = "EXCLUDE (reference genome)"
        elif genotyped or inspected:
            proposed = "DEVELOPMENT (genotyped or inspected; never 'held out')"
        elif loc == "ACOT" and role == "eval":
            proposed = "VALIDATION CANDIDATE (never genotyped; assemblies used in M8-B masked prior scoring)"
        else:
            proposed = "PANEL-EXPOSED, NOT GENOTYPED (training panel / prior scoring; not a clean holdout)"
        rows.append(dict(locus=loc, family=fam, donors=r["donors"], templates=r["templates"], m8a_status=r["status"],
                         m8b_role=role, genotyped=" | ".join(genotyped) or "-", failure_inspected=" | ".join(inspected) or "-",
                         parameter_selection=" | ".join(selection) or "-", proposed_role=proposed))
    with open(out, "w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    c = collections.Counter((r["locus"], r["proposed_role"].split(" (")[0]) for r in rows)
    for k in sorted(c):
        print(f"{k[0]}\t{k[1]}\t{c[k]}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
