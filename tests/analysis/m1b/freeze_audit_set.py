#!/usr/bin/env python3
"""M1-B: freeze the structural audit set (at most four existing cells), before any execution.

  * C4: the LOO sample with the largest RETAINED sv:7 excess at 30x in the frozen C4 package
    (genotype_c4_production_11fe681_dosage_repair/final_26401b4/package/blocks.tsv): experiment == loo, block_id == sv:7,
    coverage == 30, recommended == True; sorted by excess descending, then subject ascending (the documented
    tie-break). The same subject is taken at 10x and 30x.
  * CYP2D6: loo.NA19650 at 5x (the known sv:9 error) and its 30x twin, from the sealed pilot (run set a5).

Inputs are verified (existence and SHA-256) and recorded with their provenance; a missing input is a reported gap and is
never substituted or rebuilt. Reads only tables and input files; runs nothing.

  freeze_audit_set.py OUT_DIR                          (writes SELECTION.tsv and AUDIT_SET.tsv)
  freeze_audit_set.py OUT_DIR FROZEN_DIR PANVAR BASELINE_MANIFEST
      (verification unit: rewrite into OUT_DIR, require byte equality with the frozen files, and require the panvar binary
       SHA-256 to equal the M1-A baseline manifest's binary_panvar_sha256)
"""
import csv
import hashlib
import os
import sys

C4 = "/Users/davide.bolognini/panvar/tests/results/genotype_c4_production_11fe681"
C4PKG = "/Users/davide.bolognini/panvar/tests/results/genotype_c4_production_11fe681_dosage_repair/final_26401b4/package"
CYP = "/Users/davide.bolognini/panvar/tests/results/cyp2d6_shortread_pilot_v1"
CYP_CURATION = "/Users/davide.bolognini/panvar/tests/genotype_cyp2d6_curation.tsv"
CYP_CURATION_OF_RECORD = "/Users/davide.bolognini/panvar-worktrees/genotype-cyp2d6-pilot/tests/genotype_cyp2d6_curation.tsv"


def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def select_c4():
    rows = [r for r in csv.DictReader(open(f"{C4PKG}/blocks.tsv"), delimiter="\t")
            if r["experiment"] == "loo" and r["block_id"] == "sv:7" and r["coverage"] == "30" and r["recommended"] == "True"]
    rows.sort(key=lambda r: (-float(r["excess"]), r["subject"]))
    assert rows, "no retained LOO sv:7 row at 30x"
    top = rows[0]
    tie = len(rows) > 1 and float(rows[1]["excess"]) == float(top["excess"])
    return rows, top, tie


def main(out):
    os.makedirs(out, exist_ok=True)
    rows, top, tie = select_c4()
    with open(f"{out}/SELECTION.tsv", "w") as f:
        f.write("rank\tjob\tsubject\tcoverage\tblock_id\tdistance\tfloor\texcess\tcalled_allele1\tcalled_allele2\t"
                "call_class\tdosage_status\tfilter\tquality_filter\trecommended\n")
        for i, r in enumerate(rows[:10], 1):
            f.write("\t".join([str(i)] + [r[k] for k in ("job", "subject", "coverage", "block_id", "distance", "floor", "excess",
                                                         "called_allele1", "called_allele2", "call_class", "dosage_status",
                                                         "filter", "quality_filter", "recommended")]) + "\n")
        f.write(f"# rule: experiment=loo, block_id=sv:7, coverage=30, recommended=True; sort excess desc, subject asc\n"
                f"# rows matching: {len(rows)}; tie at the top: {'YES' if tie else 'no'}\n")
    assert not tie, "a tie at the top needs the documented tie-break reported before freezing"
    s = top["subject"]
    cases = []
    for cov in (10, 30):
        rd = f"{C4}/reads/{s}/cov{cov}"
        cases.append(dict(case=f"c4_{cov}", locus="C4", scenario="loo", sample=s, coverage=cov,
                          role="largest retained sv:7 excess at 30x (deterministic rule)" if cov == 30 else "the same sample at 10x",
                          target_block="sv:7", model_graph=f"{C4}/graphs/{s}.loo.gfa", truth_graph=f"{C4}/graphs/c4.gfa",
                          curation="none", reference="HG00096#1#haplotype1-0000024:31848049-32074186",
                          reads1=f"{rd}/reads_1.fq", reads2=f"{rd}/reads_2.fq", library=f"{rd}/library.txt",
                          historical_run=f"{C4}/loo/{s}/call.{cov}",
                          historical_eval=f"{C4}/loo/{s}/call.{cov}.eval",
                          historical_binary="mosaic-11fe681/build/panvar (campaign 11fe681; package final_26401b4 replay)",
                          historical_profile="/Users/davide.bolognini/panvar/profiles/genotype_c4_v6_campaign_local.tsv",
                          baseline_base_profile="/Users/davide.bolognini/panvar/tests/results/genotype_perf_bounded_20261015/inputs/profiles/c4_30.tsv"))
    for cov in (5, 30):
        rd = f"{CYP}/data/reads/NA19650/cov{cov}"
        cases.append(dict(case=f"cyp_{cov}", locus="CYP2D6", scenario="loo", sample="NA19650", coverage=cov,
                          role="known LOO sv:9 structural error (M1-A)" if cov == 5 else "the nearly correct 30x twin",
                          target_block="sv:9", model_graph=f"{CYP}/data/graphs/loo.NA19650.gfa",
                          truth_graph=f"{CYP}/data/graphs/cyp2d6.gfa", curation=CYP_CURATION,
                          reference="grch38#1#chr22:42031864-42245566",
                          reads1=f"{rd}/reads_1.fq", reads2=f"{rd}/reads_2.fq", library=f"{rd}/library.txt",
                          historical_run=f"{CYP}/data/runs_a5/loo.NA19650.{cov}.original",
                          historical_eval=f"{CYP}/data/eval_a5/loo.NA19650.{cov}.raw",
                          historical_binary="cyp2d6_shortread_pilot_v1/bin/panvar (pilot a5)",
                          historical_profile=f"{CYP}/resource_profile_a4.tsv",
                          baseline_base_profile=f"/Users/davide.bolognini/panvar/tests/results/genotype_perf_bounded_20261015/inputs/profiles/cyp_{cov}.tsv"))
    gaps = []
    for c in cases:
        for key in ("model_graph", "truth_graph", "reads1", "reads2", "library", "historical_profile", "baseline_base_profile"):
            if not os.path.exists(c[key]):
                gaps.append(f"{c['case']}: {key} missing ({c[key]})")
                c[key + "_sha256"] = "MISSING"
            else:
                c[key + "_sha256"] = sha(c[key])
        for suffix in (".mosaic_status.tsv", ".mosaic_candidates.tsv", ".mosaic_scores.tsv"):
            if not os.path.exists(c["historical_run"] + suffix):
                gaps.append(f"{c['case']}: historical output {suffix} missing")
        for suffix in (".summary.tsv", ".blocks.tsv"):
            if not os.path.exists(c["historical_eval"] + suffix):
                gaps.append(f"{c['case']}: historical evaluation {suffix} missing")
        c["curation_sha256"] = sha(c["curation"]) if c["curation"] != "none" else "none"
        c["read_provenance"] = open(c["library"]).readline().strip().lstrip("# ")
    if cases[2]["curation_sha256"] != sha(CYP_CURATION_OF_RECORD):
        gaps.append("CYP2D6 curation differs from the pilot's curation of record")
    keys = list(cases[0].keys())
    with open(f"{out}/AUDIT_SET.tsv", "w") as f:
        f.write("\t".join(keys) + "\n")
        for c in cases:
            f.write("\t".join(str(c[k]) for k in keys) + "\n")
    print(f"C4 selection: {top['job']} (excess {top['excess']}); next: "
          f"{rows[1]['job'] if len(rows) > 1 else '-'} ({rows[1]['excess'] if len(rows) > 1 else '-'})")
    print("gaps: " + ("; ".join(gaps) if gaps else "none"))
    return 1 if gaps else 0


def verify(out, frozen, panvar, manifest):
    code = main(out)
    for name in ("SELECTION.tsv", "AUDIT_SET.tsv"):
        same = open(f"{out}/{name}").read() == open(f"{frozen}/{name}").read()
        print(f"{name}: {'identical to the frozen file' if same else 'DIFFERS from the frozen file'}")
        code |= 0 if same else 1
    want = dict(line.rstrip("\n").split("\t", 1) for line in open(manifest) if "\t" in line)["binary_panvar_sha256"]
    got = sha(panvar)
    print(f"panvar binary {got}: {'equals' if got == want else 'DIFFERS from'} the baseline manifest")
    return code | (0 if got == want else 1)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]) if len(sys.argv) == 2 else verify(*sys.argv[1:5]))
