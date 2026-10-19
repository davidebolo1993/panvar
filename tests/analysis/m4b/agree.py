#!/usr/bin/env python3
"""M4-B caller/diagnostic agreement and incompatible-specification detection, through genotype_probe score-set.

agree:    score-set scores each arm's saved winning equivalence-set representatives under that arm's specification (the
          distance arm with --copying-transition distance_ref_core_v1_experimental --copying-length-bp <frozen lambda>), with
          their saved totals declared as expectations: score-set refuses unless every total reproduces (1e-6) and its vector
          digest equals the run's recorded transition_vector_sha256. Exit 0 iff both arms reproduce.
incompat: three specifications that must be REFUSED against the distance run: constant mode (mode mismatch), lambda * (1 +
          1e-6) (digest mismatch), and the distance mode with an explicit --switch-probability (conflict). Exit 0 iff all
          three are refused with the expected message.

  agree.py agree|incompat <probe> <const prefix> <dist prefix> <gfa> <reference> <curation|none> <reads dir> <lambda> <out dir>
"""
import csv
import subprocess
import sys

kind, probe, cpre, dpre, gfa, ref, cur, reads, lam, out = sys.argv[1:11]


def reps(prefix, path):
    rows = [r for r in csv.DictReader(open(prefix + ".equivalence_set.tsv"), delimiter="\t") if r["representative"] == "true"]
    with open(path, "w") as f:
        for r in rows:
            f.write(f"eq{r['rank']}_{r['sequence_candidate_id']}\t{r['hap1_alleles']}\t{r['hap2_alleles']}\t{r['score']}\n")
    return len(rows)


def score_set(prefix, cands, tag, extra):
    cmd = [probe, "score-set", "--out", f"{out}/{tag}", "--model-gfa", gfa, "--reference", ref, "--donor", "NA", "--run", prefix,
           "--reads1", f"{reads}/reads_1.fq", "--reads2", f"{reads}/reads_2.fq", "--library", f"{reads}/library.txt", "--threads", "4",
           "--candidates", cands] + ([] if cur == "none" else ["--panel-curation", cur]) + extra
    p = subprocess.run(cmd, capture_output=True, text=True)
    open(f"{out}/{tag}.log", "w").write(" ".join(cmd) + "\n" + p.stdout + p.stderr)
    return p.returncode, p.stdout + p.stderr


DIST = ["--copying-transition", "distance_ref_core_v1_experimental", "--copying-length-bp", lam]
if kind == "agree":
    nc, nd = reps(cpre, f"{out}/const.tsv"), reps(dpre, f"{out}/dist.tsv")
    rc, tc = score_set(cpre, f"{out}/const.tsv", "const", [])
    rd, td = score_set(dpre, f"{out}/dist.tsv", "dist", DIST)
    ok = rc == 0 and rd == 0 and "matches the run" in td and "expectations reproduced" in tc and "expectations reproduced" in td
    print(f"constant arm: {nc} representative(s), exit {rc}; distance arm: {nd} representative(s), exit {rd}")
    print([l for l in (tc + td).split("\n") if "copying_transition" in l or "expectations" in l or "REFUSED" in l])
    sys.exit(0 if ok else 1)
reps(dpre, f"{out}/dist.tsv")
checks = [("constant spec on the distance run", [], "the run used copying_transition distance_ref_core_v1_experimental"),
          ("lambda*(1+1e-6)", ["--copying-transition", "distance_ref_core_v1_experimental", "--copying-length-bp", repr(float(lam) * (1 + 1e-6))],
           "differs from the run's"),
          ("distance + --switch-probability", DIST + ["--switch-probability", "0.01"], "conflicts with --switch-probability")]
ok = True
for i, (name, extra, msg) in enumerate(checks):
    code, text = score_set(dpre, f"{out}/dist.tsv", f"incompat{i}", extra)
    refused = code != 0 and msg in text
    ok &= refused
    print(f"{name}: exit {code}, refused with expected message: {refused}")
sys.exit(0 if ok else 1)
