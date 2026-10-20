#!/usr/bin/env python3
"""M9-C prepare (frozen before any read generation or genotype outcome).

  1. families: M8-B's 30 ACOT outer-evaluation families (split.ACOT.tsv) minus the two genotyped in M8-C/M8-D
     (SL17 = HG03209, PEL25 = HG01975); the first six by sha256("M9-C:" + family); per family the single complete diploid
     donor (exactly #1 and #2); every one of the 30 families is single-donor complete diploid, checked here
  2. calling panels (curation on the fixed full graph acot.gfa, M8-C recipe): LOO excludes every outer-evaluation template;
     LZO (first two selected donors only) excludes every outer-evaluation template except the queried donor's two
  3. reads: simulate --donor D --depth 30 --read-length 100 --seed 20260916 --seed-salt D, subsample to 5x --name-salt D
     (M8-C recipe); the simulator lineage re-check is M8-C's (same genotype_probe binary, re-hashed here)
  4. profile: M8-C's r = 0.01 profile rule (P1 composition of acot_5.tsv, full overlay; constant_v1, r = 0.01) plus
     output_schema sequence_first_v2 and repair_joint_search off; it must differ from M8-C's profile.r0.01.tsv only in the
     profile name and those two added keys
  5. executables and inputs: SHA-256

  prepare.py <unit dir>
"""
import hashlib
import os
import subprocess
import sys

d = sys.argv[1]
R = "/Users/davide.bolognini/panvar/tests/results"
WT = "/Users/davide.bolognini/panvar-worktrees/s1-d-panel-references"
sys.path.insert(0, f"{WT}/tests/analysis/p1")
import compose_profile  # noqa: E402

PROBE, PANVAR = f"{WT}/build/genotype_probe", f"{WT}/build/panvar"
ACOT_GFA = f"{R}/genotype_perf_bounded_20261015/inputs/acot/acot.gfa"
ACOT_REF = "GRCh38#0#chr14:73416151-73668669"
SPLIT = f"{R}/m8b_copying_calibration_20261020/runs/inventory/split.ACOT.tsv"
BASE = f"{R}/genotype_perf_bounded_20261015/inputs/profiles/acot_5.tsv"
M8C_P = f"{R}/m8c_copying_challenge_20261020/runs/prepare.a1"
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
log = f"{d}/commands.log"


def run(cmd):
    with open(log, "a") as f:
        f.write("$ " + " ".join(cmd) + "\n")
        f.flush()
        return subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT).returncode


rec = {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(f"{M8C_P}/prepare.tsv")}
assert sha(PANVAR) == rec["panvar_sha256"], "panvar differs from the M8-C executable"
assert sha(PROBE) == rec["genotype_probe_sha256"], "genotype_probe differs from M8-C (its simulator lineage check would not carry over)"
assert sha(SPLIT) == rec["split_sha256"] and sha(ACOT_GFA) == rec["acot_gfa_sha256"]
rows = [l.rstrip("\n").split("\t") for l in open(SPLIT)]
eval_t = [r[0] for r in rows if r[2] == "eval"]
fam_of = {r[0]: r[1] for r in rows}
fams = {}
for t in eval_t:
    fams.setdefault(fam_of[t], []).append(t)
assert len(fams) == 30 and all(len({t.split("#")[0] for t in ts}) == 1 and sorted(t.split("#")[1] for t in ts) == ["1", "2"] for ts in fams.values())
used = {"SL17": "HG03209", "PEL25": "HG01975"}
for f, dn in used.items():
    assert {t.split("#")[0] for t in fams[f]} == {dn}, f"{f} is not {dn}"
pool = sorted((f for f in fams if f not in used), key=lambda f: hashlib.sha256(("M9-C:" + f).encode()).hexdigest())
assert len(pool) == 28
chosen = [(f, fams[f][0].split("#")[0], sorted(fams[f])) for f in pool[:6]]
out = [f"split_sha256\t{sha(SPLIT)}", f"acot_gfa_sha256\t{sha(ACOT_GFA)}", f"panvar_sha256\t{sha(PANVAR)}", f"genotype_probe_sha256\t{sha(PROBE)}",
       f"selection_pool\t28 families ({','.join(pool)})", "simulator_lineage\tM8-C prepare.a1 PASS (same genotype_probe SHA-256)"]
open(f"{d}/curation.LOO.txt", "w").write("# M9-C LOO calling panel (M8-C recipe): every outer-evaluation template excluded\n" + "\n".join(eval_t) + "\n")
assert open(f"{d}/curation.LOO.txt").read().split("\n")[1:] == open(f"{M8C_P}/curation.LOO.HG01975.txt").read().split("\n")[1:], "LOO panel differs from M8-C"
out.append(f"curation_LOO\t{sha(d + '/curation.LOO.txt')}")
for i, (f, dn, ts) in enumerate(chosen, 1):
    out.append(f"donor_{i}\t{dn}\tfamily {f}\ttemplates {','.join(ts)}")
    if i <= 2:
        open(f"{d}/curation.LZO.{dn}.txt", "w").write("# M9-C LZO calling panel: outer-evaluation templates excluded except the queried donor's two\n"
                                                        + "\n".join(t for t in eval_t if t not in ts) + "\n")
        out.append(f"curation_LZO_{dn}\t{sha(d + '/curation.LZO.' + dn + '.txt')}")
for f, dn, ts in chosen:
    s30, s5 = f"{d}/reads/{dn}/cov30", f"{d}/reads/{dn}/cov5"
    os.makedirs(s30)
    os.makedirs(s5)
    assert run([PROBE, "simulate", "--gfa", ACOT_GFA, "--reference", ACOT_REF, "--donor", dn, "--depth", "30", "--read-length", "100",
                "--seed", "20260916", "--seed-salt", dn, "--out", s30]) == 0
    assert run([PROBE, "subsample", "--reads1", f"{s30}/reads_1.fq", "--reads2", f"{s30}/reads_2.fq", "--library", f"{s30}/library.txt",
                "--source-coverage", "30", "--coverage", "5", "--name-salt", dn, "--out", s5]) == 0
    out.append(f"reads_{dn}\treads_1 {sha(s5 + '/reads_1.fq')} reads_2 {sha(s5 + '/reads_2.fq')} library {sha(s5 + '/library.txt')}")
compose_profile.compose(f"{d}/profile.base.tsv", "m8c-base-acot", BASE, "full", [])
base = open(f"{d}/profile.base.tsv").read().rstrip("\n").split("\n")
lines = [l if not l.startswith("profile\t") else "profile\tm8c-acot-r0.01" for l in base]
lines += ["# M8-C: the only scientific difference between the two arms", "copying_transition\tconstant_v1", "switch_probability\t0.01"]
assert "\n".join(lines) + "\n" == open(f"{M8C_P}/profile.r0.01.tsv").read(), "the re-derived M8-C r0.01 profile differs"
lines = [l if not l.startswith("profile\t") else "profile\tm9c-acot-v2-r0.01" for l in lines]
lines += ["# M9-C: reporting contract and experiment switch (declared)", "output_schema\tsequence_first_v2", "repair_joint_search\toff"]
open(f"{d}/profile.tsv", "w").write("\n".join(lines) + "\n")
out.append(f"profile\t{sha(d + '/profile.tsv')} (M8-C profile.r0.01 + output_schema sequence_first_v2 + repair_joint_search off)")
open(f"{d}/prepare.tsv", "w").write("\n".join(out) + "\n")
print("\n".join(out))
