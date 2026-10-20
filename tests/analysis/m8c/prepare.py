#!/usr/bin/env python3
"""M8-C prepare (frozen before any read generation or genotype outcome).

  1. donors: from M8-B's 30 ACOT outer-evaluation families (split.ACOT.tsv), the first two by sha256("M8-C:" + family)
     whose donor is a complete diploid (exactly two templates, #1 and #2) -- no outcome, score or complexity enters
  2. calling panels (curation manifests on the fixed full graph acot.gfa): LOO excludes every outer-evaluation template;
     LZO excludes every outer-evaluation template except the queried donor's two
  3. simulator lineage re-check: the current build regenerates the archived C4 NA19650 30x source reads byte for byte
  4. reads: simulate --donor D --depth 30 --seed 20260916 --seed-salt D (ACOT pilot rule), subsample to 5x --name-salt D
  5. profiles: the M4-B/P1 composition of acot_5.tsv (full overlay) plus copying_transition constant_v1 and
     switch_probability r in {0.01, 0.1}; the two arms must differ only in the profile name and switch_probability
  6. executables: SHA-256 of panvar and genotype_probe

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
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
log = f"{d}/commands.log"


def run(cmd):
    with open(log, "a") as f:
        f.write("$ " + " ".join(cmd) + "\n")
        f.flush()
        return subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT).returncode


rows = [l.rstrip("\n").split("\t") for l in open(SPLIT)]
eval_t = [r[0] for r in rows if r[2] == "eval"]
fam_of = {r[0]: r[1] for r in rows}
fams = sorted({fam_of[t] for t in eval_t}, key=lambda f: hashlib.sha256(("M8-C:" + f).encode()).hexdigest())
chosen = []
for f in fams:
    ts = [t for t in eval_t if fam_of[t] == f]
    donors = {t.split("#")[0] for t in ts}
    if len(donors) == 1 and sorted(t.split("#")[1] for t in ts) == ["1", "2"]:
        chosen.append((f, donors.pop(), ts))
    if len(chosen) == 2:
        break
assert len(chosen) == 2
out = [f"split_sha256\t{sha(SPLIT)}", f"acot_gfa_sha256\t{sha(ACOT_GFA)}", f"panvar_sha256\t{sha(PANVAR)}", f"genotype_probe_sha256\t{sha(PROBE)}"]
for i, (f, dn, ts) in enumerate(chosen, 1):
    out.append(f"donor_{i}\t{dn}\tfamily {f}\ttemplates {','.join(sorted(ts))}")
    open(f"{d}/curation.LOO.{dn}.txt", "w").write("# M8-C LOO calling panel: every outer-evaluation template excluded\n" + "\n".join(eval_t) + "\n")
    open(f"{d}/curation.LZO.{dn}.txt", "w").write("# M8-C LZO calling panel: outer-evaluation templates excluded except the queried donor's two\n"
                                                    + "\n".join(t for t in eval_t if t not in ts) + "\n")
# simulator lineage re-check (current build vs archived C4 NA19650 source)
o = f"{d}/lineage_c4_NA19650"
os.makedirs(o)
code = run([PROBE, "simulate", "--gfa", f"{R}/genotype_c4_production_11fe681/graphs/c4.gfa", "--reference", "HG00096#1#haplotype1-0000024:31848049-32074186",
            "--donor", "NA19650", "--depth", "30", "--read-length", "100", "--seed", "20260916", "--out", o])
src = f"{R}/genotype_c4_production_11fe681/reads/NA19650/source"
lineage_ok = code == 0 and sha(f"{o}/reads_1.fq") == sha(f"{src}/reads_1.fq") and sha(f"{o}/reads_2.fq") == sha(f"{src}/reads_2.fq")
out.append(f"simulator_lineage_recheck\t{'PASS' if lineage_ok else 'FAIL'} (current build regenerates archived C4 NA19650 30x source byte for byte)")
assert lineage_ok, "simulator lineage re-check failed"
for f, dn, ts in chosen:
    s30, s5 = f"{d}/reads/{dn}/cov30", f"{d}/reads/{dn}/cov5"
    os.makedirs(s30)
    assert run([PROBE, "simulate", "--gfa", ACOT_GFA, "--reference", ACOT_REF, "--donor", dn, "--depth", "30", "--read-length", "100",
                "--seed", "20260916", "--seed-salt", dn, "--out", s30]) == 0
    assert run([PROBE, "subsample", "--reads1", f"{s30}/reads_1.fq", "--reads2", f"{s30}/reads_2.fq", "--library", f"{s30}/library.txt",
                "--source-coverage", "30", "--coverage", "5", "--name-salt", dn, "--out", s5]) == 0
    out.append(f"reads_{dn}\treads_1 {sha(s5 + '/reads_1.fq')[:16]} reads_2 {sha(s5 + '/reads_2.fq')[:16]} library {sha(s5 + '/library.txt')[:16]}")
compose_profile.compose(f"{d}/profile.base.tsv", "m8c-base-acot", BASE, "full", [])
base = open(f"{d}/profile.base.tsv").read().rstrip("\n").split("\n")
for tag, r in (("r0.01", "0.01"), ("r0.1", "0.1")):
    lines = [l if not l.startswith("profile\t") else f"profile\tm8c-acot-{tag}" for l in base]
    lines += ["# M8-C: the only scientific difference between the two arms", "copying_transition\tconstant_v1", f"switch_probability\t{r}"]
    open(f"{d}/profile.{tag}.tsv", "w").write("\n".join(lines) + "\n")
a, b = open(f"{d}/profile.r0.01.tsv").read().split("\n"), open(f"{d}/profile.r0.1.tsv").read().split("\n")
diff = [(x, y) for x, y in zip(a, b) if x != y]
ok = len(a) == len(b) and sorted(x.split("\t")[0] for x, _ in diff) == ["profile", "switch_probability"]
out.append(f"profile_arm_difference\t{'PASS' if ok else 'FAIL'}: {diff}")
assert ok
open(f"{d}/prepare.tsv", "w").write("\n".join(out) + "\n")
print("\n".join(out))
