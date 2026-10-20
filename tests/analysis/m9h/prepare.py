#!/usr/bin/env python3
"""M9-H prepare (frozen before any read generation or outcome).

  1. families: M8-B's 30 ACOT outer-evaluation families (split.ACOT.tsv) minus every family already used for genotyping or
     confidence inspection: SL17 (HG03209) and PEL25 (HG01975) in M8-C/M8-D; 2427, NA18747, HG02178, PR03, NA21106, BD27 in
     M9-C..M9-G. The remainder is counted and recorded (expected 22). Four by ascending sha256("M9-H:" + family).
  2. donor: per family, among complete diploid donors (exactly #1 and #2), the smallest sha256("M9-H:donor:" + donor);
     every family's eligible set is recorded
  3. cells: LOO for all four; LZO for the first two selected families (fixed now, before outcomes)
  4. panels (M8-C recipe, fixed full-panel graph acot.gfa): LOO excludes every outer-evaluation template (asserted equal
     to M9-C's list); LZO excludes them except the queried donor's two
  5. reads (M8-C recipe): simulate 30x (seed 20260916, salt donor), subsample 5x (name-salt donor). genotype_probe changed
     in M9-G (a new subcommand), so the simulator lineage is re-checked here: regenerating M9-C's NA19909 30x and 5x reads
     must reproduce M9-C's recorded SHA-256 byte for byte
  6. profiles: off = M9-C's profile.tsv (the validated M9-G baseline), on = M9-G's profile.on.tsv; hashes recorded and the
     only difference asserted to be the profile name and repair_junction_proposals experimental_v1
  7. executables: SHA-256 asserted equal to M9-G's prepare record

  prepare.py <unit dir>
"""
import hashlib
import os
import subprocess
import sys

d = sys.argv[1]
R = "/Users/davide.bolognini/panvar/tests/results"
WT = "/Users/davide.bolognini/panvar-worktrees/s1-d-panel-references"
PROBE, PANVAR = f"{WT}/build/genotype_probe", f"{WT}/build/panvar"
ACOT_GFA = f"{R}/genotype_perf_bounded_20261015/inputs/acot/acot.gfa"
ACOT_REF = "GRCh38#0#chr14:73416151-73668669"
SPLIT = f"{R}/m8b_copying_calibration_20261020/runs/inventory/split.ACOT.tsv"
C9 = f"{R}/m9c_acot_reliability_20261020/runs/prepare"
G9 = f"{R}/m9g_junction_integration_20261020/runs/prepare"
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
H = lambda s: hashlib.sha256(s.encode()).hexdigest()
log = f"{d}/commands.log"


def run(cmd):
    with open(log, "a") as f:
        f.write("$ " + " ".join(cmd) + "\n")
        f.flush()
        return subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT).returncode


g = {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(f"{G9}/prepare.tsv") if not l.startswith("item")}
assert sha(PANVAR) == g["panvar"] and sha(PROBE) == g["genotype_probe"], "executables differ from M9-G"
out = [f"panvar_sha256\t{sha(PANVAR)}", f"genotype_probe_sha256\t{sha(PROBE)}", f"acot_gfa_sha256\t{sha(ACOT_GFA)}", f"split_sha256\t{sha(SPLIT)}"]
rows = [l.rstrip("\n").split("\t") for l in open(SPLIT)]
eval_t = [r[0] for r in rows if r[2] == "eval"]
fam_of = {r[0]: r[1] for r in rows}
fams = {}
for t in eval_t:
    fams.setdefault(fam_of[t], []).append(t)
used = {"SL17": "M8-C/D HG03209", "PEL25": "M8-C/D HG01975", "2427": "M9-C..G NA19909", "NA18747": "M9-C..G NA18747",
        "HG02178": "M9-C..G HG02178", "PR03": "M9-C..G HG00639", "NA21106": "M9-C..G NA21106", "BD27": "M9-C..G HG03831"}
assert set(used) <= set(fams)
remaining = sorted((f for f in fams if f not in used), key=lambda f: H("M9-H:" + f))
out.append(f"families_total\t{len(fams)}")
out.append(f"families_excluded\t{len(used)} ({'; '.join(f'{k}={v}' for k, v in used.items())})")
out.append(f"families_remaining\t{len(remaining)}")
assert len(remaining) == 22
chosen = []
for f in remaining[:4]:
    donors = {}
    for t in fams[f]:
        donors.setdefault(t.split("#")[0], []).append(t)
    eligible = sorted((dn for dn, ts in donors.items() if sorted(t.split("#")[1] for t in ts) == ["1", "2"]), key=lambda dn: H("M9-H:donor:" + dn))
    assert eligible, f"family {f} has no complete diploid donor"
    dn = eligible[0]
    chosen.append((f, dn, sorted(donors[dn])))
    out.append(f"family_{len(chosen)}\t{f}\tdonor {dn}\teligible {','.join(eligible)}\ttemplates {','.join(sorted(donors[dn]))}")
open(f"{d}/curation.LOO.txt", "w").write("# M9-H LOO calling panel (M8-C recipe): every outer-evaluation template excluded\n" + "\n".join(eval_t) + "\n")
assert open(f"{d}/curation.LOO.txt").read().split("\n")[1:] == open(f"{C9}/curation.LOO.txt").read().split("\n")[1:], "LOO list differs from M9-C"
out.append(f"curation_LOO\t{sha(d + '/curation.LOO.txt')}")
for f, dn, ts in chosen[:2]:
    open(f"{d}/curation.LZO.{dn}.txt", "w").write("# M9-H LZO calling panel: outer-evaluation templates excluded except the queried donor's two\n"
                                                    + "\n".join(t for t in eval_t if t not in ts) + "\n")
    out.append(f"curation_LZO_{dn}\t{sha(d + '/curation.LZO.' + dn + '.txt')}")
# simulator lineage re-check against M9-C (NA19909)
lin30, lin5 = f"{d}/lineage/NA19909/cov30", f"{d}/lineage/NA19909/cov5"
os.makedirs(lin30)
os.makedirs(lin5)
assert run([PROBE, "simulate", "--gfa", ACOT_GFA, "--reference", ACOT_REF, "--donor", "NA19909", "--depth", "30", "--read-length", "100",
            "--seed", "20260916", "--seed-salt", "NA19909", "--out", lin30]) == 0
assert run([PROBE, "subsample", "--reads1", f"{lin30}/reads_1.fq", "--reads2", f"{lin30}/reads_2.fq", "--library", f"{lin30}/library.txt",
            "--source-coverage", "30", "--coverage", "5", "--name-salt", "NA19909", "--out", lin5]) == 0
c9 = {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(f"{C9}/prepare.tsv")}
want = c9["reads_NA19909"]
got = f"reads_1 {sha(lin5 + '/reads_1.fq')} reads_2 {sha(lin5 + '/reads_2.fq')} library {sha(lin5 + '/library.txt')}"
lineage = got == want
out.append(f"simulator_lineage_recheck\t{'PASS' if lineage else 'FAIL'} (NA19909 5x regenerated vs M9-C prepare)")
assert lineage, "simulator lineage re-check failed"
for f, dn, ts in chosen:
    s30, s5 = f"{d}/reads/{dn}/cov30", f"{d}/reads/{dn}/cov5"
    os.makedirs(s30)
    os.makedirs(s5)
    assert run([PROBE, "simulate", "--gfa", ACOT_GFA, "--reference", ACOT_REF, "--donor", dn, "--depth", "30", "--read-length", "100",
                "--seed", "20260916", "--seed-salt", dn, "--out", s30]) == 0
    assert run([PROBE, "subsample", "--reads1", f"{s30}/reads_1.fq", "--reads2", f"{s30}/reads_2.fq", "--library", f"{s30}/library.txt",
                "--source-coverage", "30", "--coverage", "5", "--name-salt", dn, "--out", s5]) == 0
    out.append(f"reads_{dn}\treads_1 {sha(s5 + '/reads_1.fq')} reads_2 {sha(s5 + '/reads_2.fq')} library {sha(s5 + '/library.txt')}")
off, on = open(f"{C9}/profile.tsv").read(), open(f"{G9}/profile.on.tsv").read()
open(f"{d}/profile.off.tsv", "w").write(off)
open(f"{d}/profile.on.tsv", "w").write(on)
diff = sorted(set(on.splitlines()) ^ set(off.splitlines()))
assert all(l.startswith(("profile\t", "repair_junction_proposals\t", "#")) for l in diff), diff
out += [f"profile_off\t{sha(d + '/profile.off.tsv')} (M9-C prepare/profile.tsv)", f"profile_on\t{sha(d + '/profile.on.tsv')} (M9-G profile.on.tsv)",
        f"profile_difference\t{' | '.join(diff)}"]
open(f"{d}/prepare.tsv", "w").write("\n".join(out) + "\n")
print("\n".join(out))
