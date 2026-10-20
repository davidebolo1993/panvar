#!/usr/bin/env python3
"""M7-B provenance: executable identities and byte-identical regeneration of each generator lineage.

  provenance.py <unit dir>
Writes executables.tsv, provenance.tsv (one row per check) and lineage.tsv (lineage -> simulator binary, verdict).
"""
import os
import shlex
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from m7b_common import *  # noqa: E402,F401,F403

d = sys.argv[1]
log = f"{d}/commands.log"
rows, ex = [], []
for name, path in [("locityper", LOCITYPER), ("strobealign", STROBEALIGN), ("samtools", SAMTOOLS), ("jellyfish_2.3.1", JELLYFISH),
                   ("scoring_probe", SCORING_PROBE), ("c4_sim_candidate_1", C4_SIM_CANDIDATES[0]), ("c4_sim_candidate_2", C4_SIM_CANDIDATES[1]),
                   ("cyp_sim", CYP_SIM), ("brotli", BROTLI)]:
    ver = ""
    for flag in ("--version", "-V", "version"):
        try:
            p = subprocess.run([path, flag], capture_output=True, text=True, timeout=10)
            txt = (p.stdout + p.stderr).strip().splitlines()
            if txt and "usage" not in txt[0].lower() and "missing" not in txt[0].lower():
                ver = txt[0][:120]
                break
        except Exception:
            pass
    ex.append(dict(tool=name, path=path, sha256=sha(path), version=ver or "n/a"))
ex.append(dict(tool="locityper_source_checkout", path="/Users/davide.bolognini/panvar/external/locityper", sha256="git 7a6f2993227c15920df415a19045d5baf98d90e3",
               version="Cargo 1.7.2; installed-binary-to-source identity NOT established"))
write_tsv(f"{d}/executables.tsv", ex)


def same(a, b):
    return os.path.exists(a) and os.path.exists(b) and sha(a) == sha(b)


lineage = []
# C4: first candidate build reproducing both donors' archived source reads byte for byte
c4_bin = None
for i, binary in enumerate(C4_SIM_CANDIDATES, 1):
    ok_all = True
    for donor in LINEAGE_DONORS["C4"]:
        o = f"{d}/c4_cand{i}_{donor}"
        os.makedirs(o)
        code = run([binary, "simulate", "--gfa", LOCI["C4"]["gfa"], "--reference", C4REF, "--donor", donor, "--depth", "30", "--read-length", "100",
                    "--seed", "20260916", "--out", o], log)
        src = f"{R}/genotype_c4_production_11fe681/reads/{donor}/source"
        ok = code == 0 and same(f"{o}/reads_1.fq", f"{src}/reads_1.fq") and same(f"{o}/reads_2.fq", f"{src}/reads_2.fq")
        rows.append(dict(lineage="C4", check=f"simulate {donor} 30x with candidate {i}", binary=binary, exit=code, byte_identical=ok))
        ok_all &= ok
        if not ok:
            break
    if ok_all:
        c4_bin = binary
        break
lineage.append(dict(lineage="C4", simulator=c4_bin or "UNRESOLVED", source_regenerated=bool(c4_bin)))
# CYP: the recorded pilot command, output redirected
cyp_ok = True
for donor in LINEAGE_DONORS["CYP"]:
    rec = open(f"{R}/cyp2d6_shortread_pilot_v1/data/reads/{donor}/simulate.log").readline()[2:].strip().split()
    o = f"{d}/cyp_{donor}"
    os.makedirs(o)
    cmd = rec[:rec.index("--out") + 1] + [o] + rec[rec.index("--out") + 2:]
    code = run(cmd, log)
    src = f"{R}/cyp2d6_shortread_pilot_v1/data/reads/{donor}/cov30"
    ok = code == 0 and cmd[0] == CYP_SIM and same(f"{o}/reads_1.fq", f"{src}/reads_1.fq") and same(f"{o}/reads_2.fq", f"{src}/reads_2.fq")
    rows.append(dict(lineage="CYP", check=f"simulate {donor} 30x (recorded command)", binary=cmd[0], exit=code, byte_identical=ok))
    cyp_ok &= ok
lineage.append(dict(lineage="CYP", simulator=CYP_SIM if cyp_ok else "UNRESOLVED", source_regenerated=cyp_ok))
# each cell's subsample, recorded command re-run into a new directory
cell_ok = {}
for cell, c in CELLS.items():
    logf = f"{c['reads']}/subsample.log" if c["lineage"] == "C4" else f"{os.path.dirname(c['reads'])}/subsample.{c['depth']}.log"
    rec = shlex.split(open(logf).readline()[2:].strip())
    o = f"{d}/subsample_{cell}"
    os.makedirs(o)
    cmd = rec[:rec.index("--out") + 1] + [o] + rec[rec.index("--out") + 2:]
    code = run(cmd, log)
    ok = code == 0 and same(f"{o}/reads_1.fq", f"{c['reads']}/reads_1.fq") and same(f"{o}/reads_2.fq", f"{c['reads']}/reads_2.fq")
    rows.append(dict(lineage=c["lineage"], check=f"subsample {cell} (recorded command)", binary=cmd[0], exit=code, byte_identical=ok))
    cell_ok[cell] = ok
for L in lineage:
    cells = [k for k, c in CELLS.items() if c["lineage"] == L["lineage"]]
    L["cells"] = ",".join(cells)
    L["verdict"] = "PASS" if L["source_regenerated"] and all(cell_ok[k] for k in cells) else "UNRESOLVED"
write_tsv(f"{d}/provenance.tsv", rows)
write_tsv(f"{d}/lineage.tsv", lineage)
print(open(f"{d}/provenance.tsv").read())
print(open(f"{d}/lineage.tsv").read())
sys.exit(0)   # verdicts live in lineage.tsv; downstream units of an UNRESOLVED lineage refuse to run
