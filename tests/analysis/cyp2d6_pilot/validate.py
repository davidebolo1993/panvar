#!/usr/bin/env python3
"""CYP2D6/D7 pilot: package validation (writes VALIDATION.tsv; exit 0 only if every check passes).

Usage: validate.py <pilot package> <tag> <second render dir> <primary repo> <expected primary WIP sha256>
"""
import csv
import glob
import hashlib
import os
import subprocess
import sys

PKG, TAG, RENDER2, REPO, WIP = sys.argv[1:6]
PKG = os.path.abspath(PKG)
RUNS, EVAL = os.path.join(PKG, "data", f"runs_{TAG}"), os.path.join(PKG, "data", f"eval_{TAG}")
rows = []


def check(area, ok, what, detail=""):
    rows.append((area, "PASS" if ok else "FAIL", what, str(detail)))


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def kv(path):
    out = {}
    for line in open(path).read().splitlines():
        p = line.split("\t")
        if len(p) >= 2:
            out.setdefault(p[0], p[1])
    return out


def table(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


# ---- seals and inputs
for s in ("seal_inputs", "stop_evidence", "stop2_evidence", "stop3_evidence", "seal_amendment3", "seal_calls_a3",
          "seal_amendment4", "stop4_evidence", "seal_amendment5", f"seal_calls_{TAG}"):
    r = subprocess.run(["shasum", "-a", "256", "-c", "--quiet", f"{s}.sha256"], cwd=PKG, capture_output=True, text=True)
    check("seals", r.returncode == 0, f"{s}.sha256 verifies", (r.stdout + r.stderr).strip()[:200])
check("binaries", sha(os.path.join(PKG, "bin/panvar")).startswith("c869da7d"), "bin/panvar is the frozen binary")
check("binaries", sha(os.path.join(PKG, "bin/genotype_probe")).startswith("e9bdf722"), "bin/genotype_probe is the post-fix probe")
check("repository", sha(os.path.join(REPO, "tests/genotype_probe.cpp")) == WIP, "primary WIP tests/genotype_probe.cpp unchanged")
stash = subprocess.run(["git", "-C", REPO, "stash", "list", "--format=%H"], capture_output=True, text=True).stdout
check("repository", any(h.startswith("dbb170c5") for h in stash.split()), "stash dbb170c5 kept")

# ---- cells
jobs = table(os.path.join(PKG, "jobs.tsv"))
models = {r["model"]: r for r in table(os.path.join(PKG, "data/models.tsv"))}
check("cells", len(jobs) == 36, "36 frozen jobs", len(jobs))
present = [j for j in jobs if os.path.exists(os.path.join(RUNS, j["job"] + ".run.tsv"))]
check("cells", len(present) == 36, "36 cells called", len(present))
exits, notes, settings, reads, fps, cmdline = {}, [], {}, {}, [], []
for j in jobs:
    p = os.path.join(RUNS, j["job"])
    rec = table(p + ".run.tsv")[0]
    exits[j["job"]] = rec["exit"]
    st = kv(p + ".mosaic_status.tsv")
    if "cap" in st.get("note", "") or "limit" in st.get("note", ""):
        notes.append(j["job"])
    s = {}
    for line in open(p + ".settings.tsv").read().splitlines()[1:]:
        k, v = line.split("\t")[:2]
        if k not in ("gfa", "graph", "reads1", "reads2", "library", "out", "out_prefix", "panel_curation"):
            s[k] = line.split("\t", 1)[1].replace(RUNS, "<RUNS>")
    settings[j["job"]] = s
    reads[j["job"]] = (st.get("sha256_reads1"), st.get("sha256_reads2"), st.get("sha256_library"), st.get("sha256_gfa"))
    log = open(p + ".call.log").read().splitlines()[0]
    truth = os.path.join(PKG, j["truth_gfa"])
    model = os.path.join(PKG, j["model_gfa"])
    cmdline.append(j["job"] if (truth != model and truth in log) or "resource_profile_a4.tsv" not in log
                   or "--dosage-mode repair" not in log or "--threads 4" not in log or "/bin/panvar genotype" not in log
                   or "--library " + os.path.join(PKG, j["library"]) not in log else None)
    mname = "original" if j["arm"] in ("lzo", "mosaic") else f"loo.{j['sample']}"
    for layer in ("raw", "reported"):
        fp = kv(os.path.join(EVAL, f"{j['job']}.{layer}.summary.tsv")).get("model_fingerprint")
        if fp != models[mname]["model_fingerprint"]:
            fps.append(f"{j['job']}.{layer}")
check("cells", set(exits.values()) <= {"0", "1"}, "every cell exited 0 (call) or 1 (no-call)", sorted(set(exits.values())))
check("cells", not notes, "no declared resource cap or limit named in any status note", notes)
check("cells", not [c for c in cmdline if c], "every call: frozen panvar, pilot profile, 4 threads, repair mode, its own "
      "library file, no truth graph on the command line", [c for c in cmdline if c])
ref = settings[jobs[0]["job"]]
diff = [k for k, s in settings.items() if s != ref]
check("cells", not diff, "effective settings and their sources identical across all 36 cells (paths excluded)", diff)
check("cells", ref.get("max_state_classes", "").startswith("600000000") and ref.get("repair_max_neighbours", "").startswith("958"),
      "declared settings in effect: max_state_classes 600 M, repair_max_neighbours 958",
      (ref.get("max_state_classes"), ref.get("repair_max_neighbours")))
pairs = [(j["sample"], j["coverage"]) for j in jobs if j["arm"] == "lzo"]
same = all(reads[f"lzo.{s}.{c}"][:3] == reads[f"loo.{s}.{c}"][:3] for s, c in pairs)
check("cells", same, "LZO and LOO of each donor and coverage read identical reads and library (status hashes)")
mos = [reads[j["job"]][3] for j in jobs if j["arm"] in ("lzo", "mosaic")]
check("cells", len(set(mos)) == 1 and mos[0] == models["original"]["gfa_sha256"],
      "LZO and mosaic cells call on the same frozen model graph", mos[0][:16])
check("evaluation", not fps, "every evaluation's model fingerprint equals the call model's frozen fingerprint", fps)
check("evaluation", len(glob.glob(os.path.join(EVAL, "*.blocks.tsv"))) == 72, "72 evaluations (raw and reported)")
seal_time = os.path.getmtime(os.path.join(PKG, f"seal_calls_{TAG}.sha256"))
first_eval = min(os.path.getmtime(f) for f in glob.glob(os.path.join(EVAL, "*")))
check("evaluation", seal_time < first_eval, "call outputs sealed before the first evaluation file was written")

# ---- classification and tables
cls = {r["job"]: r for r in table(os.path.join(PKG, f"call_classification_{TAG}.tsv"))}
lab = [r for r in cls.values() if r["classification"] == "SEQUENCE_IDENTICAL_BLOCK_LABEL_AMBIGUITY"]
check("classification", all(r["scores_identical"] == "true" and r["unordered_sequence_identical"] == "true"
                            and len(r["sequence_candidate_ids"].split(";")) == 1
                            and r["differences_within_degeneracy_interval"] == "true" and r["caller_status"] for r in lab),
      "every label-ambiguity cell meets the five criteria; caller status kept", len(lab))
T = os.path.join(PKG, "tables")
cal = table(os.path.join(T, "callability.tsv"))
check("tables", sum(int(r["cells"]) for r in cal) == 36 and all(
    int(r["strict"]) + int(r["genotype_equivalent"]) + int(r["label_ambiguity_no_call"]) + int(r["distinct_sequence_tie_no_call"])
    + int(r["other_no_call"]) + int(r["not_completed"]) == int(r["cells"]) for r in cal),
      "callability accounts for every cell exactly once")
obs = table(os.path.join(T, "block_observations.tsv"))
nocall = {j for j, e in exits.items() if e != "0"}
check("tables", not [o for o in obs if o["run_id"] in nocall], "no block observation exists for any no-call")
caller = {}
for j, e in exits.items():
    if e == "0":
        for r in table(os.path.join(RUNS, j + ".block_genotypes.tsv")):
            caller[(j, r["block_id"])] = r["genotype_filter"]
check("tables", all(o["genotype_filter"] == caller[(o["run_id"], o["block_id"])] for o in obs if o["layer"] == "reported"),
      "retained membership is the caller's own genotype_filter (read from the call outputs)")
und = table(os.path.join(T, "truth_undefined_blocks.tsv"))
check("tables", all(r["called1_bp"] != "NA" for r in und), "truth-undefined blocks tabulated with their called genotype, not valued",
      len(und))
fc = table(os.path.join(T, "figure_check.tsv"))
check("figures", all(r["points"] == r["expected"] for r in fc), "every block genotype drawn exactly once (A raw, B retained)",
      [(r["figure"], r["points"], r["expected"]) for r in fc])
check("figures", all(r["clipped_texts"] == "0" and r["overlapping_text_pairs"] == "0" for r in fc), "no clipped or overlapping text")
diffs = []
for f in sorted(glob.glob(os.path.join(PKG, "figures", "*")) + glob.glob(os.path.join(PKG, "tables", "*"))):
    other = os.path.join(RENDER2, "figures" if "/figures/" in f else "data", os.path.basename(f))
    if not os.path.exists(other) or sha(f) != sha(other):
        diffs.append(os.path.basename(f))
check("figures", not diffs, "an independent second render is byte-identical (figures and tables)", diffs)

with open(os.path.join(PKG, "VALIDATION.tsv"), "w") as f:
    f.write("area\tstatus\tcheck\tdetail\n" + "".join("\t".join(r) + "\n" for r in rows))
failed = [r for r in rows if r[1] != "PASS"]
for r in rows:
    print("\t".join(r)[:220])
print(f"{len(rows)} checks, {len(failed)} failed")
sys.exit(1 if failed else 0)
