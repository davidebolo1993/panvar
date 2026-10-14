#!/usr/bin/env python3
"""CYP2D6/D7 short-read pilot: inputs, reads, jobs, calls and evaluation (one fixed build for every arm).

  pilot.py models <pkg> <build>       decompress the graph; LOO graphs (exactly the donor's two P lines removed);
                                      preflight every model with the curation -> data/models.tsv
  pilot.py jobs <pkg>                 jobs.tsv (36 cells) and samples.tsv from SPEC constants and mosaic_plans.tsv
  pilot.py reads <pkg> <build>        30x once per truth sample, 5x nested by subsample; checks -> data/reads.tsv
  pilot.py run <pkg> <build> <job>... calls (<= 2 at a time, 4 threads each) then evaluates the raw (.original) and
                                      reported layers; per-cell wall time and peak RSS -> data/runs/<job>.run.tsv
                                      (the original design, profile resource_profile.tsv; used by the STOP_3 sentinels)
  pilot.py call <pkg> <build> <profile> <tag> <job>...
                                      AMENDMENT_3: calls only, with the named profile, into data/runs_<tag>/<job>;
                                      no evaluation (truth is not read until the call outputs are sealed)
  pilot.py evaluate <pkg> <build> <tag> <job>...
                                      raw (.original) and reported evaluations of data/runs_<tag> -> data/eval_<tag>

No truth file is on a caller command line. The LZO and LOO arms of one donor share the same read files.
"""
import concurrent.futures
import csv
import gzip
import hashlib
import os
import re
import shutil
import subprocess
import sys
import time

REF = "grch38#1#chr22:42031864-42245566"
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
GRAPH_GZ = os.path.join(REPO, "tests/real_data/cyp2d6.gfa.gz")
CURATION = os.path.join(REPO, "tests/genotype_cyp2d6_curation.tsv")
DONORS = ["HG02769", "HG03456", "HG02953", "NA19650", "NA19705", "NA19983"]
COVERAGES = [5, 30]
SOURCE_COVERAGE = 30
SIM = {"seed": "20260921", "read-length": "100", "insert-mean": "350", "insert-sd": "40", "sigmas": "4",
       "error": "0.002", "outlier": "0.001"}
THREADS = 4
PARALLEL = 2
BACKSTOP_SECONDS = 1500          # the caller's own time_limit_seconds (1200) is the limit; this only reaps a hang
COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def run(cmd, log):
    with open(log, "w") as f:
        f.write("$ " + " ".join(map(str, cmd)) + "\n")
        f.flush()
        return subprocess.run([str(c) for c in cmd], stdout=f, stderr=subprocess.STDOUT).returncode


def kv(path):
    out = {}
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2 and p[0] != "key":
            out[p[0]] = p[1]
    return out


def read_gfa(path):
    op = gzip.open if path.endswith(".gz") else open
    seq, paths = {}, {}
    with op(path, "rt") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[0] == "S":
                seq[p[1]] = p[2]
            elif p[0] == "P":
                paths[p[1]] = p[2].split(",")
    return seq, paths


def spell(seq, steps):
    return "".join(seq[s[:-1]].translate(COMP)[::-1] if s.endswith("-") else seq[s[:-1]] for s in steps)


def fail(msg):
    print("FAIL: " + msg, file=sys.stderr)
    sys.exit(1)


# ---- models -------------------------------------------------------------------------------------------

def models(pkg, build):
    probe = os.path.join(build, "genotype_probe")
    gdir = os.path.join(pkg, "data", "graphs")
    os.makedirs(gdir, exist_ok=True)
    graph = os.path.join(gdir, "cyp2d6.gfa")
    if not os.path.exists(graph):
        with gzip.open(GRAPH_GZ, "rb") as src, open(graph, "wb") as dst:
            shutil.copyfileobj(src, dst)
    orig_lines = open(graph).read().splitlines()
    rows = []
    entries = [("original", graph)]
    for donor in DONORS:
        loo = os.path.join(gdir, f"loo.{donor}.gfa")
        if run([probe, "loo-gfa", "--gfa", graph, "--donor", donor, "--out", loo], loo + ".log") != 0:
            fail(f"loo-gfa {donor}")
        loo_lines = open(loo).read().splitlines()
        loo_set = set(loo_lines)
        removed = [l for l in orig_lines if l not in loo_set]
        donor_lines = [l for l in orig_lines if l.startswith("P\t") and l.split("\t")[1].startswith(donor + "#")]
        ok = (len(removed) == 2 and all(l.startswith("P\t") and l.split("\t")[1].startswith(donor + "#") for l in removed)
              and loo_lines == [l for l in orig_lines if l not in removed] and removed == donor_lines)
        if not ok:
            fail(f"LOO graph {donor} is not the original minus exactly the donor's two P lines")
        rows.append({"model": f"loo.{donor}", "removed_paths": ",".join(l.split("\t")[1] for l in removed)})
        entries.append((f"loo.{donor}", loo))
    out = []
    for name, path in entries:
        pf = os.path.join(gdir, f"{name}.preflight.tsv")
        if run([probe, "preflight", "--gfa", path, "--reference", REF, "--panel-curation", CURATION, "--out", pf],
               pf + ".log") != 0:
            fail(f"preflight {name}")
        p = kv(pf)
        rec = {"model": name, "gfa": os.path.relpath(path, pkg), "gfa_sha256": sha256(path),
               "removed_paths": next((r["removed_paths"] for r in rows if r["model"] == name), "none"),
               "curation": os.path.relpath(CURATION, REPO), "templates": p["templates"],
               "templates_curated_out": p["templates_curated_out"], "blocks": p["blocks"],
               "variant_blocks": p["variant_blocks"], "context_blocks": p["context_blocks"],
               "panel_round_trip": p["panel_round_trip"], "representation_independent": p["representation_independent"],
               "model_runnable": p["model_runnable"], "mosaic_arm_runnable": p["mosaic_arm_runnable"],
               "model_fingerprint": p["model_fingerprint"], "predicted_peak_rss_MB": p["predicted_peak_rss_MB"]}
        expect_templates = "126" if name == "original" else "124"
        if not (rec["templates"] == expect_templates and rec["templates_curated_out"] == "1" and rec["blocks"] == "19"
                and rec["panel_round_trip"] == "ok" and rec["model_runnable"] == "true"):
            fail(f"preflight {name}: {rec}")
        if name == "original" and rec["model_fingerprint"] != "2a5b720a1c5bf669":
            fail("original model fingerprint changed")
        out.append(rec)
    with open(os.path.join(pkg, "data", "models.tsv"), "w") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    print(f"{len(out)} models prepared and preflighted")


# ---- jobs -----------------------------------------------------------------------------------------------

def plans(pkg):
    with open(os.path.join(pkg, "mosaic_plans.tsv")) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def jobs(pkg):
    samples, rows = [], []
    for d in DONORS:
        samples.append({"sample": d, "kind": "donor", "truth_gfa": "data/graphs/cyp2d6.gfa", "simulate_donor": d,
                        "seed_salt": d, "arms": "lzo,loo"})
    for p in plans(pkg):
        s = f"mosaic_{p['plan']}"
        samples.append({"sample": s, "kind": "mosaic", "truth_gfa": f"mosaic/plans/{p['plan']}/truth.gfa",
                        "simulate_donor": "MOSAIC", "seed_salt": s, "arms": "mosaic"})
    for arm in ("lzo", "loo", "mosaic"):
        for s in samples:
            if arm not in s["arms"].split(","):
                continue
            for cov in COVERAGES:
                job = f"{arm}.{s['sample']}.{cov}"
                rows.append({
                    "job": job, "arm": arm, "sample": s["sample"], "coverage": cov, "scheme": arm,
                    "eval_donor": s["simulate_donor"],
                    "model_gfa": f"data/graphs/loo.{s['sample']}.gfa" if arm == "loo" else "data/graphs/cyp2d6.gfa",
                    "curation": os.path.relpath(CURATION, REPO), "truth_gfa": s["truth_gfa"],
                    "reads1": f"data/reads/{s['sample']}/cov{cov}/reads_1.fq",
                    "reads2": f"data/reads/{s['sample']}/cov{cov}/reads_2.fq",
                    "library": f"data/reads/{s['sample']}/cov{cov}/library.txt",
                    "out_prefix": f"data/runs/{job}"})
    for name, table in (("samples.tsv", samples), ("jobs.tsv", rows)):
        with open(os.path.join(pkg, name), "w") as f:
            w = csv.DictWriter(f, fieldnames=list(table[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(table)
    print(f"{len(samples)} truth samples, {len(rows)} jobs")


def load_jobs(pkg):
    with open(os.path.join(pkg, "jobs.tsv")) as f:
        return {r["job"]: r for r in csv.DictReader(f, delimiter="\t")}


# ---- reads ------------------------------------------------------------------------------------------------

def fastq_names(path):
    with open(path) as f:
        return [line.split()[0][1:].rsplit("/", 1)[0] for i, line in enumerate(f) if i % 4 == 0]


def reads(pkg, build):
    probe = os.path.join(build, "genotype_probe")
    with open(os.path.join(pkg, "samples.tsv")) as f:
        samples = list(csv.DictReader(f, delimiter="\t"))
    oseq, opaths = read_gfa(os.path.join(pkg, "data/graphs/cyp2d6.gfa"))
    out = []
    for s in samples:
        root = os.path.join(pkg, "data", "reads", s["sample"])
        src = os.path.join(root, f"cov{SOURCE_COVERAGE}")
        os.makedirs(src, exist_ok=True)
        cmd = [probe, "simulate", "--gfa", os.path.join(pkg, s["truth_gfa"]), "--reference", REF,
               "--donor", s["simulate_donor"], "--depth", str(SOURCE_COVERAGE)]
        for k, v in SIM.items():
            cmd += [f"--{k}", v]
        cmd += ["--seed-salt", s["seed_salt"], "--out", src]
        if run(cmd, os.path.join(root, "simulate.log")) != 0:
            fail(f"simulate {s['sample']}")
        os.remove(os.path.join(src, "truth.sorted.gfa"))
        # the simulated haplotypes are exactly the truth paths, spelled independently here
        if s["kind"] == "donor":
            tseq, tpaths = oseq, opaths
        else:
            tseq, tpaths = read_gfa(os.path.join(pkg, s["truth_gfa"]))
        want = {n: spell(tseq, st) for n, st in tpaths.items() if n.startswith(s["simulate_donor"] + "#")}
        got, name = {}, None
        for line in open(os.path.join(src, "truth.fa")):
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line[1:]
            else:
                got[name] = line
        if got != want or len(got) != 2:
            fail(f"{s['sample']}: simulated haplotypes differ from the truth paths")
        for cov in COVERAGES:
            if cov == SOURCE_COVERAGE:
                continue
            d = os.path.join(root, f"cov{cov}")
            os.makedirs(d, exist_ok=True)
            if run([probe, "subsample", "--reads1", os.path.join(src, "reads_1.fq"),
                    "--reads2", os.path.join(src, "reads_2.fq"), "--library", os.path.join(src, "library.txt"),
                    "--source-coverage", SOURCE_COVERAGE, "--coverage", cov, "--name-salt", s["seed_salt"],
                    "--out", d], os.path.join(root, f"subsample.{cov}.log")) != 0:
                fail(f"subsample {s['sample']} {cov}")
        names = {}
        for cov in COVERAGES:
            d = os.path.join(root, f"cov{cov}")
            n1, n2 = fastq_names(os.path.join(d, "reads_1.fq")), fastq_names(os.path.join(d, "reads_2.fq"))
            if n1 != n2:
                fail(f"{s['sample']} {cov}x: mates out of step")
            names[cov] = n1
        nested = set(names[min(COVERAGES)]) <= set(names[SOURCE_COVERAGE])
        if not nested:
            fail(f"{s['sample']}: {min(COVERAGES)}x is not a subset of {SOURCE_COVERAGE}x")
        for cov in COVERAGES:
            d = os.path.join(root, f"cov{cov}")
            out.append({"sample": s["sample"], "coverage": cov, "fragments": len(names[cov]),
                        "nested_in_30x": "true", "haplotypes_equal_truth_paths": "true",
                        "reads_1_sha256": sha256(os.path.join(d, "reads_1.fq")),
                        "reads_2_sha256": sha256(os.path.join(d, "reads_2.fq")),
                        "library_sha256": sha256(os.path.join(d, "library.txt")),
                        "fragment_start_intensity": kv_space(os.path.join(d, "library.txt"))["fragment_start_intensity"]})
        print(f"{s['sample']}: {len(names[SOURCE_COVERAGE])} fragments at 30x, {len(names[5])} at 5x")
    with open(os.path.join(pkg, "data", "reads.tsv"), "w") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)


def kv_space(path):
    out = {}
    for line in open(path):
        line = line.split("#", 1)[0].strip()
        if line:
            k, v = line.split(None, 1)
            out[k] = v
    return out


# ---- calls and evaluation -------------------------------------------------------------------------------

def tagged(job, tag):
    """The job with its outputs under data/runs_<tag>/ (tag None: the original data/runs/)."""
    if tag is None:
        return job
    job = dict(job)
    job["out_prefix"] = job["out_prefix"].replace("data/runs/", f"data/runs_{tag}/", 1)
    return job


def call_cmd(job, pkg, build, profile=None):
    j = lambda k: os.path.join(pkg, job[k])   # noqa: E731
    cmd = [os.path.join(build, "panvar"), "genotype", "--preset", "short-read",
           "--profile", profile or os.path.join(pkg, "resource_profile.tsv"), "--threads", str(THREADS),
           "--graph", j("model_gfa"), "--reference-path", REF, "--panel-curation", os.path.join(REPO, job["curation"]),
           "--reads1", j("reads1"), "--reads2", j("reads2"), "--library", j("library"),
           "--dosage-mode", "repair", "--out", j("out_prefix"), "--quiet"]
    if job["truth_gfa"] != job["model_gfa"] and os.path.join(pkg, job["truth_gfa"]) in cmd:
        fail(f"{job['job']}: truth graph on the call's command line")
    return cmd


def evaluate(job, pkg, build, layer, eval_dir="eval"):
    prefix = os.path.join(pkg, job["out_prefix"]) + (".original" if layer == "raw" else "")
    os.makedirs(os.path.join(pkg, "data", eval_dir), exist_ok=True)
    out = os.path.join(pkg, "data", eval_dir, job["job"] + (".raw" if layer == "raw" else ".reported"))
    cmd = [os.path.join(build, "genotype_probe"), "evaluate", "--scheme", job["scheme"], "--donor", job["eval_donor"],
           "--reference", REF, "--run-prefix", prefix, "--model-gfa", os.path.join(pkg, job["model_gfa"]),
           "--truth-gfa", os.path.join(pkg, job["truth_gfa"]),
           "--panel-curation", os.path.join(REPO, job["curation"]), "--out", out]
    code = run(cmd, out + ".evaluate.log")
    for staged in (out + ".model.sorted.gfa", out + ".truth.sorted.gfa"):
        if os.path.exists(staged):
            os.remove(staged)
    return code


def one(job, pkg, build, profile=None, evaluate_after=True):
    prefix = os.path.join(pkg, job["out_prefix"])
    if os.path.exists(prefix + ".run.tsv"):
        fail(f"{job['job']}: {prefix}.run.tsv exists; a finished cell is never rerun into the same place")
    os.makedirs(os.path.dirname(prefix), exist_ok=True)
    cmd = ["/usr/bin/time", "-l"] + call_cmd(job, pkg, build, profile)
    t0 = time.monotonic()
    with open(prefix + ".call.log", "w") as log:
        log.write("$ " + " ".join(cmd) + "\n")
        log.flush()
        try:
            code = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, timeout=BACKSTOP_SECONDS).returncode
        except subprocess.TimeoutExpired:
            code = "backstop_timeout"
    wall = time.monotonic() - t0
    text = open(prefix + ".call.log").read()
    m = re.search(r"(\d+)\s+maximum resident set size", text)
    rss = int(m.group(1)) if m else -1
    m = re.search(r"(\d+)\s+peak memory footprint", text)
    footprint = int(m.group(1)) if m else -1
    ev = {}
    if code in (0, 1) and evaluate_after:
        for layer in ("raw", "reported"):
            ev[layer] = evaluate(job, pkg, build, layer)
    rec = {"job": job["job"], "exit": code, "wall_seconds": f"{wall:.1f}", "peak_rss_bytes": rss,
           "peak_footprint_bytes": footprint, "eval_raw_exit": ev.get("raw", "NA"),
           "eval_reported_exit": ev.get("reported", "NA")}
    with open(prefix + ".run.tsv", "w") as f:
        f.write("\t".join(rec) + "\n" + "\t".join(map(str, rec.values())) + "\n")
    print("\t".join(f"{k}={v}" for k, v in rec.items()), flush=True)
    return rec


def run_jobs(pkg, build, names, profile=None, tag=None, evaluate_after=True):
    table = load_jobs(pkg)
    missing = [n for n in names if n not in table]
    if missing:
        fail(f"unknown jobs {missing}")
    with concurrent.futures.ThreadPoolExecutor(max_workers=PARALLEL) as pool:
        recs = list(pool.map(lambda n: one(tagged(table[n], tag), pkg, build, profile, evaluate_after), names))
    if evaluate_after:
        bad = [r for r in recs if r["exit"] not in (0, 1) or r["eval_raw_exit"] != 0 or r["eval_reported_exit"] != 0]
        print(f"{len(recs)} cells, {len(bad)} not completed and evaluated")
    else:
        bad = [r for r in recs if r["exit"] not in (0, 1)]
        print(f"{len(recs)} cells, {len(bad)} not completed (exit other than 0 or 1)")
    return 1 if bad else 0


def evaluate_jobs(pkg, build, tag, names):
    table = load_jobs(pkg)
    bad = []
    for n in names:
        job = tagged(table[n], tag)
        rec = tsv_row(os.path.join(pkg, job["out_prefix"] + ".run.tsv"))
        if rec["exit"] not in ("0", "1"):
            print(f"{n}: exit {rec['exit']}, not evaluated")
            continue
        for layer in ("raw", "reported"):
            code = evaluate(job, pkg, build, layer, f"eval_{tag}")
            if code != 0:
                bad.append((n, layer))
        print(f"{n}: evaluated", flush=True)
    print(f"{len(names)} cells, {len(bad)} evaluation failures {bad}")
    return 1 if bad else 0


def tsv_row(path):
    lines = open(path).read().splitlines()
    return dict(zip(lines[0].split("\t"), lines[1].split("\t")))


if __name__ == "__main__":
    cmd, pkg = sys.argv[1], os.path.abspath(sys.argv[2])
    if cmd == "models":
        models(pkg, os.path.abspath(sys.argv[3]))
    elif cmd == "jobs":
        jobs(pkg)
    elif cmd == "reads":
        reads(pkg, os.path.abspath(sys.argv[3]))
    elif cmd == "run":
        sys.exit(run_jobs(pkg, os.path.abspath(sys.argv[3]), sys.argv[4:]))
    elif cmd == "call":
        sys.exit(run_jobs(pkg, os.path.abspath(sys.argv[3]), sys.argv[6:], profile=os.path.abspath(sys.argv[4]),
                          tag=sys.argv[5], evaluate_after=False))
    elif cmd == "evaluate":
        sys.exit(evaluate_jobs(pkg, os.path.abspath(sys.argv[3]), sys.argv[4], sys.argv[5:]))
    else:
        fail(f"unknown command {cmd}")
