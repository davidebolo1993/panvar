#!/usr/bin/env python3
"""The final C4 production campaigns: resumable launch, classification, completeness and evaluation.

    scripts/c4_production_campaign.py launch --build <dir> --out <plan dir> --jobs <J> --threads <T>
    scripts/c4_production_campaign.py prepare --build <dir> --out <plan dir> [--only mosaic]
    scripts/c4_production_campaign.py status --out <plan dir>

The plan (jobs.tsv, settings.tsv) is written by tests/genotype_production_c4_campaign.sh plan; this runner
never re-selects a donor, a plan, a coverage or a setting. Every call is the one production command,

    panvar genotype --profile profiles/genotype_c4_v6_campaign_local.tsv --threads <T> ...

(the v6 campaign-local profile: the science of profiles/genotype_c4_v4.tsv -- tie branching with cap 16, the
representation guard, 32 repair rounds -- plus the reviewed campaign-local resource and exact-acceleration
settings), and truth enters only the evaluation that follows it. The production binary, the profile, the plan and
the scheduler (concurrent jobs, threads per call) are the plan's: launch refuses anything else.

Outcomes. A finished job is exactly one of
    strict_call                 COMPLETE / strict_call: genotypes, block genotypes, the haplotype FASTA and the
                                repair rounds exist
    genotype_equivalent_call    CONDITIONAL_GENOTYPE_EQUIVALENT / genotype_equivalent_call: genotypes and the
                                equivalence set, no FASTA
    no_call:<reason>            a scientific refusal the caller declared (no_call_reason): an official no-call,
                                a COMPLETED job, never missing, never retried, never a stop. None may leave a
                                genotype table, a block-genotype table or a FASTA. Each keeps its own evidence:
                                  tie_branching_incomplete  REPAIR_INCOMPLETE, tie_branched: the repair rounds and
                                                            both tie tables (a branch cap, any cap or failure while
                                                            branching)
                                  divergent_tie             RESCORE_INCOMPLETE, tie_branched: the repair rounds, both
                                                            tie tables and the equivalence set
                                  time_limit, footprint_limit   the watchdog stopped the process (not evaluated)
                                  improving_move_tie        v1 only (tie branching off): the repair rounds. Under
                                                            v6 an improving-move tie is branched, so it cannot occur
and an unfinished one is pending, or failed:<why> after an external failure (a crash, a signal, a missing or
RUNNING status, a failed evaluation). Only external failures are retried; a job that fails externally on
every attempt stops the campaign, as does any harness or provenance defect (an outcome whose files break its
contract, an input whose content moved, a binary or profile that is not the plan's or changed between
launches). A legitimate no-call never stops anything.

Inputs are bound by SHA-256 the moment they are created (inputs.tsv) and re-verified on every reuse. Jobs are
isolated: each writes only under its own prefix. completeness.tsv is rewritten after every job.
"""
import argparse
import concurrent.futures
import csv
import datetime
import hashlib
import os
import pathlib
import shutil
import subprocess
import sys
import threading
import time

REPO = pathlib.Path(__file__).resolve().parent.parent
PROFILE = REPO / "profiles" / "genotype_c4_v6_campaign_local.tsv"
REFERENCE = "HG00096#1#haplotype1-0000024:31848049-32074186"
GRAPH_GZ = REPO / "tests" / "real_data" / "c4.gfa.gz"
GRAPH_SHA = "b13af390c28dd8b43a913c43f123bb1fcf639e0569f99783fa12a778820c8733"
MOSAIC_SIGNATURES = REPO / "tests" / "results" / "genotype_c4_mosaic_c8fe2a4_raw" / "signatures"
MOSAIC_HARNESS = REPO / "tests" / "genotype_c4_heldout_mosaic.sh"
SOURCE_COVERAGE = "30"
MAX_ATTEMPTS = 3            # one run plus two retries, for external failures only
PROCESS_BACKSTOP_SECONDS = 4000   # the profile's own watchdog stops a call at 3600 s; this catches a hang
EXTERNAL_EXIT_CODES_OK = {0, 1, 124, 125}


class Defect(RuntimeError):
    """A harness or provenance defect: the campaign stops."""


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def key_values(path):
    out = {}
    try:
        with open(path) as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t", 1)
                if len(parts) == 2 and parts[0] not in out:
                    out[parts[0]] = parts[1]
    except FileNotFoundError:
        return None
    return out


def read_jobs(out):
    with open(out / "jobs.tsv") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def now():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


class Inputs:
    """inputs.tsv: every generated input, bound by content when it is created and checked on every reuse."""

    def __init__(self, out):
        self.path = out / "inputs.tsv"
        self.lock = threading.Lock()
        self.rows = {}
        if self.path.exists():
            with open(self.path) as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    self.rows[row["path"]] = row

    def bind(self, path, origin):
        path = str(path)
        with self.lock:
            digest = sha256(path)
            if path in self.rows:
                if self.rows[path]["sha256"] != digest:
                    raise Defect(f"input {path} changed content: bound {self.rows[path]['sha256']}, now {digest}")
                return digest
            self.rows[path] = {"path": path, "sha256": digest, "bytes": str(os.path.getsize(path)),
                               "origin": origin, "bound_utc": now()}
            self._write()
            return digest

    def verify(self, path):
        path = str(path)
        with self.lock:
            if path not in self.rows:
                raise Defect(f"input {path} exists but was never bound; refusing to use it")
            digest = sha256(path)
            if digest != self.rows[path]["sha256"]:
                raise Defect(f"input {path} changed content since it was bound")
            return digest

    def _write(self):
        tmp = str(self.path) + ".writing"
        with open(tmp, "w") as handle:
            writer = csv.DictWriter(handle, fieldnames=["path", "sha256", "bytes", "origin", "bound_utc"],
                                    delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in sorted(self.rows.values(), key=lambda r: r["path"]):
                writer.writerow(row)
        os.replace(tmp, self.path)


def run(cmd, log, cwd=None):
    with open(log, "a") as handle:
        handle.write(f"$ {' '.join(map(str, cmd))}\n")
        handle.flush()
        return subprocess.run([str(c) for c in cmd], stdout=handle, stderr=subprocess.STDOUT, cwd=cwd).returncode


# ---- input preparation (serial, in the dispatcher) --------------------------------------------------

def ensure_graph(out, inputs):
    graph = out / "graphs" / "c4.gfa"
    graph.parent.mkdir(parents=True, exist_ok=True)
    if not graph.exists():
        with open(graph, "wb") as target:
            subprocess.run(["gzip", "-dc", str(GRAPH_GZ)], stdout=target, check=True)
    if sha256(graph) != GRAPH_SHA:
        raise Defect("the staged C4 graph does not match the committed graph")
    inputs.bind(graph, "gzip -dc tests/real_data/c4.gfa.gz")
    return graph


def ensure_inputs(job, out, build, inputs):
    """Creates or re-verifies this job's model graph and reads. Returns nothing; raises Defect."""
    probe = build / "genotype_probe"
    campaign, subject, cov = job["campaign"], job["subject"], job["coverage"]
    model = pathlib.Path(job["model_gfa"])
    if campaign == "loo":
        if not model.exists():
            if run([probe, "loo-gfa", "--gfa", out / "graphs" / "c4.gfa", "--donor", subject, "--out", model],
                   str(model) + ".log") != 0:
                raise Defect(f"loo-gfa failed for {subject}")
            inputs.bind(model, f"genotype_probe loo-gfa --donor {subject}")
        else:
            inputs.verify(model)
    if campaign == "mosaic":
        for path in (model, job["reads1"], job["reads2"], job["library"], job["curation"]):
            inputs.verify(path)
        return
    reads = pathlib.Path(job["reads1"]).parent
    source = out / "reads" / subject / "source"
    if not pathlib.Path(job["reads1"]).exists():
        source.mkdir(parents=True, exist_ok=True)
        for mate, gz in (("reads_1.fq", job["source_reads1"]), ("reads_2.fq", job["source_reads2"])):
            if not (source / mate).exists():
                with open(source / mate, "wb") as target:
                    subprocess.run(["gzip", "-dc", gz], stdout=target, check=True)
                inputs.bind(source / mate, f"gzip -dc {gz}")
            else:
                inputs.verify(source / mate)
        reads.mkdir(parents=True, exist_ok=True)
        salt = [] if job["subsample_salt"] == "none" else ["--name-salt", job["subsample_salt"]]
        if run([probe, "subsample", "--reads1", source / "reads_1.fq", "--reads2", source / "reads_2.fq",
                "--library", job["source_library"], "--source-coverage", SOURCE_COVERAGE, "--coverage", cov]
               + salt + ["--out", reads], reads / "subsample.log") != 0:
            raise Defect(f"subsample failed for {subject} {cov}x")
        recorded = key_values(reads / "subsample.tsv") or {}
        if recorded.get("name_salt") != job["subsample_salt"]:
            raise Defect(f"{subject} {cov}x was thinned with salt {recorded.get('name_salt')}, not the plan's")
        kept = recorded.get("kept_fragments")
        for path in (job["reads1"], job["reads2"], job["library"]):
            inputs.bind(path, f"genotype_probe subsample {subject} {cov}x (unsalted, kept {kept}; the retired "
                              f"abb35a0 selector kept {job['package_kept_fragments']}) of the frozen held-out source")
    else:
        for path in (job["reads1"], job["reads2"], job["library"]):
            inputs.verify(path)


def prepare_mosaic(out, build, inputs, log):
    """The held-out synthetic-mosaic reads are genuinely missing (the c8fe2a4 package kept only their
    hashes). They are regenerated by the committed harness and checked against those recorded hashes."""
    root = out / "mosaic_inputs"
    done = root / "prepared.tsv"
    if done.exists():
        return
    env = dict(os.environ, STAGES="manifest stage verify reads nested", BUILD=str(build))
    with open(log, "a") as handle:
        code = subprocess.run(["bash", str(MOSAIC_HARNESS), str(root)], stdout=handle, stderr=subprocess.STDOUT,
                              env=env).returncode
    if code != 0:
        raise Defect(f"the held-out mosaic harness failed (exit {code}); see {log}")
    recorded = {}
    for sig in MOSAIC_SIGNATURES.glob("*.tsv"):
        kv = key_values(sig)
        if kv and "plan" in kv:
            recorded[(kv["plan"], kv["coverage"])] = kv
    rows = []
    for plan_dir in sorted(p for p in root.iterdir() if p.is_dir() and (p / "staged.gfa").exists()):
        plan = plan_dir.name
        for cov in ("2", "5", "10", "30"):
            reads = plan_dir / f"cov{cov}"
            got = {"reads_1": sha256(reads / "reads_1.fq"), "reads_2": sha256(reads / "reads_2.fq"),
                   "library": sha256(reads / "library.txt"), "staged_graph": sha256(plan_dir / "staged.gfa")}
            want = recorded.get((plan, cov), {})
            match = all(want.get(k) == v for k, v in got.items()) if want else None
            for name, path in (("reads_1", reads / "reads_1.fq"), ("reads_2", reads / "reads_2.fq"),
                               ("library", reads / "library.txt")):
                inputs.bind(path, f"tests/genotype_c4_heldout_mosaic.sh reads/nested {plan} {cov}x")
            rows.append((plan, cov, got["reads_1"], got["reads_2"], got["library"], got["staged_graph"],
                         "no_record" if match is None else ("identical" if match else "DIFFERS")))
        inputs.bind(plan_dir / "staged.gfa", f"tests/genotype_c4_heldout_mosaic.sh stage {plan}")
    inputs.bind(root / "mosaic_exclude.tsv", "tests/genotype_c4_heldout_mosaic.sh")
    differs = [r for r in rows if r[6] == "DIFFERS"]
    tmp = str(done) + ".writing"
    with open(tmp, "w") as handle:
        handle.write("plan\tcoverage\treads_1_sha256\treads_2_sha256\tlibrary_sha256\tstaged_graph_sha256\t"
                     "vs_c8fe2a4_record\n")
        for row in rows:
            handle.write("\t".join(row) + "\n")
    os.replace(tmp, done)
    if differs:
        raise Defect(f"{len(differs)} regenerated held-out mosaic inputs differ from the c8fe2a4 record; "
                     f"see {done}")


# ---- one job: the call, then its evaluation ---------------------------------------------------------

def call_command(job, build, threads):
    cmd = [build / "panvar", "genotype", "--profile", PROFILE, "--threads", str(threads), "--gfa", job["model_gfa"],
           "--reference-path", REFERENCE, "--reads1", job["reads1"], "--reads2", job["reads2"],
           "--library", job["library"], "--out-prefix", job["out_prefix"], "--quiet"]
    if job["curation"] != "none":
        cmd += ["--panel-curation", job["curation"]]
    if job["campaign"] == "loo" and job["truth_gfa"] in map(str, cmd):
        raise Defect(f"{job['job']}: the truth graph is on the call's command line")
    return cmd


def evaluate(job, evaluator, fallback=None):
    prefix = job["out_prefix"]
    # Which evaluator measured this job: a later evaluator re-measures it before consolidation.
    with open(prefix + (".fallback" if fallback else ".eval") + ".evaluator", "w") as stamp:
        stamp.write(sha256(evaluator) + "\n")
    out = prefix + (".fallback" if fallback else ".eval")
    cmd = [evaluator, "evaluate", "--scheme", job["scheme"], "--donor", job["donor_label"],
           "--reference", REFERENCE, "--run-prefix", prefix, "--model-gfa", job["model_gfa"],
           "--truth-gfa", job["truth_gfa"], "--out", out]
    if job["curation"] != "none":
        cmd += ["--panel-curation", job["curation"]]
    if fallback:
        cmd += ["--fallback-alleles", fallback]
    code = run(cmd, prefix + ".evaluate.log")
    # The evaluator stages two sorted graphs per call; they are intermediates, ~7 MB a job.
    for staged in (out + ".model.sorted.gfa", out + ".truth.sorted.gfa"):
        if os.path.exists(staged):
            os.remove(staged)
    return code


def classify(job):
    """(outcome, detail) from what the job left on disk, enforcing each outcome's file contract."""
    prefix = job["out_prefix"]
    status = key_values(prefix + ".mosaic_status.tsv")
    exists = lambda suffix: os.path.exists(prefix + suffix)
    if status is None:
        return "pending", "no status"
    state = status.get("status", "")
    if state == "RUNNING":
        return "failed", "status still RUNNING"
    # The final status and call class decide; repair_converged describes the production search alone and stays
    # false for a tie-branched call (search_conditional tie_branched).
    if state == "COMPLETE":
        if status.get("call_class") != "strict_call":
            return "defect", "COMPLETE without call_class strict_call"
        missing = [s for s in (".genotypes.tsv", ".block_genotypes.tsv", ".mosaic_haplotypes.fa", ".repair_rounds.tsv")
                   if not exists(s)]
        if missing or exists(".repair_rounds.partial.tsv"):
            return "defect", f"strict call breaks its file contract ({missing})"
        return "strict_call", status.get("phase_support", "NA")
    if state == "CONDITIONAL_GENOTYPE_EQUIVALENT":
        if exists(".mosaic_haplotypes.fa") or not exists(".genotypes.tsv") or not exists(".equivalence_set.tsv"):
            return "defect", "genotype-equivalent call breaks its file contract"
        return "genotype_equivalent_call", status.get("phase_support", "NA")
    reason = status.get("no_call_reason", "NA")
    # A no-call is its own state: a quality filter belongs to calls only and never relabels one.
    if status.get("locus_filter", "NA") not in ("NA", ""):
        return "defect", f"{state} carries locus_filter {status.get('locus_filter')}"
    if exists(".genotypes.tsv") or exists(".block_genotypes.tsv") or exists(".mosaic_haplotypes.fa"):
        return "defect", f"{state} left a genotype table or FASTA"
    if state in ("REPAIR_INCOMPLETE", "RESCORE_INCOMPLETE", "TIME_LIMIT", "RESOURCE_REFUSAL",
                 "MODEL_INCOMPLETE", "CANDIDATE_INCOMPLETE", "SPELLING_INCOMPLETE"):
        if reason == "improving_move_tie" and not (exists(".repair_rounds.tsv") or exists(".repair_rounds.partial.tsv")):
            return "defect", "improving-move tie without its repair rounds"
        # The two tie-branched no-calls keep the evidence they were decided on.
        if reason in ("tie_branching_incomplete", "divergent_tie"):
            want = "REPAIR_INCOMPLETE" if reason == "tie_branching_incomplete" else "RESCORE_INCOMPLETE"
            if state != want or status.get("search_conditional") != "tie_branched":
                return "defect", f"{reason} under {state} / search_conditional {status.get('search_conditional')}"
            missing = [s for s in (".repair_rounds.tsv", ".tie_branches.tsv", ".tie_terminal_ties.tsv") if not exists(s)]
            if reason == "divergent_tie" and not exists(".equivalence_set.tsv"):
                missing.append(".equivalence_set.tsv")
            if missing:
                return "defect", f"{reason} without its evidence ({missing})"
        if state == "TIME_LIMIT" and reason == "time_limit":
            return "no_call", "time_limit"
        if state == "RESOURCE_REFUSAL" and reason == "footprint_limit":
            return "no_call", "footprint_limit"
        return "no_call", reason if reason != "NA" else state.lower()
    return "failed", f"status {state}"


def check_provenance(job, status, profile_sha, threads, inputs):
    if status.get("profile_sha256") != profile_sha:
        raise Defect(f"{job['job']}: ran under profile {status.get('profile_sha256')}, not {profile_sha}")
    if f"threads={threads}" not in status.get("settings_overridden_after_profile", ""):
        raise Defect(f"{job['job']}: the scheduler's thread count is not recorded as an override")
    for name, path in (("gfa", job["model_gfa"]), ("reads1", job["reads1"]), ("reads2", job["reads2"]),
                       ("library", job["library"])):
        recorded = status.get(f"sha256_{name}")
        if recorded is not None and recorded != inputs.rows.get(str(path), {}).get("sha256"):
            raise Defect(f"{job['job']}: the call read a {name} that is not the bound input")


def time_fields(path):
    fields = {"wall_seconds": "NA", "user_seconds": "NA", "sys_seconds": "NA", "peak_rss_bytes": "NA",
              "peak_footprint_bytes": "NA"}
    try:
        with open(path) as handle:
            for line in handle:
                parts = line.split()
                if len(parts) >= 6 and parts[1] == "real":
                    fields["wall_seconds"], fields["user_seconds"], fields["sys_seconds"] = parts[0], parts[2], parts[4]
                elif line.strip().endswith("maximum resident set size"):
                    fields["peak_rss_bytes"] = parts[0]
                elif line.strip().endswith("peak memory footprint"):
                    fields["peak_footprint_bytes"] = parts[0]
    except FileNotFoundError:
        pass
    return fields


def run_job(job, build, evaluator, threads):
    """Call, evaluation and (for an improving-move tie) the labelled union fallback. Returns (outcome, detail)."""
    prefix = pathlib.Path(job["out_prefix"])
    prefix.parent.mkdir(parents=True, exist_ok=True)
    cmd = call_command(job, build, threads)
    with open(str(prefix) + ".command.txt", "w") as handle:
        handle.write(" ".join(map(str, cmd)) + "\n")
    with open(str(prefix) + ".log", "w") as log, open(str(prefix) + ".time", "w") as timing:
        try:
            code = subprocess.run(["/usr/bin/time", "-l"] + [str(c) for c in cmd], stdout=log, stderr=timing,
                                  timeout=PROCESS_BACKSTOP_SECONDS).returncode
        except subprocess.TimeoutExpired:
            return "failed", "process backstop timeout"
    if code not in EXTERNAL_EXIT_CODES_OK:
        return "failed", f"exit {code}"
    outcome, detail = classify(job)
    if outcome in ("failed", "pending", "defect"):
        return outcome, detail
    if outcome == "no_call" and detail in ("footprint_limit", "time_limit"):
        return outcome, detail   # the process was stopped: there is no model or call to evaluate
    if evaluate(job, evaluator) != 0 or not os.path.exists(str(prefix) + ".eval.summary.tsv"):
        return "failed", "evaluation failed"
    if outcome == "no_call" and detail == "improving_move_tie":
        status = key_values(str(prefix) + ".mosaic_status.tsv")
        fallback = f"{status.get('union_top_hap1')}|{status.get('union_top_hap2')}"
        if "NA" in fallback or evaluate(job, evaluator, fallback) != 0:
            return "failed", "fallback evaluation failed"
    return outcome, detail


# ---- the campaign -----------------------------------------------------------------------------------

class Ledger:
    """completeness.tsv and attempts.tsv, rewritten whole after every job."""

    FIELDS = ["job", "campaign", "subject", "coverage", "outcome", "detail", "attempts", "status", "call_class",
              "no_call_reason", "phase_support", "wall_seconds", "user_seconds", "sys_seconds", "peak_rss_bytes",
              "peak_footprint_bytes", "finished_utc"]

    def __init__(self, out, jobs):
        self.out = out
        self.lock = threading.Lock()
        self.rows = {j["job"]: {"job": j["job"], "campaign": j["campaign"], "subject": j["subject"],
                                "coverage": j["coverage"], "outcome": "pending", "detail": "", "attempts": "0"}
                     for j in jobs}
        path = out / "completeness.tsv"
        if path.exists():
            with open(path) as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    if row["job"] in self.rows:
                        self.rows[row["job"]].update(row)

    def record(self, job, outcome, detail, attempts):
        status = key_values(job["out_prefix"] + ".mosaic_status.tsv") or {}
        row = self.rows[job["job"]]
        row.update({"outcome": outcome, "detail": detail, "attempts": str(attempts), "status": status.get("status", ""),
                    "call_class": status.get("call_class", ""), "no_call_reason": status.get("no_call_reason", ""),
                    "phase_support": status.get("phase_support", ""), "finished_utc": now()})
        row.update(time_fields(job["out_prefix"] + ".time"))
        self.write()

    def write(self):
        with self.lock:
            tmp = self.out / "completeness.tsv.writing"
            with open(tmp, "w") as handle:
                writer = csv.DictWriter(handle, fieldnames=self.FIELDS, delimiter="\t", lineterminator="\n",
                                        extrasaction="ignore")
                writer.writeheader()
                for name in sorted(self.rows):
                    writer.writerow({k: self.rows[name].get(k, "") for k in self.FIELDS})
            os.replace(tmp, self.out / "completeness.tsv")

    def finished(self, name):
        return self.rows[name]["outcome"] in ("strict_call", "genotype_equivalent_call", "no_call")


def launch(args):
    out = pathlib.Path(args.out).resolve()
    build = pathlib.Path(args.build).resolve()
    evaluator = pathlib.Path(args.evaluator).resolve() if args.evaluator else build / "genotype_probe"
    jobs = read_jobs(out)
    settings = key_values(out / "settings.tsv")
    if settings is None or settings.get("jobs_sha256") != sha256(out / "jobs.tsv"):
        raise Defect("jobs.tsv does not match its plan's settings.tsv")
    profile_sha = sha256(PROFILE)
    binary = {"panvar": sha256(build / "panvar"), "genotype_probe": sha256(evaluator)}
    # The plan names the production profile, binary and scheduler; nothing else may run it.
    bound = {"production_profile": str(PROFILE.relative_to(REPO)), "production_profile_sha256": profile_sha,
             "production_binary_sha256": binary["panvar"], "concurrent_jobs": str(args.jobs),
             "threads_per_call": str(args.threads)}
    for name, value in bound.items():
        if settings.get(name) != value:
            raise Defect(f"launch refused: the plan's {name} is {settings.get(name)}, this launch has {value}")
    launch_file = out / "launch.tsv"
    # The production binary, the profile, the plan and the scheduler are fixed for the whole campaign. The
    # EVALUATOR (genotype_probe) may be corrected between launches -- it never touches a call -- and every
    # version is recorded; jobs it failed are retried, and jobs it measured are re-measured before the end.
    fixed = {"panvar_sha256": binary["panvar"], "profile": str(PROFILE.relative_to(REPO)),
             "profile_sha256": profile_sha, "jobs": str(args.jobs), "threads": str(args.threads),
             "jobs_sha256": settings["jobs_sha256"]}
    previous = key_values(launch_file)
    evaluator_changed = False
    if previous is not None:
        for key, value in fixed.items():
            if previous.get(key) != value:
                raise Defect(f"resume refused: {key} was {previous.get(key)} at launch, now {value}")
        history = [line.rstrip("\n").split("\t")[1] for line in open(launch_file)
                   if line.startswith("evaluator_sha256\t")]
        if binary["genotype_probe"] not in history:
            evaluator_changed = True
            with open(launch_file, "a") as handle:
                handle.write(f"evaluator_sha256\t{binary['genotype_probe']}\nresumed_utc\t{now()}\n")
        else:
            with open(launch_file, "a") as handle:
                handle.write(f"resumed_utc\t{now()}\n")
    else:
        commit = subprocess.run(["git", "-C", str(build.parent), "rev-parse", "HEAD"], capture_output=True,
                                text=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", str(build.parent), "status", "--porcelain"], capture_output=True,
                               text=True).stdout.strip()
        with open(launch_file, "w") as handle:
            handle.write("key\tvalue\n")
            for key, value in fixed.items():
                handle.write(f"{key}\t{value}\n")
            handle.write(f"evaluator_sha256\t{binary['genotype_probe']}\n")
            handle.write(f"build\t{build}\nbuild_commit\t{commit}\nbuild_tree_dirty\t{'true' if dirty else 'false'}\n"
                         f"launched_utc\t{now()}\nhost_physical_memory_bytes\t"
                         f"{subprocess.run(['sysctl', '-n', 'hw.memsize'], capture_output=True, text=True).stdout.strip()}\n")
    inputs = Inputs(out)
    ledger = Ledger(out, jobs)
    if evaluator_changed:
        for row in ledger.rows.values():
            if row["outcome"] == "failed" and "evaluation" in row.get("detail", ""):
                row.update({"outcome": "pending", "attempts": "0",
                            "detail": f"reset: failed only at evaluation, evaluator replaced {now()}"})
        ledger.write()
    ensure_graph(out, inputs)
    if any(j["campaign"] == "mosaic" for j in jobs):
        prepare_mosaic(out, build, inputs, out / "mosaic_prepare.log")
    todo = [j for j in jobs if not ledger.finished(j["job"])]
    print(f"[{now()}] {len(jobs) - len(todo)} of {len(jobs)} jobs already finished; {len(todo)} to run "
          f"with {args.jobs} job(s) x {args.threads} thread(s)", flush=True)
    stop = threading.Event()
    stop_reason = []

    def work(job):
        attempts = int(ledger.rows[job["job"]].get("attempts") or 0)
        while attempts < MAX_ATTEMPTS and not stop.is_set():
            attempts += 1
            outcome, detail = run_job(job, build, evaluator, args.threads)
            if outcome == "defect":
                ledger.record(job, "defect", detail, attempts)
                stop_reason.append(f"{job['job']}: {detail}")
                stop.set()
                return
            if outcome != "failed":
                status = key_values(job["out_prefix"] + ".mosaic_status.tsv")
                try:
                    check_provenance(job, status, profile_sha, args.threads, inputs)
                except Defect as defect:
                    ledger.record(job, "defect", str(defect), attempts)
                    stop_reason.append(str(defect))
                    stop.set()
                    return
                ledger.record(job, outcome, detail, attempts)
                print(f"[{now()}] {job['job']}: {outcome} {detail}", flush=True)
                return
            ledger.record(job, "failed", detail, attempts)
            print(f"[{now()}] {job['job']}: external failure ({detail}), attempt {attempts}", flush=True)
        if not stop.is_set():
            stop_reason.append(f"{job['job']}: repeated external failure")
            stop.set()

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        pending = set()
        for job in todo:
            if stop.is_set():
                break
            # Inputs are prepared here, one job at a time, before the job is handed to a worker.
            try:
                ensure_inputs(job, out, build, inputs)
            except Defect as defect:
                stop_reason.append(str(defect))
                stop.set()
                break
            while len(pending) >= args.jobs:
                done, pending = concurrent.futures.wait(pending, return_when=concurrent.futures.FIRST_COMPLETED)
            pending.add(pool.submit(work, job))
        concurrent.futures.wait(pending)
    ledger.write()
    if stop.is_set():
        print(f"[{now()}] STOPPED: {'; '.join(stop_reason)}", flush=True)
        return 1
    remaining = [n for n in ledger.rows if not ledger.finished(n)]
    print(f"[{now()}] all {len(jobs) - len(remaining)} of {len(jobs)} jobs accounted for", flush=True)
    if remaining:
        return 1
    current = sha256(evaluator)
    for job in jobs:
        row = ledger.rows[job["job"]]
        for kind, fallback in (("eval", None), ("fallback", "union")):
            stamp = pathlib.Path(job["out_prefix"] + f".{kind}.evaluator")
            summary = pathlib.Path(job["out_prefix"] + f".{kind}.summary.tsv")
            if not summary.exists():
                continue   # never evaluated this way (a watchdog stop; or no fallback)
            if stamp.exists() and stamp.read_text().strip() == current:
                continue   # already measured by the final evaluator
            if fallback:
                status = key_values(job["out_prefix"] + ".mosaic_status.tsv")
                fallback = f"{status.get('union_top_hap1')}|{status.get('union_top_hap2')}"
            if evaluate(job, evaluator, fallback) != 0:
                raise Defect(f"{job['job']}: re-evaluation with the final evaluator failed")
            print(f"[{now()}] {job['job']}: re-evaluated ({kind}) with the final evaluator", flush=True)
    if args.consolidate:
        return subprocess.run([sys.executable, str(REPO / "scripts" / "consolidate_c4_production.py"),
                               "--out", str(out)]).returncode
    return 0


def prepare(args):
    out = pathlib.Path(args.out).resolve()
    build = pathlib.Path(args.build).resolve()
    inputs = Inputs(out)
    ensure_graph(out, inputs)
    prepare_mosaic(out, build, inputs, out / "mosaic_prepare.log")
    print((out / "mosaic_inputs" / "prepared.tsv").read_text())
    return 0


def status(args):
    out = pathlib.Path(args.out).resolve()
    ledger = Ledger(out, read_jobs(out))
    counts = {}
    for row in ledger.rows.values():
        label = row["outcome"] + (f":{row['detail']}" if row["outcome"] in ("no_call", "failed", "defect") else "")
        counts[label] = counts.get(label, 0) + 1
    for label in sorted(counts):
        print(f"{label}\t{counts[label]}")
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="mode", required=True)
    p = sub.add_parser("launch")
    p.add_argument("--build", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--jobs", type=int, required=True)
    p.add_argument("--threads", type=int, required=True)
    p.add_argument("--no-consolidate", dest="consolidate", action="store_false")
    # The evaluator binary, when it is not the build's own genotype_probe (it may be corrected mid-campaign;
    # the production binary may not).
    p.add_argument("--evaluator")
    p = sub.add_parser("prepare")
    p.add_argument("--build", required=True)
    p.add_argument("--out", required=True)
    p = sub.add_parser("status")
    p.add_argument("--out", required=True)
    args = parser.parse_args()
    try:
        return {"launch": launch, "prepare": prepare, "status": status}[args.mode](args)
    except Defect as defect:
        print(f"[{now()}] DEFECT: {defect}", flush=True)
        return 2


if __name__ == "__main__":
    sys.exit(main())
