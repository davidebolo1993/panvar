#!/usr/bin/env python3
"""M15-C: the documented examples of docs/guides/genotype_short_read.md, run and checked (development controls only; no
reserved family). Each example runs its command(s), asserts the expected exit status and outputs, and exits 0 only if every
expectation holds.  Usage: examples.py <example> <unit dir>"""
import csv
import hashlib
import os
import shutil
import subprocess
import sys

csv.field_size_limit(sys.maxsize)
WT = "/Users/davide.bolognini/panvar-worktrees/s1-d-panel-references"
R = "/Users/davide.bolognini/panvar/tests/results"
PANVAR = f"{WT}/build/panvar"
M11D = f"{R}/m11d_dev_pilot_20261020/runs"
CELLS = {
    "c4": dict(graph=f"{R}/genotype_c4_production_11fe681/graphs/NA19705.loo.gfa", ref="HG00096#1#haplotype1-0000024:31848049-32074186",
               cur=None, reads=f"{R}/genotype_c4_production_11fe681/reads/NA19705/cov5", archived=f"{M11D}/call.c4.NA19705/call",
               profile=f"{WT}/profiles/genotype_experimental_sr_v1_c4.tsv"),
    "cyp2d6": dict(graph=f"{R}/cyp2d6_shortread_pilot_v1/data/graphs/loo.HG02769.gfa", ref="grch38#1#chr22:42031864-42245566",
                   cur=f"{WT}/tests/genotype_cyp2d6_curation.tsv", reads=f"{R}/cyp2d6_shortread_pilot_v1/data/reads/HG02769/cov5",
                   archived=f"{M11D}/call.cyp.HG02769/call", profile=f"{WT}/profiles/genotype_experimental_sr_v1_cyp2d6.tsv"),
    "acot": dict(graph=f"{R}/genotype_perf_bounded_20261015/inputs/acot/acot.gfa", ref="GRCh38#0#chr14:73416151-73668669",
                 cur=f"{R}/m9h_jp11_validation_20261020/runs/prepare/curation.LOO.txt",
                 reads=f"{R}/m8c_copying_challenge_20261020/runs/prepare.a1/reads/HG03209/cov5", archived=f"{M11D}/call.acot.HG03209/call",
                 profile=f"{WT}/profiles/genotype_experimental_sr_v1_acot.tsv"),
}


def T(p):
    return list(csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t")) if os.path.exists(p) else []


def KV(p):
    return dict(l.rstrip("\n").split("\t", 1) for l in open(p) if "\t" in l and not l.startswith("#")) if os.path.exists(p) else {}


def base_cmd(cell, out, profile=None):
    c = CELLS[cell]
    cmd = [PANVAR, "genotype", "--profile", profile or c["profile"], "--threads", "4", "--graph", c["graph"], "--reference-path", c["ref"]]
    if c["cur"]:
        cmd += ["--panel-curation", c["cur"]]
    return cmd + ["--reads1", f"{c['reads']}/reads_1.fq", "--reads2", f"{c['reads']}/reads_2.fq", "--library", f"{c['reads']}/library.txt",
                  "--out", out, "--quiet"]


def run(cmd):
    print("$", " ".join(cmd).replace(R, "<results>").replace(WT, "<worktree>"))
    p = subprocess.run(cmd, capture_output=True, text=True)
    print(f"exit {p.returncode}")
    if p.stderr.strip():
        print("stderr:", p.stderr.strip()[-600:])
    return p


def prediction_key(prefix):
    win = next(r for r in T(prefix + ".equivalence_set.tsv") if r["rank"] == "1" and r["representative"] == "true")
    hs = {r["haplotype"]: (r["sha256"], r["origin"]) for r in T(prefix + ".sequence_haplotypes.tsv")}
    st = KV(prefix + ".mosaic_status.tsv")
    return (win["score"], win["hap1_alleles"], win["hap2_alleles"], tuple(sorted(hs.items())), st.get("status"), st.get("model_fingerprint"))


def settings(prefix):
    return {r["key"]: r["value"] for r in T(prefix + ".settings.tsv")}


CHECKS = []


def check(ok, what):
    CHECKS.append((ok, what))
    print(("PASS " if ok else "FAIL ") + what)


def standard(cell, d):
    out = f"{d}/{cell}"
    p = run(base_cmd(cell, out))
    check(p.returncode == 0, f"{cell}: exit 0 (PREDICTION)")
    sl = KV(out + ".sequence_locus.tsv")
    check(sl.get("prediction") == "PREDICTION" and sl.get("completion") == "completed_ranking", f"{cell}: completed prediction")
    a = CELLS[cell]["archived"]
    check(prediction_key(out) == prediction_key(a), f"{cell}: prediction identical to the archived development control (M11-D)")
    new, old = settings(out), settings(a)
    expected_diff = {"time_limit_seconds", "block_diagnostics"}   # the archived runs declared a watchdog guard and diagnostics
    diff = sorted(k for k in set(new) | set(old) if new.get(k) != old.get(k) and k not in expected_diff)
    check(not diff, f"{cell}: every other effective setting equals the archived run's ({diff or 'no difference'})")
    print(f"time_limit_seconds new {new.get('time_limit_seconds')} archived {old.get('time_limit_seconds')}; "
          f"dosage {KV(out + '.dosage_status.tsv').get('dosage_mode_source', KV(out + '.dosage_status.tsv').get('dosage_audit'))}")
    return out


def ex_standard_c4(d):
    out = standard("c4", d)
    # the off-panel FASTA case: hap2 OFF_PANEL -> one FASTA with both records, origins annotated, bytes = the table's SHA-256
    hs = {r["haplotype"]: r for r in T(out + ".sequence_haplotypes.tsv")}
    fa = out + ".sequence_haplotypes.fa"
    check(any(r["origin"] == "OFF_PANEL" for r in hs.values()) and os.path.exists(fa), "c4: an off-panel haplotype writes sequence_haplotypes.fa")
    recs, name = {}, None
    for line in open(fa) if os.path.exists(fa) else []:
        line = line.rstrip("\n")
        if line.startswith(">"):
            name = line[1:].split()[0]
            recs[name] = [line, ""]
        elif name:
            recs[name][1] += line
    ok = len(recs) == 2 and all(hashlib.sha256(v[1].encode()).hexdigest() in {r["sha256"] for r in hs.values()} for v in recs.values())
    check(ok, "c4: both FASTA records match the SHA-256 in sequence_haplotypes.tsv")
    check(not os.path.exists(out + ".mosaic_haplotypes.fa"), "c4: the legacy FASTA is suppressed under v2")
    return 0


def ex_standard_cyp2d6(d):
    standard("cyp2d6", d)
    return 0


def ex_standard_acot(d):
    out = standard("acot", d)
    hs = T(out + ".sequence_haplotypes.tsv")
    check(all(r["origin"] == "PANEL" for r in hs) and not os.path.exists(out + ".sequence_haplotypes.fa"),
          "acot: both haplotypes on the panel -> no FASTA (PANEL_REPRESENTED_PAIR)")
    p = run([sys.executable, f"{WT}/tests/analysis/s1d/extract_panel_haplotype.py", "--check", CELLS["acot"]["graph"], out + ".sequence_haplotypes.tsv"])
    print(p.stdout[-800:])
    check(p.returncode == 0, "acot: panel sequences recovered from the graph by extraction reference and checked (hash, length)")
    return 0


def ex_resource_adjusted(d):
    prof = f"{d}/genotype_experimental_sr_v1_c4.resource.tsv"
    keep = {"time_limit_seconds": "120", "footprint_limit_bytes": "4294967296", "emission_cache_bytes": "268435456"}
    lines = []
    for line in open(CELLS["c4"]["profile"]):
        k = line.split("\t")[0]
        if k == "profile":
            line = "profile\tgenotype-experimental-sr-v1-c4-resource-adjusted\n"
        elif k in keep:
            line = f"{k}\t{keep[k]}\n"
        lines.append(line)
    open(prof, "w").writelines(lines)
    out = f"{d}/c4.resource"
    p = run(base_cmd("c4", out, prof))
    check(p.returncode == 0, "resource-adjusted: exit 0")
    check(prediction_key(out) == prediction_key(CELLS["c4"]["archived"]), "resource-adjusted: prediction identical (limits did not bind)")
    s = settings(out)
    check(all(s.get(k) == v for k, v in keep.items()), f"resource-adjusted: effective {sorted(keep)} recorded")
    return 0


def ex_time_limit_binding(d):
    prof = f"{d}/acot.tight.tsv"
    lines = [l if l.split("\t")[0] not in ("profile", "time_limit_seconds") else
             ("profile\tgenotype-experimental-sr-v1-acot-tight\n" if l.startswith("profile\t") else "time_limit_seconds\t8\n")
             for l in open(CELLS["acot"]["profile"])]
    open(prof, "w").writelines(lines)
    out = f"{d}/acot.tight"
    p = run(base_cmd("acot", out, prof))
    sl, st = KV(out + ".sequence_locus.tsv"), KV(out + ".mosaic_status.tsv")
    print(f"prediction {sl.get('prediction')} status {st.get('status')} stop {st.get('stop_reason', st.get('resource_stop', 'NA'))}")
    check(p.returncode in (3, 124), f"binding time limit: exit {p.returncode} (3 PROVISIONAL after a checkpoint, or 124 before one)")
    check((p.returncode == 3 and sl.get("prediction") == "PROVISIONAL") or (p.returncode == 124 and sl.get("prediction") != "PREDICTION"),
          "binding time limit: never reported as a completed prediction")
    return 0


def ex_invalid(d):
    c4 = CELLS["c4"]["profile"]

    def refused(label, cmd, needle, prefix):
        p = run(cmd)
        msg = (p.stderr + p.stdout)
        left = [f for f in os.listdir(d) if f.startswith(os.path.basename(prefix) + ".")]
        check(p.returncode == 2 and needle in msg and not left, f"refused ({label}): exit 2, message names it, nothing written at the prefix")

    refused("dosage repair under v2", base_cmd("c4", f"{d}/r1") + ["--dosage-mode", "repair"], "sequence_first_v2", f"{d}/r1")
    legacy = f"{d}/legacy_early.tsv"
    open(legacy, "w").write("profile\tbad-early-legacy\noutput_schema\tlegacy\nearly_incumbent\texperimental_v1\n")
    refused("early_incumbent without sequence_first_v2", base_cmd("c4", f"{d}/r2", legacy), "early_incumbent", f"{d}/r2")
    refused("adaptive repair off with rounds", base_cmd("c4", f"{d}/r3") + ["--adaptive-repair", "off", "--repair-rounds", "3"], "repair", f"{d}/r3")
    dup = f"{d}/dup.tsv"
    open(dup, "w").write("profile\tdup\nmismatch_band\t6\nmismatch_band\t6\n")
    refused("key set twice in one profile", base_cmd("c4", f"{d}/r4", dup), "set twice", f"{d}/r4")
    unk = f"{d}/unknown.tsv"
    open(unk, "w").write("profile\tunknown\nmemory_per_thread_bytes\t1000\n")
    refused("unknown profile key (no memory-per-thread budget exists)", base_cmd("c4", f"{d}/r5", unk), "memory_per_thread_bytes", f"{d}/r5")
    noname = f"{d}/noname.tsv"
    open(noname, "w").write("mismatch_band\t6\n")
    refused("profile that does not name itself", base_cmd("c4", f"{d}/r6", noname), "does not name itself", f"{d}/r6")
    refused("removed --max-candidates", base_cmd("c4", f"{d}/r7") + ["--max-candidates", "5"], "--max-mosaic-candidates", f"{d}/r7")
    return 0


if __name__ == "__main__":
    name, d = sys.argv[1], sys.argv[2]
    os.makedirs(d, exist_ok=True)
    globals()["ex_" + name](d)
    bad = [w for ok, w in CHECKS if not ok]
    print(f"{len(CHECKS) - len(bad)} of {len(CHECKS)} checks passed" + (f"; FAILED: {bad}" if bad else ""))
    sys.exit(1 if bad else 0)
