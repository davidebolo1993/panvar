#!/usr/bin/env python3
"""M8-D continuation analysis (frozen with DESIGN.md).

  analyse.py hashcheck <out>   the frozen M8-C inputs (executables, graph, reads, library, curation, profiles) by SHA-256
  analyse.py freeze <out>      after both calls: prediction outcome per run (COMPLETE with locus genotype, else INCOMPLETE)
                               and SHA-256 of its authoritative outputs, before any truth is read
  analyse.py contrast <out>    if both runs are COMPLETE and their unordered sequence pairs differ: search reach and
                               fragment / exposure / prior of both winners under both runs (existing score-set)
  analyse.py summary <out>     corrected paired semantics (both complete: compare unordered sequence pairs; either
                               incomplete: same_call = UNKNOWN), with both-missing and one-missing fixtures; the six M8-C
                               completed cells plus these two; the M8-C partial repair traces (search moves) for the record
"""
import csv
import hashlib
import os
import subprocess
import sys

mode, out = sys.argv[1], sys.argv[2]
R = "/Users/davide.bolognini/panvar/tests/results"
WT = "/Users/davide.bolognini/panvar-worktrees/s1-d-panel-references"
C = f"{R}/m8c_copying_challenge_20261020"
N = f"{R}/m8d_hg01975_loo_20261020/runs"
P = f"{C}/runs/prepare.a1"
DN = "HG01975"
ACOT_GFA = f"{R}/genotype_perf_bounded_20261015/inputs/acot/acot.gfa"
ACOT_REF = "GRCh38#0#chr14:73416151-73668669"
PROBE, PANVAR = f"{WT}/build/genotype_probe", f"{WT}/build/panvar"
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
T = lambda p: list(csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t")) if os.path.exists(p) else []
KV = lambda p: {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(p) if "\t" in l and not l.startswith("#")} if os.path.exists(p) else {}
os.makedirs(out, exist_ok=True)


def write(name, rows):
    if rows:
        with open(f"{out}/{name}.tsv", "w") as f:
            w = csv.DictWriter(f, fieldnames=list(dict.fromkeys(k for r in rows for k in r)), delimiter="\t", lineterminator="\n", restval="")
            w.writeheader()
            w.writerows(rows)


if mode == "hashcheck":
    rec = {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(f"{P}/prepare.tsv")}
    rows = [dict(item="panvar", recorded=rec["panvar_sha256"], now=sha(PANVAR)), dict(item="genotype_probe", recorded=rec["genotype_probe_sha256"], now=sha(PROBE)),
            dict(item="acot.gfa", recorded=rec["acot_gfa_sha256"], now=sha(ACOT_GFA))]
    reads = rec[f"reads_{DN}"].split()
    for k, f in (("reads_1", "reads_1.fq"), ("reads_2", "reads_2.fq"), ("library", "library.txt")):
        rows.append(dict(item=f"{DN} cov5 {f}", recorded=reads[reads.index(k) + 1] + "...", now=sha(f"{P}/reads/{DN}/cov5/{f}")))
    # curation and profiles were not hashed by M8-C's prepare; re-derive them with its frozen rules and require byte equality
    sys.path.insert(0, f"{WT}/tests/analysis/p1")
    import compose_profile  # noqa: E402
    split = [l.rstrip("\n").split("\t") for l in open(f"{R}/m8b_copying_calibration_20261020/runs/inventory/split.ACOT.tsv")]
    rows.append(dict(item="split.ACOT.tsv", recorded=rec["split_sha256"], now=sha(f"{R}/m8b_copying_calibration_20261020/runs/inventory/split.ACOT.tsv")))
    eval_t = [x[0] for x in split if x[2] == "eval"]
    want = {f"curation.LOO.{DN}.txt": "# M8-C LOO calling panel: every outer-evaluation template excluded\n" + "\n".join(eval_t) + "\n"}
    compose_profile.compose(f"{out}/profile.base.rederived.tsv", "m8c-base-acot", f"{R}/genotype_perf_bounded_20261015/inputs/profiles/acot_5.tsv", "full", [])
    base = open(f"{out}/profile.base.rederived.tsv").read().rstrip("\n").split("\n")
    for tag, r in (("r0.01", "0.01"), ("r0.1", "0.1")):
        lines = [l if not l.startswith("profile\t") else f"profile\tm8c-acot-{tag}" for l in base]
        want[f"profile.{tag}.tsv"] = "\n".join(lines + ["# M8-C: the only scientific difference between the two arms", "copying_transition\tconstant_v1", f"switch_probability\t{r}"]) + "\n"
    for f, text in want.items():
        rows.append(dict(item=f, recorded=hashlib.sha256(text.encode()).hexdigest(), now=sha(f"{P}/{f}"), basis="re-derived by M8-C prepare rule"))
    # the original commands must reference exactly these inputs (only --out-prefix differs)
    for r in ("r0.01", "r0.1"):
        cmd = open(f"{C}/runs/call.{DN}.LOO.{r}/run.log").readline()
        rows.append(dict(item=f"M8-C command {r}", recorded="profile+curation+reads paths", now="match" if all(x in cmd for x in (
            f"{P}/profile.{r}.tsv", f"{P}/curation.LOO.{DN}.txt", f"{P}/reads/{DN}/cov5/reads_1.fq", "--dosage-mode repair", "--threads 4")) else "MISMATCH"))
    ok = True
    for r in rows:
        r["ok"] = r["now"] == "match" if r["item"].startswith("M8-C command") else r["now"].startswith(r["recorded"].rstrip("."))
        ok &= r["ok"]
    write("hashcheck", rows)
    print(open(f"{out}/hashcheck.tsv").read())
    sys.exit(0 if ok else 1)

runs = {r: f"{N}/call.{DN}.LOO.{r}" for r in ("r0.01", "r0.1")}

if mode == "freeze":
    rows = []
    for r, d in runs.items():
        lg = T(f"{d}/call.locus_genotype.tsv")
        complete = bool(lg) and lg[0]["run_status"] == "COMPLETE"
        row = dict(run=f"{DN}.LOO.{r}", outcome="COMPLETE" if complete else "INCOMPLETE (no completed genotype; never evaluated)",
                   run_status=lg[0]["run_status"] if lg else "no locus genotype")
        row["switch_probability_settings"] = next((l.split("\t")[1] for l in open(f"{d}/call.settings.tsv") if l.startswith("switch_probability\t")), "absent") \
            if os.path.exists(f"{d}/call.settings.tsv") else "absent"
        for f in ("call.locus_genotype.tsv", "call.equivalence_set.tsv", "call.mosaic_haplotypes.fa"):
            row[f] = sha(f"{d}/{f}")[:16] if os.path.exists(f"{d}/{f}") else "absent"
        rows.append(row)
    write("freeze", rows)
    print(open(f"{out}/freeze.tsv").read())
    sys.exit(0)

if mode == "contrast":
    fz = {r["run"]: r for r in T(f"{N}/freeze/freeze.tsv")}
    rows = []
    calls = {r: T(f"{d}/call.locus_genotype.tsv")[0] for r, d in runs.items() if fz.get(f"{DN}.LOO.{r}", {}).get("outcome") == "COMPLETE"}
    if len(calls) < 2:
        rows.append(dict(status="not assessable (a run is incomplete)"))
    elif frozenset([calls["r0.01"]["hap1_sha256"], calls["r0.01"]["hap2_sha256"]]) == frozenset([calls["r0.1"]["hap1_sha256"], calls["r0.1"]["hap2_sha256"]]):
        rows.append(dict(status="UNCHANGED (same unordered sequence pair)"))
    else:
        reads = f"{P}/reads/{DN}/cov5"
        for r, d in runs.items():
            cand = {c["candidate_id"]: c for c in T(f"{d}/call.mosaic_candidates.tsv")}
            scored = set()
            for s in T(f"{d}/call.mosaic_scores.tsv"):
                c = cand.get(s["representative_beam_candidate"])
                if c:
                    scored.add(frozenset([c["hap1_alleles"], c["hap2_alleles"]]))
            for x in T(f"{d}/call.repair_rounds.tsv"):   # states the bounded repair moved to
                scored.add(frozenset([x["hap1_alleles"], x["hap2_alleles"]]))
            cf = f"{out}/{r}.candidates.tsv"
            with open(cf, "w") as f:
                for wr in ("r0.01", "r0.1"):
                    c = calls[wr]
                    f.write(f"winner_{wr}\t{c['hap1_block_path']}\t{c['hap2_block_path']}\t{c['total_score'] if wr == r else 'NA'}\n")
            code = subprocess.run([PROBE, "score-set", "--out", f"{out}/{r}.fs", "--model-gfa", ACOT_GFA, "--reference", ACOT_REF, "--donor", DN, "--run", f"{d}/call",
                                   "--reads1", f"{reads}/reads_1.fq", "--reads2", f"{reads}/reads_2.fq", "--library", f"{reads}/library.txt", "--threads", "4",
                                   "--switch-probability", r[1:], "--candidates", cf, "--panel-curation", f"{P}/curation.LOO.{DN}.txt"], capture_output=True, text=True)
            fs = {x["label"]: x for x in T(f"{out}/{r}.fs.scores_truth_free.tsv")}
            for wr in ("r0.01", "r0.1"):
                c = calls[wr]
                x = fs.get(f"winner_{wr}", {})
                rows.append(dict(status="CHANGED", scored_in_run=r, candidate=f"winner_{wr}",
                                 reached_by_this_run=frozenset([c["hap1_block_path"], c["hap2_block_path"]]) in scored,
                                 fragment_term=x.get("fragment_term", ""), exposure_term=x.get("exposure_term", ""), prior=x.get("prior", ""), total=x.get("total", ""),
                                 score_set_exit=code.returncode, note=(code.stdout + code.stderr).strip().splitlines()[-1][:160] if (code.stdout + code.stderr).strip() else ""))
    write("contrast", rows)
    print(open(f"{out}/contrast.tsv").read())
    sys.exit(0)

if mode == "summary":
    def same_call(a, b):
        """Corrected semantics: both complete -> compare unordered sequence pairs; either incomplete -> UNKNOWN."""
        if not (a and b and a.get("complete") and b.get("complete")):
            return "UNKNOWN"
        return frozenset([a["h1"], a["h2"]]) == frozenset([b["h1"], b["h2"]])
    fixtures = [("both missing", same_call(None, None), "UNKNOWN"), ("one missing", same_call({"complete": True, "h1": "x", "h2": "y"}, None), "UNKNOWN"),
                ("one incomplete", same_call({"complete": True, "h1": "x", "h2": "y"}, {"complete": False}), "UNKNOWN"),
                ("both complete, same pair swapped", same_call({"complete": True, "h1": "x", "h2": "y"}, {"complete": True, "h1": "y", "h2": "x"}), True),
                ("both complete, different", same_call({"complete": True, "h1": "x", "h2": "y"}, {"complete": True, "h1": "x", "h2": "z"}), False)]
    fx_ok = all(got == want for _, got, want in fixtures)
    write("fixtures", [dict(case=c, got=g, expected=w, result="PASS" if g == w else "FAIL") for c, g, w in fixtures])
    old = T(f"{C}/runs/summary.a2/cells.tsv")
    led = {r["unit"]: r for r in T(f"{N}/../ledger.tsv")}
    fz = {r["run"]: r for r in T(f"{N}/freeze/freeze.tsv")}
    cells = {}
    for row in old:
        if row["cell"].startswith(f"{DN}.LOO."):
            row = dict(row, note="M8-C TIMEOUT at the 120-s unit bound (historical; superseded by M8-D for the comparison)")
            cells[row["cell"] + " (M8-C)"] = row
        else:
            cells[row["cell"]] = row
    for r, d in runs.items():
        lg = T(f"{d}/call.locus_genotype.tsv")
        ev = KV(f"{N}/eval.{r}/eval.summary.tsv")
        st = KV(f"{d}/call.mosaic_status.tsv")
        row = dict(cell=f"{DN}.LOO.{r}", note="M8-D continuation (930-s watchdog, caller guard 900 s)", call_exit=led.get(f"call.{r}", {}).get("exit", "not run"),
                   call_seconds=led.get(f"call.{r}", {}).get("wall_seconds", ""), call_peak_rss_gib=led.get(f"call.{r}", {}).get("peak_rss_gib", ""),
                   eval_exit=led.get(f"eval.{r}", {}).get("exit", "not run"), outcome=fz.get(f"{DN}.LOO.{r}", {}).get("outcome", ""))
        if lg:
            c = lg[0]
            row.update(run_status=c["run_status"], call_class=c["call_class"], phase_support=c["phase_support"], hap1_sha=c["hap1_sha256"][:16], hap2_sha=c["hap2_sha256"][:16],
                       total_score=c["total_score"], fragment_term=c["fragment_term"], exposure_term=c["exposure_term"], prior=c["sequence_pair_prior"],
                       switches=f"{c.get('hap1_switches')}+{c.get('hap2_switches')}")
        row.update({k: ev.get(k, "") for k in ("status", "availability_floor_total", "unordered_distance_total", "excess_over_floor_total", "phased_distance_total",
                                                "phase_switch_blocks", "blocks_exact", "blocks", "sv_blocks_exact", "undefined_truth_blocks", "truth_candidate_source",
                                                "floor_attained_in_candidates", "winner_excess_over_floor")})
        row.update({k: st.get(k, "") for k in ("sequence_candidates", "mosaic_candidates", "proposal_states_pruned", "seconds_total", "repair_rounds_run",
                                                "numeric_separation_gap", "winner_candidate_set_probability")})
        cells[row["cell"]] = row
    write("cells", list(cells.values()))
    # corrected paired table over completed outcomes (M8-C pairs + the M8-D pair)
    paired = []
    hap = lambda x: {"complete": x.get("run_status") == "COMPLETE", "h1": x.get("hap1_sha") or x.get("hap1"), "h2": x.get("hap2_sha") or x.get("hap2")}
    for dn, sch in (("HG03209", "LZO"), ("HG03209", "LOO"), ("HG01975", "LZO"), ("HG01975", "LOO")):
        a, b = cells.get(f"{dn}.{sch}.r0.01", {}), cells.get(f"{dn}.{sch}.r0.1", {})
        sc = same_call(hap(a) if a else None, hap(b) if b else None)
        num = lambda x: float(x["excess_over_floor_total"]) if x.get("excess_over_floor_total") not in ("", None, "NA") else None
        ea, eb = num(a), num(b)
        verdict = ("not assessable" if sc == "UNKNOWN" or ea is None or eb is None else "both exact at floor (non-discriminating)" if ea == 0 and eb == 0
                   else "r=0.1 better" if eb < ea else "r=0.1 worse" if eb > ea else "equal excess")
        paired.append(dict(donor=dn, scheme=sch, same_call=sc, excess_r001=ea, excess_r01=eb, floor=a.get("availability_floor_total", ""), verdict=verdict))
    write("paired_corrected", paired)
    # M8-C partial search traces (record only)
    tr = []
    for r in ("r0.01", "r0.1"):
        moves = [x["repaired_block"] for x in T(f"{C}/runs/call.{DN}.LOO.{r}/call.repair_rounds.partial.tsv") if x.get("moved") == "true"]
        tr.append(dict(run=f"M8-C {DN}.LOO.{r} (partial)", moved_blocks_in_order=",".join(moves)))
    for r, d in runs.items():
        moves = [x["repaired_block"] for x in T(f"{d}/call.repair_rounds.tsv") if x.get("moved") == "true"]
        tr.append(dict(run=f"M8-D {DN}.LOO.{r}", moved_blocks_in_order=",".join(moves)))
    write("search_traces", tr)
    for n in ("fixtures", "paired_corrected", "search_traces"):
        print(open(f"{out}/{n}.tsv").read())
    sys.exit(0 if fx_ok else 1)
