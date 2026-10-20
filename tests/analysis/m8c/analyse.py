#!/usr/bin/env python3
"""M8-C analysis (frozen with DESIGN.md).

  analyse.py checks <runs dir> <out dir>     every run: every recorded switch_probability equals its arm's r; copying
                                             transition rows/digests (if any) agree within a run; byte-identical reads and
                                             library per donor across its four runs (from the run commands)
  analyse.py contrast <runs dir> <out dir>   per donor x scheme whose reported unordered sequence pair changes between arms:
                                             whether each arm's winner was scored by the other run (search reach), and
                                             fragment / exposure / prior of both winners under both runs (existing score-set)
  analyse.py summary <runs dir> <out dir>    per cell: evaluation metrics, call class, sequence-based panel membership and
                                             truth identity of each called haplotype, search work, runtime, memory; paired
                                             r = 0.1 vs r = 0.01 per donor x scheme
"""
import csv
import glob
import hashlib
import os
import re
import subprocess
import sys

mode, runs, out = sys.argv[1:4]
R = "/Users/davide.bolognini/panvar/tests/results"
WT = "/Users/davide.bolognini/panvar-worktrees/s1-d-panel-references"
PROBE = f"{WT}/build/genotype_probe"
ACOT_GFA = f"{R}/genotype_perf_bounded_20261015/inputs/acot/acot.gfa"
ACOT_REF = "GRCh38#0#chr14:73416151-73668669"
T = lambda p: list(csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t"))
KV = lambda p: {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(p) if "\t" in l and not l.startswith("#")} if os.path.exists(p) else {}
os.makedirs(out, exist_ok=True)
prep = f"{runs}/prepare"
donors = [l.split("\t")[1] for l in open(f"{prep}/prepare.tsv") if l.startswith("donor_")]
CELLS = [(dn, sch, r) for dn in donors for sch in ("LZO", "LOO") for r in ("r0.01", "r0.1")]
cell_dir = lambda dn, sch, r: f"{runs}/call.{dn}.{sch}.{r}"


def write(name, rows):
    if rows:
        with open(f"{out}/{name}.tsv", "w") as f:
            w = csv.DictWriter(f, fieldnames=list(dict.fromkeys(k for r in rows for k in r)), delimiter="\t", lineterminator="\n", restval="")
            w.writeheader()
            w.writerows(rows)


if mode == "checks":
    rows, ok = [], True
    for dn, sch, r in CELLS:
        d = cell_dir(dn, sch, r)
        want = float(r[1:])
        vals, digests = [], set()
        for f in glob.glob(f"{d}/call*.tsv"):
            for line in open(f, errors="ignore"):
                k = line.split("\t")[0]
                if k == "switch_probability":
                    vals.append((os.path.basename(f), line.rstrip("\n").split("\t")[1]))
                if "transition" in k and "digest" in k:
                    digests.add((k, line.rstrip("\n").split("\t")[1]))
        cmd = open(f"{d}/run.log").readline()
        reads = re.findall(r"--(reads1|reads2|library) (\S+)", cmd)
        good = bool(vals) and all(abs(float(v) - want) < 1e-12 for _, v in vals)
        ok &= good
        rows.append(dict(cell=f"{dn}.{sch}.{r}", switch_probability_records=";".join(f"{a}={b}" for a, b in vals), all_equal_arm_r=good,
                         transition_digests=";".join(f"{a}={b[:16]}" for a, b in sorted(digests)) or "none recorded",
                         inputs=";".join(f"{k}:{hashlib.sha256(open(p, 'rb').read()).hexdigest()[:12]}" for k, p in reads)))
    by_donor = {}
    for row in rows:
        by_donor.setdefault(row["cell"].split(".")[0], set()).add(row["inputs"])
    same = all(len(v) == 1 for v in by_donor.values())
    write("checks", rows + [dict(cell="ALL", all_equal_arm_r=ok, inputs=f"byte-identical reads+library across each donor's 4 runs: {same}")])
    print(open(f"{out}/checks.tsv").read())
    sys.exit(0 if ok and same else 1)

if mode == "contrast":
    rows = []
    for dn in donors:
        for sch in ("LZO", "LOO"):
            calls = {}
            for r in ("r0.01", "r0.1"):
                lg = T(f"{cell_dir(dn, sch, r)}/call.locus_genotype.tsv")
                if lg:
                    calls[r] = lg[0]
            if len(calls) < 2:
                rows.append(dict(donor=dn, scheme=sch, status="incomplete (a run missing)"))
                continue
            pair = lambda c: frozenset([c["hap1_sha256"], c["hap2_sha256"]])
            if pair(calls["r0.01"]) == pair(calls["r0.1"]):
                rows.append(dict(donor=dn, scheme=sch, status="UNCHANGED (same unordered sequence pair)"))
                continue
            cur = f"{prep}/curation.{sch}.{dn}.txt"
            reads = f"{prep}/reads/{dn}/cov5"
            for r in ("r0.01", "r0.1"):
                d = cell_dir(dn, sch, r)
                cand = {c["candidate_id"]: c for c in T(f"{d}/call.mosaic_candidates.tsv")}
                scored = set()
                for s in T(f"{d}/call.mosaic_scores.tsv"):
                    c = cand.get(s["representative_beam_candidate"])
                    if c:
                        scored.add(frozenset([c["hap1_alleles"], c["hap2_alleles"]]))
                cf = f"{out}/{dn}.{sch}.{r}.candidates.tsv"
                with open(cf, "w") as f:
                    for wr in ("r0.01", "r0.1"):
                        c = calls[wr]
                        f.write(f"winner_{wr}\t{c['hap1_block_path']}\t{c['hap2_block_path']}\t{c['total_score'] if wr == r else 'NA'}\n")
                code = subprocess.run([PROBE, "score-set", "--out", f"{out}/{dn}.{sch}.{r}.fs", "--model-gfa", ACOT_GFA, "--reference", ACOT_REF, "--donor", dn,
                                       "--run", f"{d}/call", "--reads1", f"{reads}/reads_1.fq", "--reads2", f"{reads}/reads_2.fq", "--library", f"{reads}/library.txt",
                                       "--threads", "4", "--candidates", cf, "--panel-curation", cur], capture_output=True, text=True)
                fs = {x["label"]: x for x in T(f"{out}/{dn}.{sch}.{r}.fs.scores_truth_free.tsv")} if code.returncode == 0 else {}
                for wr in ("r0.01", "r0.1"):
                    c = calls[wr]
                    x = fs.get(f"winner_{wr}", {})
                    rows.append(dict(donor=dn, scheme=sch, status="CHANGED", scored_in_run=r, candidate=f"winner_{wr}",
                                     scored_by_this_run_search=frozenset([c["hap1_block_path"], c["hap2_block_path"]]) in scored,
                                     fragment_term=x.get("fragment_term", ""), exposure_term=x.get("exposure_term", ""), prior=x.get("prior", ""), total=x.get("total", ""),
                                     score_set_exit=code.returncode, score_set_note=(code.stdout + code.stderr).strip().splitlines()[-1][:160] if (code.stdout + code.stderr).strip() else ""))
    write("contrast", rows)
    print(open(f"{out}/contrast.tsv").read() if rows else "no contrasts")

if mode == "summary":
    def spell(path_name, segs, paths):
        out_ = []
        for st in paths[path_name]:
            s = segs[st[:-1]]
            out_.append(s if st[-1] == "+" else s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1])
        return "".join(out_).upper()
    segs, paths = {}, {}
    for line in open(ACOT_GFA):
        if line.startswith("S\t"):
            f = line.rstrip("\n").split("\t")
            segs[f[1]] = f[2]
        elif line.startswith("P\t"):
            f = line.rstrip("\n").split("\t")
            paths[f[1]] = f[2].split(",")
    rc = lambda s: s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]
    shaset = {}
    for p in paths:
        s = spell(p, segs, paths)
        shaset[hashlib.sha256(s.encode()).hexdigest()] = p
        shaset[hashlib.sha256(rc(s).encode()).hexdigest()] = p
    led = {r["unit"]: r for r in T(f"{runs}/../ledger.tsv")}
    idx = {dn: f"d{k + 1}" for k, dn in enumerate(donors)}   # ledger units use d1/d2
    rows = []
    for dn, sch, r in CELLS:
        d, e = cell_dir(dn, sch, r), f"{runs}/eval.{dn}.{sch}.{r}"
        excluded = {l.strip() for l in open(f"{prep}/curation.{sch}.{dn}.txt") if l.strip() and not l.startswith("#")}
        lg = T(f"{d}/call.locus_genotype.tsv")
        ev = KV(f"{e}/eval.summary.tsv")
        st = KV(f"{d}/call.mosaic_status.tsv")
        row = dict(cell=f"{dn}.{sch}.{r}", call_exit=led.get(f"call.{idx[dn]}.{sch}.{r}", {}).get("exit", "not run"),
                   call_seconds=led.get(f"call.{idx[dn]}.{sch}.{r}", {}).get("wall_seconds", ""), call_peak_rss_gib=led.get(f"call.{idx[dn]}.{sch}.{r}", {}).get("peak_rss_gib", ""),
                   eval_exit=led.get(f"eval.{idx[dn]}.{sch}.{r}", {}).get("exit", "not run"))
        if lg:
            c = lg[0]
            mem = []
            for k in ("hap1_sha256", "hap2_sha256"):
                p = shaset.get(c[k])
                mem.append("off-panel (mosaic/reconstruction)" if p is None else (f"TRUTH template {p.split('#')[0]}#{p.split('#')[1]}" if p.startswith(dn + "#")
                           else ("panel member " + p.split("#")[0] + "#" + p.split("#")[1] + ("" if p not in excluded else " (EXCLUDED from calling panel?!)"))))
            row.update(run_status=c["run_status"], call_class=c["call_class"], phase_support=c["phase_support"], hap1=mem[0], hap2=mem[1],
                       total_score=c["total_score"], fragment_term=c["fragment_term"], exposure_term=c["exposure_term"], prior=c["sequence_pair_prior"],
                       switches=f"{c.get('hap1_switches')}+{c.get('hap2_switches')}")
        row.update({k: ev.get(k, "") for k in ("status", "availability_floor_total", "unordered_distance_total", "excess_over_floor_total", "phased_distance_total",
                                                "phase_switch_blocks", "blocks_exact", "blocks", "sv_blocks_exact", "undefined_truth_blocks", "truth_hap1_spellable",
                                                "truth_hap2_spellable", "truth_candidate_source", "floor_attained_in_candidates", "winner_excess_over_floor")})
        row.update({k: st.get(k, "") for k in ("sequence_candidates", "mosaic_candidates", "proposal_states_pruned", "seconds_total", "projection_ambiguous")})
        rows.append(row)
    write("cells", rows)
    paired = []
    for dn in donors:
        for sch in ("LZO", "LOO"):
            a = next(x for x in rows if x["cell"] == f"{dn}.{sch}.r0.01")
            b = next(x for x in rows if x["cell"] == f"{dn}.{sch}.r0.1")
            num = lambda x, k: float(x[k]) if x.get(k) not in ("", None, "NA") else None
            ea, eb = num(a, "excess_over_floor_total"), num(b, "excess_over_floor_total")
            paired.append(dict(donor=dn, scheme=sch, excess_r001=ea, excess_r01=eb, distance_r001=a.get("unordered_distance_total"), distance_r01=b.get("unordered_distance_total"),
                               floor=a.get("availability_floor_total"), same_call=a.get("hap1") == b.get("hap1") and a.get("hap2") == b.get("hap2") and a.get("total_score") != "" and
                               a.get("fragment_term") == b.get("fragment_term"),
                               verdict=("not assessable" if ea is None or eb is None else "both exact at floor (non-discriminating)" if ea == 0 and eb == 0
                                        else "r=0.1 better" if eb < ea else "r=0.1 worse" if eb > ea else "equal excess")))
    write("paired", paired)
    print(open(f"{out}/cells.tsv").read())
    print(open(f"{out}/paired.tsv").read())
