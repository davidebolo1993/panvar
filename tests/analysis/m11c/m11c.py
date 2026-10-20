#!/usr/bin/env python3
"""M11-C helpers (frozen with DESIGN.md). Nothing here genotypes, fits, ranks or pools.

  m11c.py compare <off_dir> <on_dir> <out>   off/on scientific identity of every pre-existing output (block_diagnostics on
                                             adds only <prefix>.block_diagnostics.tsv and the declared settings row)
  m11c.py diag_check <out>                   the smoke runs' block_diagnostics tables: block association, values, failure
                                             handling, determinism
  m11c.py truth <out>                        whole-model truth block sequences and allele sequences for every evaluated
                                             current-schema development run (genotype_probe block-sequences, the
                                             evaluator's own model / curation / truth arguments)
  m11c.py qvfrac <out>                       offline availability-relative QV fraction per block, beside the existing
                                             excess and normalized excess

QVfrac (Davide; COSIGT interpretation), offline evaluation only:
  QV(x, t) = -10 log10(max(0.5, d) / L), d = exact global edit distance, L = alignment columns (panvar::nw_edit_distance;
             the benchmark's qv_from; L = 0 gives 60)
  predicted: the reported block alleles (p1, p2) against the truth block haplotypes (t1, t2), paired by the existing minimum
             summed-edit criterion (the evaluator's unordered distance); a tie selects (p1-t1, p2-t2)
  best:      per truth haplotype, the evaluator's recorded closest allele over the active model's block alleles (its
             allowed allele-pair space; LZO: the truth itself; tie rule: increasing length gap, then allele order, strict
             improvement only); floor = closest1 + closest2 (verified)
  QVfrac = (QV(p, t1) + QV(p', t2)) / (QV(b1, t1) + QV(b2, t2))   -- near 1: near the best available
  normalized_excess (existing, scripts/genotype_report/load.py) -- near 0: near the best available
  Not the QV / perfect-match-ceiling ratio. The pairing that maximizes summed QV is computed and any disagreement with the
  minimum-edit pairing is reported, never substituted.
"""
import csv
import hashlib
import importlib.util
import math
import os
import re
import shlex
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
WT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
R = "/Users/davide.bolognini/panvar/tests/results"
OUT = f"{R}/m11c_confidence_inputs_20261020"
M11B = f"{R}/m11a_confidence_design_20261020/m11b/runs/manifest.a1/manifest.tsv"
ALIGN = f"{WT}/build/m10c_align_eval"


def T(p):
    return list(csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t")) if os.path.exists(p) else []


def KV(p):
    return dict(l.rstrip("\n").split("\t", 1) for l in open(p) if "\t" in l and not l.startswith("#")) if os.path.exists(p) else {}


def write(out, name, rows):
    with open(f"{out}/{name}.tsv", "w") as f:
        if not rows:
            f.write("none\n")
            return
        w = csv.DictWriter(f, fieldnames=list(dict.fromkeys(k for r in rows for k in r)), delimiter="\t", lineterminator="\n", restval="")
        w.writeheader()
        w.writerows(rows)


def qv(d, L):
    return 60.0 if L == 0 else -10.0 * math.log10(max(0.5, float(d)) / float(L))


# ------------------------------------------------------------------------------------------------ off/on identity
def compare(off, on, out):
    spec = importlib.util.spec_from_file_location("m10a_compare", os.path.join(HERE, "..", "m10a", "m10a_compare.py"))
    mc = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mc)
    m = mc.load_v32("operational,profile,checkpoints")
    m.SCIENTIFIC_SKIP = re.compile(r"^(?:" + m.SCIENTIFIC_SKIP.pattern + r"|call\.block_diagnostics\.tsv)$")
    argv = sys.argv
    sys.argv = ["v32", off, on, f"{out}/compare.tsv", "--categories", "operational,profile,checkpoints", "--declared-settings", "block_diagnostics"]
    try:
        code = m.main()
    finally:
        sys.argv = argv
    new_only = os.path.exists(f"{on}/call.block_diagnostics.tsv") and not os.path.exists(f"{off}/call.block_diagnostics.tsv")
    with open(f"{out}/compare.tsv", "a") as f:
        f.write(f"block_diagnostics_only_in_on\t{new_only}\nverdict\t{'PASS' if code == 0 and new_only else 'FAIL'}\n")
    print(open(f"{out}/compare.tsv").read())
    return 0 if code == 0 and new_only else 1


# ------------------------------------------------------------------------------------------------ diagnostics checks
def diag_check(out):
    runs = f"{OUT}/runs"
    rows, ok = [], True

    def table(d):
        return T(f"{runs}/{d}/call.block_diagnostics.tsv")

    def check(name, cond, detail=""):
        nonlocal ok
        ok = ok and bool(cond)
        rows.append(dict(check=name, ok=bool(cond), detail=detail))
        print(("PASS " if cond else "FAIL ") + name + (f" [{detail}]" if detail else ""))

    for d in ("smoke.c4.on", "smoke.acot.on", "smoke.cyp.on", "smoke.c4.dosageoff.on"):
        t = table(d)
        seq = [r["block_id"] for r in T(f"{runs}/{d}/call.sequence_blocks.tsv")]
        check(f"{d}: one row per model block, in model order, indices 0..n-1",
              t and [r["block_id"] for r in t] == seq and [int(r["block_index"]) for r in t] == list(range(len(t))), f"{len(t)} rows")
        legal_route = {"repeated_unit", "pair", "not_sensitive"}
        check(f"{d}: every census value is a route or an explicit NOT_EVALUATED / UNAVAILABLE",
              all(r["census_route"] in legal_route or r["census_route"].startswith(("NOT_EVALUATED:", "UNAVAILABLE:")) for r in t))
        check(f"{d}: every JP value is TIED_ELIGIBLE / NOT_ELIGIBLE or an explicit NOT_EVALUATED",
              all(r["jp_emission_tie"] in ("TIED_ELIGIBLE", "NOT_ELIGIBLE") or r["jp_emission_tie"].startswith("NOT_EVALUATED:") for r in t))
    c4, ac, cy, off = (table(d) for d in ("smoke.c4.on", "smoke.acot.on", "smoke.cyp.on", "smoke.c4.dosageoff.on"))
    route = lambda t, b: next((r["census_route"] for r in t if r["block_id"] == b), None)  # noqa: E731
    check("C4: sv:7 is repeated_unit (as in the saved dosage tables)", route(c4, "sv:7") == "repeated_unit")
    for d, t in (("smoke.c4.on", c4), ("smoke.acot.on", ac)):
        st = KV(f"{runs}/{d}/call.mosaic_status.tsv")
        tied = [r["block_id"] for r in t if r["jp_emission_tie"] == "TIED_ELIGIBLE"]
        check(f"{d}: TIED_ELIGIBLE count equals the status junction_eligible_blocks", str(len(tied)) == st.get("junction_eligible_blocks"),
              f"{tied} vs {st.get('junction_eligible_blocks')}")
    check("ACOT: the eligible block is sv:5 (M9-E)", [r["block_id"] for r in ac if r["jp_emission_tie"] == "TIED_ELIGIBLE"] == ["sv:5"])
    cst = KV(f"{runs}/smoke.cyp.on/call.dosage_status.tsv")
    check("CYP2D6 NA19650 5x: the dosage audit FAILED after the census, and the census routes are retained",
          cst.get("dosage_audit_status") == "FAILED" and all(r["census_route"] in ("repeated_unit", "pair", "not_sensitive") for r in cy),
          cst.get("dosage_audit_note", "")[:90])
    check("CYP2D6 NA19650 5x: sv:9 is repeated_unit", route(cy, "sv:9") == "repeated_unit")
    check("CYP2D6 (JP off): every JP value NOT_EVALUATED:junction_proposals_off",
          all(r["jp_emission_tie"] == "NOT_EVALUATED:junction_proposals_off" for r in cy))
    check("dosage mode off: every census value NOT_EVALUATED:dosage_audit_off (no route inferred)",
          all(r["census_route"] == "NOT_EVALUATED:dosage_audit_off" for r in off))
    a, b = (open(f"{runs}/{d}/call.block_diagnostics.tsv").read() for d in ("smoke.c4.on", "smoke.c4.on.r"))
    check("determinism: the C4 table is byte-identical on a repeat run", a == b)
    write(out, "diag_checks", rows)
    return 0 if ok else 1


# ------------------------------------------------------------------------------------------------ truth and QVfrac
def eval_args(eval_dir):
    first = open(f"{eval_dir}/run.log").readline()
    toks = shlex.split(first[first.index("genotype_probe"):])
    a = {toks[i][2:]: toks[i + 1] for i in range(len(toks) - 1) if toks[i].startswith("--")}
    return a


def v2_runs():
    rows = [r for r in T(M11B) if r["status"] == "INCLUDED" and r["stratum"].startswith("v2:")]
    for r in rows:
        r["eval_dir"] = os.path.dirname(r["summary"])
        r["args"] = eval_args(r["eval_dir"])
        cur = r["args"].get("panel-curation", "")
        r["key"] = hashlib.sha256("|".join([r["args"]["model-gfa"], r["args"]["reference"], r["args"]["truth-gfa"], r["args"]["donor"],
                                            hashlib.sha256(open(cur, "rb").read()).hexdigest() if cur else "none"]).encode()).hexdigest()[:16]
    return rows


def truth(out):
    done, rows = set(), []
    for r in v2_runs():
        if r["key"] in done:
            continue
        done.add(r["key"])
        a = r["args"]
        blocks = [x["block_id"] for x in T(r["run"] + ".sequence_blocks.tsv")]
        d = f"{out}/{r['key']}"
        os.makedirs(d, exist_ok=True)
        cmd = [f"{WT}/build/genotype_probe", "block-sequences", "--out", f"{d}/bs", "--reference", a["reference"], "--model-gfa", a["model-gfa"],
               "--blocks", ",".join(blocks), "--truth-gfa", a["truth-gfa"], "--donor", a["donor"]] + \
              (["--panel-curation", a["panel-curation"]] if a.get("panel-curation") else [])
        p = subprocess.run(cmd, capture_output=True, text=True)
        fp_ok = f"model fingerprint {r['model_fingerprint']}" in p.stdout
        rows.append(dict(key=r["key"], donor=a["donor"], model_gfa=a["model-gfa"], curation=a.get("panel-curation", "none"),
                         model_fingerprint=r["model_fingerprint"], fingerprint_matches=fp_ok, exit=p.returncode))
        print(r["key"], a["donor"], p.returncode, fp_ok)
    write(out, "truth_units", rows)
    return 0 if all(x["exit"] == 0 and x["fingerprint_matches"] for x in rows) else 1


def fasta(p):
    out, name = {}, None
    for line in open(p):
        line = line.rstrip("\n")
        if line.startswith(">"):
            name = line[1:].split()[0]
            out[name] = []
        elif name:
            out[name].append(line)
    return {k: "".join(v) for k, v in out.items()}


def qvfrac(out):
    tdir = f"{OUT}/runs/truth"
    pairs, meta = [], []
    for r in v2_runs():
        d = f"{tdir}/{r['key']}"
        for e in T(f"{r['eval_dir']}/eval.blocks.tsv"):
            bid = e["block_id"]
            m = dict(stratum=r["stratum"], donor=r["donor"], run=r["run"], block_id=bid, kind=e["kind"])
            if e["truth1_bp"] == "NA" or e["truth2_bp"] == "NA":
                meta.append(dict(m, status="TRUTH_UNDEFINED"))
                continue
            al, tr = fasta(f"{d}/bs.{bid}.alleles.fa"), fasta(f"{d}/bs.{bid}.truth.fa")
            p1, p2 = al[f"allele_{e['called_allele1']}"], al[f"allele_{e['called_allele2']}"]
            b1, b2 = al[f"allele_{e['truth1_closest_allele']}"], al[f"allele_{e['truth2_closest_allele']}"]
            t1, t2 = tr["truth1"], tr["truth2"]
            pid = hashlib.sha256((r["run"] + bid).encode()).hexdigest()[:16]
            for tag, x, y in (("p1t1", p1, t1), ("p2t2", p2, t2), ("p1t2", p1, t2), ("p2t1", p2, t1), ("b1t1", b1, t1), ("b2t2", b2, t2)):
                pairs.append((f"{pid}.{tag}", x, y))
            meta.append(dict(m, status="DEFINED", pid=pid, distance=float(e["unordered_distance"]), floor=float(e["availability_floor"]),
                             excess=float(e["excess"]), c1=float(e["truth1_closest_distance"]), c2=float(e["truth2_closest_distance"]),
                             den=max(len(t1) + len(t2), len(p1) + len(p2))))
    os.makedirs(f"{out}/seq", exist_ok=True)
    man = f"{out}/pairs.manifest.tsv"
    with open(man, "w") as f:
        for pid, x, y in pairs:
            open(f"{out}/seq/{pid}.A.txt", "w").write(x)
            open(f"{out}/seq/{pid}.B.txt", "w").write(y)
            f.write(f"{pid}\t{out}/seq/{pid}.A.txt\t{out}/seq/{pid}.B.txt\t-\n")
    subprocess.run([ALIGN, "--distance", man, f"{out}/pairs.distance.tsv"], check=True)
    D = {r["id"]: (int(r["edits"]), int(r["aln_len"])) for r in T(f"{out}/pairs.distance.tsv")}
    rows, consistency = [], {"pred_distance": 0, "floor": 0, "n": 0, "qv_pairing_differs": 0}
    for m in meta:
        if m["status"] != "DEFINED":
            rows.append(dict((k, m[k]) for k in ("stratum", "donor", "block_id", "kind", "status")))
            continue
        g = {t: D[f"{m['pid']}.{t}"] for t in ("p1t1", "p2t2", "p1t2", "p2t1", "b1t1", "b2t2")}
        e1, e2 = g["p1t1"][0] + g["p2t2"][0], g["p1t2"][0] + g["p2t1"][0]
        sel = ("p1t1", "p2t2") if e1 <= e2 else ("p1t2", "p2t1")
        other = ("p1t2", "p2t1") if sel == ("p1t1", "p2t2") else ("p1t1", "p2t2")
        qvp = qv(*g[sel[0]]) + qv(*g[sel[1]])
        qvo = qv(*g[other[0]]) + qv(*g[other[1]])
        qvb = qv(*g["b1t1"]) + qv(*g["b2t2"])
        consistency["n"] += 1
        consistency["pred_distance"] += min(e1, e2) == m["distance"]
        consistency["floor"] += (g["b1t1"][0] == m["c1"] and g["b2t2"][0] == m["c2"])
        consistency["qv_pairing_differs"] += qvo > qvp
        rows.append(dict(stratum=m["stratum"], donor=m["donor"], block_id=m["block_id"], kind=m["kind"], status="DEFINED",
                         excess=m["excess"], normalized_excess=m["excess"] / m["den"] if m["den"] else None, floor=m["floor"],
                         pred_edits=min(e1, e2), pred_pairing="p1-t1,p2-t2" if sel[0] == "p1t1" else "p1-t2,p2-t1",
                         qv_pred_1=qv(*g[sel[0]]), qv_pred_2=qv(*g[sel[1]]), L_pred=f"{g[sel[0]][1]}+{g[sel[1]][1]}",
                         qv_best_1=qv(*g["b1t1"]), qv_best_2=qv(*g["b2t2"]), L_best=f"{g['b1t1'][1]}+{g['b2t2'][1]}",
                         QVfrac=qvp / qvb if qvb else None, qv_max_pairing_differs=qvo > qvp, run=m["run"][-60:]))
    write(out, "qvfrac_blocks", rows)
    c = consistency
    summary = [dict(blocks_defined=c["n"], predicted_distance_matches_evaluator=c["pred_distance"], best_pair_edits_match_evaluator=c["floor"],
                    qv_maximizing_pairing_differs_from_min_edit_pairing=c["qv_pairing_differs"],
                    truth_undefined=sum(1 for m in meta if m["status"] != "DEFINED"))]
    write(out, "qvfrac_summary", summary)
    print(summary)
    return 0 if c["pred_distance"] == c["n"] and c["floor"] == c["n"] else 1


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "compare":
        os.makedirs(sys.argv[4], exist_ok=True)
        sys.exit(compare(sys.argv[2], sys.argv[3], sys.argv[4]))
    out = sys.argv[2]
    os.makedirs(out, exist_ok=True)
    sys.exit({"diag_check": diag_check, "truth": truth, "qvfrac": qvfrac}[mode](out))
