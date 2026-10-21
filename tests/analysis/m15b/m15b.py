#!/usr/bin/env python3
"""M15-B: final preflight and reserved ACOT evaluation (descriptive; the release candidate as a named experimental profile).

Modes:
  prepare <out>   frozen manifest of the 18 reserved families (donor, graph paths, LOO-curation exclusion check); simulator
                  lineage check (regenerate M9-H's HG02723 30x/5x and match its recorded SHA-256); reads for the 18 donors
                  (M8-C recipe: simulate 30x seed 20260916 salt donor; subsample 5x name-salt donor); every input hash
  freeze <out>    SHA-256 of every prediction output of the 18 calls, written before any truth evaluation
  truth <donor> <out>  block-sequences projection of the donor's truth (for the representation interval only)
  analyse <out>   the frozen analysis (DESIGN.md §4)
"""
import collections
import csv
import hashlib
import json
import os
import random
import statistics
import subprocess
import sys

csv.field_size_limit(sys.maxsize)
R = "/Users/davide.bolognini/panvar/tests/results"
WT = "/Users/davide.bolognini/panvar-worktrees/s1-d-panel-references"
OUT = f"{R}/m15b_reserved_acot_20261021"
PROBE, PANVAR = f"{WT}/build/genotype_probe", f"{WT}/build/panvar"
ACOT_GFA = f"{R}/genotype_perf_bounded_20261015/inputs/acot/acot.gfa"
ACOT_REF = "GRCh38#0#chr14:73416151-73668669"
CURATION = f"{R}/m9h_jp11_validation_20261020/runs/prepare/curation.LOO.txt"
FAMILY_TABLE = f"{R}/m11a_confidence_design_20261020/family_exposure.tsv"
M9H_PREP = f"{R}/m9h_jp11_validation_20261020/runs/prepare/prepare.tsv"
INTERVAL = ("sv:5", "ctx:6")       # ACOT degeneracy interval of this model (M9-I inventory; same fingerprint)
LARGE = 1000.0
BOOT, SEED = 2000, 11


def T(p):
    return list(csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t")) if os.path.exists(p) else []


def KV(p):
    return dict(l.rstrip("\n").split("\t", 1) for l in open(p) if "\t" in l and not l.startswith("#")) if os.path.exists(p) else {}


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def families():
    rows = [x for x in T(FAMILY_TABLE) if x["locus"] == "ACOT" and x["proposed_role"].startswith("VALIDATION")]
    return sorted(((x["family"], x["donors"]) for x in rows), key=lambda t: t[1])


def reads_dir(donor):
    return f"{OUT}/runs/prepare/reads/{donor}/cov5"


def run(cmd):
    p = subprocess.run(cmd, capture_output=True, text=True)
    print(" ".join(cmd[:4]), "->", p.returncode)
    if p.returncode != 0:
        print(p.stdout[-2000:], p.stderr[-2000:])
    return p.returncode


def simulate(donor, base):
    s30, s5 = f"{base}/{donor}/cov30", f"{base}/{donor}/cov5"
    os.makedirs(s30)
    os.makedirs(s5)
    assert run([PROBE, "simulate", "--gfa", ACOT_GFA, "--reference", ACOT_REF, "--donor", donor, "--depth", "30", "--read-length", "100",
                "--seed", "20260916", "--seed-salt", donor, "--out", s30]) == 0
    assert run([PROBE, "subsample", "--reads1", f"{s30}/reads_1.fq", "--reads2", f"{s30}/reads_2.fq", "--library", f"{s30}/library.txt",
                "--source-coverage", "30", "--coverage", "5", "--name-salt", donor, "--out", s5]) == 0
    return s5


def prepare(out):
    lines = [f"panvar_sha256\t{sha(PANVAR)}", f"genotype_probe_sha256\t{sha(PROBE)}", f"graph\t{ACOT_GFA}\t{sha(ACOT_GFA)}",
             f"curation\t{CURATION}\t{sha(CURATION)}"]
    cur = {l.strip() for l in open(CURATION) if l.strip() and not l.startswith("#")}
    paths = [l.split("\t")[1] for l in open(ACOT_GFA) if l.startswith("P\t")]
    fams = families()
    assert len(fams) == 18, fams
    for f, d in fams:
        ps = [p for p in paths if p.startswith(d + "#")]
        assert len(ps) == 2 and all(p in cur for p in ps), (f, d, ps)   # exactly two truth paths, both outside the active panel
        lines.append(f"family\t{f}\t{d}\t{ps[0]}\t{ps[1]}\texcluded_from_active_panel")
    # simulator lineage: regenerating M9-H's HG02723 reads must reproduce the recorded hashes byte for byte
    rec = next(l for l in open(M9H_PREP) if l.startswith("reads_HG02723\t")).split()
    want = dict(zip(rec[1::2], rec[2::2]))
    s5 = simulate("HG02723", f"{out}/lineage")
    got = {"reads_1": sha(f"{s5}/reads_1.fq"), "reads_2": sha(f"{s5}/reads_2.fq"), "library": sha(f"{s5}/library.txt")}
    ok = all(want[k] == got[k] for k in want)
    lines.append(f"lineage_HG02723\t{'REPRODUCED' if ok else 'DIFFERS'}\t{json.dumps(got)}")
    if not ok:
        open(f"{out}/prepare.tsv", "w").write("\n".join(lines) + "\n")
        print("simulator lineage check FAILED")
        return 1
    for f, d in fams:
        s5 = simulate(d, f"{out}/reads")
        lines.append(f"reads\t{d}\treads_1 {sha(s5 + '/reads_1.fq')}\treads_2 {sha(s5 + '/reads_2.fq')}\tlibrary {sha(s5 + '/library.txt')}")
    open(f"{out}/prepare.tsv", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


def freeze(out):
    rows = []
    for f, d in families():
        base = f"{OUT}/runs/call.{d}"
        if not os.path.isdir(base):
            rows.append(f"{d}\t-\tNO_CALL_DIRECTORY\n")
            continue
        for name in sorted(os.listdir(base)):
            p = f"{base}/{name}"
            if os.path.isfile(p) and name != "run.log":
                rows.append(f"{d}\t{name}\t{sha(p)}\n")
    open(f"{out}/prediction_outputs.sha256.tsv", "w").writelines(["donor\tfile\tsha256\n"] + rows)
    print(f"{len(rows)} prediction outputs hashed before truth evaluation")
    return 0


def truth(donor, out):
    blocks = [r["block_id"] for r in T(f"{OUT}/runs/call.{donor}/call.sequence_blocks.tsv")]
    if not blocks:
        print("no sequence_blocks: nothing to project")
        return 0
    cmd = [PROBE, "block-sequences", "--out", f"{out}/bs", "--reference", ACOT_REF, "--model-gfa", ACOT_GFA, "--blocks", ",".join(blocks),
           "--truth-gfa", ACOT_GFA, "--donor", donor, "--panel-curation", CURATION]
    p = subprocess.run(cmd, capture_output=True, text=True)
    print(p.stdout + p.stderr)
    fp = KV(f"{OUT}/runs/call.{donor}/call.mosaic_status.tsv").get("model_fingerprint")
    return 0 if p.returncode == 0 and f"model fingerprint {fp}" in p.stdout else 1


def fasta(p):
    out, name = {}, None
    for l in open(p):
        l = l.rstrip("\n")
        if l.startswith(">"):
            name = l[1:].split()[0]
            out[name] = []
        elif name:
            out[name].append(l)
    return {k: "".join(v) for k, v in out.items()}


def edit_distances(pairs, out, tag):
    os.makedirs(f"{out}/seq", exist_ok=True)
    man = f"{out}/{tag}.manifest.tsv"
    with open(man, "w") as f:
        for pid, a, b in pairs:
            open(f"{out}/seq/{pid}.A.txt", "w").write(a)
            open(f"{out}/seq/{pid}.B.txt", "w").write(b)
            f.write(f"{pid}\t{out}/seq/{pid}.A.txt\t{out}/seq/{pid}.B.txt\t-\n")
    subprocess.run([f"{WT}/build/m10c_align_eval", "--distance", man, f"{out}/{tag}.distance.tsv"], check=True, capture_output=True)
    return {r["id"]: int(r["edits"]) for r in T(f"{out}/{tag}.distance.tsv")}


def repeated_unit(route):
    if route == "repeated_unit":
        return "YES"
    if route in ("not_sensitive", "pair"):
        return "NO"
    return "UNKNOWN"


def ledger():
    p = f"{OUT}/ledger.tsv"
    if not os.path.exists(p):
        return {}
    lines = [l.rstrip("\n").split("\t") for l in open(p)]
    return {r[0]: dict(zip(lines[0], r)) for r in lines[1:]}


def analyse(out):
    led = ledger()
    runs, blocks = [], []
    for fam, d in families():
        call = f"{OUT}/runs/call.{d}/call"
        st, sl, lg = KV(call + ".mosaic_status.tsv"), KV(call + ".sequence_locus.tsv"), (T(call + ".locus_genotype.tsv") or [{}])[0]
        ds = KV(call + ".dosage_status.tsv")
        es = KV(f"{OUT}/runs/eval.{d}/eval.summary.tsv")
        cl = led.get(f"call.{d}", {})
        valid = cl.get("exit") == "0" and sl.get("prediction") == "PREDICTION" and sl.get("completion") == "completed_ranking" \
            and sl.get("provisional", "no") == "no" and st.get("representation_guard_complete", "true") == "true"
        runs.append(dict(family=fam, donor=d, call_exit=cl.get("exit", "NA"), call_watchdog=cl.get("watchdog", "NA"), wall_s=cl.get("wall_seconds", "NA"),
                         peak_gib=cl.get("peak_gib", "NA"), status=st.get("status", "NA"), prediction=sl.get("prediction", "NA"),
                         sequence_status=sl.get("sequence_status", "NA"), completion=sl.get("completion", "NA"), run_status=lg.get("run_status", "NA"),
                         call_class=lg.get("call_class", "NA"), no_call_reason=lg.get("no_call_reason", "NA"), completed_prediction=valid,
                         dosage_audit=ds.get("dosage_audit_status", "NA"), dosage_note=ds.get("dosage_audit_note", "NA"),
                         eval_status=es.get("status", "NOT_EVALUATED"), seconds_total=st.get("seconds_total", "NA"),
                         sequence_ambiguous=sl.get("sequence_status", "NA") not in ("RESOLVED", "NA")))
        ev = T(f"{OUT}/runs/eval.{d}/eval.blocks.tsv")
        bg = {r["block_id"]: r for r in T(call + ".block_genotypes.tsv")}
        sb = {r["block_id"]: r for r in T(call + ".sequence_blocks.tsv")}
        dg = {r["block_id"]: r for r in T(call + ".block_diagnostics.tsv")}
        db = {r["block_id"]: r for r in T(call + ".dosage_blocks.tsv")}
        # representation interval sv:5..ctx:6: summed block distances against the concatenated interval distance
        rep = {}
        order = [e["block_id"] for e in ev]
        if ev and all(b in order for b in INTERVAL):
            blk = order[order.index(INTERVAL[0]):order.index(INTERVAL[1]) + 1]
            es_ = {e["block_id"]: e for e in ev}
            ok = all(es_[b]["excess"] not in ("NA", "") and es_[b]["truth1_bp"] != "NA" for b in blk)
            tdir = f"{OUT}/runs/truth.{d}"
            if ok and all(os.path.exists(f"{tdir}/bs.{b}.truth.fa") for b in blk):
                cat = {"p1": "", "p2": "", "t1": "", "t2": ""}
                for b in blk:
                    al, tr = fasta(f"{tdir}/bs.{b}.alleles.fa"), fasta(f"{tdir}/bs.{b}.truth.fa")
                    cat["p1"] += al[f"allele_{es_[b]['called_allele1']}"]
                    cat["p2"] += al[f"allele_{es_[b]['called_allele2']}"]
                    cat["t1"] += tr["truth1"]
                    cat["t2"] += tr["truth2"]
                D = edit_distances([(f"{d}.{t}", cat[x], cat[y]) for t, x, y in (("p1t1", "p1", "t1"), ("p2t2", "p2", "t2"), ("p1t2", "p1", "t2"),
                                                                                  ("p2t1", "p2", "t1"))], out, d)
                iv_d = min(D[f"{d}.p1t1"] + D[f"{d}.p2t2"], D[f"{d}.p1t2"] + D[f"{d}.p2t1"])
                blk_d = sum(float(es_[b]["unordered_distance"]) for b in blk)
                if blk_d - iv_d > 0:
                    for b in blk:
                        rep[b] = ("REPRESENTATION_INTERVAL_EXACT" if iv_d == 0 else "REPRESENTATION_UNRESOLVED", blk_d - iv_d, iv_d)
        for e in ev:
            b = e["block_id"]
            if e["truth1_bp"] in ("NA", "") or e["truth2_bp"] in ("NA", ""):
                cat = "UNDEFINED_TRUTH"
            elif e["excess"] in ("NA", ""):
                cat = "NO_PREDICTION"
            elif all(e[k] == "0" for k in ("truth1_bp", "truth2_bp", "called1_bp", "called2_bp")):
                cat = "ALL_EMPTY"
            elif b in rep:
                cat = rep[b][0]
            else:
                cat = "BLOCK_OUTCOME"
            ex = float(e["excess"]) if e["excess"] not in ("NA", "") else None
            den = max(float(e["truth1_bp"]) + float(e["truth2_bp"]), float(e["called1_bp"]) + float(e["called2_bp"])) \
                if cat not in ("UNDEFINED_TRUTH", "NO_PREDICTION") else 0
            dstat = db.get(b, {}).get("dosage_status", "")
            audit = KV(call + ".dosage_status.tsv").get("dosage_audit_status", "NA")
            blocks.append(dict(family=fam, donor=d, block_id=b, kind=e["kind"], category=cat, distance=e["unordered_distance"],
                               floor=e["availability_floor"], excess=ex, normalized_excess=(ex / den) if (ex is not None and den) else None,
                               catastrophic=bool(ex is not None and ex >= LARGE), floor_attained=bool(ex == 0),
                               truth_bp=f"{e['truth1_bp']}+{e['truth2_bp']}", called_bp=f"{e['called1_bp']}+{e['called2_bp']}",
                               representation_artefact_bp=rep.get(b, ("", "", ""))[1], interval_concatenated_distance=rep.get(b, ("", "", ""))[2],
                               genotype_support=bg.get(b, {}).get("genotype_support", "NA"), block_status=sb.get(b, {}).get("block_status", "NA"),
                               genotype_filter=bg.get(b, {}).get("genotype_filter", "NA"), dosage_filter=bg.get(b, {}).get("dosage_filter", "NA"),
                               error_mass_upper=sb.get(b, {}).get("error_mass_upper", "NA"),
                               legacy_GQ=sb.get(b, {}).get("legacy_representative_GQ", "NA"),
                               census_route=dg.get(b, {}).get("census_route", "NA"),
                               REPEATED_UNIT_STRUCTURE=repeated_unit(dg.get(b, {}).get("census_route", "")),
                               DOSAGE_AUDIT_DISAGREES=dstat in ("CONTENT_REPAIRED", "TRIGGERED_UNREPAIRED", "SUSPECT"),
                               DOSAGE_AUDIT_UNAVAILABLE=(audit != "COMPLETE") or dstat in ("UNDERPOWERED", "NOT_EVALUATED"),
                               PRIOR_ONLY_SUPPORT=bg.get(b, {}).get("genotype_support", "") == "PRIOR_ONLY",
                               BLOCK_REPRESENTATION_AMBIGUOUS=sb.get(b, {}).get("block_status", "") == "REPRESENTATION_AMBIGUOUS",
                               BLOCK_GENOTYPE_UNRESOLVED=sb.get(b, {}).get("block_status", "") == "GENOTYPE_UNRESOLVED"))
    write(out, "runs", runs)
    write(out, "blocks", blocks)
    inf = [b for b in blocks if b["category"] == "BLOCK_OUTCOME"]
    cats = collections.Counter(b["category"] for b in blocks)
    # per family
    fam_rows = []
    for fam, d in families():
        S = [b for b in inf if b["donor"] == d]
        ex = [b["excess"] for b in S]
        fam_rows.append(dict(family=fam, donor=d, informative_blocks=len(S), excess_sum=sum(ex) if ex else "NA",
                             excess_median=statistics.median(ex) if ex else "NA", excess_max=max(ex) if ex else "NA",
                             floor_attained=sum(b["floor_attained"] for b in S), catastrophic=sum(b["catastrophic"] for b in S),
                             blocks_total=sum(1 for b in blocks if b["donor"] == d)))
    write(out, "families", fam_rows)
    # family bootstrap
    by_fam = collections.defaultdict(list)
    for b in inf:
        by_fam[b["donor"]].append(b)
    fams = sorted(by_fam)
    rng = random.Random(SEED)

    def stats(sample):
        bl = [b for f in sample for b in by_fam[f]]
        if not bl:
            return None
        return (statistics.mean(b["excess"] for b in bl), statistics.mean(b["normalized_excess"] or 0 for b in bl),
                sum(b["catastrophic"] for b in bl) / len(bl), sum(b["floor_attained"] for b in bl) / len(bl),
                sum(any(b["catastrophic"] for b in by_fam[f]) for f in sample) / len(sample))
    obs = stats(fams)
    boots = [stats([rng.choice(fams) for _ in fams]) for _ in range(BOOT)] if fams else []
    names = ("mean_excess_per_block", "mean_normalized_excess", "catastrophic_block_rate", "floor_attainment_rate", "families_with_catastrophe")
    summ = []
    for i, n in enumerate(names):
        v = sorted(x[i] for x in boots if x)
        summ.append(dict(measure=n, observed=obs[i] if obs else "NA", ci95_low=v[int(0.025 * len(v))] if v else "NA",
                         ci95_high=v[int(0.975 * len(v)) - 1] if v else "NA"))
    ex_all = sorted(b["excess"] for b in inf)
    nx_all = sorted(b["normalized_excess"] for b in inf if b["normalized_excess"] is not None)
    q = lambda v, p: v[min(len(v) - 1, int(p * (len(v) - 1)))] if v else "NA"  # noqa: E731
    summ += [dict(measure=f"excess_q{int(p * 100)}", observed=q(ex_all, p), ci95_low="", ci95_high="") for p in (0.5, 0.9, 0.99, 1.0)]
    summ += [dict(measure=f"normalized_excess_q{int(p * 100)}", observed=q(nx_all, p), ci95_low="", ci95_high="") for p in (0.5, 0.9, 0.99, 1.0)]
    summ += [dict(measure=f"blocks_{k}", observed=v, ci95_low="", ci95_high="") for k, v in sorted(cats.items())]
    summ += [dict(measure="informative_families", observed=len(fams), ci95_low="", ci95_high=""),
             dict(measure="runs_completed_prediction", observed=sum(r["completed_prediction"] for r in runs), ci95_low="", ci95_high=""),
             dict(measure="runs_sequence_ambiguous", observed=sum(r["sequence_ambiguous"] for r in runs), ci95_low="", ci95_high="")]
    write(out, "summary", summ)
    write(out, "catastrophic", [b for b in blocks if b["catastrophic"]])
    # descriptive annotation associations over informative blocks
    assoc = []
    for ann in ("REPEATED_UNIT_STRUCTURE", "DOSAGE_AUDIT_DISAGREES", "DOSAGE_AUDIT_UNAVAILABLE", "PRIOR_ONLY_SUPPORT", "BLOCK_REPRESENTATION_AMBIGUOUS",
                "BLOCK_GENOTYPE_UNRESOLVED", "genotype_filter", "dosage_filter"):
        for val, n in sorted(collections.Counter(str(b[ann]) for b in inf).items()):
            S = [b for b in inf if str(b[ann]) == val]
            assoc.append(dict(annotation=ann, value=val, blocks=n, families=len({b["donor"] for b in S}),
                              catastrophic=sum(b["catastrophic"] for b in S), floor_attained=sum(b["floor_attained"] for b in S),
                              mean_excess=statistics.mean(b["excess"] for b in S)))
    write(out, "annotations", assoc)
    for s in summ:
        print(s)
    return 0


def write(out, name, rows):
    with open(f"{out}/{name}.tsv", "w") as f:
        if not rows:
            f.write("none\n")
            return
        w = csv.DictWriter(f, fieldnames=list(dict.fromkeys(k for r in rows for k in r)), delimiter="\t", lineterminator="\n", restval="")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "truth":
        os.makedirs(sys.argv[3], exist_ok=True)
        sys.exit(truth(sys.argv[2], sys.argv[3]))
    out = sys.argv[2]
    os.makedirs(out, exist_ok=True)
    sys.exit({"prepare": prepare, "freeze": freeze, "analyse": analyse}[mode](out))
