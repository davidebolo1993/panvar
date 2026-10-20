#!/usr/bin/env python3
"""M11-D: bounded unified-configuration development pilot (frozen with DESIGN.md). Development estimates only.

  m11d.py selftest <out>     fitting / grouping gate on synthetic data (no genotyping)
  m11d.py inventory <out>    degeneracy inventory of every pilot LOO model (existing tools; model-only)
  m11d.py truth <cell> <out> whole-model truth and allele block sequences (genotype_probe block-sequences, the evaluation's args)
  m11d.py outcomes <out>     per-block outcomes: excess, normalized excess, floor attainment, catastrophic, QVfrac (M11-C
                             definition), synonymous-boundary artefact flags; per-block inputs (error mass / GQ, read support,
                             block diagnostics)
  m11d.py fit <out>          the frozen family-grouped development comparison of B0 / C1 / C2
"""
import collections
import csv
import hashlib
import math
import os
import random
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "m11c"))
import m11c  # noqa: E402  (qv, fasta, T, KV, write)

WT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
R = "/Users/davide.bolognini/panvar/tests/results"
OUT = f"{R}/m11d_dev_pilot_20261020"
T, KV, write, qv, fasta = m11c.T, m11c.KV, m11c.write, m11c.qv, m11c.fasta
TOP = 0.8                 # the historical COSIGT top availability-relative category (an evaluation label, not a threshold)
CATASTROPHIC_BP = 1000.0
RIDGE = 1e-4              # frozen: L2 on standardized slopes only, for numerical stability under separation
BOOT, SEED = 2000, 11

C4R, CYPR, ACOTR = "HG00096#1#haplotype1-0000024:31848049-32074186", "grch38#1#chr22:42031864-42245566", "GRCh38#0#chr14:73416151-73668669"
C4D, CYPD, ACOTG = f"{R}/genotype_c4_production_11fe681", f"{R}/cyp2d6_shortread_pilot_v1/data", f"{R}/genotype_perf_bounded_20261015/inputs/acot/acot.gfa"
CYPCUR = "/Users/davide.bolognini/panvar/tests/genotype_cyp2d6_curation.tsv"
ACOTCUR = f"{R}/m9h_jp11_validation_20261020/runs/prepare/curation.LOO.txt"   # every ACOT outer-evaluation template excluded
ACOT_READS = {"HG03209": f"{R}/m8c_copying_challenge_20261020/runs/prepare.a1/reads/HG03209/cov5",
              "HG01975": f"{R}/m8c_copying_challenge_20261020/runs/prepare.a1/reads/HG01975/cov5",
              **{d: f"{R}/m9c_acot_reliability_20261020/runs/prepare/reads/{d}/cov5" for d in ("NA19909", "NA18747", "HG02178", "HG00639", "NA21106", "HG03831")},
              **{d: f"{R}/m9h_jp11_validation_20261020/runs/prepare/reads/{d}/cov5" for d in ("HG02723", "HG02074", "NA18971", "HG02083")}}
# frozen cells (DESIGN.md §2): CYP2D6 = the 6 pilot donors; C4 = the first 12 eligible development families by
# sha256("M11-D:C4:" + family) (donor by sha256("M11-D:C4:donor:" + donor)); ACOT = the 12 previously genotyped families
FAMILY = {"HG02769": "GB57", "HG02953": "NG16", "HG03456": "SL51", "NA19650": "m001", "NA19705": "2368", "NA19983": "2436",
          "NA19317": "NA19317", "HG03065": "SL05", "HG02018": "VN047", "HG00513": "SH032", "NA12329": "1328", "HG00358": "HG00358",
          "NA19331": "LWK001", "HG01457": "CLM51", "HG03732": "IT003", "NA24385": "NA24385", "NA19434": "LWK003",
          "HG03209": "SL17", "HG01975": "PEL25", "NA19909": "2427", "NA18747": "NA18747", "HG02178": "HG02178", "HG00639": "PR03",
          "NA21106": "NA21106", "HG03831": "BD27", "HG02723": "GB52", "HG02074": "VN066", "NA18971": "NA18971", "HG02083": "VN069"}
CELLS = ([dict(id=f"cyp.{d}", locus="CYP2D6", donor=d) for d in ("HG02769", "HG02953", "HG03456", "NA19650", "NA19705", "NA19983")] +
         [dict(id=f"c4.{d}", locus="C4", donor=d) for d in ("NA19705", "NA19317", "HG03065", "HG02018", "HG00513", "NA12329", "HG00358",
                                                             "NA19331", "HG01457", "HG03732", "NA24385", "NA19434")] +
         [dict(id=f"acot.{d}", locus="ACOT", donor=d) for d in ACOT_READS])
WATCHDOG = {"C4": 150, "CYP2D6": 150, "ACOT": 300}   # time_limit_seconds = watchdog - 20 (declared operational guard)
ARTEFACT_INTERVAL = {"ACOT": ("sv:5", "ctx:6"), "CYP2D6": ("sv:3", "ctx:4", "sv:5")}


def cell(cid):
    return next(c for c in CELLS if c["id"] == cid)


def inputs(c):
    d = c["donor"]
    if c["locus"] == "C4":
        rd = f"{C4D}/reads/{d}/cov5"
        return dict(model=f"{C4D}/graphs/{d}.loo.gfa", ref=C4R, cur=None, truth=f"{C4D}/graphs/c4.gfa", reads=rd, base="c4_5")
    if c["locus"] == "CYP2D6":
        return dict(model=f"{CYPD}/graphs/loo.{d}.gfa", ref=CYPR, cur=CYPCUR, truth=f"{CYPD}/graphs/cyp2d6.gfa", reads=f"{CYPD}/reads/{d}/cov5",
                    base="cyp_5")
    return dict(model=ACOTG, ref=ACOTR, cur=ACOTCUR, truth=ACOTG, reads=ACOT_READS[d], base="acot_5")


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


# ------------------------------------------------------------------------------------------------ fitting
def fit_logistic(X, y):
    """Logistic regression by Newton iterations; slopes standardized and ridge-stabilized (RIDGE); intercept free."""
    X = np.asarray(X, float)
    y = np.asarray(y, float)
    mu = X.mean(0) if X.shape[1] else np.zeros(0)
    sd = X.std(0) if X.shape[1] else np.zeros(0)
    sd[sd == 0] = 1.0
    Z = np.column_stack([np.ones(len(y)), (X - mu) / sd]) if X.shape[1] else np.ones((len(y), 1))
    w = np.zeros(Z.shape[1])
    pen = np.full(Z.shape[1], RIDGE)
    pen[0] = 0.0
    for _ in range(100):
        p = 1.0 / (1.0 + np.exp(-np.clip(Z @ w, -35, 35)))
        g = Z.T @ (p - y) + pen * w
        H = (Z * (p * (1 - p))[:, None]).T @ Z + np.diag(pen) + 1e-12 * np.eye(len(w))
        step = np.linalg.solve(H, g)
        w -= step
        if np.max(np.abs(step)) < 1e-10:
            break
    return dict(w=w, mu=mu, sd=sd)


def predict(model, X):
    X = np.asarray(X, float)
    Z = np.column_stack([np.ones(len(X)), (X - model["mu"]) / model["sd"]]) if X.shape[1] else np.ones((len(X), 1))
    return 1.0 / (1.0 + np.exp(-np.clip(Z @ model["w"], -35, 35)))


def logloss(p, y):
    p = np.clip(np.asarray(p, float), 1e-12, 1 - 1e-12)
    y = np.asarray(y, float)
    return -(y * np.log(p) + (1 - y) * np.log(1 - p))


def auc(p, y):
    pos = [a for a, b in zip(p, y) if b == 1]
    neg = [a for a, b in zip(p, y) if b == 0]
    if not pos or not neg:
        return None
    s = sum((1.0 if a > b else 0.5 if a == b else 0.0) for a in pos for b in neg)
    return s / (len(pos) * len(neg))


def lofo(rows, features):
    """Leave-one-family-out predictions (families kept together across loci)."""
    fams = sorted({r["family"] for r in rows})
    out = {}
    for f in fams:
        train = [r for r in rows if r["family"] != f]
        test = [r for r in rows if r["family"] == f]
        m = fit_logistic([[r[k] for k in features] for r in train], [r["y"] for r in train])
        for r, p in zip(test, predict(m, [[r[k] for k in features] for r in test])):
            out[id(r)] = float(p)
    return out


def selftest(out):
    rng = np.random.default_rng(5)
    n = 4000
    x = rng.normal(size=n)
    z = rng.integers(0, 2, size=n)
    p = 1 / (1 + np.exp(-(0.5 + 1.5 * x - 1.0 * z)))
    y = (rng.random(n) < p).astype(float)
    m = fit_logistic(np.column_stack([x, z]), y)
    b = m["w"][1:] / m["sd"]
    ok1 = abs(b[0] - 1.5) < 0.15 and abs(b[1] + 1.0) < 0.15
    m0 = fit_logistic(np.zeros((n, 0)), y)
    ok2 = abs(predict(m0, np.zeros((3, 0)))[0] - y.mean()) < 1e-9
    rows = [dict(family=f"F{i % 7}", y=float(y[i]), x=float(x[i]), z=float(z[i])) for i in range(700)]
    pr = lofo(rows, ["x"])
    ok3 = len(pr) == len(rows)
    sep = fit_logistic(np.array([[0.0], [1.0], [2.0], [3.0]]), np.array([0, 0, 1, 1.0]))
    ok4 = np.all(np.isfinite(sep["w"]))
    res = [dict(check="recovers simulated slopes (1.5, -1.0) within 0.15", ok=ok1, detail=f"{b}"),
           dict(check="constant baseline equals the training prevalence", ok=ok2),
           dict(check="leave-one-family-out predicts every row once", ok=ok3),
           dict(check="perfect separation stays finite (ridge)", ok=bool(ok4))]
    write(out, "selftest", res)
    print(res)
    return 0 if all(r["ok"] for r in res) else 1


# ------------------------------------------------------------------------------------------------ inventory and truth
def inventory(out):
    specs, rows = [], []
    for c in CELLS:
        if c["locus"] != "C4":
            continue
        i = inputs(c)
        code = subprocess.run([f"{WT}/build/cyp2d6_block_alleles", i["model"], i["ref"], "none", f"{out}/{c['id']}"]).returncode
        if code:
            return code
        specs.append(f"{c['id']}={out}/{c['id']}")
    code = subprocess.run([sys.executable, f"{WT}/tests/analysis/cyp2d6_pilot/degeneracy_inventory.py", f"{out}/degeneracy.tsv",
                           f"{out}/degeneracy_summary.tsv", *specs]).returncode
    print(open(f"{out}/degeneracy_summary.tsv").read() if os.path.exists(f"{out}/degeneracy_summary.tsv") else "no summary")
    return code


def truth(cid, out):
    c = cell(cid)
    i = inputs(c)
    blocks = [r["block_id"] for r in T(f"{OUT}/runs/call.{cid}/call.sequence_blocks.tsv")]
    if not blocks:
        print("no sequence_blocks: nothing to project")
        return 1
    cmd = [f"{WT}/build/genotype_probe", "block-sequences", "--out", f"{out}/bs", "--reference", i["ref"], "--model-gfa", i["model"],
           "--blocks", ",".join(blocks), "--truth-gfa", i["truth"], "--donor", c["donor"]] + (["--panel-curation", i["cur"]] if i["cur"] else [])
    p = subprocess.run(cmd, capture_output=True, text=True)
    print(p.stdout + p.stderr)
    fp = KV(f"{OUT}/runs/call.{cid}/call.mosaic_status.tsv").get("model_fingerprint")
    return 0 if p.returncode == 0 and f"model fingerprint {fp}" in p.stdout else 1


# ------------------------------------------------------------------------------------------------ outcomes
def degenerate_intervals():
    out = {}
    s = f"{OUT}/runs/inventory/degeneracy_summary.tsv"
    for r in T(s):
        iv = r.get("coalescence_intervals", "none")
        out[r["model"]] = [] if iv in ("", "none", "NA") else [tuple(x.split("..")) for x in iv.split(",")]
    return out


def align(pairs, d, tag):
    os.makedirs(f"{d}/seq", exist_ok=True)
    man = f"{d}/{tag}.manifest.tsv"
    with open(man, "w") as f:
        for pid, a, b in pairs:
            open(f"{d}/seq/{pid}.A.txt", "w").write(a)
            open(f"{d}/seq/{pid}.B.txt", "w").write(b)
            f.write(f"{pid}\t{d}/seq/{pid}.A.txt\t{d}/seq/{pid}.B.txt\t-\n")
    subprocess.run([f"{WT}/build/m10c_align_eval", "--distance", man, f"{d}/{tag}.distance.tsv"], check=True)
    return {r["id"]: (int(r["edits"]), int(r["aln_len"])) for r in T(f"{d}/{tag}.distance.tsv")}


def gq_of(sb):
    """C1 input on a stable log scale: -10 log10(error_mass_upper); a printed 0 uses legacy_representative_GQ only where
    lower == upper (the correspondence verified in this run), else NA."""
    try:
        m = float(sb["error_mass_upper"])
    except (KeyError, ValueError):
        return None, "missing"
    if m > 0:
        return -10.0 * math.log10(m), "error_mass_upper"
    if sb.get("error_mass_lower") == sb.get("error_mass_upper"):
        try:
            return float(sb["legacy_representative_GQ"]), "legacy_GQ_for_printed_zero"
        except ValueError:
            return None, "printed_zero_no_GQ"
    return None, "printed_zero_bounds_differ"


def outcomes(out):
    degen = degenerate_intervals()
    rows, runs, corr = [], [], []
    for c in CELLS:
        cid = c["id"]
        call = f"{OUT}/runs/call.{cid}/call"
        ev = f"{OUT}/runs/eval.{cid}/eval"
        st, sl = KV(call + ".mosaic_status.tsv"), KV(call + ".sequence_locus.tsv")
        es = KV(ev + ".summary.tsv")
        assessable = sl.get("completion") == "completed_ranking" and sl.get("provisional", "no") == "no" and st.get("representation_guard_complete", "true") == "true"
        runs.append(dict(cell=cid, locus=c["locus"], donor=c["donor"], family=FAMILY[c["donor"]], prediction=sl.get("prediction", "MISSING"),
                         sequence_status=sl.get("sequence_status", "NA"), completion=sl.get("completion", "NA"), legacy_status=st.get("status", "NA"),
                         assessable=assessable, eval_status=es.get("status", "NOT_EVALUATED"), model_fingerprint=st.get("model_fingerprint", "NA"),
                         decision_set_probability=sl.get("decision_set_probability", "NA"), seconds_total=st.get("seconds_total", "NA")))
        if not os.path.exists(ev + ".blocks.tsv"):
            continue
        sb = {r["block_id"]: r for r in T(call + ".sequence_blocks.tsv")}
        bg = {r["block_id"]: r for r in T(call + ".block_genotypes.tsv")}
        dg = {r["block_id"]: r for r in T(call + ".block_diagnostics.tsv")}
        tdir = f"{OUT}/runs/truth.{cid}"
        evb = T(ev + ".blocks.tsv")
        pairs, meta = [], []
        for e in evb:
            bid = e["block_id"]
            if e["truth1_bp"] == "NA" or e["truth2_bp"] == "NA" or not os.path.exists(f"{tdir}/bs.{bid}.truth.fa"):
                meta.append((e, None))
                continue
            if "NA" in (e["called_allele1"], e["called_allele2"], e["excess"]):
                # amendment A1: a block without a called genotype (a caller no-call) is NO_PREDICTION, never an outcome
                meta.append((e, "NO_PREDICTION"))
                continue
            al, tr = fasta(f"{tdir}/bs.{bid}.alleles.fa"), fasta(f"{tdir}/bs.{bid}.truth.fa")
            seqs = dict(p1=al[f"allele_{e['called_allele1']}"], p2=al[f"allele_{e['called_allele2']}"], b1=al[f"allele_{e['truth1_closest_allele']}"],
                        b2=al[f"allele_{e['truth2_closest_allele']}"], t1=tr["truth1"], t2=tr["truth2"])
            pid = f"{cid}.{bid}".replace(":", "_")
            for tag, x, y in (("p1t1", "p1", "t1"), ("p2t2", "p2", "t2"), ("p1t2", "p1", "t2"), ("p2t1", "p2", "t1"), ("b1t1", "b1", "t1"), ("b2t2", "b2", "t2")):
                pairs.append((f"{pid}.{tag}", seqs[x], seqs[y]))
            meta.append((e, (pid, seqs)))
        # synonymous-boundary artefact intervals: from the model's degeneracy inventory (C4: this unit's inventory; ACOT / CYP2D6:
        # the M9-I inventory of the same LOO models, sv:5..ctx:6 / sv:3..sv:5)
        ivs = degen.get(cid) if c["locus"] == "C4" else [ARTEFACT_INTERVAL[c["locus"]][::max(1, len(ARTEFACT_INTERVAL[c["locus"]]) - 1)]]
        iv_pairs = []
        for k, iv in enumerate(ivs or []):
            order = [e["block_id"] for e in evb]
            if iv[0] not in order or iv[1] not in order:
                continue
            blk = order[order.index(iv[0]):order.index(iv[1]) + 1]
            parts = [m for e, m in meta if e["block_id"] in blk]
            if len(parts) != len(blk) or any(not isinstance(m, tuple) for m in parts):
                continue
            cat = {k2: "".join(m[1][k2] for m in parts) for k2 in ("p1", "p2", "t1", "t2")}
            iv_pairs.append((k, blk, cat))
            for tag, x, y in (("p1t1", "p1", "t1"), ("p2t2", "p2", "t2"), ("p1t2", "p1", "t2"), ("p2t1", "p2", "t1")):
                pairs.append((f"{cid}.iv{k}.{tag}", cat[x], cat[y]))
        D = align(pairs, out, cid) if pairs else {}
        artefact_blocks = {}
        for k, blk, cat in iv_pairs:
            g = {t: D[f"{cid}.iv{k}.{t}"][0] for t in ("p1t1", "p2t2", "p1t2", "p2t1")}
            iv_d = min(g["p1t1"] + g["p2t2"], g["p1t2"] + g["p2t1"])
            blk_d = sum(float(e["unordered_distance"]) for e in evb if e["block_id"] in blk)
            for b in blk:
                artefact_blocks[b] = max(0.0, blk_d - iv_d)
        for e, m in meta:
            bid = e["block_id"]
            base = dict(cell=cid, locus=c["locus"], donor=c["donor"], family=FAMILY[c["donor"]], block_id=bid, kind=e["kind"], assessable=assessable)
            gq, src = gq_of(sb.get(bid, {}))
            if sb.get(bid, {}).get("error_mass_lower") == sb.get(bid, {}).get("error_mass_upper") and gq is not None and src == "error_mass_upper":
                try:
                    corr.append(abs(gq - float(sb[bid]["legacy_representative_GQ"])))
                except (KeyError, ValueError):
                    pass
            inp = dict(gq=gq, gq_source=src, error_mass_upper=sb.get(bid, {}).get("error_mass_upper", "NA"),
                       read_support=bg.get(bid, {}).get("genotype_support", "NA"), census_route=dg.get(bid, {}).get("census_route", "NA"),
                       jp_emission_tie=dg.get(bid, {}).get("jp_emission_tie", "NA"), block_status=sb.get(bid, {}).get("block_status", "NA"))
            if m is None:
                rows.append(dict(base, truth_defined=False, predicted=False, **inp))
                continue
            if m == "NO_PREDICTION":
                rows.append(dict(base, truth_defined=True, predicted=False, **inp))
                continue
            pid, seqs = m
            g = {t: D[f"{pid}.{t}"] for t in ("p1t1", "p2t2", "p1t2", "p2t1", "b1t1", "b2t2")}
            e1, e2 = g["p1t1"][0] + g["p2t2"][0], g["p1t2"][0] + g["p2t1"][0]
            sel = ("p1t1", "p2t2") if e1 <= e2 else ("p1t2", "p2t1")
            qvp = qv(*g[sel[0]]) + qv(*g[sel[1]])
            qvb = qv(*g["b1t1"]) + qv(*g["b2t2"])
            ex = float(e["excess"])
            den = max(len(seqs["t1"]) + len(seqs["t2"]), len(seqs["p1"]) + len(seqs["p2"]))
            qf = qvp / qvb if qvb else None
            rows.append(dict(base, truth_defined=True, predicted=True, excess=ex, normalized_excess=ex / den if den else None, floor_attained=ex == 0,
                             catastrophic=ex >= CATASTROPHIC_BP, QVfrac=qf, QVfrac_above_1=bool(qf is not None and qf > 1),
                             top_category=bool(qf is not None and qf >= TOP), quintile=None if qf is None else min(4, int(qf * 5)) if qf < 1 else 4,
                             distance_matches_evaluator=min(e1, e2) == float(e["unordered_distance"]),
                             floor_matches_evaluator=(g["b1t1"][0] == float(e["truth1_closest_distance"]) and g["b2t2"][0] == float(e["truth2_closest_distance"])),
                             synonymous_artefact_bp=artefact_blocks.get(bid), synonymous_artefact=bool(artefact_blocks.get(bid)) and float(e["unordered_distance"]) > 0,
                             **inp))
    write(out, "blocks", rows)
    write(out, "runs", runs)
    write(out, "gq_correspondence", [dict(blocks_checked=len(corr), max_abs_difference=max(corr) if corr else None,
                                          within_0_01=sum(x < 0.01 for x in corr))])
    bad = [r for r in rows if r.get("predicted") and not (r["distance_matches_evaluator"] and r["floor_matches_evaluator"])]
    print(f"blocks {len(rows)}; evaluator consistency failures {len(bad)}; GQ correspondence {len(corr)} blocks, max |diff| {max(corr) if corr else None}")
    return 0 if not bad else 1


# ------------------------------------------------------------------------------------------------ the frozen comparison
def fit(out):
    # amendment A1: read the outcomes unit with the same suffix (fit.a1 <- outcomes.a1); a missing table is a failure, not an empty fit
    src = os.path.join(os.path.dirname(out), "outcomes" + os.path.basename(out)[len("fit"):], "blocks.tsv")
    if not os.path.exists(src):
        print(f"no outcomes table at {src}")
        return 1
    rows = T(src)
    usable = []
    for r in rows:
        if r["truth_defined"] != "True" or r.get("predicted") != "True" or r["assessable"] != "True" or r["gq"] in ("", "None"):
            continue
        usable.append(dict(r, y=1.0 if r["top_category"] == "True" else 0.0, x=float(r["gq"]),
                           prior=1.0 if r["read_support"] == "PRIOR_ONLY" else 0.0))
    preds = {name: lofo(usable, feats) for name, feats in (("B0", []), ("C1", ["x"]), ("C2", ["x", "prior"]))}
    for r in usable:
        for n in preds:
            r[f"p_{n}"] = preds[n][id(r)]
    perf = []
    rng = random.Random(SEED)
    for locus in ("ALL", "C4", "CYP2D6", "ACOT"):
        S = [r for r in usable if locus == "ALL" or r["locus"] == locus]
        if not S:
            continue
        y = np.array([r["y"] for r in S])
        row = dict(locus=locus, blocks=len(S), families=len({r["family"] for r in S}), donors=len({r["donor"] for r in S}),
                   top_rate=f"{y.mean():.4f}")
        ll = {n: logloss([r[f"p_{n}"] for r in S], y) for n in ("B0", "C1", "C2")}
        for n in ("B0", "C1", "C2"):
            p = np.array([r[f"p_{n}"] for r in S])
            row[f"logloss_{n}"] = f"{ll[n].mean():.4f}"
            row[f"brier_{n}"] = f"{np.mean((p - y) ** 2):.4f}"
            a = auc(list(p), list(y))
            row[f"auc_{n}"] = "NA" if a is None else f"{a:.3f}"
        fams = sorted({r["family"] for r in S})
        idx = {f: [i for i, r in enumerate(S) if r["family"] == f] for f in fams}
        for a_, b_ in (("C1", "B0"), ("C2", "B0"), ("C2", "C1")):
            diff = ll[a_] - ll[b_]
            obs = float(diff.mean())
            bs = []
            for _ in range(BOOT):
                pick = [rng.choice(fams) for _ in fams]
                ii = [i for f in pick for i in idx[f]]
                bs.append(float(diff[ii].mean()))
            bs.sort()
            row[f"dlogloss_{a_}_minus_{b_}"] = f"{obs:.4f} [{bs[int(0.025 * BOOT)]:.4f}, {bs[int(0.975 * BOOT) - 1]:.4f}]"
        perf.append(row)
    write(out, "performance", perf)
    write(out, "predictions", [{k: r[k] for k in ("cell", "locus", "family", "block_id", "kind", "y", "x", "prior", "p_B0", "p_C1", "p_C2",
                                                  "excess", "QVfrac", "read_support", "census_route", "jp_emission_tie", "synonymous_artefact")}
                               for r in usable])
    for p in perf:
        print(p)
    return 0


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "truth":
        os.makedirs(sys.argv[3], exist_ok=True)
        sys.exit(truth(sys.argv[2], sys.argv[3]))
    out = sys.argv[2]
    os.makedirs(out, exist_ok=True)
    sys.exit({"selftest": selftest, "inventory": inventory, "outcomes": outcomes, "fit": fit}[mode](out))
