#!/usr/bin/env python3
"""M12-C: bounded completion and interpretation of M12-B (M12-B's frozen rule and results preserved).

Modes:
  bind <out>      verify the continuation's bindings (executables, frozen candidate lists, plan, selection, inputs) against
                  M12-B's records, and link M12-B's frozen select/count/score outputs into this package (read-only reuse)
  evaluate <out>  M12-B's own evaluate() over the augmented set (its 12 scored blocks + the 6 completed here), unchanged code;
                  a reproduction check of M12-B's 12 rows; the K5 record comparison; severity transitions; the combined
                  coverage ledger; the frozen decision rule applied to the augmented results (M12-B's decision kept apart)
"""
import collections
import csv
import hashlib
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "m12b"))
import m12b  # noqa: E402

csv.field_size_limit(sys.maxsize)
T, KV, write = m12b.T, m12b.KV, m12b.write
M12B = m12b.OUT
OUT = f"{m12b.m11d.R}/m12c_dosage_first_completion_20261021"
CONTINUE = ["B10", "B12", "B13", "B14", "B16", "B17"]


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def bind(out):
    rows, ok = [], True
    man = open(f"{M12B}/MANIFEST.tsv").read()
    for exe in ("genotype_probe", "m10c_align_eval"):
        want = re.search(rf"# executable\t{exe}\t([0-9a-f]+)", man).group(1)
        got = sha(f"{m12b.WT}/build/{exe}")
        rows.append(dict(item=f"executable {exe}", expected=want, observed=got, ok=want == got))
    sel = {r["bid"]: r for r in T(f"{M12B}/runs/select/selection.tsv")}
    plan = {r["bid"]: r for r in T(f"{M12B}/runs/select/score_plan.tsv")}
    led = {r["unit"]: r for r in m12b.m12b_ledger()}
    for b in CONTINUE:
        p = f"{M12B}/runs/select/{b}.candidates.tsv"
        n = sum(1 for _ in open(p)) - 1
        rows.append(dict(item=f"{b} candidate list", expected=f"{sel[b]['candidates_in_selected']} candidates + incumbent",
                         observed=f"{n} candidates; sha256 {sha(p)[:16]}", ok=str(n) == sel[b]["candidates_in_selected"]))
        rows.append(dict(item=f"{b} M12-B plan and status", expected="admitted by bound; NOT RUN on budget",
                         observed=f"estimate {plan[b]['estimate_s']} s, timeout {plan[b]['timeout_s']}; ledger {led[f'score.{b}']['note'][:40]}",
                         ok=sel[b]["admitted"] == "True" and led[f"score.{b}"]["exit"] == "not_run"))
        cell = m12b.m11d.cell(sel[b]["cell"])
        reads = m12b.m11d.inputs(cell)["reads"] + "/reads_1.fq"
        want = re.search(rf"# block\t{b}\t[^\n]*reads1_sha256 ([0-9a-f]+)", man).group(1)
        rows.append(dict(item=f"{b} reads", expected=want, observed=sha(reads)[:16], ok=sha(reads)[:16] == want))
    ok = all(r["ok"] for r in rows)
    write(out, "bindings", rows)
    # read-only reuse: link M12-B's frozen outputs (selection, counts, the 12 scored blocks) into this package's runs/
    runs = f"{OUT}/runs"
    links = ["select"] + [f"count.{b['bid']}" for b in m12b.blocks_manifest()] + \
            [u for u in (f"score.{b['bid']}" for b in m12b.blocks_manifest()) if os.path.isdir(f"{M12B}/runs/{u}")]
    for name in links:
        dst = f"{runs}/{name}"
        if not os.path.lexists(dst):
            os.symlink(f"{M12B}/runs/{name}", dst)
    print(f"bindings {sum(r['ok'] for r in rows)}/{len(rows)} ok; linked {len(links)} M12-B outputs read-only")
    return 0 if ok else 1


def severity(x):
    x = float(x)
    return "floor" if x == 0 else "near_floor(<=10)" if x <= 10 else "moderate(<1000)" if x < 1000 else "large(>=1000)"


def evaluate(out):
    m12b.OUT = OUT                      # M12-B's unchanged evaluate(), reading this package (links + 6 new score units)
    rc = m12b.evaluate(out)
    rows = {r["bid"]: r for r in T(f"{out}/evaluation.tsv")}
    old = {r["bid"]: r for r in T(f"{M12B}/runs/evaluate.a3/evaluation.tsv")}
    # 1. M12-B's 12 scored rows reproduce exactly (outcome fields)
    keys = ("outcome", "selected_alleles", "selected_state", "selected_excess", "selected_distance", "d_total_vs_incumbent", "reproduced")
    repro = [dict(bid=b, identical=all(rows[b].get(k) == old[b].get(k) for k in keys)) for b in old if old[b].get("reproduced")]
    # 2. severity transitions over every scored block
    sev = []
    for b, r in rows.items():
        if "selected_excess_max" not in r or r["selected_excess_max"] == "":
            continue
        o, s = float(r["original_excess"]), float(r["selected_excess_max"])
        sev.append(dict(bid=b, cell=r["cell"], original_excess=o, selected_excess=s, delta_excess=s - o,
                        original_severity=severity(o), selected_severity=severity(s),
                        change=("unchanged" if s == o else "improved" if s < o else "worsened"),
                        worsened_while_large=o >= 1000 and s > o, source="M12-C" if b in CONTINUE else "M12-B"))
    write(out, "severity", sev)
    # 3. combined coverage ledger: every original block
    ledger = []
    for b in m12b.blocks_manifest():
        r = rows[b["bid"]]
        src = "M12-C (completed here)" if b["bid"] in CONTINUE else ("M12-B" if r.get("reproduced") else "not scored")
        ledger.append(dict(bid=b["bid"], cell=b["cell"], block=b["block"], scored_in=src, outcome=r["outcome"],
                           original_excess=r["original_excess"], selected_excess=r.get("selected_excess", ""),
                           original_state=r.get("original_state", ""), selected_state=r.get("selected_state", r.get("selected_states", "")),
                           truth_total_state=r.get("truth_total_state", ""), total_state_recovered=r.get("total_copy_recovered", ""),
                           allocation_matches_truth=r.get("allocation_matches_truth", ""), limitation=r.get("limitation", "")))
    write(out, "coverage_ledger", ledger)
    # 4. the K5 record comparison (existing M12-B records only; no new scoring)
    k5 = k5_comparison(out)
    summ = T(f"{out}/summary.tsv")[0]
    summ.update(m12b_historical_decision=T(f"{M12B}/runs/evaluate.a3/summary.tsv")[0]["decision"],
                m12b_rows_reproduced=f"{sum(x['identical'] for x in repro)}/{len(repro)}",
                worsened_while_large=sum(x["worsened_while_large"] for x in sev),
                worsened_any=sum(x["change"] == "worsened" for x in sev))
    write(out, "summary", [summ])
    write(out, "m12b_reproduction", repro)
    print(json.dumps(summ))
    for x in k5:
        print(x)
    return rc


def k5_comparison(out):
    bid, cell, blk = "B06", "cyp.NA19983", "sv:9"
    feats = [list(map(int, r["copies_by_allele"].split(","))) for r in T(f"{M12B}/runs/count.{bid}/ca.features.tsv")]
    kmers = [r["kmers"] for r in T(f"{M12B}/runs/count.{bid}/ca.features.tsv")]
    loss = {(int(r["allele_h1"]), int(r["allele_h2"])): r["loss"] for r in T(f"{M12B}/runs/count.{bid}/ca.catalogue.tsv") if r["profile"] == "real"}
    sc = {r["label"]: r for r in T(f"{M12B}/runs/score.{bid}/ss.scores_truth_free.tsv")}
    by_pair = {tuple(map(int, l.split("_")[1:])): r for l, r in sc.items() if l != "incumbent"}
    ev = next(x for x in m12b.m11e.eval_rows(cell) if x["block_id"] == blk)
    al = m12b.fasta(f"{m12b.M11D}/runs/truth.{cell}/bs.{blk}.alleles.fa")
    tr = m12b.fasta(f"{m12b.M11D}/runs/truth.{cell}/bs.{blk}.truth.fa")
    floor = float(ev["availability_floor"])
    t1s, _ = m12b.truth_state(bid, tr["truth1"])
    t2s, _ = m12b.truth_state(bid, tr["truth2"])
    rows = []
    inc = sc["incumbent"]
    for label, pair in (("selected", (104, 97)), ("evaluator_best_A", (95, 10)), ("evaluator_best_B", (10, 95)),
                        ("incumbent", (int(ev["called_allele1"]), int(ev["called_allele2"])))):
        r = inc if label == "incumbent" else by_pair[pair]
        st = [tuple(f[a] for f in feats) for a in pair]
        d = m12b.distances(f"k5.{pair[0]}_{pair[1]}", (al[f"allele_{pair[0]}"], al[f"allele_{pair[1]}"]), tr["truth1"], tr["truth2"], out)
        rows.append(dict(candidate=label, pair=f"{pair[0]},{pair[1]}", hap_states=f"{st[0]}|{st[1]}",
                         total_state=str(tuple(x + y for x, y in zip(*st))), feature0_total=st[0][0] + st[1][0],
                         count_loss=loss.get(pair, "NA"), fragment=r["fragment_term"], exposure=r["exposure_term"], prior=r["prior"],
                         total=r["total"], d_total_vs_selected=None, distance=d, excess=d - floor,
                         block_bp=f"{len(al[f'allele_{pair[0]}'])}+{len(al[f'allele_{pair[1]}'])}"))
    sel_total = float(rows[0]["total"])
    for x in rows:
        x["d_total_vs_selected"] = float(x["total"]) - sel_total
    rows.append(dict(candidate="truth (diagnostic)", pair="-", hap_states=f"{t1s}|{t2s}", total_state=str(tuple(a + b for a, b in zip(t1s, t2s))),
                     feature0_total=t1s[0] + t2s[0], block_bp=f"{len(tr['truth1'])}+{len(tr['truth2'])}", excess="-"))
    write(out, "k5_records", rows)
    write(out, "k5_features", [dict(feature=i, kmers=k) for i, k in enumerate(kmers)])
    return rows


if __name__ == "__main__":
    mode, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    sys.exit({"bind": bind, "evaluate": evaluate}[mode](out))
