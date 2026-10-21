#!/usr/bin/env python3
"""M14-A: bounded representative-panel experiment (approximate, default-off route; diagnostic).

Mechanism (one): whole panel HAPLOTYPES are removed from the active panel with the caller's existing --panel-curation. The
graph frame and block structure are unchanged (curation never changes blocks); alleles carried only by removed haplotypes
become uncallable; mosaic generation, copying-prior donors and placement run over the retained (representative) haplotypes,
each counted once.

Modes: cluster <out> | analyse <out>
  cluster: per cell, Mash distances from the probe's panel sketches (panel sequences only), deterministic farthest-point
           selection of K = ceil(N/2) and ceil(N/4) representatives seeded at the reference template, nearest-representative
           assignment, curation files (base curation + every non-representative), member->representative tables
  analyse: work and time, per-block excess against BOTH floors, locus-level whole-region distance, transitions vs full panel
"""
import csv
import math
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "m11d"))
import m11d  # noqa: E402

csv.field_size_limit(sys.maxsize)
T, KV, write = m11d.T, m11d.KV, m11d.write
R_, WT = m11d.R, m11d.WT
OUT = f"{R_}/m14a_representative_panel_20261021"
K_MER, SKETCH = 21, 2000
LARGE, NEAR = 1000.0, 10.0
SETTINGS = ("full", "r2", "r4")      # full panel; ~2x; ~4x reduction
REF = {"C4": m11d.C4R, "CYP2D6": m11d.CYPR, "ACOT": m11d.ACOTR}

# frozen cells (DESIGN.md §3): the M11-D development manifest prefix per locus -- all 6 CYP2D6, the first 6 C4, the first 3 ACOT
CELLS = [c for c in m11d.CELLS if c["locus"] == "CYP2D6"] + [c for c in m11d.CELLS if c["locus"] == "C4"][:6] + \
        [c for c in m11d.CELLS if c["locus"] == "ACOT"][:3]


def mash(a, b):
    """Mash distance from two sorted bottom-s sketches."""
    i = j = shared = taken = 0
    while taken < SKETCH and i < len(a) and j < len(b):
        if a[i] == b[j]:
            shared += 1
            i += 1
            j += 1
        elif a[i] < b[j]:
            i += 1
        else:
            j += 1
        taken += 1
    while taken < SKETCH and (i < len(a) or j < len(b)):
        taken += 1
        if i < len(a):
            i += 1
        else:
            j += 1
    jac = shared / taken if taken else 0.0
    return 1.0 if jac == 0 else max(0.0, -math.log(2 * jac / (1 + jac)) / K_MER)


def cluster(out):
    rows, maps = [], []
    for c in CELLS:
        sk = T(f"{OUT}/runs/sketch.{c['id']}/panel.sketch.tsv")
        names = [r["template"] for r in sk]
        sketch = [sorted(int(x) for x in r["sketch"].split(",")) for r in sk]
        n = len(names)
        D = [[0.0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                D[i][j] = D[j][i] = mash(sketch[i], sketch[j])
        ref = REF[c["locus"]]
        seed = names.index(ref) if ref in names else min(range(n), key=lambda i: names[i])
        base = m11d.inputs(c)["cur"]
        base_lines = [l for l in open(base).read().splitlines() if l.strip() and not l.startswith("#")] if base else []
        for setting, frac in (("r2", 2), ("r4", 4)):
            k = math.ceil(n / frac)
            reps = [seed]
            near = list(D[seed])
            while len(reps) < k:   # farthest point; ties by template name
                cand = max((i for i in range(n) if i not in reps), key=lambda i: (near[i], [-ord(ch) for ch in names[i]]))
                reps.append(cand)
                near = [min(near[x], D[cand][x]) for x in range(n)]
            rep_set = set(reps)
            assign = {}
            for i in range(n):
                best = min(reps, key=lambda r: (D[i][r], names[r]))
                assign[i] = best
                maps.append(dict(cell=c["id"], setting=setting, member=names[i], representative=names[best], distance=D[i][best],
                                 is_representative=i in rep_set))
            removed = [names[i] for i in range(n) if i not in rep_set]
            with open(f"{out}/{c['id']}.{setting}.curation.txt", "w") as f:
                f.write(f"# M14-A {setting}: base curation ({base or 'none'}) + {len(removed)} non-representative templates removed "
                        f"(farthest-point K={k} of {n}, Mash k={K_MER} s={SKETCH}, seed {names[seed]})\n")
                for l in base_lines + removed:
                    f.write(l + "\n")
            dist = sorted(D[i][assign[i]] for i in range(n) if i not in rep_set)
            q = lambda p: dist[min(len(dist) - 1, int(p * (len(dist) - 1)))] if dist else 0.0  # noqa: E731
            sizes = sorted((sum(1 for i in range(n) if assign[i] == r) for r in reps), reverse=True)
            rows.append(dict(cell=c["id"], locus=c["locus"], setting=setting, active_templates=n, representatives=k, removed=len(removed),
                             seed=names[seed], member_distance_median=round(q(0.5), 6), member_distance_p90=round(q(0.9), 6),
                             member_distance_max=round(q(1.0), 6), largest_cluster=sizes[0], singleton_clusters=sum(s == 1 for s in sizes)))
    write(out, "clusters", rows)
    write(out, "membership", maps)
    for r in rows:
        print(r)
    return 0


def curation_of(c, setting):
    return m11d.inputs(c)["cur"] if setting == "full" else f"{OUT}/runs/cluster/{c['id']}.{setting}.curation.txt"


# ------------------------------------------------------------------------------------------------ analyse (truth opened here)
STATUS_KEYS = ("status", "seconds_total", "seconds_model", "seconds_markers", "seconds_proposal", "seconds_rescoring", "panel_diplotypes",
               "mosaic_candidates_found", "sequence_candidates", "unique_haplotypes", "physical_origins", "accepted_mate_placements",
               "resources_phys_footprint_peak_bytes")


def provenance_truth(c):
    d = c["donor"]
    if c["locus"] == "C4":
        sim = f"{R_}/genotype_coverage_c4_held_out_abb35a0/{d}/source"
    elif c["locus"] == "CYP2D6":
        sim = f"{R_}/cyp2d6_shortread_pilot_v1/data/reads/{d}/cov30"
    else:
        sim = m11d.ACOT_READS[d].replace("/cov5", "/cov30")
    s = KV(f"{sim}/simulation.tsv")
    t = m11d.fasta(f"{sim}/truth.fa")
    return t[s["haplotype_0"]], t[s["haplotype_1"]]


def locus_distance(c, setting, out):
    fa = f"{OUT}/runs/call.{c['id']}.{setting}/call.sequence_haplotypes.fa"
    if not os.path.exists(fa):
        return None
    h = m11d.fasta(fa)
    hs = [h[k] for k in sorted(h)]
    if len(hs) != 2:
        return None
    t1, t2 = provenance_truth(c)
    tag = f"{c['id']}.{setting}"
    os.makedirs(f"{out}/seq", exist_ok=True)
    for nm, s in (("h1", hs[0]), ("h2", hs[1]), ("t1", t1), ("t2", t2)):
        open(f"{out}/seq/{tag}.{nm}.txt", "w").write(s)
    man = f"{out}/{tag}.manifest.tsv"
    with open(man, "w") as f:
        for x, y in (("h1", "t1"), ("h2", "t2"), ("h1", "t2"), ("h2", "t1")):
            f.write(f"{x}{y}\t{out}/seq/{tag}.{x}.txt\t{out}/seq/{tag}.{y}.txt\t-\n")
    subprocess.run([f"{WT}/build/m10c_align_eval", "--distance", man, f"{out}/{tag}.distance.tsv"], check=True, capture_output=True)
    d = {r["id"]: int(r["edits"]) for r in T(f"{out}/{tag}.distance.tsv")}
    for nm in ("h1", "h2", "t1", "t2"):
        os.remove(f"{out}/seq/{tag}.{nm}.txt")
    return min(d["h1t1"] + d["h2t2"], d["h1t2"] + d["h2t1"])


def ledger():
    p = f"{OUT}/ledger.tsv"
    if not os.path.exists(p):
        return {}
    lines = [l.rstrip("\n").split("\t") for l in open(p)]
    return {r[0]: dict(zip(lines[0], r)) for r in lines[1:]}


def analyse(out):
    led = ledger()
    runs, blocks = [], []
    for c in CELLS:
        full_ev = {r["block_id"]: r for r in T(f"{OUT}/runs/eval.{c['id']}.full/eval.blocks.tsv")}
        for setting in SETTINGS:
            st = KV(f"{OUT}/runs/call.{c['id']}.{setting}/call.mosaic_status.tsv")
            sl = KV(f"{OUT}/runs/call.{c['id']}.{setting}/call.sequence_locus.tsv")
            es = KV(f"{OUT}/runs/eval.{c['id']}.{setting}/eval.summary.tsv")
            ev = T(f"{OUT}/runs/eval.{c['id']}.{setting}/eval.blocks.tsv")
            call_led = led.get(f"call.{c['id']}.{setting}", {})
            row = dict(cell=c["id"], locus=c["locus"], setting=setting, call_exit=call_led.get("exit", "NA"), call_wall=call_led.get("wall_seconds", "NA"),
                       call_peak_gib=call_led.get("peak_gib", "NA"), sequence_status=sl.get("sequence_status", "NA"),
                       completion=sl.get("completion", "NA"), eval_status=es.get("status", "NOT_EVALUATED"))
            for k in STATUS_KEYS:
                row[k] = st.get(k, "NA")
            ld = locus_distance(c, setting, out)
            row["locus_distance"] = ld if ld is not None else "NA"
            tot_full_excess = tot_red_excess = 0.0
            nblk = 0
            for e in ev:
                if e["excess"] in ("NA", "") or e["block_id"] not in full_ev:
                    continue
                f = full_ev[e["block_id"]]
                d = float(e["unordered_distance"])
                floor_full = float(f["availability_floor"])
                floor_red = float(e["availability_floor"])
                ex_full = d - floor_full
                den = max(float(e["truth1_bp"]) + float(e["truth2_bp"]), float(e["called1_bp"]) + float(e["called2_bp"]))
                base_full_excess = float(f["excess"]) if f["excess"] not in ("NA", "") else None
                blocks.append(dict(cell=c["id"], locus=c["locus"], setting=setting, block_id=e["block_id"], kind=e["kind"], distance=d,
                                   floor_full=floor_full, floor_reduced=floor_red, floor_increase=floor_red - floor_full,
                                   excess_vs_full_floor=ex_full, normalized_excess_vs_full_floor=ex_full / den if den else None,
                                   excess_vs_reduced_floor=d - floor_red, full_panel_excess=base_full_excess,
                                   delta_excess_vs_full_panel=(ex_full - base_full_excess) if base_full_excess is not None else None))
                tot_full_excess += ex_full
                tot_red_excess += d - floor_red
                nblk += 1
            row.update(blocks_evaluated=nblk, excess_sum_vs_full_floor=tot_full_excess, excess_sum_vs_reduced_floor=tot_red_excess)
            runs.append(row)
    write(out, "runs", runs)
    write(out, "blocks", blocks)
    # transitions against the matched full-panel run (excess measured against the ORIGINAL full-panel floor)
    trans = []
    for b in blocks:
        if b["setting"] == "full" or b["full_panel_excess"] is None:
            continue
        o, s = b["full_panel_excess"], b["excess_vs_full_floor"]
        cls = ("catastrophe_introduced" if o < LARGE <= s else "catastrophe_removed" if s < LARGE <= o else
               "catastrophe_worsened" if o >= LARGE and s > o else "catastrophe_reduced" if o >= LARGE and s < o else
               "catastrophe_unchanged" if o >= LARGE else "control_damaged" if o <= NEAR < s else
               "unchanged" if s == o else "worse" if s > o else "better")
        trans.append(dict(b, transition=cls))
    write(out, "transitions", trans)
    summ = []
    for setting in ("r2", "r4"):
        for locus in ("ALL", "C4", "CYP2D6", "ACOT"):
            S = [t for t in trans if t["setting"] == setting and (locus == "ALL" or t["locus"] == locus)]
            rf = [r for r in runs if r["setting"] == "full" and (locus == "ALL" or r["locus"] == locus)]
            rs = [r for r in runs if r["setting"] == setting and (locus == "ALL" or r["locus"] == locus)]
            wall = lambda R: sum(float(r["call_wall"]) for r in R if r["call_wall"] not in ("NA", ""))  # noqa: E731
            cnt = lambda k: sum(t["transition"] == k for t in S)  # noqa: E731
            summ.append(dict(setting=setting, locus=locus, blocks=len(S), call_wall_full=round(wall(rf), 1), call_wall_reduced=round(wall(rs), 1),
                             speedup=round(wall(rf) / wall(rs), 3) if wall(rs) else "NA",
                             excess_sum_full=sum(float(r["excess_sum_vs_full_floor"]) for r in rf),
                             excess_sum_reduced_vs_full_floor=sum(float(r["excess_sum_vs_full_floor"]) for r in rs),
                             catastrophe_introduced=cnt("catastrophe_introduced"), catastrophe_worsened=cnt("catastrophe_worsened"),
                             catastrophe_reduced=cnt("catastrophe_reduced"), catastrophe_removed=cnt("catastrophe_removed"),
                             control_damaged=cnt("control_damaged"), worse=cnt("worse"), better=cnt("better"), unchanged=cnt("unchanged"),
                             no_call_or_incomplete_full=sum(r["eval_status"] not in ("COMPLETE", "CONDITIONAL_GENOTYPE_EQUIVALENT") for r in rf),
                             no_call_or_incomplete_reduced=sum(r["eval_status"] not in ("COMPLETE", "CONDITIONAL_GENOTYPE_EQUIVALENT") for r in rs)))
    write(out, "summary", summ)
    for s in summ:
        print(s)
    return 0


if __name__ == "__main__":
    mode, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    sys.exit({"cluster": cluster, "analyse": analyse}[mode](out))
