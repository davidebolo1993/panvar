#!/usr/bin/env python3
"""M9-G gates and analysis (frozen with DESIGN.md).

  m9g.py parity <out>     C++ JP1.1 (probe junction-proposals) vs the M9-F Python reference on the eight M9-F contexts:
                          per-pair scores; top-32 membership vs the JP1.1 reference (M9-F scores, content tie-break) and the
                          disclosed membership changes vs frozen JP1 (index tie-break)
  m9g.py refusal <out>    the CLI refuses experimental_v1 with repair_rounds 0 before writing any output
  m9g.py offcheck <out>   "off" arm vs the archived M9-C outputs: scientific files byte-identical (declared allowlists for
                          the rounds and branch tables; timing, memory and cache-work columns excluded)
  m9g.py freeze <out>     both arms' predictions, before any truth-based evaluation
  m9g.py analyse <out>    paired off/on comparison
"""
import csv
import hashlib
import os
import subprocess
import sys

mode, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
R = "/Users/davide.bolognini/panvar/tests/results"
C = f"{R}/m9c_acot_reliability_20261020/runs"
D = f"{R}/m9d_attribution_20261020/runs"
F = f"{R}/m9f_junction_proposals_20261020/runs"
OUT = f"{R}/m9g_junction_integration_20261020/runs"
WT = "/Users/davide.bolognini/panvar-worktrees/s1-d-panel-references"
T = lambda p: list(csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t")) if os.path.exists(p) else []
KV = lambda p: {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(p) if "\t" in l and not l.startswith("#")} if os.path.exists(p) else {}
F8 = ["LOO.d1", "LOO.d3", "LOO.d6", "LZO.d1", "LZO.d2", "LOO.d2", "LOO.d4", "LOO.d5"]
PAIRED = ["LOO.d1", "LOO.d3", "LOO.d6", "LZO.d1", "LZO.d2"]
LED = {r["unit"]: r for r in T(f"{OUT}/../ledger.tsv")}


def write(name, rows):
    with open(f"{out}/{name}.tsv", "w") as f:
        if not rows:
            f.write("none\n")
            return
        w = csv.DictWriter(f, fieldnames=list(dict.fromkeys(k for r in rows for k in r)), delimiter="\t", lineterminator="\n", restval="")
        w.writeheader()
        w.writerows(rows)


def dump(name):
    blocks = []
    for r in T(f"{D}/models/{name}.block_alleles.tsv"):
        b = int(r["block_index"])
        while len(blocks) <= b:
            blocks.append([])
        blocks[b].append(r["sequence"])
    return blocks


if mode == "prepare":
    P = f"{C}/prepare"
    base = open(f"{P}/profile.tsv").read()
    open(f"{out}/profile.on.tsv", "w").write(base.replace("profile\tm9c-acot-v2-r0.01", "profile\tm9g-acot-v2-r0.01-jp11")
                                            + "# M9-G: the only difference from the off arm (M9-C's own profile, used unchanged)\nrepair_junction_proposals\texperimental_v1\n")
    sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
    rows = [dict(item="panvar", sha256=sha(f"{WT}/build/panvar")), dict(item="genotype_probe", sha256=sha(f"{WT}/build/genotype_probe")),
            dict(item="profile.off (M9-C prepare/profile.tsv)", sha256=sha(f"{P}/profile.tsv")), dict(item="profile.on", sha256=sha(f"{out}/profile.on.tsv"))]
    diff = [l for l in open(f"{out}/profile.on.tsv").read().splitlines() if l not in base.splitlines()]
    rows.append(dict(item="profile.on minus profile.off", sha256=" | ".join(diff)))
    write("prepare", rows)
    print(open(f"{out}/prepare.tsv").read())
    sys.exit(0)

if mode == "parity":
    rows = []
    for cell in F8:
        cpp = T(f"{OUT}/parity.{cell}/jp.junction_ranks.tsv")
        py = T(f"{F}/propose.{cell}/ranks.tsv")
        if not cpp or not py:
            rows.append(dict(cell=cell, status="missing output"))
            continue
        b5 = [r["block_id"] for r in T(f"{C}/call.{cell}/call.sequence_blocks.tsv")].index("sv:5")
        al = dump("LOO" if cell.startswith("LOO") else cell)[b5]
        P = {(int(r["x"]), int(r["y"])): r for r in py}
        worst = max(max(abs(float(c[k]) - float(P[(int(c["x"]), int(c["y"]))][p])) for k, p in (("score_ab", "score_ab"), ("score_ba", "score_ba"), ("best", "best")))
                    for c in cpp)
        cpp_top = [(int(r["x"]), int(r["y"])) for r in cpp[:32]]
        # JP1.1 reference from the Python scores: Python tie classes, content order inside each class
        ref = sorted(py, key=lambda r: (int(r["tie_class"]), min(al[int(r["x"])], al[int(r["y"])]), max(al[int(r["x"])], al[int(r["y"])]), int(r["x"]), int(r["y"])))
        ref_top = [(int(r["x"]), int(r["y"])) for r in ref[:32]]
        jp1_top = [(int(r["x"]), int(r["y"])) for r in sorted(py, key=lambda r: int(r["rank"]))[:32]]
        rows.append(dict(cell=cell, pairs=len(cpp), max_abs_score_diff=worst, distinct_classes_cpp=len({r["tie_class"] for r in cpp}),
                         distinct_classes_python=len({r["tie_class"] for r in py}), top32_equals_jp11_reference=set(cpp_top) == set(ref_top),
                         top32_order_equals_jp11_reference=cpp_top == ref_top, top32_membership_changes_vs_frozen_jp1=len(set(cpp_top) ^ set(jp1_top)) // 2,
                         ok=worst <= 1e-6 and set(cpp_top) == set(ref_top)))
    write("parity", rows)
    print(open(f"{out}/parity.tsv").read())
    sys.exit(0 if all(r.get("ok") for r in rows) else 1)

if mode == "refusal":
    P = f"{C}/prepare"
    prof = f"{out}/profile.on.tsv"
    lines = open(f"{OUT}/prepare/profile.on.tsv").read()
    if len(sys.argv) > 3 and sys.argv[3] == "legacy":   # Amendment 2: without v2, whose own repair-rounds refusal fires first
        lines = "".join(l + "\n" for l in lines.splitlines() if not l.startswith("output_schema\t"))
    open(prof, "w").write(lines)
    reads = f"{P}/reads/NA19909/cov5"
    cmd = [f"{WT}/build/panvar", "genotype", "--profile", prof, "--threads", "4", "--gfa", f"{R}/genotype_perf_bounded_20261015/inputs/acot/acot.gfa",
           "--reference-path", "GRCh38#0#chr14:73416151-73668669", "--reads1", f"{reads}/reads_1.fq", "--reads2", f"{reads}/reads_2.fq",
           "--library", f"{reads}/library.txt", "--out-prefix", f"{out}/refused", "--quiet", "--panel-curation", f"{P}/curation.LOO.txt",
           "--dosage-mode", "audit", "--repair-rounds", "0"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    files = [f for f in os.listdir(out) if f.startswith("refused")]
    msg = (r.stdout + r.stderr).strip().splitlines()[-1] if (r.stdout + r.stderr).strip() else ""
    ok = r.returncode != 0 and "repair_junction_proposals" in msg and not files
    write("refusal", [dict(exit=r.returncode, message=msg[:300], outputs_written=",".join(files) or "none", ok=ok)])
    print(open(f"{out}/refusal.tsv").read())
    sys.exit(0 if ok else 1)

ROUND_COLS = ["round", "start_total", "moved", "repaired_block", "local_pair_rank", "assignment", "allele_x", "allele_y", "hap1_alleles", "hap2_alleles",
              "fragment_term", "exposure_term", "prior", "total", "delta_total", "status", "unique_best", "separation", "neighbours_generated",
              "neighbours_scored", "duplicate", "suppressed_visited", "self_skipped", "rejected_infinite_prior", "unique_haplotypes", "scored_pairs",
              "order_invariant", "swap_invariant"]
BRANCH_COLS = ["branch", "parent_tie", "start_hap1", "start_hap2", "start_total", "rounds", "moves", "stop", "outcome", "terminal_hap1", "terminal_hap2",
               "terminal_genotype", "terminal_total"]
WHOLE = ["sequence_locus.tsv", "sequence_blocks.tsv", "sequence_haplotypes.tsv", "sequence_haplotypes.fa", "sequence_alternatives.tsv",
         "representation_intervals.tsv", "locus_genotype.tsv", "block_genotypes.tsv", "genotypes.tsv", "mosaic_scores.tsv", "mosaic_candidates.tsv",
         "equivalence_set.tsv", "dosage_status.tsv", "settings.tsv"]
STATUS_KEYS = ["status", "model_fingerprint", "sequence_candidates", "mosaic_candidates", "proposal_states_pruned", "numeric_separation_gap",
               "winner_candidate_set_probability", "repair_rounds_run", "phase_support", "tie_branches_run", "tie_distinct_terminal_genotypes",
               "sequence_prediction", "sequence_status", "sequence_fasta", "marker_depth_rate_per_copy", "marker_background"]


def projected(path, cols):
    rows = T(path)
    return "\n".join("\t".join(r.get(c, "") for c in cols) for r in rows)


if mode == "offcheck":
    rows = []
    for cell in PAIRED:
        a, b = f"{C}/call.{cell}/call", f"{OUT}/call.off.{cell}/call"
        for f in WHOLE:
            pa, pb = f"{a}.{f}", f"{b}.{f}"
            if not os.path.exists(pa) and not os.path.exists(pb):
                continue
            # each arm's own output prefix (written into locus tables and file links) is normalised; nothing else
            norm = lambda p, pre: open(p).read().replace(pre, "<PREFIX>") if os.path.exists(p) else None
            same = norm(pa, a) is not None and norm(pa, a) == norm(pb, b)
            rows.append(dict(cell=cell, file=f, identical=same))
        rows.append(dict(cell=cell, file="repair_rounds.tsv (scientific columns)", identical=projected(f"{a}.repair_rounds.tsv", ROUND_COLS) == projected(f"{b}.repair_rounds.tsv", ROUND_COLS)))
        if os.path.exists(f"{a}.tie_branches.tsv") or os.path.exists(f"{b}.tie_branches.tsv"):
            rows.append(dict(cell=cell, file="tie_branches.tsv (scientific columns)", identical=projected(f"{a}.tie_branches.tsv", BRANCH_COLS) == projected(f"{b}.tie_branches.tsv", BRANCH_COLS)))
        sa, sb = KV(f"{a}.mosaic_status.tsv"), KV(f"{b}.mosaic_status.tsv")
        diff = [k for k in STATUS_KEYS if sa.get(k) != sb.get(k)]
        rows.append(dict(cell=cell, file="mosaic_status.tsv (scientific keys)", identical=not diff, detail=",".join(diff)))
        rows.append(dict(cell=cell, file="mosaic_status.tsv has no junction rows", identical=not any(k.startswith("junction") or k == "repair_junction_proposals" for k in sb)))
    write("offcheck", rows)
    bad = [r for r in rows if not r["identical"]]
    print(open(f"{out}/offcheck.tsv").read())
    print(f"identical {len(rows) - len(bad)}/{len(rows)}")
    sys.exit(0 if not bad else 1)

if mode == "freeze":
    rows = []
    for arm in ("off", "on"):
        for cell in PAIRED:
            d = f"{OUT}/call.{arm}.{cell}/call"
            sl = KV(f"{d}.sequence_locus.tsv")
            led = LED.get(f"call.{arm}.{cell}", {})
            rows.append(dict(arm=arm, cell=cell, exit=led.get("exit", "not run"), timed_out=led.get("timed_out", ""), prediction=sl.get("prediction", "absent"),
                             sequence_status=sl.get("sequence_status", "NA"), provisional=sl.get("provisional", "NA"), legacy_run_status=sl.get("legacy_run_status", "NA"),
                             legacy_no_call_reason=sl.get("legacy_no_call_reason", "NA"), hap1=sl.get("representative_hap1_sha256", "NA")[:16],
                             hap2=sl.get("representative_hap2_sha256", "NA")[:16], decision_set_groups=sl.get("decision_set_groups", "NA"),
                             representative_score=sl.get("representative_score", "NA"),
                             sha_sequence_blocks=hashlib.sha256(open(f"{d}.sequence_blocks.tsv", "rb").read()).hexdigest()[:16] if os.path.exists(f"{d}.sequence_blocks.tsv") else "absent"))
    write("freeze", rows)
    print(open(f"{out}/freeze.tsv").read())
    sys.exit(0)

if mode == "analyse":
    res, blocks = [], []
    for cell in PAIRED:
        arms = {}
        for arm in ("off", "on"):
            d = f"{OUT}/call.{arm}.{cell}/call"
            sl, st = KV(f"{d}.sequence_locus.tsv"), KV(f"{d}.mosaic_status.tsv")
            lg = (T(f"{d}.locus_genotype.tsv") or [{}])[0]
            ev = KV(f"{OUT}/eval.{arm}.{cell}/eval.summary.tsv")
            eb = {r["block_id"]: r for r in T(f"{OUT}/eval.{arm}.{cell}/eval.blocks.tsv")}
            sb = {r["block_id"]: r for r in T(f"{d}.sequence_blocks.tsv")}
            moves = [f"{r['repaired_block']}:{r['allele_x']}/{r['allele_y']}" for r in T(f"{d}.repair_rounds.tsv") if r.get("moved") == "true"]
            arms[arm] = dict(sl=sl, st=st, lg=lg, ev=ev, eb=eb, sb=sb, moves=moves)
        o, n = arms["off"], arms["on"]
        exact = lambda a, b: (lambda e, s: "NA" if not e or not s or s["block_status"] != "IDENTIFIED" else
                              ("unattainable" if not (e["truth1_closest_distance"] == "0" and e["truth2_closest_distance"] == "0") else
                               sorted(map(int, s["genotype"].split("/"))) == sorted((int(e["truth1_closest_allele"]), int(e["truth2_closest_allele"])))))(a["eb"].get(b), a["sb"].get(b))
        changed = [b for b in o["sb"] if o["sb"][b].get("genotype") != n["sb"].get(b, {}).get("genotype")]
        for b in o["sb"]:
            if b in changed or b == "sv:5":
                blocks.append(dict(cell=cell, block=b, off_genotype=o["sb"][b].get("genotype"), on_genotype=n["sb"].get(b, {}).get("genotype"),
                                   off_status=o["sb"][b].get("block_status"), on_status=n["sb"].get(b, {}).get("block_status"),
                                   off_error_mass_upper=o["sb"][b].get("error_mass_upper"), on_error_mass_upper=n["sb"].get(b, {}).get("error_mass_upper"),
                                   off_exact=exact(o, b), on_exact=exact(n, b),
                                   off_legacy_distance=o["eb"].get(b, {}).get("unordered_distance"), on_legacy_distance=n["eb"].get(b, {}).get("unordered_distance")))
        g = lambda a, k: a["lg"].get(k, "NA")
        res.append(dict(cell=cell, off_prediction=f"{o['sl'].get('prediction')}/{o['sl'].get('sequence_status')}", on_prediction=f"{n['sl'].get('prediction')}/{n['sl'].get('sequence_status')}",
                        same_sequence_pair=(o["sl"].get("representative_hap1_sha256"), o["sl"].get("representative_hap2_sha256")) == (n["sl"].get("representative_hap1_sha256"), n["sl"].get("representative_hap2_sha256")),
                        blocks_changed=",".join(changed) or "none",
                        off_excess=o["ev"].get("excess_over_floor_total"), on_excess=n["ev"].get("excess_over_floor_total"),
                        off_distance=o["ev"].get("unordered_distance_total"), on_distance=n["ev"].get("unordered_distance_total"), floor=o["ev"].get("availability_floor_total"),
                        off_total=g(o, "total_score"), on_total=g(n, "total_score"), off_fragment=g(o, "fragment_term"), on_fragment=g(n, "fragment_term"),
                        off_exposure=g(o, "exposure_term"), on_exposure=g(n, "exposure_term"), off_prior=g(o, "sequence_pair_prior"), on_prior=g(n, "sequence_pair_prior"),
                        off_moves=";".join(o["moves"]), on_moves=";".join(n["moves"]),
                        off_tie_branches=o["st"].get("tie_branches_run", "0"), on_tie_branches=n["st"].get("tie_branches_run", "0"),
                        on_junction_eligible=n["st"].get("junction_eligible_blocks"), on_junction_invocations=n["st"].get("junction_invocations"),
                        on_junction_cache_hits=n["st"].get("junction_cache_hits"), on_junction_seconds=n["st"].get("junction_seconds"),
                        on_junction_cutoff_tie_crossings=n["st"].get("junction_cutoff_tie_crossings"), on_read_table_bytes=n["st"].get("junction_read_table_bytes"),
                        on_read_table_seconds=n["st"].get("junction_read_table_seconds"),
                        off_seconds=LED.get(f"call.off.{cell}", {}).get("wall_seconds"), on_seconds=LED.get(f"call.on.{cell}", {}).get("wall_seconds"),
                        off_peak_gib=LED.get(f"call.off.{cell}", {}).get("peak_rss_gib"), on_peak_gib=LED.get(f"call.on.{cell}", {}).get("peak_rss_gib")))
    write("paired", res)
    write("blocks", blocks)
    for r in res:
        print(r)
    print(open(f"{out}/blocks.tsv").read())
    sys.exit(0)
