#!/usr/bin/env python3
"""M10-B driver (frozen with DESIGN.md; no reads, no genotyping; saved models and sequences only).

  m10b.py dump <name> <out>   m10b_walk_dump on one frozen model (MODELS); its fingerprint must equal the frozen value
  m10b.py coverage <out>      what the model supplies, per locus (walk multiplicity and orientation, duplicate-sequence
                              nodes, labelling invariance of raw walks / long-node skeletons / parse-set projections over
                              every degeneracy-inventory equivalence)
  m10b.py cases <out>         the frozen saved-case set (CASES): projections and diploid comparisons
"""
import collections
import csv
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import project  # noqa: E402

WT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
R = "/Users/davide.bolognini/panvar/tests/results"
OUT = f"{R}/m10b_structural_projection_20261020"
C4REF, CYPREF, ACOTREF = "HG00096#1#haplotype1-0000024:31848049-32074186", "grch38#1#chr22:42031864-42245566", "GRCh38#0#chr14:73416151-73668669"
CYPCUR = "/Users/davide.bolognini/panvar/tests/genotype_cyp2d6_curation.tsv"
MODELS = {   # name: graph, reference, curation, frozen fingerprint (from the saved runs that used the model)
    "c4dev": (f"{R}/genotype_c4_production_11fe681/graphs/c4.gfa", C4REF, "none", "d2b8c6c2e0b4079d"),
    "c4loo": (f"{R}/genotype_c4_production_11fe681/graphs/NA19650.loo.gfa", C4REF, "none", None),
    "cypdev": (f"{R}/cyp2d6_shortread_pilot_v1/data/graphs/cyp2d6.gfa", CYPREF, CYPCUR, "2a5b720a1c5bf669"),
    "cyploo": (f"{R}/cyp2d6_shortread_pilot_v1/data/graphs/loo.NA19650.gfa", CYPREF, CYPCUR, None),
    "acotdev": (f"{R}/genotype_perf_bounded_20261015/inputs/acot/acot.gfa", ACOTREF, "none", "d622742f9cf8e7a5"),
    "acotloo": (f"{R}/genotype_perf_bounded_20261015/inputs/acot/acot.gfa", ACOTREF,
                f"{R}/m9h_jp11_validation_20261020/runs/prepare/curation.LOO.txt", "12cbe6beb1ad28d1"),
}
M1B, M9H, M9I, M10A = f"{R}/m1b_structural_20261019/runs", f"{R}/m9h_jp11_validation_20261020/runs", \
    f"{R}/m9i_compound_targets_20261020/runs", f"{R}/m10a_dev_configuration_20261020/runs"
DEGEN = f"{M10A}/declare/degeneracy.tsv"   # dev-model inventory (c4 / cyp / acot)
KV = lambda p: dict(l.rstrip("\n").split("\t", 1) for l in open(p) if "\t" in l and not l.startswith("#"))  # noqa: E731


def write(out, name, rows):
    with open(f"{out}/{name}.tsv", "w") as f:
        if not rows:
            f.write("none\n")
            return
        w = csv.DictWriter(f, fieldnames=list(dict.fromkeys(k for r in rows for k in r)), delimiter="\t", lineterminator="\n", restval="")
        w.writeheader()
        w.writerows(rows)


def model(name):
    return project.Model(f"{OUT}/runs/dump.{name}/{name}")


def dump(name, out):
    g, ref, cur, fp = MODELS[name]
    code = subprocess.run([f"{WT}/build/m10b_walk_dump", g, ref, cur, f"{out}/{name}"]).returncode
    if code:
        return code
    got = KV(f"{out}/{name}.summary.tsv")["model_fingerprint"]
    print(f"{name}: fingerprint {got} (frozen {fp or 'checked against the saved run in cases'})")
    return 0 if fp is None or got == fp else 1


def coverage(out):
    rows, inv = [], []
    for name in MODELS:
        m = model(name)
        mult_blocks, orient_blocks, maxmult = set(), set(), 0
        for (b, a), ws in m.walks.items():
            for w in ws:
                c = collections.Counter(s[:-1] for s in w if len(m.nodes[s[:-1]]) >= project.ANCHOR_BP)
                if c and max(c.values()) > 1:
                    mult_blocks.add(m.block_ids[b])
                    maxmult = max(maxmult, max(c.values()))
                if any(s.endswith("-") and len(m.nodes[s[:-1]]) >= project.ANCHOR_BP for s in w):
                    orient_blocks.add(m.block_ids[b])
        by_seq = collections.Counter(min(s, project.rc(s)) for s in m.nodes.values())
        multi_walk = sum(1 for ws in m.walks.values() if len({tuple(w) for w in ws}) > 1)
        rows.append(dict(model=name, fingerprint=m.fingerprint, blocks=len(m.block_ids), alleles=len(m.seq),
                         blocks_with_repeated_long_node=";".join(sorted(mult_blocks)) or "none", max_long_node_multiplicity=maxmult,
                         blocks_with_reverse_long_node=";".join(sorted(orient_blocks)) or "none",
                         duplicate_sequence_node_classes=sum(1 for n in by_seq.values() if n > 1),
                         alleles_with_differing_walks=multi_walk))
        locus = {"c4dev": "c4", "cypdev": "cyp", "acotdev": "acot"}.get(name)
        if not locus:
            continue
        eq = [r for r in project.T(DEGEN) if r["model"] == locus]
        raw = skel = resolved = ambiguous = 0
        for r in eq:
            f, l = int(r["first_block_index"]), int(r["last_block_index"])
            ta, tb = [list(map(int, r[k].split(","))) for k in ("tuple_a", "tuple_b")]
            wa = [s for i, a in enumerate(ta) for s in m.walks[(f + i, a)][0]]
            wb = [s for i, a in enumerate(tb) for s in m.walks[(f + i, a)][0]]
            raw += wa == wb
            sa = [c for c in (m.step_class(s) for s in wa) if c]
            sb = [c for c in (m.step_class(s) for s in wb) if c]
            skel += sa == sb
            pa, pb = m.project(m.spell(f, ta), f, l), m.project(m.spell(f, tb), f, l)
            assert pa[1] == pb[1], "parse-set projection must be labelling-invariant by construction"
            resolved += pa[0] == "RESOLVED"
            ambiguous += pa[0] == "AMBIGUOUS"
        inv.append(dict(model=name, equivalences=len(eq), identical_raw_walks=raw, identical_long_node_skeletons=skel,
                        parse_set_identical=len(eq), parse_set_resolved=resolved, parse_set_ambiguous=ambiguous))
    write(out, "model_supply", rows)
    write(out, "labelling_invariance", inv)
    print(open(f"{out}/model_supply.tsv").read() + "\n" + open(f"{out}/labelling_invariance.tsv").read())
    return 0


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


def labelling(run):
    sl = KV(f"{run}/call.sequence_locus.tsv")
    return [list(map(int, sl[f"representative_labelling_hap{h}"].split(","))) for h in (1, 2)]


def alts(path):
    out = {}
    for line in open(path):
        if line.startswith("#") or "\t" not in line:
            continue
        f = line.rstrip("\n").split("\t")
        out[f[0]] = [list(map(int, f[1].split(","))), list(map(int, f[2].split(",")))]
    return out


def truth(dirpath, prefix, blocks):
    hs = ["", ""]
    for b in blocks:
        recs = fasta(f"{dirpath}/{prefix}.{b}.truth.fa")
        hs = [hs[0] + recs["truth1"], hs[1] + recs["truth2"]]
    return hs


CASES = [   # frozen: id, model, interval (first..last block id), haplotype pairs (label -> source), comparisons, purpose
    dict(id="C4_NA19650_allocation", model="c4loo", run=f"{M1B}/call.c4_10", interval=("sv:7", "sv:7"),
         pairs={"reported_3+1": ("alts", f"{M1B}/alts.c4_10/c4_10.candidates.tsv", "reported"),
                "saved_2+2": ("alts", f"{M1B}/alts.c4_10/c4_10.candidates.tsv", "saved_2+2"),
                "saved_1+2": ("alts", f"{M1B}/alts.c4_10/c4_10.candidates.tsv", "saved_1+2"),
                "truth": ("truth", f"{M1B}/seq.c4", "seq")},
         compare=[("reported_3+1", "saved_2+2"), ("reported_3+1", "truth"), ("saved_2+2", "truth"), ("reported_3+1", "saved_1+2")],
         purpose="3+1 vs 2+2 at equal total dosage (M1-B: truth 2+2, reported 3+1)"),
    dict(id="C4_HG00096_control", model="c4dev", run=f"{M10A}/D.c4_5", interval=("ctx:0", "ctx:10"),
         pairs={"call": ("call", f"{M10A}/D.c4_5"), "truth_panel": ("panel", "HG00096#")},
         compare=[("call", "truth_panel")], purpose="correct control (development cell; truth = the donor's intact panel templates)"),
    dict(id="CYP_NA19650_dosage", model="cyploo", run=f"{M1B}/call.cyp_5", interval=("sv:9", "sv:9"),
         pairs={"reported_5x_2+2": ("call", f"{M1B}/call.cyp_5"), "twin_30x_3+2": ("call", f"{M1B}/call.cyp_30"),
                "saved_2+3": ("alts", f"{M1B}/alts.cyp_5/cyp_5.candidates.tsv", "saved_2+3")},
         compare=[("reported_5x_2+2", "twin_30x_3+2"), ("reported_5x_2+2", "saved_2+3"), ("twin_30x_3+2", "saved_2+3")],
         purpose="dosage failure (M1-B: truth 3+2, 5x reported 2+2); no saved truth block sequence at sv:9"),
    dict(id="CYP_NA19650_synonymous", model="cyploo", run=f"{M1B}/call.cyp_30", interval=("sv:3", "sv:5"),
         pairs={"call_30x": ("call", f"{M1B}/call.cyp_30"), "truth": ("truth", f"{M9I}/truth.CYP.loo.NA19650", "bs")},
         compare=[("call_30x", "truth")], purpose="block labels wrong by 118+118 bp, compound target exact (M9-I)"),
    dict(id="ACOT_L4_synonymous", model="acotloo", run=f"{M9H}/call.on.L4", interval=("sv:5", "ctx:6"),
         pairs={"jp_on": ("call", f"{M9H}/call.on.L4"), "jp_off": ("call", f"{M9H}/call.off.L4"),
                "truth": ("truth", f"{M9I}/truth.ACOT.LOO.HG02083", "bs")},
         compare=[("jp_on", "truth"), ("jp_off", "truth"), ("jp_on", "jp_off")],
         purpose="JP on: interval exact (labels 3+3 bp); JP off: interval distance 14 (M9-H)"),
    dict(id="ACOT_L2_few_bases", model="acotloo", run=f"{M9H}/call.on.L2", interval=("sv:5", "ctx:6"),
         pairs={"jp_on": ("call", f"{M9H}/call.on.L2"), "jp_off": ("call", f"{M9H}/call.off.L2"),
                "truth": ("truth", f"{M9I}/truth.ACOT.LOO.HG02074", "bs")},
         compare=[("jp_on", "truth"), ("jp_off", "truth"), ("jp_on", "jp_off")],
         purpose="remaining difference a few bases: interval distance 3 (on) / 5 (off) (M9-H)"),
]


def cases(out):
    rows, comps, ok = [], [], True
    cache = {}
    for c in CASES:
        m = cache.setdefault(c["model"], model(c["model"]))
        run_fp = KV(f"{c['run']}/call.mosaic_status.tsv").get("model_fingerprint")
        fp_ok = m.fingerprint == run_fp
        ok = ok and fp_ok
        f, l = m.block(c["interval"][0]), m.block(c["interval"][1])
        seqs, labels = {}, {}
        for label, src in c["pairs"].items():
            if src[0] == "truth":
                seqs[label] = truth(src[1], src[2], m.block_ids[f:l + 1])
                labels[label] = ["truth sequence", "truth sequence"]
            else:
                tuples = labelling(src[1]) if src[0] == "call" else alts(src[1])[src[2]] if src[0] == "alts" else \
                    [v for k, v in sorted(m.panel.items()) if k.startswith(src[1])]
                assert len(tuples) == 2, (c["id"], label, len(tuples))
                seqs[label] = [m.spell(f, t[f:l + 1]) for t in tuples]
                labels[label] = [",".join(map(str, t[f:l + 1])) for t in tuples]
        proj = {}
        for label, hs in seqs.items():
            proj[label] = [m.project(h, f, l) for h in hs]
            for i, p in enumerate(proj[label]):
                prof = [project.profile(s) for s in sorted(p[1])] if p[0] in ("RESOLVED", "AMBIGUOUS") else []
                rows.append(dict(case=c["id"], model=c["model"], model_fingerprint_matches_run=fp_ok, interval="..".join(c["interval"]),
                                 haplotype=f"{label}.h{i + 1}", block_labels=labels[label][i], bp=len(seqs[label][i]),
                                 status=p[0], parses=p[2], skeletons=len(p[1]) if p[0] in ("RESOLVED", "AMBIGUOUS") else 0,
                                 skeleton_steps=";".join(str(x["steps"]) for x in prof) or "-",
                                 max_multiplicity=";".join(str(x["max_multiplicity"]) for x in prof) or "-",
                                 multiplicity_signature=";".join(x["multiplicity_signature"] for x in prof) or "-",
                                 reason=p[1] if p[0] in ("UNKNOWN", "OVERFLOW") else ""))
        for a, b in c["compare"]:
            verdict, detail = project.compare(proj[a], proj[b])
            same_bytes = sorted(seqs[a]) == sorted(seqs[b])
            same_labels = sorted(labels[a]) == sorted(labels[b])
            comps.append(dict(case=c["id"], a=a, b=b, verdict=verdict, detail=detail, exact_sequence_pair_equal=same_bytes,
                              block_labels_equal=same_labels if "truth" not in (a, b) else "n/a (truth has no labels)",
                              purpose=c["purpose"]))
    write(out, "projections", rows)
    write(out, "comparisons", comps)
    print(open(f"{out}/comparisons.tsv").read())
    return 0 if ok else 1


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "dump":
        sys.exit(dump(sys.argv[2], sys.argv[3]))
    out = sys.argv[2]
    os.makedirs(out, exist_ok=True)
    sys.exit({"coverage": coverage, "cases": cases}[mode](out))
