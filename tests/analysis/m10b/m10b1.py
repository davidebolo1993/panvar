#!/usr/bin/env python3
"""M10-B.1 driver (frozen with DESIGN_B1.md; no reads, no genotyping; saved models and sequences only).

  m10b1.py dump <name> <out>   the extended m10b_walk_dump (adds blocks / sources / segments); fingerprint must equal FROZEN_FP;
                               the four M10-B files must be byte-identical to M10-B's dump of the same model
  m10b1.py audit <out>         boundary-transformation audit of all six models: per-block relation, whole-template respelling,
                               raw per-block respelling, source kinds, case-scope support, degeneracy projections under v2
  m10b1.py cases <out>         replay of M10-B's frozen case set with project_v2, beside M10-B's version-1 verdicts
"""
import collections
import csv
import filecmp
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import m10b  # noqa: E402  (frozen M10-B definitions: models, cases, truth readers)
import project_v2 as p2  # noqa: E402

WT, R = m10b.WT, m10b.R
OUT = f"{R}/m10b1_projection_hardening_20261020"
M10B = m10b.OUT
FROZEN_FP = {"c4dev": "d2b8c6c2e0b4079d", "c4loo": "14851130f453a071", "cypdev": "2a5b720a1c5bf669", "cyploo": "e41f77fda080c143",
             "acotdev": "d622742f9cf8e7a5", "acotloo": "12cbe6beb1ad28d1"}
KV = m10b.KV


def write(out, name, rows):
    with open(f"{out}/{name}.tsv", "w") as f:
        if not rows:
            f.write("none\n")
            return
        w = csv.DictWriter(f, fieldnames=list(dict.fromkeys(k for r in rows for k in r)), delimiter="\t", lineterminator="\n", restval="")
        w.writeheader()
        w.writerows(rows)


def model(name):
    return p2.Model(f"{OUT}/runs/dump.{name}/{name}")


def dump(name, out):
    g, ref, cur, _ = m10b.MODELS[name]
    code = subprocess.run([f"{WT}/build/m10b_walk_dump", g, ref, cur, f"{out}/{name}"]).returncode
    if code:
        return code
    got = KV(f"{out}/{name}.summary.tsv")["model_fingerprint"]
    same = all(filecmp.cmp(f"{out}/{name}.{x}.tsv", f"{M10B}/runs/dump.{name}/{name}.{x}.tsv", shallow=False)
               for x in ("block_alleles", "walks", "nodes", "panel", "summary"))
    print(f"{name}: fingerprint {got} (frozen {FROZEN_FP[name]}); M10-B dump files identical: {same}")
    return 0 if got == FROZEN_FP[name] and same else 1


def audit(out):
    blocks, models, scopes, degen = [], [], [], []
    for name in m10b.MODELS:
        m = model(name)
        aud = m.audit()
        for b, r in aud.items():
            raw = sum(1 for a, s in m.alleles[b] for w in m.allele_walks(b, a)[1] if "".join(m.step(x) for x in w) == s)
            blocks.append(dict(model=name, raw_walk_spells_normalized=raw, **r))
        tmpl_ok = 0
        for t, al in m.panel.items():
            walk = [st for b, a in enumerate(al) for w in m.allele_walks(b, a)[1][:1] for st in w]
            tmpl_ok += "".join(m.step(x) for x in walk) == m.spell(0, al)
        kinds = collections.Counter()
        missing = 0
        for (b, a) in m.seq:
            src = m.sources.get((b, a), {})
            for k, v in src.items():
                kinds[k] += v
            missing += m.allele_walks(b, a)[0] == "MISSING"
        segs = [int(r["length"]) for r in p2.T(f"{OUT}/runs/dump.{name}/{name}.segments.tsv")]
        models.append(dict(model=name, blocks=len(m.block_ids), alleles=len(m.seq), blocks_verified=sum(r["verified"] for r in aud.values()),
                           blocks_with_trim=sum(1 for r in aud.values() if r["trimmed_prefix"] or r["trimmed_suffix"]),
                           templates_respelled_whole=f"{tmpl_ok}/{len(m.panel)}", nonempty_invariant_segments=sum(1 for x in segs if x),
                           alleles_missing_metadata=missing, **{k: v for k, v in sorted(kinds.items())}))
        locus = {"c4dev": "c4", "cypdev": "cyp", "acotdev": "acot"}.get(name)
        if locus:
            st = collections.Counter()
            for r in p2.T(m10b.DEGEN):
                if r["model"] != locus:
                    continue
                f, l = int(r["first_block_index"]), int(r["last_block_index"])
                ta = list(map(int, r["tuple_a"].split(",")))
                st[m.project(m.spell(f, ta), f, l)[0]] += 1
            degen.append(dict(model=name, **{k: v for k, v in sorted(st.items())}))
    for c in m10b.CASES:
        m = model(c["model"])
        f, l = m.block(c["interval"][0]), m.block(c["interval"][1])
        ok, cl, cr, why = m.scope(f, l)
        scopes.append(dict(case=c["id"], model=c["model"], interval="..".join(c["interval"]),
                           kind="whole haplotype" if f == 0 and l == len(m.block_ids) - 1 else "single block" if f == l else "compound interval",
                           supported=ok, clip_left_bp=cl, clip_right_bp=cr, reason=why))
    write(out, "block_audit", blocks)
    write(out, "model_audit", models)
    write(out, "case_scopes", scopes)
    write(out, "degeneracy_v2", degen)
    for n in ("model_audit", "case_scopes", "degeneracy_v2"):
        print(f"## {n}\n" + open(f"{out}/{n}.tsv").read())
    return 0 if all(r["blocks_verified"] == r["blocks"] for r in models) else 1


def cases(out):
    v1 = {(r["case"], r["a"], r["b"]): r for r in p2.T(f"{M10B}/runs/cases/comparisons.tsv")}
    rows, comps = [], []
    for c in m10b.CASES:
        m = model(c["model"])
        f, l = m.block(c["interval"][0]), m.block(c["interval"][1])
        seqs = {}
        for label, src in c["pairs"].items():
            if src[0] == "truth":
                seqs[label] = m10b.truth(src[1], src[2], m.block_ids[f:l + 1])
            else:
                tuples = m10b.labelling(src[1]) if src[0] == "call" else m10b.alts(src[1])[src[2]] if src[0] == "alts" else \
                    [v for k, v in sorted(m.panel.items()) if k.startswith(src[1])]
                seqs[label] = [m.spell(f, t[f:l + 1]) for t in tuples]
        proj = {lab: [m.project(h, f, l) for h in hs] for lab, hs in seqs.items()}
        for lab, ps in proj.items():
            for i, p in enumerate(ps):
                prof = [p2.profile(s) for s in sorted(p[1])] if p[0] in ("RESOLVED", "AMBIGUOUS") else []
                rows.append(dict(case=c["id"], haplotype=f"{lab}.h{i + 1}", status=p[0], parses=p[2],
                                 skeletons=len(p[1]) if prof else 0, max_multiplicity=";".join(str(x["max_multiplicity"]) for x in prof) or "-",
                                 multiplicity_signature=";".join(x["multiplicity_signature"] for x in prof) or "-",
                                 reason=p[1] if p[0] in p2.FAILING else p[3]))
        for a, b in c["compare"]:
            verdict, detail = p2.compare(proj[a], proj[b])
            old = v1.get((c["id"], a, b), {})
            comps.append(dict(case=c["id"], a=a, b=b, v2_category=verdict, v2_detail=detail[:200], v1_verdict=old.get("verdict", "NA"),
                              v1_detail=old.get("detail", "")[:80], exact_sequence_pair_equal=sorted(seqs[a]) == sorted(seqs[b]),
                              changed=("yes" if {"EQUAL": "UNIQUE_AGREEMENT", "EQUAL_SET_VALUED": "SET_IDENTICAL", "DISJOINT": "DISAGREEMENT",
                                                 "COMPATIBLE": "COMPATIBLE", "UNKNOWN": "UNSUPPORTED"}.get(old.get("verdict")) != verdict else "no")))
    write(out, "projections_v2", rows)
    write(out, "comparisons_v2", comps)
    print(open(f"{out}/comparisons_v2.tsv").read())
    return 0


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "dump":
        sys.exit(dump(sys.argv[2], sys.argv[3]))
    out = sys.argv[2]
    os.makedirs(out, exist_ok=True)
    sys.exit({"audit": audit, "cases": cases}[mode](out))
