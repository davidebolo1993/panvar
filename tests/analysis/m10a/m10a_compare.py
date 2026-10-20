#!/usr/bin/env python3
"""M10-A comparator (frozen with DESIGN.md). Every mode is an existing reviewed comparator plus EXPLICIT, named allowances;
nothing scientific is relaxed. Modes and what they establish:

  backend   REF vs O (5x)   p1_compare.py 'backend' (the P1 strict + qualified chain, unchanged) AND byte identity of every
                            sequence_first_v2 output (sequence_*, representation_intervals, sequence_compound), which P1's
                            SCIENCE list predates                                     -> computational-backend equivalence
  route     O vs D0         s2l_compare_v3.2, categories operational,profile,route,admission,checkpoints + 'incumbent'
                            (declared settings: the three retained-route keys)      -> retained-route equivalence
  jp        D0 vs D         s2l_compare_v3.2, categories operational,profile + 'junction' (declared setting
                            repair_junction_proposals). IDENTITY IS DEMANDED ONLY IF D reports junction_eligible_blocks 0;
                            otherwise the verdict is JP_ATTRIBUTABLE and the differences are reported, not failed
  prepared  D vs P          s2l_compare_v3.2, categories operational,profile + 'prepared'; no declared setting
                                                                                    -> persisted-index route equivalence
  repeat    D vs D (again)  s2l_compare_v3.2, categories operational,profile,checkpoints  -> determinism

Added categories (exact names, from saved S2-A off/on and P1 cold/loaded pairs):
  incumbent  one-sided early_incumbent* and checkpoint_* rows; the extra work counters of the early evaluation
  junction   one-sided junction_* rows (JP1.1 status)
  prepared   one-sided prepared_locus* rows and the model-preparation timing rows that a loaded run replaces

  m10a_compare.py MODE A B OUT.tsv
Exit: 0 PASS (or JP_ATTRIBUTABLE with valid contracts), 1 FAIL, 2 UNRESOLVED (an incomplete side).
"""
import importlib.util
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ANALYSIS = os.path.join(HERE, "..")
ROUTE_KEYS = "early_incumbent,union_fragment_blocks,repair_round1_origin_admission"
INCUMBENT_ONE_SIDED = re.compile(r"^(early_incumbent\w*|checkpoint_\w+)$")
INCUMBENT_WORK = {"combination_hypothetical_scalar_logs", "combination_hypothetical_scalar_validations", "combination_net_logs_avoided",
                  "combination_net_validations_avoided", "combination_prepared_calls", "combination_prepared_logs",
                  "combination_prepared_terms", "combination_prepared_validations", "emission_cache_executed_dps",
                  "emission_cache_executed_prefilters",
                  # the early pair is scored through the same indexed backend, so on the ordinary route its work adds here
                  "indexed_union_executed_dps", "indexed_union_executed_prefilters", "indexed_union_emission_cache_component_bytes",
                  "indexed_union_worker_scratch_bytes"}
JUNCTION_ONE_SIDED = re.compile(r"^junction_\w+$")
PREPARED_ONE_SIDED = re.compile(r"^(prepared_locus\w*|seconds_detail_model_\w+)$")
SEQUENCE_FILES = re.compile(r"^call\.(sequence_\w+\.(tsv|fa)|representation_intervals\.tsv)$")
SKIP_EXTRA = re.compile(r"^(args\.txt|driver\..*)$")


def load_v32(extra):
    spec = importlib.util.spec_from_file_location("s2l_v32", os.path.join(ANALYSIS, "s2l", "s2l_compare_v32.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    base = m.classify

    def classify(k, in_a, in_b, cats, settings_ok):
        got = base(k, in_a, in_b, cats, settings_ok)
        if got:
            return got
        one = in_a != in_b
        if "incumbent" in cats and ((one and INCUMBENT_ONE_SIDED.match(k)) or (not one and k in INCUMBENT_WORK)):
            return "incumbent"
        if "junction" in cats and one and JUNCTION_ONE_SIDED.match(k):
            return "junction"
        if "prepared" in cats and one and PREPARED_ONE_SIDED.match(k):
            return "prepared"
        return None
    m.classify = classify
    m.SCIENTIFIC_SKIP = re.compile(r"^(?:" + m.SCIENTIFIC_SKIP.pattern + "|" + SKIP_EXTRA.pattern + r")$")
    return m


def run_v32(a, b, out, cats, declared):
    m = load_v32(cats)
    argv = sys.argv
    sys.argv = ["s2l_compare_v32", a, b, out, "--categories", cats, "--declared-settings", declared]
    try:
        return m.main()
    finally:
        sys.argv = argv


def kv(p):
    return dict(l.split("\t", 1) for l in open(p).read().split("\n") if "\t" in l) if os.path.exists(p) else {}


def outcome(d):
    spec = importlib.util.spec_from_file_location("p1c", os.path.join(ANALYSIS, "p1", "p1_compare.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m.outcome(d), m.comparable(m.outcome(d))


def append(out, rows):
    with open(out, "a") as f:
        for k, v in rows:
            f.write(f"{k}\t{v}\n")


def main():
    mode, a, b, out = sys.argv[1], os.path.abspath(sys.argv[2]), os.path.abspath(sys.argv[3]), sys.argv[4]
    (oa, ca), (ob, cb) = outcome(a), outcome(b)
    if not (ca and cb):
        open(out, "w").write(f"key\tvalue\nmode\t{mode}\nverdict\tUNRESOLVED\noutcome_a\t{oa}\noutcome_b\t{ob}\n")
        print(open(out).read(), end="")
        return 2
    if mode == "backend":
        r = subprocess.run([sys.executable, os.path.join(ANALYSIS, "p1", "p1_compare.py"), "backend", a, b], capture_output=True, text=True)
        seq_diff = [f for f in sorted(set(os.listdir(a)) | set(os.listdir(b))) if SEQUENCE_FILES.match(f) and
                    (not os.path.exists(os.path.join(a, f)) or not os.path.exists(os.path.join(b, f)) or
                     open(os.path.join(a, f)).read().replace(a, "<D>") != open(os.path.join(b, f)).read().replace(b, "<D>"))]
        ok = r.returncode == 0 and not seq_diff
        open(out, "w").write(f"key\tvalue\nmode\tbackend\np1_compare_exit\t{r.returncode}\np1_compare_row\t{r.stdout.strip()}\n"
                             f"sequence_outputs_differing\t{','.join(seq_diff) or 'none'}\nverdict\t{'PASS' if ok else 'FAIL'}\n")
        print(open(out).read(), end="")
        return 0 if ok else 1
    if mode == "route":
        code = run_v32(a, b, out, "operational,profile,route,admission,checkpoints,incumbent", ROUTE_KEYS)
    elif mode == "prepared":
        code = run_v32(a, b, out, "operational,profile,prepared", "")
    elif mode == "repeat":
        code = run_v32(a, b, out, "operational,profile,checkpoints", "")
    elif mode == "jp":
        eligible = kv(os.path.join(b, "call.mosaic_status.tsv")).get("junction_eligible_blocks", "missing")
        code = run_v32(a, b, out, "operational,profile,junction", "repair_junction_proposals")
        if eligible == "0":
            append(out, [("junction_eligible_blocks", "0"), ("identity_demanded", "yes"), ("verdict", "PASS" if code == 0 else "FAIL")])
            print(f"jp: eligible 0, identity demanded: {'PASS' if code == 0 else 'FAIL'}")
            return code
        la, lb = kv(os.path.join(a, "call.sequence_locus.tsv")), kv(os.path.join(b, "call.sequence_locus.tsv"))
        rows = [("junction_eligible_blocks", eligible), ("identity_demanded", "no (JP1.1 changes candidate generation)"),
                ("identity_observed", "yes" if code == 0 else "no")]
        for k in ("prediction", "sequence_status", "representative_hap1_sha256", "representative_hap2_sha256"):
            rows.append((f"{k}_D0|D", f"{la.get(k, 'NA')}|{lb.get(k, 'NA')}"))
        contract = eligible != "missing" and all(os.path.exists(os.path.join(d, f"call.{f}")) for d in (a, b)
                                                 for f in ("sequence_locus.tsv", "sequence_blocks.tsv", "sequence_compound.tsv"))
        rows.append(("output_contract_valid", "yes" if contract else "NO"))
        rows.append(("verdict", "JP_ATTRIBUTABLE" if contract else "FAIL"))
        append(out, rows)
        print("\n".join(f"{k}\t{v}" for k, v in rows))
        return 0 if contract else 1
    else:
        raise SystemExit(f"unknown mode {mode}")
    append(out, [("verdict", "PASS" if code == 0 else "FAIL")])
    return code


if __name__ == "__main__":
    sys.exit(main())
