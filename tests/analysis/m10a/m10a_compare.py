#!/usr/bin/env python3
"""M10-A comparator (frozen with DESIGN.md). Every mode is an existing reviewed comparator plus EXPLICIT, named allowances;
nothing scientific is relaxed. Modes and what they establish:

  backend   REF vs O (5x)   p1_compare.py 'backend' (the P1 strict + qualified chain, unchanged) AND byte identity of every
                            sequence_first_v2 output (sequence_*, representation_intervals, sequence_compound), which P1's
                            SCIENCE list predates                                     -> computational-backend equivalence
  route     O vs D0         s2l_compare_v3.2, categories operational,profile,route,admission,checkpoints + 'incumbent'
                            (declared settings: the three retained-route keys)      -> retained-route equivalence
  jp        D0 vs D         Amendment 1 (prospective; Sol's reproducer): invariants that JP eligibility NEVER waives, then
                            identity or attribution. Always, in this order (any violation -> FAIL):
                              J1 junction_eligible_blocks in D is a nonnegative integer; D0 has no junction_* rows and JP
                                 off; D has JP experimental_v1
                              J2 required outputs present and well-formed in both arms (REQUIRED, plus sequence_compound.tsv
                                 when compound_annotation is contract_v2); sequence_locus prediction in PREDICTIONS
                              J3 effective settings: equal key sets and values except repair_junction_proposals
                              J4 protected status keys equal: input/model/library identity (IDENTITY_STATUS) and configured
                                 caps and budgets (CONFIG_STATUS)
                            eligible 0: full s2l_compare_v3.2 identity (categories operational,profile,junction)
                            eligible > 0: J5 every hypothesis scored in both arms -- matched by its labelling (the unordered
                                 pair of allele tuples), never by candidate or group id -- has EXACTLY equal printed objective
                                 components (OBJECTIVE, and the labelling's exact_path_prior); at least one common
                                 hypothesis. candidate_set_probability is excluded (it depends on pool membership).
                                 Only then JP_ATTRIBUTABLE; v3.2 identity is run and reported (not failed)
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
# Amendment 1 (JP invariants)
REQUIRED = ("settings.tsv", "mosaic_status.tsv", "mosaic_scores.tsv", "mosaic_candidates.tsv", "repair_rounds.tsv", "sequence_locus.tsv",
            "sequence_blocks.tsv", "sequence_alternatives.tsv", "sequence_haplotypes.tsv")
PREDICTIONS = {"PREDICTION", "PROVISIONAL", "NO_PREDICTION"}
IDENTITY_STATUS = re.compile(r"^(model_fingerprint|sha256_\w+|gfa|reference_path|locus_geometry_\w+|transition_vector_sha256|"
                             r"transition_boundaries|transition_zero_distance_\w+|library_provenance|curation_digest|"
                             r"boundary_bases_moved_to_context)$")
CONFIG_STATUS = re.compile(r"^(max_\w+|\w+_limit_\w+|\w+_budget_(bytes|basis|scope)|\w+_budgeted_bytes|proposal_state_limit|"
                           r"representation_guard_bound_nats|union_cap_\w+|union_caps_scope|threads_requested)$")
# Amendment 3: the plan declares different per-arm timeouts for D0 and D (ACOT 5x, every 30x cell); the wall-clock watchdog
# never changes a completed result (compose_profile OPERATIONAL_OVERRIDES), and both sides must have completed (outcome check)
JP_GUARDS = {"time_limit_seconds"}
OBJECTIVE = ("score", "fragment_term", "exposure_term", "sequence_pair_prior", "numeric_error_bound", "fragments_with_mass",
             "fragments_total", "mismatch_band")


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
        # Amendment 2: the annotation's own timer is timing (inspected: 0.00014-0.00016 s)
        if "operational" in cats and not one and k == "compound_seconds":
            return "operational"
        # Amendment 2: the status echo of the declared JP setting (its value is checked in the settings table)
        if "junction" in cats and one and k == "repair_junction_proposals":
            return "junction"
        if "junction" in cats and not one and k in JP_GUARDS:   # Amendment 3 (jp mode only)
            return "guard"
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


def table(p):
    import csv
    return list(csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t"))


def settings(d):
    return {k: v.split("\t")[0] for k, v in kv(os.path.join(d, "call.settings.tsv")).items() if k != "key"}


def jp_invariants(a, b):
    """J1-J4: (problems, eligible). Problems are 'J<n> reason' strings; an empty list means every invariant holds."""
    bad = []
    sta, stb = kv(os.path.join(a, "call.mosaic_status.tsv")), kv(os.path.join(b, "call.mosaic_status.tsv"))
    sa, sb = settings(a), settings(b)
    eligible = stb.get("junction_eligible_blocks", "missing")
    if not re.fullmatch(r"[0-9]+", eligible):
        bad.append(f"J1 junction_eligible_blocks is not a nonnegative integer: {eligible!r}")
    if any(k.startswith("junction_") for k in sta) or sa.get("repair_junction_proposals", "off") != "off":
        bad.append("J1 the JP-off arm carries JP rows or settings")
    if sb.get("repair_junction_proposals") != "experimental_v1":
        bad.append("J1 the JP-on arm does not run repair_junction_proposals experimental_v1")
    for d, st in ((a, sa), (b, sb)):
        need = list(REQUIRED) + (["sequence_compound.tsv"] if st.get("compound_annotation") == "contract_v2" else [])
        for f in need:
            p = os.path.join(d, "call." + f)
            if not os.path.exists(p) or os.path.getsize(p) == 0:
                bad.append(f"J2 missing or empty required output {os.path.basename(d)}/call.{f}")
            elif f not in ("settings.tsv", "mosaic_status.tsv", "sequence_locus.tsv") and not any(
                    "\t" in l and not l.startswith("#") for l in open(p, errors="replace").read().split("\n")[:50]):
                bad.append(f"J2 malformed required output {os.path.basename(d)}/call.{f} (no TSV header)")
        pred = kv(os.path.join(d, "call.sequence_locus.tsv")).get("prediction")
        if os.path.exists(os.path.join(d, "call.sequence_locus.tsv")) and pred not in PREDICTIONS:
            bad.append(f"J2 sequence_locus prediction {pred!r} in {os.path.basename(d)}")
    for k in sorted(set(sa) | set(sb)):
        if k != "repair_junction_proposals" and k not in JP_GUARDS and sa.get(k) != sb.get(k):
            bad.append(f"J3 undeclared effective setting {k}: {sa.get(k)} vs {sb.get(k)}")
    for k in sorted(set(sta) | set(stb)):
        if (IDENTITY_STATUS.match(k) or CONFIG_STATUS.match(k)) and k not in JP_GUARDS and sta.get(k) != stb.get(k):
            kind = "identity" if IDENTITY_STATUS.match(k) else "configured cap/budget"
            bad.append(f"J4 {kind} status key {k}: {sta.get(k)} vs {stb.get(k)}")
    return bad, eligible


def jp_objective(a, b):
    """J5: (problems, common hypotheses). Matched by labelling, never by id; exact printed equality."""
    def load(d):
        scores = {r["sequence_candidate_id"]: r for r in table(os.path.join(d, "call.mosaic_scores.tsv"))}
        labs = {}
        for r in table(os.path.join(d, "call.mosaic_candidates.tsv")):
            labs.setdefault(tuple(sorted((r["hap1_alleles"], r["hap2_alleles"]))), []).append(r)
        return scores, labs
    (sa, la), (sb, lb) = load(a), load(b)
    common = set(la) & set(lb)
    bad = []
    for k in sorted(common):
        prior_a = {r["exact_path_prior"] for r in la[k]}
        prior_b = {r["exact_path_prior"] for r in lb[k]}
        if prior_a != prior_b:
            bad.append(f"J5 exact_path_prior of labelling {k[0][:40]}|{k[1][:40]}: {sorted(prior_a)} vs {sorted(prior_b)}")
        rows_a = [sa[r["sequence_candidate_id"]] for r in la[k] if r["sequence_candidate_id"] in sa]
        rows_b = [sb[r["sequence_candidate_id"]] for r in lb[k] if r["sequence_candidate_id"] in sb]
        if not rows_a or not rows_b:
            bad.append(f"J5 labelling {k[0][:40]}|{k[1][:40]} has no score row in one arm")
            continue
        va = {tuple(r[c] for c in OBJECTIVE) for r in rows_a}
        vb = {tuple(r[c] for c in OBJECTIVE) for r in rows_b}
        if va != vb or len(va) != 1:
            bad.append(f"J5 objective of labelling {k[0][:40]}|{k[1][:40]}: {sorted(va)[:1]} vs {sorted(vb)[:1]}")
        if len(bad) > 20:
            bad.append("J5 ... (truncated after 20)")
            break
    if not common:
        bad.append("J5 no hypothesis is scored in both arms: the objective check would be vacuous")
    return bad, len(common)


def jp(a, b, out):
    bad, eligible = jp_invariants(a, b)
    tl = (kv(os.path.join(a, "call.mosaic_status.tsv")).get("time_limit_seconds"), kv(os.path.join(b, "call.mosaic_status.tsv")).get("time_limit_seconds"))
    rows = [("mode", "jp"), ("junction_eligible_blocks", eligible), ("time_limit_seconds_D0|D", f"{tl[0]}|{tl[1]}"),
            ("invariants_J1_J4", "pass" if not bad else "FAIL")]
    if not bad and eligible == "0":
        code = run_v32(a, b, out, "operational,profile,junction,checkpoints", "repair_junction_proposals,time_limit_seconds")   # Amendment 2: checkpoints
        append(out, rows + [("identity_demanded", "yes (eligible 0)"), ("verdict", "PASS" if code == 0 else "FAIL")])
        print(open(out).read(), end="")
        return code
    if not bad:
        obj, common = jp_objective(a, b)
        rows += [("common_hypotheses", common), ("objective_J5", "pass" if not obj else "FAIL")]
        bad += obj
    if bad:
        open(out, "w").write("key\tvalue\n")
        append(out, rows + [("problem", p) for p in bad] + [("verdict", "FAIL")])
        print(open(out).read(), end="")
        return 1
    code = run_v32(a, b, out, "operational,profile,junction,checkpoints", "repair_junction_proposals,time_limit_seconds")
    la, lb = kv(os.path.join(a, "call.sequence_locus.tsv")), kv(os.path.join(b, "call.sequence_locus.tsv"))
    rows += [("identity_demanded", "no (JP1.1 changes candidate generation; invariants J1-J5 passed)"),
             ("identity_observed", "yes" if code == 0 else "no")]
    for k in ("prediction", "sequence_status", "representative_hap1_sha256", "representative_hap2_sha256"):
        rows.append((f"{k}_D0|D", f"{la.get(k, 'NA')}|{lb.get(k, 'NA')}"))
    rows.append(("verdict", "JP_ATTRIBUTABLE"))
    append(out, rows)
    print(open(out).read(), end="")
    return 0


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
        code = run_v32(a, b, out, "operational,profile,prepared,checkpoints", "")   # Amendment 2 (prospective): checkpoints
    elif mode == "repeat":
        code = run_v32(a, b, out, "operational,profile,checkpoints", "")
    elif mode == "jp":
        return jp(a, b, out)
    else:
        raise SystemExit(f"unknown mode {mode}")
    append(out, [("verdict", "PASS" if code == 0 else "FAIL")])
    return code


if __name__ == "__main__":
    sys.exit(main())
