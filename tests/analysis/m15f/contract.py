#!/usr/bin/env python3
"""M15-F: the simplified filter contract, as an executable specification (read-only; no production code).

  current(facts)  -- a faithful port of src/genotype_headline.cpp (headline_filters_v2), verified against the values
                     recorded in saved runs (`fidelity`)
  proposed(facts) -- the proposed primary view (CONTRACT.md section 3)
  facts(prefix)   -- the facts both read, from a saved run's own outputs (never modified)

Usage: contract.py fidelity|synthetic|saved <out dir>"""
import csv
import hashlib
import os
import sys

csv.field_size_limit(sys.maxsize)
R = "/Users/davide.bolognini/panvar/tests/results"
SEV = {"NO_CALL": 4, "NOT_EVALUATED": 3, "LOW_QUALITY": 2, "PASS": 1, "NOT_APPLICABLE": 0}


# ---------------------------------------------------------------- facts ------------------------------------------------
def kv(p):
    if not os.path.exists(p):
        return {}
    out = {}
    for line in open(p):
        if "\t" in line:
            k, v = line.rstrip("\n").split("\t", 1)
            out.setdefault(k, v)
    return out


def rows(p):
    return list(csv.DictReader(open(p), delimiter="\t")) if os.path.exists(p) else []


def facts(prefix):
    """Everything the current headline and the proposed view read, from a run prefix. Mirrors apply_headline's reading
    (call class derived from status when absent; declared library -> NOT_APPLICABLE; block facts from the block table,
    dosage facts falling back to .dosage_blocks.tsv) and adds the v2 sequence report and genotype_support."""
    st, sl = kv(prefix + ".mosaic_status.tsv"), kv(prefix + ".sequence_locus.tsv")
    locus_row = (rows(prefix + ".locus_genotype.tsv") or [{}])[0]
    run_status = st.get("status", "NA")
    call_class = st.get("call_class") or ("strict_call" if run_status == "COMPLETE" else
                                          "genotype_equivalent_call" if run_status == "CONDITIONAL_GENOTYPE_EQUIVALENT" else "no_call")
    no_call_reason = st.get("no_call_reason", "NA")
    if call_class == "no_call" and not st.get("no_call_reason") and run_status != "NA":
        no_call_reason = "status_" + run_status
    if st.get("library_provenance", "").startswith("declared"):
        lib = "NOT_APPLICABLE"
    else:
        lib = st.get("library_sensitivity_status", st.get("library_sensitivity", "NOT_EVALUATED"))
    ids = lambda key: set() if st.get(key, "none") in ("none", "") else set(st[key].split(","))
    dosage_mode = locus_row.get("dosage_mode", sl.get("dosage_mode", "off"))
    if dosage_mode in ("NA", "", None, "not_run"):
        dosage_mode = "off"
    ds = kv(prefix + ".dosage_status.tsv")
    audit = locus_row.get("dosage_audit_status") or ds.get("dosage_audit_status", "NA")
    f = dict(prefix=prefix, call_class=call_class, no_call_reason=no_call_reason, phase_support=st.get("phase_support", "NA"),
             library_sensitivity=lib, dosage_mode=dosage_mode, dosage_audit_status=audit if dosage_mode != "off" else "NA",
             dosage_audit_note=ds.get("dosage_audit_note", "none"), library_incomplete=st.get("library_sensitivity_incomplete", ""),
             prediction=sl.get("prediction", "ABSENT"), completion=sl.get("completion", "NA"),
             sequence_status=sl.get("sequence_status", "NA"), schema=sl.get("schema_version", "legacy"), blocks=[])
    seq_blocks = {r["block_id"]: r for r in rows(prefix + ".sequence_blocks.tsv")}
    f["seq_blocks"] = seq_blocks
    called = call_class in ("strict_call", "genotype_equivalent_call")
    dosage_rows = {r["block_id"]: r for r in rows(prefix + ".dosage_blocks.tsv")}
    sens, unaud = ids("library_sensitive_blocks"), ids("library_sensitivity_not_evaluated_blocks")
    table = rows(prefix + ".block_genotypes.tsv") if called else []
    for r in table:
        d = dosage_rows.get(r["block_id"], {})
        pick = lambda k: r.get(k) or d.get(k, "NA")
        b = dict(block_id=r["block_id"], zygosity=r.get("zygosity", "NA"), representation_filter=r.get("filter", "NA"),
                 library_sensitive=r["block_id"] in sens, library_not_evaluated=r["block_id"] in unaud,
                 dosage_status=pick("dosage_status"), repair_class=pick("repair_class"), dosage_content_scope=pick("dosage_content_scope"),
                 dosage_reason=pick("dosage_reason"), dosage_repair_applied=r.get("dosage_repair_applied", "NA") == "true",
                 repair_not_applied_reason="not_recorded", genotype_support=r.get("genotype_support", "NA"),
                 recorded=(r.get("genotype_filter"), r.get("genotype_filter_reasons"), r.get("dosage_filter"), r.get("dosage_filter_reasons")))
        # the audit-mode not-applied record is not persisted (only through dosage_filter_reasons): recovered from there
        rec = r.get("dosage_filter_reasons", "")
        if "COMPLETE_REPAIR_AVAILABLE_NOT_APPLIED" in rec:
            b["repair_not_applied_reason"] = "dosage_mode_audit"
        elif rec.startswith("COMPLETE_REPAIR_NOT_APPLIED:"):
            b["repair_not_applied_reason"] = rec.split(":", 1)[1]
        elif r.get("dosage_stage_filter") == "DOSAGE_REPAIR_NOT_APPLIED":
            b["repair_not_applied_reason"] = r.get("dosage_stage_filter_reason", "not_recorded")
        f["blocks"].append(b)
    f["block_ids"] = [b["block_id"] for b in f["blocks"]] or list(seq_blocks) or ["locus"]
    f["recorded"] = {k: locus_row.get(k) for k in ("genotype_filter", "genotype_filter_reasons", "dosage_filter",
                                                   "dosage_filter_reasons", "phase_filter", "phase_filter_reasons", "call_summary")}
    return f


# ---------------------------------------------------------------- current (port of headline_filters_v2) ---------------
def v(status, reason=None):
    return (status, [reason] if reason else [])


def combine(parts):
    if not parts:
        return ("PASS", [])
    s, rs = parts[0][0], []
    for p in parts:
        if SEV[p[0]] > SEV[s]:
            s = p[0]
        for r in p[1]:
            if r not in rs:
                rs.append(r)
    return (s, rs)


def no_call(L):
    return L["call_class"] not in ("strict_call", "genotype_equivalent_call")


def dosage_ran(L):
    return L["dosage_mode"] != "off" and L["dosage_audit_status"] == "COMPLETE"


def cur_genotype(L, b):
    if no_call(L):
        return v("NO_CALL", "RUN_NO_CALL:" + L["no_call_reason"])
    parts, rep = [], b["representation_filter"]
    if rep != "PASS" and not rep.startswith("LOWQUAL"):
        parts.append(v("NOT_EVALUATED", "REPRESENTATION_NOT_EVALUATED"))
    if L["library_sensitivity"] == "NOT_EVALUATED" or b["library_not_evaluated"]:
        parts.append(v("NOT_EVALUATED", "LIBRARY_SENSITIVITY_NOT_EVALUATED"))
    if rep == "LOWQUAL_REPRESENTATION_CONFLICT":
        parts.append(v("LOW_QUALITY", "REPRESENTATION_CONFLICT"))
    if L["library_sensitivity"] == "SENSITIVE" and b["library_sensitive"]:
        parts.append(v("LOW_QUALITY", "LIBRARY_SENSITIVE"))
    if dosage_ran(L) and (b["dosage_status"] in ("CONTENT_REPAIRED", "TRIGGERED_UNREPAIRED") or
                          (b["dosage_status"] == "REPAIRED" and not b["dosage_repair_applied"])):
        parts.append(v("LOW_QUALITY", "EXACT_DIPLOTYPE_UNRESOLVED"))
    return combine(parts)


def cur_dosage(L, b):
    if no_call(L):
        return v("NO_CALL", "RUN_NO_CALL")
    if L["dosage_mode"] == "off":
        return v("NOT_EVALUATED", "DOSAGE_NOT_EVALUATED:dosage_mode_off")
    if L["dosage_audit_status"] != "COMPLETE":
        return v("NOT_EVALUATED", "DOSAGE_NOT_EVALUATED:audit_" + L["dosage_audit_status"])
    s = b["dosage_status"]
    if s == "CONSISTENT":
        return v("PASS", "DOSAGE_CONSISTENT")
    if s == "REPAIRED":
        complete = b["repair_class"] == "COMPLETE_DIPLOTYPE_REPAIR"
        if complete and b["dosage_repair_applied"]:
            return v("PASS", "COMPLETE_REPAIR_APPLIED")
        if complete and L["dosage_mode"] == "audit" and b["repair_not_applied_reason"] == "dosage_mode_audit":
            return v("PASS", "COMPLETE_REPAIR_AVAILABLE_NOT_APPLIED")
        return v("LOW_QUALITY", "COMPLETE_REPAIR_NOT_APPLIED:" + b["repair_not_applied_reason"])
    if s == "CONTENT_REPAIRED" and b["dosage_content_scope"] == "COMPLETE":
        return v("PASS", "CONTENT_CORRECTED_EXACT_DIPLOTYPE_UNRESOLVED")
    if s == "CONTENT_REPAIRED":
        return v("LOW_QUALITY", "CONTENT_CORRECTION_CONDITIONAL")
    if s == "SUSPECT":
        return v("LOW_QUALITY", "DOSAGE_SUSPECT")
    if s == "TRIGGERED_UNREPAIRED":
        return v("LOW_QUALITY", "TRIGGERED_UNREPAIRED")
    if s == "UNDERPOWERED":
        return v("NOT_EVALUATED", "DOSAGE_UNDERPOWERED")
    if s == "NOT_SENSITIVE":
        return v("NOT_APPLICABLE", "NOT_DOSAGE_SENSITIVE")
    return v("NOT_EVALUATED", "DOSAGE_NOT_EVALUATED")


def cur_phase(L, a, b):
    if no_call(L):
        return v("NO_CALL", "RUN_NO_CALL")
    if a["zygosity"] != "HET" or b["zygosity"] != "HET":
        return v("NOT_APPLICABLE", "NO_PHASE_INFORMATIVE_BOUNDARY")
    ps = L["phase_support"]
    if L["call_class"] == "genotype_equivalent_call" or ps == "UNRESOLVED":
        base = v("NO_CALL", "PHASE_UNRESOLVED")
    elif ps in ("READS_AUDITED", "READS_WITHIN_CANDIDATES"):
        base = v("PASS")
    elif ps == "PRIOR_ONLY":
        base = v("LOW_QUALITY", "MODEL_ONLY_PHASE")
    else:
        base = v("NOT_EVALUATED", "PHASE_ASSESSMENT_INCOMPLETE")
    parts = [base]
    if a["dosage_repair_applied"] or b["dosage_repair_applied"]:
        parts.append(v("NOT_EVALUATED", "PHASE_NOT_ASSESSED_AFTER_DOSAGE_REPAIR"))
    return combine(parts)


def aggregate(units):
    counts = {s: 0 for s in SEV}
    status, any_app = "NOT_APPLICABLE", False
    for _, (s, _) in units:
        counts[s] += 1
        if s == "NOT_APPLICABLE":
            continue
        if not any_app or SEV[s] > SEV[status]:
            status = s
        any_app = True
    by = []
    for uid, (s, rs) in units:
        if s == "NOT_APPLICABLE" and any_app:
            continue
        for r in rs:
            e = next((x for x in by if x[0] == r), None)
            if e is None:
                e = (r, [])
                by.append(e)
            e[1].append(uid)
    reasons = ";".join(f"{r}@{','.join(i)}" for r, i in by) or "none"
    return status, reasons


NONE_BLOCK = dict(zygosity="NA", representation_filter="NA", library_sensitive=False, library_not_evaluated=False,
                  dosage_status="NA", repair_class="NA", dosage_content_scope="NA", dosage_repair_applied=False,
                  repair_not_applied_reason="not_recorded", genotype_support="NA", dosage_reason="NA")


def current(L):
    blocks = L["blocks"]
    if not blocks:
        ids = [i for i in L["block_ids"] if i != "locus"]
        g = [(i, cur_genotype(L, NONE_BLOCK)) for i in ids]
        d = [(i, cur_dosage(L, NONE_BLOCK)) for i in ids]
        p = [(f"{ids[k]}>{ids[k + 1]}", cur_phase(L, NONE_BLOCK, NONE_BLOCK)) for k in range(len(ids) - 1)]
        if no_call(L):
            g = g or [("locus", cur_genotype(L, NONE_BLOCK))]
            d = d or [("locus", cur_dosage(L, NONE_BLOCK))]
            p = p or [("locus", cur_phase(L, NONE_BLOCK, NONE_BLOCK))]
    else:
        g = [(b["block_id"], cur_genotype(L, b)) for b in blocks]
        d = [(b["block_id"], cur_dosage(L, b)) for b in blocks]
        p = [(f"{blocks[k]['block_id']}>{blocks[k + 1]['block_id']}", cur_phase(L, blocks[k], blocks[k + 1])) for k in range(len(blocks) - 1)]
    out = {}
    for name, units in (("genotype", g), ("dosage", d), ("phase", p)):
        out[name + "_filter"], out[name + "_filter_reasons"] = aggregate(units)
    out["call_summary"] = f"GENOTYPE_{out['genotype_filter']};DOSAGE_{out['dosage_filter']};PHASE_{out['phase_filter']}"
    return out


# ---------------------------------------------------------------- proposed (assessment_v1) ----------------------------
def prop_prediction(L):
    """The primary outcome is the v2 sequence report's; the legacy classification only where no v2 report exists."""
    if L["prediction"] in ("PREDICTION", "PROVISIONAL", "NO_PREDICTION"):
        return L["prediction"], L["completion"]
    if L["call_class"] in ("strict_call", "genotype_equivalent_call"):
        return "PREDICTION", "legacy_" + L["call_class"]
    return "NO_PREDICTION", L["no_call_reason"]


def prop_block_genotype(L, b):
    """Detected genotype conflicts on the REPORTED block (dosage is assessed separately, never folded in here)."""
    found, missing = [], []
    rep = b["representation_filter"]
    if rep == "LOWQUAL_REPRESENTATION_CONFLICT":
        found.append("REPRESENTATION_CONFLICT")
    elif rep != "PASS":
        missing.append("REPRESENTATION_GUARD:" + (rep or "NA"))
    if L["library_sensitivity"] == "SENSITIVE" and b["library_sensitive"]:
        found.append("LIBRARY_SENSITIVE")
    elif L["library_sensitivity"] == "NOT_EVALUATED" or b["library_not_evaluated"]:
        missing.append("LIBRARY_SENSITIVITY:" + (L.get("library_incomplete") or "not_evaluated"))
    value = "CONFLICT" if found else ("NONE_DETECTED" if rep in ("PASS", "LOWQUAL_REPRESENTATION_CONFLICT") else "NOT_ASSESSED")
    return value, found, missing


def prop_block_dosage(L, b):
    """(dosage_check, reason, dosage_alternative): the check concerns the REPORTED alleles; any correction is an
    alternative reported separately."""
    mode, audit = L["dosage_mode"], L["dosage_audit_status"]
    if mode == "off":
        return "NOT_REQUESTED", "dosage_check_not_requested", "NA"
    if audit != "COMPLETE":
        why = {"FAILED": "AUDIT_FAILED", "NOT_RUN": "AUDIT_NOT_RUN", "REFUSED_LIBRARY": "AUDIT_UNSUPPORTED_LIBRARY"}.get(audit, "AUDIT_" + audit)
        return "UNAVAILABLE", why, "NA"
    s = b["dosage_status"]
    if s == "NOT_SENSITIVE":
        return "NOT_APPLICABLE", "NOT_DOSAGE_SENSITIVE", "NONE"
    if s == "CONSISTENT":
        return "AGREES", "CONSISTENT", "NONE"
    if s == "SUSPECT":
        return "SUSPECT", "SUSPECT", "NONE"
    if s == "UNDERPOWERED":
        return "UNAVAILABLE", "UNDERPOWERED", "NONE"
    if s == "REPAIRED" and b["dosage_repair_applied"]:
        return "AGREES", "REPORTED_IS_DOSAGE_REPAIR", "APPLIED"   # legacy --dosage-mode repair only
    if s == "REPAIRED":
        alt = ("COMPLETE_REPAIR_AVAILABLE" if b["repair_not_applied_reason"] == "dosage_mode_audit"
               else "COMPLETE_REPAIR_REFUSED:" + b["repair_not_applied_reason"])
        return "DISAGREES", "TRIGGERED", alt
    if s == "CONTENT_REPAIRED":
        return "DISAGREES", "TRIGGERED", "CONTENT_CORRECTION_AVAILABLE"
    if s in ("TRIGGERED_UNREPAIRED", "AUTO_REPAIR"):
        return "DISAGREES", "TRIGGERED", "REPAIR_NOT_FOUND:" + (b.get("dosage_reason") or "none")
    return "UNAVAILABLE", "NOT_EVALUATED:" + (b.get("dosage_reason") or "NA"), "NONE"


def prop_phase(L, pred, het_boundaries):
    """Locus-level: phase support is computed once per locus (decision-level), never per boundary."""
    if pred == "NO_PREDICTION":
        return "NA", "no_prediction"
    if pred == "PROVISIONAL":
        return "NOT_ASSESSED", "provisional_checkpoint"
    if L["sequence_status"] in ("GENOTYPE_RESOLVED_SEQUENCE_SET", "SEQUENCE_TIED") or L["call_class"] == "genotype_equivalent_call":
        return "UNRESOLVED", "several_sequence_pairs"
    if het_boundaries == 0:
        return "NOT_APPLICABLE", "no_phase_informative_boundary"
    ps = L["phase_support"]
    if ps in ("READS_AUDITED", "READS_WITHIN_CANDIDATES"):
        return "READ_SUPPORTED", ps.lower()
    if ps == "PRIOR_ONLY":
        return "PRIOR_ONLY", "prior_only"
    if ps == "UNRESOLVED" and L["sequence_status"] == "RESOLVED":
        # one sequence pair (phase fixed by the haplotypes) but its read support was not computed: the legacy decision
        # tied over labellings of that same pair
        return "NOT_ASSESSED", "label_tie_support_not_computed"
    return "NOT_ASSESSED", "phase_support_" + ps


def proposed(L):
    pred, detail = prop_prediction(L)
    out = {"prediction": pred, "prediction_detail": detail, "sequence_resolution": L["sequence_status"] if pred != "NO_PREDICTION" else "NA"}
    findings, missing = [], []

    def note(lst, code, bid):
        e = next((x for x in lst if x[0] == code), None)
        if e is None:
            e = (code, [])
            lst.append(e)
        if bid not in e[1]:
            e[1].append(bid)

    gvals, dvals, svals, alts = [], [], [], []
    for b in (L["blocks"] if pred != "NO_PREDICTION" else []):
        bid = b["block_id"]
        gv, gf, gm = prop_block_genotype(L, b)
        gvals.append(gv)
        for c in gf:
            note(findings, c, bid)
        for c in gm:
            note(missing, c, bid)
        dv, dr, alt = prop_block_dosage(L, b)
        dvals.append(dv)
        if dv in ("SUSPECT", "DISAGREES"):
            note(findings, "DOSAGE_" + dv, bid)
        if dv == "UNAVAILABLE":
            if dr.startswith("AUDIT_"):   # a run-level failure: recorded once, at the locus
                note(missing, "DOSAGE:" + dr + (":" + L["dosage_audit_note"] if dr == "AUDIT_FAILED" else ""), "locus")
            else:
                note(missing, "DOSAGE:" + dr, bid)
        if alt not in ("NA", "NONE"):
            alts.append(f"{alt}@{bid}")
        svals.append(b["genotype_support"])
        if b["genotype_support"] == "PRIOR_ONLY":
            note(findings, "PRIOR_DEPENDENT", bid)
    for bid, sb in (L["seq_blocks"].items() if pred == "PREDICTION" else []):
        if sb.get("block_status") == "REPRESENTATION_AMBIGUOUS":
            note(findings, "LABELS_AMBIGUOUS", bid)
        if sb.get("block_status") == "GENOTYPE_UNRESOLVED":
            note(findings, "GENOTYPE_UNRESOLVED", bid)
    if pred == "NO_PREDICTION":
        out.update(genotype_check="NA", dosage_check="NA", read_support="NA")
    elif not L["blocks"]:
        # a v2 prediction whose legacy block tables were not written (legacy no-call beside a v2 PREDICTION/PROVISIONAL)
        out.update(genotype_check="NOT_ASSESSED", read_support="NOT_ASSESSED",
                   dosage_check="NOT_REQUESTED" if L["dosage_mode"] == "off" else "UNAVAILABLE")
        note(missing, "BLOCK_CHECKS:" + ("provisional_checkpoint" if pred == "PROVISIONAL" else "legacy_decision_" + L["no_call_reason"]), "locus")
        if L["dosage_mode"] != "off":
            note(missing, "DOSAGE:" + {"NOT_RUN": "AUDIT_NOT_RUN", "FAILED": "AUDIT_FAILED"}.get(L["dosage_audit_status"], "AUDIT_" + L["dosage_audit_status"]), "locus")
    else:
        out["genotype_check"] = "CONFLICT" if "CONFLICT" in gvals else "NONE_DETECTED" if "NONE_DETECTED" in gvals else "NOT_ASSESSED"
        order = ["DISAGREES", "SUSPECT", "AGREES", "UNAVAILABLE", "NOT_REQUESTED", "NOT_APPLICABLE"]
        out["dosage_check"] = next((x for x in order if x in dvals), "NA")
        out["read_support"] = "PRIOR_DEPENDENT" if "PRIOR_ONLY" in svals else "READS" if "READS" in svals else "NOT_ASSESSED"
    het = sum(1 for k in range(len(L["blocks"]) - 1) if L["blocks"][k]["zygosity"] == "HET" and L["blocks"][k + 1]["zygosity"] == "HET")
    if not L["blocks"] and pred != "NO_PREDICTION":
        het = 1   # unknown zygosity: do not claim NOT_APPLICABLE
    out["phase"], out["phase_reason"] = prop_phase(L, pred, het)
    if any(b["dosage_repair_applied"] for b in L["blocks"]) and out["phase"] not in ("NA", "NOT_APPLICABLE"):
        out["phase"], out["phase_reason"] = "NOT_ASSESSED", "after_dosage_repair"   # legacy --dosage-mode repair only
    if out["phase"] in ("PRIOR_ONLY", "UNRESOLVED"):
        note(findings, "PHASE_" + out["phase"], "locus")
    if out["phase"] == "NOT_ASSESSED":
        note(missing, "PHASE:" + out["phase_reason"], "locus")
    fmt = lambda lst: ";".join(f"{c}@{','.join(i)}" for c, i in lst) or "none"
    out["findings"], out["not_assessed"] = fmt(findings), fmt(missing)
    out["dosage_alternatives"] = ";".join(alts) or "none"
    out["summary"] = (f"{pred};SEQUENCE_{out['sequence_resolution']};GENOTYPE_{out['genotype_check']};DOSAGE_{out['dosage_check']};"
                      f"SUPPORT_{out['read_support']};PHASE_{out['phase']}")
    return out


# ---------------------------------------------------------------- units ------------------------------------------------
CURRENT_KEYS = ["genotype_filter", "genotype_filter_reasons", "dosage_filter", "dosage_filter_reasons", "phase_filter",
                "phase_filter_reasons", "call_summary"]
PROPOSED_KEYS = ["prediction", "prediction_detail", "sequence_resolution", "genotype_check", "dosage_check", "read_support",
                 "phase", "findings", "not_assessed", "dosage_alternatives", "summary"]


def saved_prefixes():
    out = []
    for root, dirs, files in os.walk(R):
        if "m15b_reserved" in root or "/index" in root:
            dirs[:] = []
            continue
        for f in files:
            if f.endswith(".locus_genotype.tsv"):
                out.append(os.path.join(root, f[: -len(".locus_genotype.tsv")]))
    return sorted(out)


def unit_fidelity(d):
    """The Python port reproduces every recorded headline value (locus filters, reasons, call summary)."""
    ok = bad = skipped = 0
    lines = ["prefix\tresult\tfirst_difference"]
    for p in saved_prefixes():
        L = facts(p)
        if not L["recorded"].get("genotype_filter"):
            skipped += 1
            continue
        c = current(L)
        diff = next((k for k in CURRENT_KEYS if c[k] != L["recorded"].get(k)), None)
        ok += diff is None
        bad += diff is not None
        lines.append(f"{p.replace(R + '/', '')}\t{'match' if diff is None else 'DIFF'}\t"
                     f"{'' if diff is None else diff + ': port=' + c[diff][:120] + ' recorded=' + str(L['recorded'].get(diff))[:120]}")
    open(f"{d}/fidelity.tsv", "w").write("\n".join(lines) + "\n")
    print(f"fidelity: {ok} match, {bad} differ, {skipped} without a recorded headline (pre-headline runs)")
    return 0 if bad == 0 else 1


def block(bid="sv:1", **k):
    b = dict(NONE_BLOCK, block_id=bid, zygosity="HET", representation_filter="PASS", genotype_support="READS",
             dosage_status="CONSISTENT", dosage_reason="none")
    b.update(k)
    return b


def locus(**k):
    L = dict(prefix="synthetic", call_class="strict_call", no_call_reason="none", phase_support="READS_WITHIN_CANDIDATES",
             library_sensitivity="NOT_APPLICABLE", dosage_mode="audit", dosage_audit_status="COMPLETE", dosage_audit_note="none",
             library_incomplete="", prediction="PREDICTION", completion="completed_ranking", sequence_status="RESOLVED",
             schema="sequence_first_v2", seq_blocks={}, blocks=[block("sv:1"), block("sv:3")])
    L.update(k)
    L["block_ids"] = [b["block_id"] for b in L["blocks"]]
    return L


SYNTHETIC = [
    ("S1 no dosage requested", locus(dosage_mode="off", dosage_audit_status="NA")),
    ("S2 requested audit failed", locus(dosage_audit_status="FAILED", dosage_audit_note="unit library: alleles disagree in sv:3")),
    ("S3 audit underpowered at one block", locus(blocks=[block("sv:1"), block("sv:3", dosage_status="UNDERPOWERED")])),
    ("S4 dosage consistent", locus()),
    ("S5 complete correction available, not applied (audit)", locus(blocks=[block("sv:1"), block(
        "sv:3", dosage_status="REPAIRED", repair_class="COMPLETE_DIPLOTYPE_REPAIR", dosage_content_scope="COMPLETE",
        repair_not_applied_reason="dosage_mode_audit")])),
    ("S6 complete correction refused by a guard", locus(blocks=[block("sv:1"), block(
        "sv:3", dosage_status="REPAIRED", repair_class="COMPLETE_DIPLOTYPE_REPAIR", dosage_content_scope="COMPLETE",
        repair_not_applied_reason="reported_haplotypes_not_reproduced")])),
    ("S7 conditional content correction (exact diplotype unresolved)", locus(blocks=[block("sv:1"), block(
        "sv:3", dosage_status="CONTENT_REPAIRED", repair_class="CONDITIONAL_DOSAGE_CONTENT", dosage_content_scope="COMPLETE")])),
    ("S8 triggered, repair not run (cap)", locus(blocks=[block("sv:1"), block("sv:3", dosage_status="TRIGGERED_UNREPAIRED",
                                                                                dosage_reason="cap_refused")])),
    ("S9 dosage suspect", locus(blocks=[block("sv:1"), block("sv:3", dosage_status="SUSPECT")])),
    ("S10 representation conflict + library audit incomplete", locus(
        library_sensitivity="NOT_EVALUATED", library_incomplete="primary_no_call",
        blocks=[block("sv:1", representation_filter="LOWQUAL_REPRESENTATION_CONFLICT"), block("sv:3")])),
    ("S11 representation conflict + dosage failed", locus(dosage_audit_status="FAILED", dosage_audit_note="exception",
                                                          blocks=[block("sv:1", representation_filter="LOWQUAL_REPRESENTATION_CONFLICT"), block("sv:3")])),
    ("S12 prior-dependent, otherwise unflagged", locus(phase_support="PRIOR_ONLY", blocks=[block("sv:1"), block("sv:3", genotype_support="PRIOR_ONLY")])),
    ("S13 one sequence pair, ambiguous labels (legacy tie)", locus(call_class="no_call", no_call_reason="distinct_tie", phase_support="UNRESOLVED",
                                                                   dosage_audit_status="NOT_RUN", blocks=[],
                                                                   seq_blocks={"sv:1": {"block_status": "IDENTIFIED"}, "sv:3": {"block_status": "REPRESENTATION_AMBIGUOUS"}})),
    ("S14 genotype available, phase unresolved (genotype-equivalent)", locus(call_class="genotype_equivalent_call", phase_support="UNRESOLVED",
                                                                             sequence_status="GENOTYPE_RESOLVED_SEQUENCE_SET")),
    ("S15 provisional (resource stop after a checkpoint)", locus(call_class="no_call", no_call_reason="time_limit", prediction="PROVISIONAL",
                                                                 completion="budget_stop", sequence_status="NOT_DETERMINED", blocks=[],
                                                                 dosage_mode="off", dosage_audit_status="NA")),
    ("S16 no prediction (stop before a checkpoint)", locus(call_class="no_call", no_call_reason="time_limit", prediction="NO_PREDICTION",
                                                           completion="budget_stop_before_first_evaluated_candidate", sequence_status="NA",
                                                           blocks=[], dosage_mode="off", dosage_audit_status="NA")),
    ("S17 one problematic block among many", locus(blocks=[block(f"b{i}") for i in range(8)] + [block("b8", dosage_status="SUSPECT")] +
                                                   [block("b9", dosage_status="NOT_SENSITIVE")])),
    ("S18 library sensitive at one block", locus(library_sensitivity="SENSITIVE", blocks=[block("sv:1", library_sensitive=True), block("sv:3")])),
    ("S19 legacy repair mode, repair applied", locus(dosage_mode="repair", blocks=[block("sv:1"), block(
        "sv:3", dosage_status="REPAIRED", repair_class="COMPLETE_DIPLOTYPE_REPAIR", dosage_content_scope="COMPLETE", dosage_repair_applied=True)])),
    ("S20 homozygous everywhere", locus(blocks=[block("sv:1", zygosity="HOM"), block("sv:3", zygosity="HOM")])),
]


def unit_synthetic(d):
    with open(f"{d}/synthetic.tsv", "w") as fh:
        fh.write("case\t" + "\t".join("current_" + k for k in CURRENT_KEYS) + "\t" + "\t".join("proposed_" + k for k in PROPOSED_KEYS) + "\n")
        for name, L in SYNTHETIC:
            c, p = current(L), proposed(L)
            fh.write(name + "\t" + "\t".join(c[k] for k in CURRENT_KEYS) + "\t" + "\t".join(p[k] for k in PROPOSED_KEYS) + "\n")
            print(f"{name}\n  current : {c['call_summary']}  [{c['genotype_filter_reasons']} | {c['dosage_filter_reasons']} | {c['phase_filter_reasons']}]\n"
                  f"  proposed: {p['summary']}\n            findings={p['findings']}  not_assessed={p['not_assessed']}  alternatives={p['dosage_alternatives']}")
    return 0


SAVED = [
    ("CYP2D6 HG02769, --dosage-check, audit FAILED (M15-D.1)", "m15d1_correctness_20261021/runs/ex.cyp_audit/cyp2d6"),
    ("C4 NA19434, content correction not applied (M11-D)", "m11d_dev_pilot_20261020/runs/call.c4.NA19434/call"),
    ("C4 NA24385, content correction not applied (M11-D)", "m11d_dev_pilot_20261020/runs/call.c4.NA24385/call"),
    ("CYP2D6 30x, v2 PREDICTION RESOLVED beside legacy distinct_tie (S2-L control)", "s2l_cross_locus_20261019/runs/control.cyp_30/call"),
    ("C4 NA19705, no dosage requested (M15-D no-profile)", "m15d_basic_command_20261021/runs/ex.noprofile_c4/c4"),
    ("C4 NA19705, declared + --dosage-check (M15-D.1 not run; M15-D a2)", "m15d_basic_command_20261021/runs/a2.ex.declared_c4/c4"),
    ("C4 checkpoint stop, PROVISIONAL (M15-D.1)", "m15d1_correctness_20261021/runs/ex.term_checkpoint/c4"),
    ("C4 time stop before a checkpoint (M15-D.1)", "m15d1_correctness_20261021/runs/ex.term_time/c4"),
    ("CYP2D6 NA19983, representation conflict (M14-A full)", "m14a_representative_panel_20261021/runs/call.cyp.NA19983.full/call"),
    ("CYP2D6 NA19983, genotype-equivalent (M14-A r2)", "m14a_representative_panel_20261021/runs/call.cyp.NA19983.r2/call"),
    ("ACOT 5x, dosage suspect (S2-J)", "s2j_origin_admission_20261018/runs/caller.acot5.reuse/call"),
]


def unit_saved(d):
    hashes_before = {}
    for _, rel in SAVED:
        for f in sorted(os.listdir(os.path.dirname(f"{R}/{rel}"))):
            path = os.path.join(os.path.dirname(f"{R}/{rel}"), f)
            if os.path.isfile(path):
                hashes_before[path] = hashlib.sha256(open(path, "rb").read()).hexdigest()
    with open(f"{d}/saved.tsv", "w") as fh:
        fh.write("example\tprefix\t" + "\t".join("current_" + k for k in CURRENT_KEYS) + "\t" + "\t".join("proposed_" + k for k in PROPOSED_KEYS) + "\n")
        for name, rel in SAVED:
            L = facts(f"{R}/{rel}")
            c, p = current(L), proposed(L)
            rec = L["recorded"].get("call_summary") or "(no locus table: minimal status only)"
            fh.write(f"{name}\t{rel}\t" + "\t".join(c[k] for k in CURRENT_KEYS) + "\t" + "\t".join(p[k] for k in PROPOSED_KEYS) + "\n")
            print(f"{name}\n  recorded: {rec}\n  current : {c['call_summary']}  [{c['genotype_filter_reasons'][:110]} | {c['dosage_filter_reasons'][:110]}]\n"
                  f"  proposed: {p['summary']}\n            findings={p['findings'][:200]}\n            not_assessed={p['not_assessed'][:200]}\n"
                  f"            alternatives={p['dosage_alternatives']}")
    changed = [p for p, h in hashes_before.items() if hashlib.sha256(open(p, "rb").read()).hexdigest() != h]
    print(f"originals unchanged: {len(hashes_before) - len(changed)} of {len(hashes_before)} files")
    return 0 if not changed else 1


if __name__ == "__main__":
    unit, d = sys.argv[1], sys.argv[2]
    os.makedirs(d, exist_ok=True)
    sys.exit(globals()["unit_" + unit](d))
