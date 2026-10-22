#!/usr/bin/env python3
"""M15-G: the corrected assessment contract (assessment_v1, as amended), as an executable specification.

It reads a v2 run prefix's own outputs and returns the exact strings the C++ implementation must write: the locus rows
(appended to .sequence_locus.tsv and mirrored into .locus_genotype.tsv) and the block columns (appended to
.sequence_blocks.tsv and mirrored into the legacy block tables). CONTRACT.md is the prose; this is the reference.

Corrections over M15-F's draft (each with a fixture in FIXTURES):
  C1 the legacy classification never overrides v2 sequence resolution;
  C2 phase applicability is two or more heterozygous blocks anywhere, not adjacent HET/HET pairs;
  C3 unknown read support is never converted to READS (per block NOT_ASSESSED; locus coverage counts);
  C4 a provisional result never yields a decision set or completed support, whatever block rows exist;
  C5 findings and missing assessments are listed together, with compact coverage counts.
Verified context: a legacy block fact is used only for a completed v2 PREDICTION whose legacy decision exists and whose
legacy reported alleles (sorted) equal the v2 published genotype at that block; otherwise the block is NOT_ASSESSED."""
import csv
import os
import sys

csv.field_size_limit(sys.maxsize)
RULE = "assessment_v1"
LOCUS_KEYS = ["assessment_rule", "sequence_resolution", "genotype_check", "genotype_check_coverage", "dosage_check",
              "dosage_check_coverage", "read_support", "read_support_coverage", "phase", "phase_reason", "findings",
              "not_assessed", "dosage_alternatives", "assessment_summary"]
BLOCK_COLUMNS = ["block_resolution", "genotype_check", "genotype_check_reasons", "dosage_check", "dosage_check_reason",
                 "dosage_alternative", "read_support", "read_support_reason"]
FINDING_ORDER = ["REPRESENTATION_CONFLICT", "LIBRARY_SENSITIVE", "DOSAGE_DISAGREES", "DOSAGE_SUSPECT", "PRIOR_DEPENDENT",
                 "LABELS_AMBIGUOUS", "GENOTYPE_UNRESOLVED"]


def kv(p):
    out = {}
    if os.path.exists(p):
        for line in open(p):
            if "\t" in line:
                k, v = line.rstrip("\n").split("\t", 1)
                out.setdefault(k, v)
    return out


def rows(p):
    return list(csv.DictReader(open(p), delimiter="\t")) if os.path.exists(p) else []


def clean(text):
    """A free-text note made safe inside a CODE(note)@ids list."""
    return "".join("_" if c in ";@\t\n()" else c for c in text)


def read(prefix):
    """The facts, from the prefix's own files only."""
    st, sl = kv(prefix + ".mosaic_status.tsv"), kv(prefix + ".sequence_locus.tsv")
    run_status = st.get("status", "NA")
    call_class = st.get("call_class") or ("strict_call" if run_status == "COMPLETE" else
                                          "genotype_equivalent_call" if run_status == "CONDITIONAL_GENOTYPE_EQUIVALENT" else "no_call")
    if call_class == "no_call" and sl.get("legacy_call_class"):
        call_class = sl["legacy_call_class"]
    no_call_reason = st.get("no_call_reason") or sl.get("legacy_no_call_reason") or ("status_" + run_status)
    declared = st.get("library_provenance", "").startswith("declared")
    lib = "NOT_APPLICABLE" if declared else st.get("library_sensitivity_status", st.get("library_sensitivity", "NOT_EVALUATED"))
    ids = lambda key: set() if st.get(key, "none") in ("none", "") else set(st[key].split(","))
    mode = sl.get("dosage_mode", "off")
    if mode in ("", "NA", "not_run"):
        mode = "off"
    ds = kv(prefix + ".dosage_status.tsv")
    return dict(
        schema=sl.get("schema_version", ""), prediction=sl.get("prediction", "NO_PREDICTION"), completion=sl.get("completion", "NA"),
        sequence_status=sl.get("sequence_status", "NA"), call_class=call_class, no_call_reason=no_call_reason,
        phase_support=st.get("phase_support", sl.get("legacy_phase_support", "NA")), library=lib,
        library_incomplete=st.get("library_sensitivity_incomplete", "not_evaluated") or "not_evaluated",
        sensitive=ids("library_sensitive_blocks"), unaudited=ids("library_sensitivity_not_evaluated_blocks"),
        dosage_mode=mode, audit=ds.get("dosage_audit_status", "NA") if mode != "off" else "NA",
        audit_note=ds.get("dosage_audit_note", "none"),
        v2_blocks=rows(prefix + ".sequence_blocks.tsv"),
        legacy={r["block_id"]: r for r in rows(prefix + ".block_genotypes.tsv")},
        dosage={r["block_id"]: r for r in rows(prefix + ".dosage_blocks.tsv")})


def sorted_pair(a, b):
    try:
        x, y = sorted((int(a), int(b)))
        return f"{x}/{y}"
    except ValueError:
        return "NA"


def zygosity(genotype):
    if genotype in ("", "NA") or "/" not in genotype:
        return "UNKNOWN"
    a, b = genotype.split("/", 1)
    return "HOM" if a == b else "HET"


def unverified_reason(F, v2, legacy):
    """Why a legacy block fact may not be used for the published prediction ('' when it may)."""
    if F["prediction"] == "PROVISIONAL":
        return "provisional"
    if F["prediction"] != "PREDICTION":
        return "no_prediction"
    if v2.get("block_status") != "IDENTIFIED" or zygosity(v2.get("genotype", "NA")) == "UNKNOWN":
        return "genotype_not_identified"
    if F["call_class"] not in ("strict_call", "genotype_equivalent_call"):
        return "legacy_decision_" + F["no_call_reason"]
    if legacy is None:
        return "no_legacy_block"
    if sorted_pair(legacy.get("allele1", ""), legacy.get("allele2", "")) != v2.get("genotype"):
        return "legacy_block_differs"
    return ""


def repair_not_applied(legacy):
    """The audit-mode record of a complete repair not applied, from this run's own headline reasons."""
    rec = legacy.get("dosage_filter_reasons", "") if legacy else ""
    if "COMPLETE_REPAIR_AVAILABLE_NOT_APPLIED" in rec:
        return "dosage_mode_audit"
    if rec.startswith("COMPLETE_REPAIR_NOT_APPLIED:"):
        return rec.split(":", 1)[1]
    if legacy and legacy.get("dosage_stage_filter") == "DOSAGE_REPAIR_NOT_APPLIED":
        return legacy.get("dosage_stage_filter_reason", "not_recorded")
    return "not_recorded"


def block_dosage(F, why, legacy, drow):
    if F["dosage_mode"] == "off":
        return "NOT_REQUESTED", "dosage_check_not_requested", "NA"
    if F["audit"] != "COMPLETE":
        return "UNAVAILABLE", {"FAILED": "AUDIT_FAILED", "NOT_RUN": "AUDIT_NOT_RUN",
                               "REFUSED_LIBRARY": "AUDIT_UNSUPPORTED_LIBRARY"}.get(F["audit"], "AUDIT_" + F["audit"]), "NA"
    if why:
        return "UNAVAILABLE", "NOT_VERIFIED:" + why, "NA"
    pick = lambda k: (legacy.get(k) if legacy and legacy.get(k) else drow.get(k, "NA"))
    s = pick("dosage_status")
    applied = legacy is not None and legacy.get("dosage_repair_applied", "NA") == "true"
    if s == "NOT_SENSITIVE":
        return "NOT_APPLICABLE", "NOT_DOSAGE_SENSITIVE", "NONE"
    if s == "CONSISTENT":
        return "AGREES", "CONSISTENT", "NONE"
    if s == "SUSPECT":
        return "SUSPECT", "SUSPECT", "NONE"
    if s == "UNDERPOWERED":
        return "UNAVAILABLE", "UNDERPOWERED", "NONE"
    if s == "REPAIRED" and applied:
        return "AGREES", "REPORTED_IS_DOSAGE_REPAIR", "APPLIED"
    if s == "REPAIRED":
        r = repair_not_applied(legacy)
        return "DISAGREES", "TRIGGERED", "COMPLETE_REPAIR_AVAILABLE" if r == "dosage_mode_audit" else "COMPLETE_REPAIR_REFUSED:" + r
    if s == "CONTENT_REPAIRED":
        return "DISAGREES", "TRIGGERED", "CONTENT_CORRECTION_AVAILABLE"
    if s in ("TRIGGERED_UNREPAIRED", "AUTO_REPAIR"):
        return "DISAGREES", "TRIGGERED", "REPAIR_NOT_FOUND:" + (pick("dosage_reason") or "none")
    return "UNAVAILABLE", "NOT_EVALUATED:" + (pick("dosage_reason") or "NA"), "NONE"


def assess(F):
    """(locus: dict LOCUS_KEYS -> str, blocks: list of (block_id, dict BLOCK_COLUMNS -> str))"""
    pred = F["prediction"]
    blocks_out, findings, missing = [], [], []
    all_ids = [r["block_id"] for r in F["v2_blocks"]]

    def note(lst, code, bid):
        e = next((x for x in lst if x[0] == code), None)
        if e is None:
            e = [code, []]
            lst.append(e)
        if bid not in e[1]:
            e[1].append(bid)

    gvals, dvals, svals, alts, zyg, verified_all = [], [], [], [], [], True
    for v2 in (F["v2_blocks"] if pred != "NO_PREDICTION" else []):
        bid = v2["block_id"]
        legacy = F["legacy"].get(bid)
        why = unverified_reason(F, v2, legacy)
        verified_all = verified_all and not why
        zyg.append(zygosity(v2.get("genotype", "NA")))
        status = v2.get("block_status", "NA")
        res = {"REPRESENTATION_AMBIGUOUS": "LABELS_AMBIGUOUS"}.get(status, status)
        codes, reasons = [], []
        # genotype check
        if why:
            g = "NOT_ASSESSED"
            reasons.append(why)
            note(missing, "GENOTYPE_CHECK:" + why, bid)
        else:
            rep = legacy.get("filter", "NA")
            found = []
            if rep == "LOWQUAL_REPRESENTATION_CONFLICT":
                found.append("REPRESENTATION_CONFLICT")
            if F["library"] == "SENSITIVE" and bid in F["sensitive"]:
                found.append("LIBRARY_SENSITIVE")
            if rep not in ("PASS", "LOWQUAL_REPRESENTATION_CONFLICT"):
                reasons.append("representation_guard_" + rep)
                note(missing, "GENOTYPE_CHECK:representation_guard_" + rep, bid)
            if F["library"] == "NOT_EVALUATED" or bid in F["unaudited"]:
                if not (F["library"] == "SENSITIVE" and bid in F["sensitive"]):
                    note(missing, "LIBRARY_SENSITIVITY:" + F["library_incomplete"], bid)
            reasons = found + reasons
            codes += found
            g = "CONFLICT" if found else ("NONE_DETECTED" if rep in ("PASS", "LOWQUAL_REPRESENTATION_CONFLICT") else "NOT_ASSESSED")
        gvals.append(g)
        # dosage check (the reported alleles); the alternative is diagnostic
        d, dr, alt = block_dosage(F, why, legacy, F["dosage"].get(bid, {}))
        dvals.append(d)
        if d in ("DISAGREES", "SUSPECT"):
            codes.append("DOSAGE_" + d)
        if d == "UNAVAILABLE":
            if dr.startswith("AUDIT_"):
                note(missing, "DOSAGE:" + dr + ("(" + clean(F["audit_note"]) + ")" if dr == "AUDIT_FAILED" else ""), "locus")
            else:
                note(missing, "DOSAGE:" + dr, bid)
        if alt not in ("NA", "NONE"):
            alts.append(f"{alt}@{bid}")
        # read support
        support = legacy.get("genotype_support", "NA") if (legacy and not why) else "NA"
        if why:
            s, sr = "NOT_ASSESSED", why
        elif support == "READS":
            s, sr = "READS", "reads"
        elif support == "PRIOR_ONLY":
            s, sr = "PRIOR_DEPENDENT", "prior_only"
            codes.append("PRIOR_DEPENDENT")
        else:
            s, sr = "NOT_ASSESSED", "support_" + support
        if s == "NOT_ASSESSED":
            note(missing, "READ_SUPPORT:" + sr, bid)
        svals.append(s)
        if pred == "PREDICTION" and res == "LABELS_AMBIGUOUS":
            codes.append("LABELS_AMBIGUOUS")
        if pred == "PREDICTION" and res == "GENOTYPE_UNRESOLVED":
            codes.append("GENOTYPE_UNRESOLVED")
        for c in FINDING_ORDER:
            if c in codes:
                note(findings, c, bid)
        blocks_out.append((bid, {"block_resolution": res, "genotype_check": g, "genotype_check_reasons": ",".join(reasons) or "none",
                                 "dosage_check": d, "dosage_check_reason": dr, "dosage_alternative": alt,
                                 "read_support": s, "read_support_reason": sr}))
    # locus aggregation
    n = len(blocks_out)
    if pred == "NO_PREDICTION":
        L = {"sequence_resolution": "NA", "genotype_check": "NA", "genotype_check_coverage": "NA", "dosage_check": "NA",
             "dosage_check_coverage": "NA", "read_support": "NA", "read_support_coverage": "NA"}
    else:
        ga = sum(1 for x in gvals if x in ("CONFLICT", "NONE_DETECTED"))
        g = "CONFLICT" if "CONFLICT" in gvals else "NONE_DETECTED" if ga else "NOT_ASSESSED"
        if F["dosage_mode"] == "off":
            d, dcov = "NOT_REQUESTED", "NA"
        else:
            applicable = [x for x in dvals if x != "NOT_APPLICABLE"]
            da = sum(1 for x in applicable if x in ("AGREES", "SUSPECT", "DISAGREES"))
            d = next((x for x in ("DISAGREES", "SUSPECT", "AGREES") if x in dvals), None) or \
                ("NOT_APPLICABLE" if dvals and not applicable else "UNAVAILABLE")
            dcov = f"{da}/{len(applicable)}"
        sa = sum(1 for x in svals if x in ("READS", "PRIOR_DEPENDENT"))
        s = "PRIOR_DEPENDENT" if "PRIOR_DEPENDENT" in svals else "READS" if sa else "NOT_ASSESSED"
        L = {"sequence_resolution": F["sequence_status"] if pred == "PREDICTION" else "NOT_DETERMINED",
             "genotype_check": g, "genotype_check_coverage": f"{ga}/{n}", "dosage_check": d, "dosage_check_coverage": dcov,
             "read_support": s, "read_support_coverage": f"{sa}/{n}"}
    # phase: locus-level by construction
    het, unknown = zyg.count("HET"), zyg.count("UNKNOWN")
    if pred == "NO_PREDICTION":
        p, pr = "NA", "no_prediction"
    elif pred == "PROVISIONAL":
        p, pr = "NOT_ASSESSED", "provisional"
    elif F["sequence_status"] in ("GENOTYPE_RESOLVED_SEQUENCE_SET", "SEQUENCE_TIED"):
        p, pr = "UNRESOLVED", "several_sequence_pairs"
    elif unknown == 0 and het < 2:
        p, pr = "NOT_APPLICABLE", "fewer_than_two_heterozygous_blocks"
    elif F["call_class"] != "strict_call":
        p, pr = "NOT_ASSESSED", "legacy_decision_" + ("genotype_equivalent" if F["call_class"] == "genotype_equivalent_call" else F["no_call_reason"])
    elif not verified_all:
        p, pr = "NOT_ASSESSED", "blocks_not_verified"
    elif F["phase_support"] in ("READS_AUDITED", "READS_WITHIN_CANDIDATES"):
        p, pr = "READ_SUPPORTED", F["phase_support"].lower()
    elif F["phase_support"] == "PRIOR_ONLY":
        p, pr = "PRIOR_ONLY", "prior_only"
    else:
        p, pr = "NOT_ASSESSED", "phase_support_" + F["phase_support"]
    if p in ("PRIOR_ONLY", "UNRESOLVED"):
        note(findings, "PHASE_" + p, "locus")
    if p == "NOT_ASSESSED":
        note(missing, "PHASE:" + pr, "locus")

    def fmt(lst):
        parts = []
        for code, ids in lst:
            idtext = "all" if (len(all_ids) > 1 and ids == all_ids) else ",".join(ids)
            parts.append(f"{code}@{idtext}")
        return ";".join(parts) or "none"

    cov = lambda key: (f"[{L[key]}]" if L[key] not in ("NA",) and "/" in L[key] and
                       0 < int(L[key].split("/")[0]) < int(L[key].split("/")[1]) else "")
    L.update({"assessment_rule": RULE, "phase": p, "phase_reason": pr, "findings": fmt(findings), "not_assessed": fmt(missing),
              "dosage_alternatives": ";".join(alts) or "none"})
    L["assessment_summary"] = (f"{pred};SEQUENCE_{L['sequence_resolution']};GENOTYPE_{L['genotype_check']}{cov('genotype_check_coverage')};"
                               f"DOSAGE_{L['dosage_check']}{cov('dosage_check_coverage')};SUPPORT_{L['read_support']}{cov('read_support_coverage')};"
                               f"PHASE_{p}")
    return {k: L[k] for k in LOCUS_KEYS}, blocks_out


# ---------------------------------------------------------------- fixtures (as files) -----------------------------------
def write_case(d, name, *, prediction="PREDICTION", completion="completed_ranking", sequence_status="RESOLVED",
               status="COMPLETE", call_class="strict_call", no_call_reason="none", phase_support="READS_WITHIN_CANDIDATES",
               library="declared", sensitive=(), unaudited=(), incomplete="", dosage_mode="audit", audit="COMPLETE",
               audit_note="none", blocks=(), legacy=True):
    """blocks: (block_id, block_status, genotype, legacy dict or None, dosage dict or None). Legacy dict keys: allele1,
    allele2, filter, genotype_support, dosage_filter_reasons, dosage_repair_applied."""
    p = os.path.join(d, name)
    st = [("status", status), ("call_class", call_class), ("no_call_reason", no_call_reason), ("phase_support", phase_support)]
    if library == "declared":
        st.append(("library_provenance", "declared: library.txt"))
    else:
        st += [("library_provenance", "auto: test"), ("library_sensitivity_status", library),
               ("library_sensitive_blocks", ",".join(sensitive) or "none"),
               ("library_sensitivity_not_evaluated_blocks", ",".join(unaudited) or "none")]
        if incomplete:
            st.append(("library_sensitivity_incomplete", incomplete))
    open(p + ".mosaic_status.tsv", "w").write("key\tvalue\n" + "".join(f"{k}\t{v}\n" for k, v in st))
    sl = [("schema_version", "sequence_first_v2"), ("prediction", prediction), ("sequence_status", sequence_status),
          ("completion", completion), ("dosage_mode", dosage_mode), ("legacy_call_class", call_class),
          ("legacy_no_call_reason", no_call_reason), ("exit_code", {"PREDICTION": "0", "PROVISIONAL": "3"}.get(prediction, "1"))]
    open(p + ".sequence_locus.tsv", "w").write("key\tvalue\n" + "".join(f"{k}\t{v}\n" for k, v in sl))
    open(p + ".sequence_blocks.tsv", "w").write("block_id\tblock_status\tgenotype\n" + "".join(f"{b[0]}\t{b[1]}\t{b[2]}\n" for b in blocks))
    lg = [b for b in blocks if b[3] is not None]
    if legacy and lg:
        cols = ["block_id", "allele1", "allele2", "zygosity", "genotype_support", "filter", "genotype_filter", "genotype_filter_reasons",
                "dosage_filter", "dosage_filter_reasons", "dosage_repair_applied"]
        with open(p + ".block_genotypes.tsv", "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for b in lg:
                r = dict(allele1="0", allele2="1", filter="PASS", genotype_support="READS", dosage_filter_reasons="none",
                         dosage_repair_applied="NA", genotype_filter="PASS", genotype_filter_reasons="none", dosage_filter="PASS")
                r.update(b[3])
                r["block_id"], r["zygosity"] = b[0], "HOM" if r["allele1"] == r["allele2"] else "HET"
                fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    if dosage_mode != "off":
        open(p + ".dosage_status.tsv", "w").write(f"key\tvalue\ndosage_audit\tshadow\ndosage_audit_status\t{audit}\ndosage_audit_note\t{audit_note}\n")
        with open(p + ".dosage_blocks.tsv", "w") as fh:
            fh.write("block_id\tdosage_status\trepair_class\tdosage_content_scope\tdosage_reason\n")
            for b in blocks:
                dd = dict(dosage_status="CONSISTENT", repair_class="NA", dosage_content_scope="NA", dosage_reason="none")
                dd.update(b[4] or {})
                fh.write(f"{b[0]}\t{dd['dosage_status']}\t{dd['repair_class']}\t{dd['dosage_content_scope']}\t{dd['dosage_reason']}\n")
    return p


def B(bid, genotype="0/1", status="IDENTIFIED", legacy="auto", dosage=None, **leg):
    """A block: v2 row + a matching legacy row (alleles from the genotype) unless legacy=None."""
    if legacy == "auto":
        a, b = (genotype.split("/") + ["0"])[:2] if "/" in genotype else ("0", "1")
        legacy = dict(allele1=b, allele2=a)   # legacy order is hap1/hap2, not sorted
        legacy.update(leg)
    return (bid, status, genotype, legacy, dosage)


# (name, kwargs). The 20 M15-F decision-table cases and the five M15-G counterexamples (G1-G5).
FIXTURES = [
    ("S01_no_dosage_requested", dict(dosage_mode="off", blocks=[B("sv:1"), B("sv:3")])),
    ("S02_requested_audit_failed", dict(audit="FAILED", audit_note="unit library: alleles disagree in sv:3", blocks=[B("sv:1"), B("sv:3")])),
    ("S03_underpowered_one_block", dict(blocks=[B("sv:1"), B("sv:3", dosage=dict(dosage_status="UNDERPOWERED"))])),
    ("S04_dosage_consistent", dict(blocks=[B("sv:1"), B("sv:3")])),
    ("S05_complete_repair_available", dict(blocks=[B("sv:1"), B("sv:3", dosage=dict(dosage_status="REPAIRED", repair_class="COMPLETE_DIPLOTYPE_REPAIR",
                                                                                     dosage_content_scope="COMPLETE"),
                                                                 dosage_filter_reasons="COMPLETE_REPAIR_AVAILABLE_NOT_APPLIED")])),
    ("S06_repair_refused_by_guard", dict(blocks=[B("sv:1"), B("sv:3", dosage=dict(dosage_status="REPAIRED", repair_class="COMPLETE_DIPLOTYPE_REPAIR"),
                                                              dosage_filter_reasons="COMPLETE_REPAIR_NOT_APPLIED:reported_haplotypes_not_reproduced")])),
    ("S07_content_correction", dict(blocks=[B("sv:1"), B("sv:3", dosage=dict(dosage_status="CONTENT_REPAIRED", repair_class="CONDITIONAL_DOSAGE_CONTENT",
                                                                             dosage_content_scope="COMPLETE"))])),
    ("S08_triggered_repair_not_run", dict(blocks=[B("sv:1"), B("sv:3", dosage=dict(dosage_status="TRIGGERED_UNREPAIRED", dosage_reason="cap_refused"))])),
    ("S09_dosage_suspect", dict(blocks=[B("sv:1"), B("sv:3", dosage=dict(dosage_status="SUSPECT"))])),
    ("S10_conflict_plus_library_incomplete", dict(library="NOT_EVALUATED", incomplete="primary_no_call",
                                                  blocks=[B("sv:1", filter="LOWQUAL_REPRESENTATION_CONFLICT"), B("sv:3")])),
    ("S11_conflict_plus_audit_failed", dict(audit="FAILED", audit_note="exception", blocks=[B("sv:1", filter="LOWQUAL_REPRESENTATION_CONFLICT"), B("sv:3")])),
    ("S12_prior_dependent", dict(phase_support="PRIOR_ONLY", blocks=[B("sv:1"), B("sv:3", genotype_support="PRIOR_ONLY")])),
    ("S13_legacy_tie_one_pair", dict(status="RESCORE_INCOMPLETE", call_class="no_call", no_call_reason="distinct_tie", phase_support="UNRESOLVED",
                                     audit="NOT_RUN", blocks=[B("sv:1", legacy=None), B("sv:3", genotype="NA", status="REPRESENTATION_AMBIGUOUS", legacy=None)])),
    ("S14_genotype_equivalent", dict(status="CONDITIONAL_GENOTYPE_EQUIVALENT", call_class="genotype_equivalent_call", phase_support="UNRESOLVED",
                                     sequence_status="GENOTYPE_RESOLVED_SEQUENCE_SET", blocks=[B("sv:1"), B("sv:3")])),
    ("S15_provisional", dict(prediction="PROVISIONAL", completion="budget_stop", sequence_status="NOT_DETERMINED", status="TIME_LIMIT",
                             call_class="no_call", no_call_reason="time_limit", dosage_mode="off",
                             blocks=[B("sv:1", status="PROVISIONAL_LABEL", legacy=None), B("sv:3", status="PROVISIONAL_LABEL", legacy=None)])),
    ("S16_no_prediction", dict(prediction="NO_PREDICTION", completion="budget_stop_before_first_evaluated_candidate", sequence_status="NA",
                               status="TIME_LIMIT", call_class="no_call", no_call_reason="time_limit", dosage_mode="off", blocks=[])),
    ("S17_one_problem_among_ten", dict(blocks=[B(f"b{i}") for i in range(8)] + [B("b8", dosage=dict(dosage_status="SUSPECT")),
                                                                               B("b9", dosage=dict(dosage_status="NOT_SENSITIVE"))])),
    ("S18_library_sensitive", dict(library="SENSITIVE", sensitive=("sv:1",), blocks=[B("sv:1"), B("sv:3")])),
    ("S19_legacy_repair_applied", dict(dosage_mode="repair", blocks=[B("sv:1"), B("sv:3", genotype="0/2", allele1="1", allele2="0",
                                                                                   dosage_repair_applied="true",
                                                                                   dosage=dict(dosage_status="REPAIRED", repair_class="COMPLETE_DIPLOTYPE_REPAIR"))])),
    ("S20_homozygous", dict(blocks=[B("sv:1", genotype="0/0"), B("sv:3", genotype="1/1")])),
    # M15-G counterexamples
    ("G1_unique_pair_beside_legacy_equivalent", dict(status="CONDITIONAL_GENOTYPE_EQUIVALENT", call_class="genotype_equivalent_call",
                                                    phase_support="UNRESOLVED", sequence_status="RESOLVED", blocks=[B("sv:1"), B("sv:3")])),
    ("G2_het_separated_by_hom", dict(blocks=[B("sv:1", genotype="0/1"), B("ctx:2", genotype="2/2"), B("sv:3", genotype="0/1")])),
    ("G3_mixed_read_support", dict(blocks=[B("sv:1"), B("sv:3", genotype_support="NA"), B("sv:5", legacy=None)])),
    ("G4_provisional_with_block_rows", dict(prediction="PROVISIONAL", completion="round_cap", sequence_status="NOT_DETERMINED",
                                            status="REPAIR_INCOMPLETE", call_class="no_call", no_call_reason="round_cap",
                                            blocks=[B("sv:1", status="PROVISIONAL_LABEL"), B("sv:3", status="PROVISIONAL_LABEL", genotype_support="PRIOR_ONLY",
                                                                                         dosage=dict(dosage_status="SUSPECT"))])),
    ("G5_agreement_beside_unavailable", dict(blocks=[B("sv:1"), B("sv:3", dosage=dict(dosage_status="UNDERPOWERED")),
                                                     B("sv:5", dosage=dict(dosage_status="NOT_EVALUATED", dosage_reason="no_tested_feature")),
                                                     B("sv:7", dosage=dict(dosage_status="NOT_SENSITIVE"))])),
    # verified context: a legacy block whose alleles differ from the published v2 genotype is not used
    ("G6_legacy_block_differs", dict(blocks=[B("sv:1"), B("sv:3", genotype="0/1", allele1="2", allele2="0", filter="LOWQUAL_REPRESENTATION_CONFLICT")])),
]

# What each counterexample must show (checked by `check_fixtures`; C1-C5 of the module docstring).
EXPECT = {
    "G1_unique_pair_beside_legacy_equivalent": lambda L, Bk: L["sequence_resolution"] == "RESOLVED" and L["phase"] == "NOT_ASSESSED"
        and L["phase_reason"] == "legacy_decision_genotype_equivalent",
    "G2_het_separated_by_hom": lambda L, Bk: L["phase"] == "READ_SUPPORTED",
    "G3_mixed_read_support": lambda L, Bk: L["read_support"] == "READS" and L["read_support_coverage"] == "1/3" and "SUPPORT_READS[1/3]" in L["assessment_summary"]
        and "READ_SUPPORT:support_NA@sv:3" in L["not_assessed"] and "READ_SUPPORT:no_legacy_block@sv:5" in L["not_assessed"],
    "G4_provisional_with_block_rows": lambda L, Bk: L["sequence_resolution"] == "NOT_DETERMINED" and L["read_support"] == "NOT_ASSESSED"
        and L["dosage_check"] == "UNAVAILABLE" and L["genotype_check"] == "NOT_ASSESSED" and L["findings"] == "none" and L["phase"] == "NOT_ASSESSED"
        and all(b["read_support"] == "NOT_ASSESSED" and b["block_resolution"] == "PROVISIONAL_LABEL" for _, b in Bk),
    "G5_agreement_beside_unavailable": lambda L, Bk: L["dosage_check"] == "AGREES" and L["dosage_check_coverage"] == "1/3"
        and "DOSAGE:UNDERPOWERED@sv:3" in L["not_assessed"] and "DOSAGE:NOT_EVALUATED:no_tested_feature@sv:5" in L["not_assessed"],
    "G6_legacy_block_differs": lambda L, Bk: "REPRESENTATION_CONFLICT" not in L["findings"] and "GENOTYPE_CHECK:legacy_block_differs@sv:3" in L["not_assessed"],
    "S10_conflict_plus_library_incomplete": lambda L, Bk: "REPRESENTATION_CONFLICT@sv:1" in L["findings"] and "LIBRARY_SENSITIVITY:primary_no_call@all" in L["not_assessed"],
    "S13_legacy_tie_one_pair": lambda L, Bk: L["sequence_resolution"] == "RESOLVED" and "LABELS_AMBIGUOUS@sv:3" in L["findings"],
    "S16_no_prediction": lambda L, Bk: L["genotype_check"] == "NA" and not Bk,
}


def write_fixtures(d):
    os.makedirs(d, exist_ok=True)
    return [(name, write_case(d, name, **kw)) for name, kw in FIXTURES]


def check_fixtures(d):
    bad = []
    for name, prefix in write_fixtures(d):
        L, Bk = assess(read(prefix))
        if name in EXPECT and not EXPECT[name](L, Bk):
            bad.append(name)
        print(f"{name}: {L['assessment_summary']}\n    findings={L['findings']}  not_assessed={L['not_assessed']}  alternatives={L['dosage_alternatives']}")
    print(f"counterexample expectations: {len([n for n in EXPECT]) - len(bad)} of {len(EXPECT)} hold" + (f"; FAILED {bad}" if bad else ""))
    return 0 if not bad else 1


def render(prefix):
    """The exact text the C++ `genotype_assessment_test show <prefix>` prints, for byte comparison."""
    L, Bk = assess(read(prefix))
    out = ["#locus"] + [f"{k}\t{L[k]}" for k in LOCUS_KEYS] + ["#blocks", "block_id\t" + "\t".join(BLOCK_COLUMNS)]
    out += [bid + "\t" + "\t".join(b[c] for c in BLOCK_COLUMNS) for bid, b in Bk]
    return "\n".join(out) + "\n"


if __name__ == "__main__":
    if sys.argv[1] == "fixtures":
        sys.exit(check_fixtures(sys.argv[2]))
    elif sys.argv[1] == "render":
        sys.stdout.write(render(sys.argv[2]))
