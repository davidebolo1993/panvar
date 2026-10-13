#!/usr/bin/env python3
"""--dosage-mode integration gate (docs/modules/genotype.md, Dosage modes), on one C4 campaign cell re-run
end to end with `panvar genotype`:

  identity   the flag absent and `--dosage-mode off` equal the sealed 11fe681 outputs (normalised as the shadow gate
             normalises them) and add only <prefix>.locus_genotype.tsv; audit leaves every pre-existing output equal to
             off and adds only <prefix>.dosage_*;
  repair     every block that is not an applied COMPLETE_DIPLOTYPE_REPAIR equals audit in every legacy column; an
             applied block differs only in its genotype columns; the search record (scores, candidates, equivalence
             set, trajectory) is unchanged; the status is audit's with dosage keys appended;
  round trip <prefix>.original.* equal audit's outputs;
  schema     the trailing block columns and the locus table's columns, in order;
  leakage    no conditional representative in a primary allele column, the locus digests or the main FASTA;
  locus      the locus table agrees with the block table, the FASTA and the equivalence set;
  evaluator  (--probe) the evaluator scores <repair>.original exactly as it scores audit, and the repaired prefix
             from the locus table only when a repair was applied;
  threads    (--threads-sweep) repair mode at 1/2/4/8 threads: every deterministic output byte-identical.

Usage: genotype_dosage_mode_gate.py --panvar BIN --jobs JOBS.tsv --cell CELL --scratch DIR [--probe BIN] [--threads-sweep]
"""
import argparse, hashlib, os, re, subprocess, sys

LEGACY = [".genotypes.tsv", ".block_genotypes.tsv", ".mosaic_haplotypes.fa", ".mosaic_status.tsv", ".mosaic_scores.tsv",
          ".mosaic_candidates.tsv", ".equivalence_set.tsv"]
SEARCH_RECORD = [".mosaic_scores.tsv", ".mosaic_candidates.tsv", ".equivalence_set.tsv", ".repair_rounds.tsv",
                 ".tie_branches.tsv", ".tie_terminal_ties.tsv"]
DOSAGE_COLS = ["dosage_status", "dosage_reason", "dosage_stage_filter", "dosage_stage_filter_reason", "repair_class", "dosage_repair_applied",
               "dosage_content_scope", "diplotype_scope", "allocation_equivalence_scope", "serialization_representative",
               "original_allele1", "original_allele2", "original_content", "corrected_content", "admissible_content_vectors",
               "allocation_status", "original_GQ", "gq_status", "repair_separation", "conditional_repair_gap", "original_filter",
               "original_quality_filter", "original_genotype_support", "filter_basis", "quality_filter_basis"]
# headline_filters_v2: trailing columns every mode writes; they legitimately differ between dosage modes, so the
# legacy comparisons below strip them, and the headline gate (tests/results/headline_filters_impl) checks them.
HEADLINE_BLOCK_COLS = ["genotype_filter", "genotype_filter_reasons", "dosage_filter", "dosage_filter_reasons"]
HEADLINE_LOCUS_COLS = ["headline_rule", "genotype_filter", "genotype_filter_reasons", "dosage_filter", "dosage_filter_reasons",
                       "phase_filter", "phase_filter_reasons", "call_summary"] + \
                      [f"{f}_filter_{s}" for f in ("genotype", "dosage", "phase")
                       for s in ("units", "applicable", "pass", "low_quality", "not_evaluated", "no_call", "not_applicable")]
def strip_headline(text, path):
    drop = HEADLINE_BLOCK_COLS if path.endswith((".genotypes.tsv", ".block_genotypes.tsv")) else \
           HEADLINE_LOCUS_COLS if path.endswith(".locus_genotype.tsv") else None
    if drop is None or not text: return text
    rows = [l.split("\t") for l in text.split("\n")]
    keep = [i for i, c in enumerate(rows[0]) if c not in drop]
    return "\n".join("\t".join(r[i] for i in keep if i < len(r)) if r != [""] else "" for r in rows)
LOCUS_COLS = ["prefix", "run_status", "call_class", "no_call_reason", "locus_filter", "quality_rule", "quality_pass_blocks",
              "quality_lowqual_blocks", "filter_basis", "phase_support", "candidate_scope", "search_conditional", "reported_genotype_source",
              "hap1_id", "hap1_sha256", "hap1_bp", "hap1_block_path", "hap2_id", "hap2_sha256", "hap2_bp", "hap2_block_path",
              "unordered_genotype", "equivalence_members", "total_score", "fragment_term", "exposure_term", "sequence_pair_prior",
              "score_basis", "caller_total_score", "separation", "separation_basis", "separation_certified", "hap1_switches",
              "hap1_switch_blocks", "hap2_switches", "hap2_switch_blocks", "switch_basis", "dosage_mode", "dosage_audit_status",
              "dosage_complete_class_blocks", "dosage_complete_applied_blocks", "dosage_complete_applied_block_ids",
              "dosage_complete_not_applied_blocks", "dosage_complete_not_applied_detail", "dosage_content_corrected_blocks",
              "dosage_content_corrected_block_ids", "dosage_triggered_unrepaired_blocks", "dosage_triggered_unrepaired_block_ids",
              "dosage_suspect_blocks", "dosage_underpowered_blocks", "dosage_orientation_tied_block_ids", "locus_repair_scope",
              "block_table", "genotypes_table", "main_fasta", "original_fasta", "original_block_table", "dosage_blocks_table",
              "dosage_representative_fasta", "dosage_equivalence_table", "call_equivalence_set"]
GENOTYPE_COLS = {"allele1", "allele2", "allele1_bp", "allele2_bp", "zygosity", "conditional_genotype_posterior", "conditional_GQ",
                 "allele1_sources", "allele2_sources"}
# the reported genotype's own annotation, recomputed on every block of a run whose genotype a repair changed (section 23)
FILTER_COLS = {"filter", "filter_reason", "representation_marker_margin", "representation_competitor_pair", "representation_fragment_diff",
               "representation_exposure_diff", "representation_likelihood_margin", "representation_conflict_statistic",
               "genotype_support", "quality_filter", "quality_reason"}
SUMMARY_KEYS = {"locus_filter", "lowqual_blocks", "not_evaluated_blocks", "representation_guard_complete", "quality_pass_blocks",
                "quality_lowqual_blocks", "quality_pass_block_ids", "quality_lowqual_block_ids", "quality_reason_representation_conflict",
                "quality_reason_representation_not_evaluated", "quality_reason_no_read_support", "blocks_genotype_read_supported"}
STATUS_VOLATILE = re.compile(r"^(seconds_|resources_|threads_|phys_footprint|rss_|emission_cache_dp_hits|emission_cache_dp_misses|"
                             r"emission_cache_executed_dps|emission_cache_executed_prefilters|scorer_emission_cache_arena_bytes|"
                             r"union_cached_prefilters_consumed|union_prefilter_cache_misses|union_prefilter_cache_decisions_not_retained|worker_seconds_|sha256_|profile|preset|"
                             r"settings_overridden|scorer_placement_capacity_bytes_peak)")
SKIP = (".eval.", ".gate_eval", ".evaluate.log", ".command.txt", ".log", ".time", ".done", ".repair_memory.tsv", ".tie_branch_memory.tsv")
REFERENCE = "HG00096#1#haplotype1-0000024:31848049-32074186"

failures = []
def check(ok, what):
    if not ok:
        failures.append(what)
        print("FAIL", what)
    return ok

def norm(path, prefix, threads_exempt=False):
    text = strip_headline(open(path).read().replace(prefix, "PREFIX"), path)
    if path.endswith(".mosaic_status.tsv"):
        lines = [l for l in text.split("\n") if not STATUS_VOLATILE.match(l) and "PREFIX" not in l]
        if threads_exempt: lines = [l for l in lines if not l.startswith("threads")]
        return "\n".join(lines)
    if path.endswith(".repair_rounds.tsv"): return "\n".join("\t".join(l.split("\t")[:28]) for l in text.split("\n"))
    if path.endswith(".tie_branches.tsv"): return "\n".join("\t".join(l.split("\t")[:14]) for l in text.split("\n"))
    return text

def table(path):
    rows = [l.split("\t") for l in strip_headline(open(path).read(), path).rstrip("\n").split("\n")]
    return rows[0], rows[1:]

def fasta(path):
    recs, cur = [], None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"): recs.append([line[1:], ""]); cur = recs[-1]
        elif cur is not None: cur[1] += line
    return recs

def files(prefix):
    d, b = os.path.dirname(prefix), os.path.basename(prefix)
    return sorted(f[len(b):] for f in os.listdir(d) if f.startswith(b + "."))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panvar", required=True); ap.add_argument("--jobs", required=True); ap.add_argument("--cell", required=True)
    ap.add_argument("--scratch", required=True); ap.add_argument("--probe"); ap.add_argument("--threads-sweep", action="store_true")
    ap.add_argument("--no-run", action="store_true", help="check the arms already under --scratch")
    ap.add_argument("--minimal", action="store_true", help="off, audit and repair only (the flag-absent arm equals off on the small cells)")
    a = ap.parse_args()
    job = next(l.rstrip("\n").split("\t") for l in open(a.jobs) if l.split("\t")[0] == a.cell)
    model, r1, r2, lib, cur, truth, sealed = job[6], job[7], job[8], job[9], job[10], job[11], job[-1]
    donor, scheme = job[5], job[4]   # the donor label, as the campaign evaluated it (MOSAIC for mosaic cells)
    # the profile the sealed run was configured by, as its status records it (fresh cells keep no command.txt)
    profile = next(l.rstrip("\n").split("\t")[1] for l in open(sealed + ".mosaic_status.tsv") if l.startswith("profile\t"))
    base = ["--profile", profile, "--gfa", model, "--reference-path", REFERENCE, "--reads1", r1, "--reads2", r2, "--library", lib, "--quiet"]
    if cur != "none": base += ["--panel-curation", cur]
    arms = {"absent": [], "off": ["--dosage-mode", "off"], "audit": ["--dosage-mode", "audit"], "repair": ["--dosage-mode", "repair"]}
    if a.minimal: del arms["absent"]
    threads = {k: 8 for k in arms}
    if a.threads_sweep:
        for t in (1, 2, 4):
            arms[f"repair_t{t}"] = ["--dosage-mode", "repair"]; threads[f"repair_t{t}"] = t
    P, rc = {}, {}
    for arm, extra in arms.items():
        d = os.path.join(a.scratch, arm); os.makedirs(d, exist_ok=True)
        P[arm] = os.path.join(d, "call")
        if a.no_run: rc[arm] = 0; continue
        with open(P[arm] + ".log", "w") as log:
            rc[arm] = subprocess.run([a.panvar, "genotype", *base, "--threads", str(threads[arm]), *extra, "--out-prefix", P[arm]],
                                     stdout=log, stderr=subprocess.STDOUT).returncode
    check(len(set(rc.values())) == 1, f"exit statuses differ: {rc}")

    # identity: absent / off against the sealed campaign, and against each other
    sealed_files = [s for s in files(sealed) if not any(x in s for x in SKIP)]
    for arm in [x for x in ("absent", "off") if x in P]:
        for s in sealed_files:
            check(os.path.exists(P[arm] + s) and norm(sealed + s, sealed) == norm(P[arm] + s, P[arm]), f"{arm}{s} differs from the campaign")
        extra = set(files(P[arm])) - set(files(sealed)) - {".log"}
        check(extra == {".locus_genotype.tsv"}, f"{arm} adds {sorted(extra)}")
    for s in (files(P["off"]) if "absent" in P else []):
        if s == ".log" or any(x in s for x in SKIP): continue
        check(norm(P["absent"] + s, P["absent"]) == norm(P["off"] + s, P["off"]), f"absent{s} != off{s}")
    # audit: every pre-existing output equal to off; only dosage files added
    for s in files(P["off"]):
        if s in (".log", ".locus_genotype.tsv") or any(x in s for x in SKIP): continue
        check(norm(P["audit"] + s, P["audit"]) == norm(P["off"] + s, P["off"]), f"audit{s} differs from off")
    extra = {s for s in set(files(P["audit"])) - set(files(P["off"])) if not any(x in s for x in SKIP)}
    check(all(s.startswith(".dosage_") for s in extra), f"audit adds non-dosage files {sorted(s for s in extra if not s.startswith('.dosage_'))}")

    # locus table: schema and agreement with the outputs, in every arm
    locus = {}
    for arm in ("off", "audit", "repair"):
        h, rows = table(P[arm] + ".locus_genotype.tsv")
        check(h == LOCUS_COLS, f"{arm} locus schema")
        check(len(rows) == 1, f"{arm} locus rows {len(rows)}")
        locus[arm] = {k: v.replace(P[arm], "PREFIX") for k, v in zip(h, rows[0])}
    L = locus["repair"]
    cls = L["call_class"]
    applied = set() if L["dosage_complete_applied_block_ids"] in ("none", "NA") else set(L["dosage_complete_applied_block_ids"].split(","))
    print(f"cell {a.cell}: {cls}, applied {sorted(applied) or 'none'}, content-corrected {L['dosage_content_corrected_block_ids']}, "
          f"triggered-unrepaired {L['dosage_triggered_unrepaired_block_ids']}, locus scope {L['locus_repair_scope']}")
    for arm in ("off", "audit", "repair"):
        l = locus[arm]
        for c in ("block_table", "genotypes_table", "main_fasta", "original_fasta", "original_block_table", "dosage_blocks_table",
                  "dosage_representative_fasta", "dosage_equivalence_table", "call_equivalence_set"):
            if l[c] != "NA": check(os.path.exists(os.path.join(os.path.dirname(P[arm]), l[c])), f"{arm} locus link {c} missing")
        if cls == "no_call":
            check(all(l[k] == "NA" for k in ("hap1_sha256", "hap2_sha256", "hap1_block_path", "main_fasta")), f"{arm} no-call carries a haplotype")
            continue
        bh, brows = table(P[arm] + ".block_genotypes.tsv")
        ug = ";".join(f"{r[bh.index('block_id')]}={min(int(r[bh.index('allele1')]), int(r[bh.index('allele2')]))}/"
                      f"{max(int(r[bh.index('allele1')]), int(r[bh.index('allele2')]))}" for r in brows)
        check(l["unordered_genotype"] == ug, f"{arm} locus unordered genotype disagrees with the block table")
        if cls == "strict_call":
            fa = fasta(P[arm] + ".mosaic_haplotypes.fa")
            for i, hap in enumerate(("hap1", "hap2")):
                check(l[f"{hap}_sha256"] == hashlib.sha256(fa[i][1].encode()).hexdigest(), f"{arm} {hap} digest != main FASTA")
                path = l[f"{hap}_block_path"].split(",")
                check([r[bh.index("allele1" if i == 0 else "allele2")] for r in brows] == path, f"{arm} {hap} path disagrees with the block table")
            if l["reported_genotype_source"] == "caller":
                eh, erows = table(P[arm] + ".equivalence_set.tsv")
                rep = next(r for r in erows if r[eh.index("rank")] == "1" and r[eh.index("representative")] == "true")
                check((l["hap1_block_path"], l["hap2_block_path"]) == (rep[eh.index("hap1_alleles")], rep[eh.index("hap2_alleles")]),
                      f"{arm} locus paths != the caller's representative")
                check(l["total_score"] == rep[eh.index("score")], f"{arm} locus total != the caller's")
        else:
            check(l["hap1_sha256"] == "NA" and l["hap1_block_path"] == "NA" and l["main_fasta"] == "NA", f"{arm} equivalence call claims phase")
    check(locus["off"]["dosage_mode"] == "off" and locus["off"]["locus_repair_scope"] == "NOT_RUN", "off locus dosage fields")
    for arm in ("off", "audit", "repair"):   # the locus filter and quality counts reconcile with the final block table
        l = locus[arm]
        if l["call_class"] == "no_call" or not os.path.exists(P[arm] + ".block_genotypes.tsv"): continue
        bh, brows = table(P[arm] + ".block_genotypes.tsv")
        f = [r[bh.index("filter")] for r in brows]; q = [r[bh.index("quality_filter")] for r in brows]
        want = "LOWQUAL_REPRESENTATION_CONFLICT" if "LOWQUAL_REPRESENTATION_CONFLICT" in f else "PASS" if all(x == "PASS" for x in f) else "NOT_EVALUATED"
        check(l["locus_filter"] == want, f"{arm} locus_filter {l['locus_filter']} does not reconcile with the block table ({want})")
        check(l["quality_pass_blocks"] == str(q.count("PASS")) and l["quality_lowqual_blocks"] == str(len(q) - q.count("PASS")),
              f"{arm} locus quality counts do not reconcile with the block table")
        if l["call_class"] == "strict_call":
            heads = [h for h, _ in fasta(P[arm] + ".mosaic_haplotypes.fa")]
            check(all(f"locus_filter={want}" in h for h in heads), f"{arm} FASTA locus_filter disagrees with the block table")
    check(L["filter_basis"] == ("caller" if not applied else L["filter_basis"]) and (L["filter_basis"] != "caller") == bool(applied),
          "repair locus filter_basis does not follow the applied repair")
    for k in LOCUS_COLS:
        if k.startswith("dosage_") or k in ("locus_repair_scope", "original_fasta", "original_block_table", "dosage_blocks_table",
                                            "dosage_representative_fasta", "dosage_equivalence_table"): continue
        check(locus["off"][k] == locus["audit"][k], f"audit locus {k} differs from off")
        if not applied and k != "filter_basis": check(locus["audit"][k] == L[k], f"repair locus {k} differs from audit without an applied repair")

    # repair: round trip, search record, status, block tables, FASTA, leakage
    for s in LEGACY:
        if os.path.exists(P["audit"] + s):
            check(norm(P["repair"] + ".original" + s, P["repair"]) == norm(P["audit"] + s, P["audit"]), f"repair.original{s} != audit{s}")
    for s in SEARCH_RECORD:
        if os.path.exists(P["audit"] + s):
            check(norm(P["repair"] + s, P["repair"]) == norm(P["audit"] + s, P["audit"]), f"repair{s} (search record) changed")
    for s in files(P["audit"]):
        if s.startswith(".dosage_") and s not in (".dosage_status.tsv", ".dosage_repair_scoring.tsv"):
            check(norm(P["repair"] + s, P["repair"]) == norm(P["audit"] + s, P["audit"]), f"repair{s} differs from audit")
    sa, sr = norm(P["audit"] + ".mosaic_status.tsv", P["audit"]).rstrip("\n"), norm(P["repair"] + ".mosaic_status.tsv", P["repair"]).rstrip("\n")
    if not applied:
        check(sr.startswith(sa) and all(x.startswith("dosage_") for x in sr[len(sa):].strip("\n").split("\n")), "repair status is not audit's plus dosage keys")
    else:
        al, rl = sa.split("\n"), sr.split("\n")
        check(len(rl) >= len(al), "repair status shorter than audit's")
        changed = {a_.split("\t")[0] for a_, r_ in zip(al, rl) if a_ != r_}
        check(changed <= SUMMARY_KEYS and all(r_.split("\t")[0] == a_.split("\t")[0] for a_, r_ in zip(al, rl)),
              f"repair status changed caller lines outside the filter summary: {sorted(changed - SUMMARY_KEYS)}")
        tail = rl[len(al):]
        check(all(x.startswith("dosage_") or x.startswith("original_") for x in tail if x), "repair status appends only dosage_ and original_ keys")
        for a_ in al:
            if a_.split("\t")[0] in SUMMARY_KEYS: check("original_" + a_ in tail, f"original value of {a_.split(chr(9))[0]} not kept")
    dosage_rows = {}
    if os.path.exists(P["repair"] + ".dosage_blocks.tsv"):
        dh, drows = table(P["repair"] + ".dosage_blocks.tsv")
        dosage_rows = {r[0]: dict(zip(dh, r)) for r in drows}
    for s in (".block_genotypes.tsv", ".genotypes.tsv"):
        if not os.path.exists(P["audit"] + s): continue
        ah, arows = table(P["audit"] + s); rh, rrows = table(P["repair"] + s)
        check(rh == ah + DOSAGE_COLS, f"repair{s} trailing schema")
        check(len(arows) == len(rrows), f"repair{s} row count")
        for ar, rr in zip(arows, rrows):
            bid = ar[0]
            lead = dict(zip(ah, rr[:len(ah)])); old = dict(zip(ah, ar)); tail = dict(zip(DOSAGE_COLS, rr[len(ah):]))
            changed = {k for k in ah if lead[k] != old[k]}
            if applied:
                check(tail["filter_basis"] != "caller" and tail["quality_filter_basis"] != "caller", f"repair{s} {bid}: filters still the caller's")
                check((tail["original_filter"], tail["original_quality_filter"]) == (old["filter"], old["quality_filter"]),
                      f"repair{s} {bid}: original filters not kept")
                if "genotype_support" in old: check(tail["original_genotype_support"] == old["genotype_support"], f"repair{s} {bid}: original support")
                if bid not in applied:
                    check(changed <= FILTER_COLS, f"repair{s} {bid}: columns beyond the filters changed {sorted(changed - FILTER_COLS)}")
            else:
                check(tail["filter_basis"] == "caller" and tail["quality_filter_basis"] == "caller" and tail["original_filter"] == old["filter"]
                      and tail["original_quality_filter"] == old["quality_filter"], f"repair{s} {bid}: filter provenance without a repair")
            if bid in applied:
                check(changed <= GENOTYPE_COLS | FILTER_COLS, f"repair{s} {bid}: non-genotype columns changed {sorted(changed - GENOTYPE_COLS - FILTER_COLS)}")
                comp = lead["representation_competitor_pair"]
                reported = "/".join(sorted([lead["allele1"], lead["allele2"]], key=int))
                check(comp != reported, f"repair{s} {bid}: the representation competitor is the reported pair itself (a stale filter)")
                check(tail["gq_status"] == "UNAVAILABLE_REPAIRED" and lead["conditional_GQ"] == "NA", f"repair{s} {bid}: GQ not withdrawn")
                check(tail["dosage_status"] == "REPAIRED" and tail["dosage_repair_applied"] == "true", f"repair{s} {bid}: status")
            elif not applied:
                check(not changed, f"repair{s} {bid}: legacy columns changed {sorted(changed)} without an applied repair")
            check((tail["original_allele1"], tail["original_allele2"]) == (old["allele1"], old["allele2"]), f"repair{s} {bid}: original alleles")
            if bid in dosage_rows:
                st = dosage_rows[bid]["dosage_status"]
                check(tail["dosage_status"] == st, f"repair{s} {bid}: dosage_status != dosage_blocks")
                want = {"CONSISTENT": "PASS", "NOT_SENSITIVE": "PASS", "CONTENT_REPAIRED": "DOSAGE_CONTENT_CORRECTED",
                        "TRIGGERED_UNREPAIRED": "DOSAGE_TRIGGERED_UNREPAIRED", "SUSPECT": "DOSAGE_SUSPECT", "UNDERPOWERED": "DOSAGE_UNDERPOWERED",
                        "NOT_EVALUATED": "DOSAGE_NOT_EVALUATED", "REPAIRED": "PASS" if bid in applied else "DOSAGE_REPAIR_NOT_APPLIED"}.get(st)
                check(tail["dosage_stage_filter"] == want, f"repair{s} {bid}: dosage_stage_filter {tail['dosage_stage_filter']} for {st}")
    if cls == "strict_call":
        fa_a, fa_r = fasta(P["audit"] + ".mosaic_haplotypes.fa"), fasta(P["repair"] + ".mosaic_haplotypes.fa")
        if not applied:
            check(norm(P["audit"] + ".mosaic_haplotypes.fa", P["audit"]) == norm(P["repair"] + ".mosaic_haplotypes.fa", P["repair"]), "repair FASTA changed without a repair")
        else:
            for (ha, sa_), (hr, sr_) in zip(fa_a, fa_r):
                ta = [t for t in ha.replace(P["audit"], "PREFIX").split() if not t.startswith(("bp=", "locus_filter="))]
                tr = [t for t in hr.replace(P["repair"], "PREFIX").split() if not t.startswith(("bp=", "locus_filter="))]
                check(tr[:len(ta)] == ta and all(t.startswith("dosage_") for t in tr[len(ta):]), "repaired FASTA header")
        reps = set()
        if os.path.exists(P["repair"] + ".dosage_repaired.fa"):
            reps = {h.split()[0] for h, _ in fasta(P["repair"] + ".dosage_repaired.fa") if "REPRESENTATIVE_ONLY_NOT_RESOLVED" in h}
        for h, seq in fa_r:
            check(hashlib.sha256(seq.encode()).hexdigest() not in reps, "a conditional representative is in the main FASTA")
        check(not ({L["hap1_sha256"], L["hap2_sha256"]} & reps), "a conditional representative is a locus digest")

    # evaluator: the originals score as audit does; the repaired prefix from the locus table only when applied
    if a.probe and truth != "none":
        ev = {}
        for arm, prefix in (("audit", P["audit"]), ("original", P["repair"] + ".original"), ("repaired", P["repair"])):
            out = prefix + ".gate_eval"
            r = subprocess.run([a.probe, "evaluate", "--scheme", scheme, "--donor", donor, "--reference", REFERENCE, "--run-prefix", prefix,
                                "--model-gfa", model, "--truth-gfa", truth, "--out", out] + (["--panel-curation", cur] if cur != "none" else []),
                               capture_output=True, text=True)
            check(r.returncode == 0, f"evaluate {arm}: {r.stderr.strip()[-200:]}")
            # the summary copies the run's wall times; those and the run path are the only lines two runs may differ in
            ev[arm] = (open(out + ".blocks.tsv").read(),
                       [l for l in open(out + ".summary.tsv") if not l.startswith("run\t") and not l.startswith("seconds_")])
        check(ev["original"] == ev["audit"], "the evaluator scores the originals differently from audit")
        has_source = any(l.startswith("call_source\t") for l in ev["repaired"][1])
        check(has_source == bool(applied), "the evaluator's call source does not follow the applied repair")
        if not applied: check(ev["repaired"] == ev["audit"], "the evaluator scores an unrepaired repair-mode run differently")

    if a.threads_sweep:
        for t in (1, 2, 4):
            arm = f"repair_t{t}"
            for s in files(P["repair"]):
                if s in (".log", ".dosage_status.tsv", ".dosage_repair_scoring.tsv", ".repair_memory.tsv", ".tie_branch_memory.tsv") or ".gate_eval" in s:
                    continue
                check(os.path.exists(P[arm] + s) and norm(P[arm] + s, P[arm], True) == norm(P["repair"] + s, P["repair"], True),
                      f"repair{s} differs at {t} threads")
            check(set(f for f in files(P[arm]) if ".gate_eval" not in f) == set(f for f in files(P["repair"]) if ".gate_eval" not in f),
                  f"file sets differ at {t} threads")
    print(("PASS " if not failures else "FAILED ") + a.cell + f" ({len(failures)} failures)")
    return 1 if failures else 0

if __name__ == "__main__":
    sys.exit(main())
