// Gates for the three headline filters (headline_filters_v2): each mapping row, the fail-closed precedence, locus
// aggregation (worst unit, NOT_APPLICABLE only when every unit is, REASON@ids) and the per-status counts, which must
// reconcile exactly with the units.

#include "panvar/genotype_headline.hpp"

#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace panvar::headline;

namespace {
void require(bool ok, const std::string& what) { if (!ok) throw std::runtime_error(what); }
LocusFacts strict(const std::string& phase = "PRIOR_ONLY", const std::string& mode = "repair") {
    LocusFacts l; l.call_class = "strict_call"; l.phase_support = phase; l.dosage_mode = mode; l.dosage_audit_status = "COMPLETE";
    l.library_sensitivity = "NOT_APPLICABLE"; return l;
}
BlockFacts block(const std::string& id, const std::string& dosage = "NOT_SENSITIVE", const std::string& zyg = "HET") {
    BlockFacts b; b.block_id = id; b.zygosity = zyg; b.representation_filter = "PASS"; b.dosage_status = dosage; return b;
}
void expect(const Verdict& got, const std::string& status, const std::string& reasons, const std::string& what) {
    require(got.status == status && got.reason_text() == reasons, what + ": got " + got.status + " / " + got.reason_text());
}
} // namespace

int main() {
    try {
        // precedence: NO_CALL > NOT_EVALUATED > LOW_QUALITY > PASS > NOT_APPLICABLE; every reason kept
        Verdict lq; lq.status = "LOW_QUALITY"; lq.reasons = {"A"};
        Verdict ne; ne.status = "NOT_EVALUATED"; ne.reasons = {"B"};
        expect(combine({lq, ne}), "NOT_EVALUATED", "A,B", "fail-closed precedence");
        // genotype mapping
        const LocusFacts L = strict();
        expect(genotype_verdict(L, block("sv:1")), "PASS", "none", "plain block");
        BlockFacts rc = block("sv:7"); rc.representation_filter = "LOWQUAL_REPRESENTATION_CONFLICT";
        expect(genotype_verdict(L, rc), "LOW_QUALITY", "REPRESENTATION_CONFLICT", "representation conflict");
        BlockFacts ne_rep = block("sv:7"); ne_rep.representation_filter = "NOT_EVALUATED";
        expect(genotype_verdict(L, ne_rep), "NOT_EVALUATED", "REPRESENTATION_NOT_EVALUATED", "representation not evaluated");
        BlockFacts content = block("sv:7", "CONTENT_REPAIRED"); content.repair_class = "CONDITIONAL_DOSAGE_CONTENT"; content.dosage_content_scope = "COMPLETE";
        expect(genotype_verdict(L, content), "LOW_QUALITY", "EXACT_DIPLOTYPE_UNRESOLVED", "content correction: genotype");
        expect(dosage_verdict(L, content), "PASS", "CONTENT_CORRECTED_EXACT_DIPLOTYPE_UNRESOLVED", "content correction: dosage");
        BlockFacts rep = block("sv:1", "REPAIRED"); rep.repair_class = "COMPLETE_DIPLOTYPE_REPAIR"; rep.dosage_repair_applied = true;
        expect(genotype_verdict(L, rep), "PASS", "none", "applied repair: genotype");
        expect(dosage_verdict(L, rep), "PASS", "COMPLETE_REPAIR_APPLIED", "applied repair: dosage");
        // repair mode, a guard refused the candidate: dosage LOW_QUALITY with the guard's reason
        BlockFacts unapplied = rep; unapplied.dosage_repair_applied = false; unapplied.repair_not_applied_reason = "not_a_strict_call";
        expect(genotype_verdict(L, unapplied), "LOW_QUALITY", "EXACT_DIPLOTYPE_UNRESOLVED", "refused complete repair: genotype");
        expect(dosage_verdict(L, unapplied), "LOW_QUALITY", "COMPLETE_REPAIR_NOT_APPLIED:not_a_strict_call", "refused complete repair: dosage");
        BlockFacts unrecorded = rep; unrecorded.dosage_repair_applied = false;
        expect(dosage_verdict(L, unrecorded), "LOW_QUALITY", "COMPLETE_REPAIR_NOT_APPLIED:not_recorded", "no recorded reason: never PASS");
        // audit mode, the repair is available and only the mode kept it out: usable dosage content, unresolved genotype
        const LocusFacts audit = strict("PRIOR_ONLY", "audit");
        BlockFacts available = unapplied; available.repair_not_applied_reason = kAuditModeNotApplied;
        expect(genotype_verdict(audit, available), "LOW_QUALITY", "EXACT_DIPLOTYPE_UNRESOLVED", "audit-available repair: genotype");
        expect(dosage_verdict(audit, available), "PASS", "COMPLETE_REPAIR_AVAILABLE_NOT_APPLIED", "audit-available repair: dosage");
        expect(dosage_verdict(L, available), "LOW_QUALITY", "COMPLETE_REPAIR_NOT_APPLIED:dosage_mode_audit",
               "the audit reason is PASS only in audit mode");
        expect(dosage_verdict(audit, unapplied), "LOW_QUALITY", "COMPLETE_REPAIR_NOT_APPLIED:not_a_strict_call", "audit mode, guard refused");
        expect(dosage_verdict(L, block("sv:7", "SUSPECT")), "LOW_QUALITY", "DOSAGE_SUSPECT", "suspect");
        expect(dosage_verdict(L, block("sv:7", "UNDERPOWERED")), "NOT_EVALUATED", "DOSAGE_UNDERPOWERED", "underpowered");
        expect(dosage_verdict(L, block("ctx:0")), "NOT_APPLICABLE", "NOT_DOSAGE_SENSITIVE", "not sensitive");
        expect(dosage_verdict(L, block("sv:7", "CONSISTENT")), "PASS", "DOSAGE_CONSISTENT", "consistent");
        BlockFacts tu = block("sv:7", "TRIGGERED_UNREPAIRED");
        expect(genotype_verdict(L, tu), "LOW_QUALITY", "EXACT_DIPLOTYPE_UNRESOLVED", "triggered unrepaired: genotype");
        expect(dosage_verdict(L, tu), "LOW_QUALITY", "TRIGGERED_UNREPAIRED", "triggered unrepaired: dosage");
        // dosage mode off: honestly not evaluated, never PASS; and a dosage status is then never used for the genotype
        const LocusFacts off = strict("PRIOR_ONLY", "off");
        expect(dosage_verdict(off, block("sv:7", "CONSISTENT")), "NOT_EVALUATED", "DOSAGE_NOT_EVALUATED:dosage_mode_off", "mode off");
        expect(genotype_verdict(off, content), "PASS", "none", "mode off: no dosage component in the genotype");
        // library sensitivity
        LocusFacts sens = strict(); sens.library_sensitivity = "SENSITIVE";
        BlockFacts ls = block("ctx:0"); ls.library_sensitive = true;
        expect(genotype_verdict(sens, ls), "LOW_QUALITY", "LIBRARY_SENSITIVE", "library sensitive block");
        expect(genotype_verdict(sens, block("ctx:2")), "PASS", "none", "an unaffected block of a sensitive locus");
        LocusFacts lne = strict(); lne.library_sensitivity = "NOT_EVALUATED";
        expect(genotype_verdict(lne, block("ctx:2")), "NOT_EVALUATED", "LIBRARY_SENSITIVITY_NOT_EVALUATED", "library not evaluated");
        // phase
        expect(phase_boundary_verdict(L, block("a"), block("b")), "LOW_QUALITY", "MODEL_ONLY_PHASE", "prior-only phase");
        expect(phase_boundary_verdict(strict("READS_WITHIN_CANDIDATES"), block("a"), block("b")), "PASS", "none", "read phase");
        expect(phase_boundary_verdict(strict("READS_AUDITED"), block("a"), block("b")), "PASS", "none", "audited phase");
        LocusFacts ge = strict("UNRESOLVED"); ge.call_class = "genotype_equivalent_call";
        expect(phase_boundary_verdict(ge, block("a"), block("b")), "NO_CALL", "PHASE_UNRESOLVED", "genotype-equivalent");
        expect(phase_boundary_verdict(ge, block("a", "NOT_SENSITIVE", "HOM"), block("b")), "NOT_APPLICABLE",
               "NO_PHASE_INFORMATIVE_BOUNDARY", "homozygous endpoint");
        expect(phase_boundary_verdict(strict("NA"), block("a"), block("b")), "NOT_EVALUATED", "PHASE_ASSESSMENT_INCOMPLETE", "incomplete");
        expect(phase_boundary_verdict(L, rep, block("b")), "NOT_EVALUATED", "MODEL_ONLY_PHASE,PHASE_NOT_ASSESSED_AFTER_DOSAGE_REPAIR",
               "next to an applied repair");
        // no-call
        LocusFacts nc; nc.no_call_reason = "tie_branching_incomplete";
        expect(genotype_verdict(nc, BlockFacts{}), "NO_CALL", "RUN_NO_CALL:tie_branching_incomplete", "no-call genotype");
        // aggregation: worst unit; NOT_APPLICABLE only when every unit is; REASON@ids; counts reconcile
        Verdict pass, na; na.status = "NOT_APPLICABLE"; na.reasons = {"NOT_DOSAGE_SENSITIVE"};
        Verdict under; under.status = "NOT_EVALUATED"; under.reasons = {"DOSAGE_UNDERPOWERED"};
        Verdict cons; cons.status = "PASS"; cons.reasons = {"DOSAGE_CONSISTENT"};
        const Aggregate a = aggregate({{"ctx:0", na}, {"sv:1", under}, {"sv:7", cons}, {"ctx:8", na}});
        require(a.status == "NOT_EVALUATED" && a.reasons == "DOSAGE_UNDERPOWERED@sv:1;DOSAGE_CONSISTENT@sv:7",
                "aggregate status and reasons: " + a.status + " / " + a.reasons);
        require(a.units == 4 && a.applicable == 2 && a.counts.at("NOT_APPLICABLE") == 2 && a.counts.at("PASS") == 1 &&
                a.counts.at("NOT_EVALUATED") == 1, "aggregate counts");
        std::size_t sum = 0;
        for (const auto& [s, n] : a.counts) sum += n;
        require(sum == a.units && a.units - a.counts.at("NOT_APPLICABLE") == a.applicable, "counts reconcile");
        const Aggregate all_na = aggregate({{"ctx:0", na}, {"ctx:2", na}});
        require(all_na.status == "NOT_APPLICABLE" && all_na.reasons == "NOT_DOSAGE_SENSITIVE@ctx:0,ctx:2", "all not applicable");
        const Aggregate lq_ne = aggregate({{"sv:1", lq}, {"sv:7", ne}});
        require(lq_ne.status == "NOT_EVALUATED", "locus precedence is fail-closed");
        // file level: a development-view run's status has no call_class; a COMPLETE run is still a strict call
        {
            const std::string dir = std::string(std::getenv("TMPDIR") ? std::getenv("TMPDIR") : "/tmp") + "/panvar_headline_test";
            std::filesystem::create_directories(dir);
            const std::string p = dir + "/run";
            std::ofstream(p + ".mosaic_status.tsv") << "key\tvalue\nstatus\tCOMPLETE\nphase_support\tPRIOR_ONLY\nlibrary_provenance\tdeclared: x\n";
            std::ofstream(p + ".block_genotypes.tsv") << "block_id\tzygosity\tfilter\nctx:0\tHET\tPASS\nsv:1\tHET\tNOT_EVALUATED\n";
            std::remove((p + ".genotypes.tsv").c_str());
            const LocusHeadline h = apply_headline(p, {"ctx:0", "sv:1"}, "off", "NA");
            require(h.genotype.status == "NOT_EVALUATED" && h.genotype.counts.at("PASS") == 1 && h.genotype.counts.at("NOT_EVALUATED") == 1 &&
                        h.phase.status == "LOW_QUALITY" && h.dosage.status == "NOT_EVALUATED",
                    "development-view run: headline " + h.call_summary);
            std::ifstream in(p + ".block_genotypes.tsv"); std::string header; std::getline(in, header);
            require(header.find("\tgenotype_filter\tgenotype_filter_reasons\tdosage_filter\tdosage_filter_reasons") != std::string::npos,
                    "block columns appended");
            std::ofstream(p + ".mosaic_status.tsv") << "key\tvalue\nstatus\tREPAIR_INCOMPLETE\ncall_class\tno_call\nno_call_reason\ttie_branching_incomplete\n";
            const LocusHeadline n = apply_headline(p, {"ctx:0", "sv:1"}, "repair", "NA");
            require(n.call_summary == "GENOTYPE_NO_CALL;DOSAGE_NO_CALL;PHASE_NO_CALL" && n.genotype.units == 2 && n.phase.units == 1,
                    "no-call run: " + n.call_summary);
        }
    } catch (const std::exception& e) {
        std::cerr << "FAIL  " << e.what() << '\n';
        return 1;
    }
    std::cout << "PASS  headline filters: every mapping row, fail-closed precedence, dosage-mode off, library, phase, no-call,"
                 " aggregation and reconciled counts\n";
    return 0;
}
