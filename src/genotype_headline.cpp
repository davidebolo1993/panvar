#include "panvar/genotype_headline.hpp"

#include "panvar/genotype_repair.hpp"   // atomic_write

#include <algorithm>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>

namespace panvar::headline {

const std::array<const char*, 5> kStatuses = {"NO_CALL", "NOT_EVALUATED", "LOW_QUALITY", "PASS", "NOT_APPLICABLE"};

int severity(const std::string& status) {
    if (status == "NO_CALL") return 4;
    if (status == "NOT_EVALUATED") return 3;
    if (status == "LOW_QUALITY") return 2;
    if (status == "PASS") return 1;
    if (status == "NOT_APPLICABLE") return 0;
    throw std::invalid_argument("unknown headline status " + status);
}

std::string Verdict::reason_text() const {
    std::string text;
    for (const std::string& r : reasons) text += (text.empty() ? "" : ",") + r;
    return text.empty() ? "none" : text;
}

Verdict combine(const std::vector<Verdict>& parts) {
    Verdict out;
    if (parts.empty()) return out;   // PASS, no reason
    out.status = parts.front().status;
    for (const Verdict& part : parts) {
        if (severity(part.status) > severity(out.status)) out.status = part.status;
        for (const std::string& r : part.reasons) {
            if (std::find(out.reasons.begin(), out.reasons.end(), r) == out.reasons.end()) out.reasons.push_back(r);
        }
    }
    return out;
}

namespace {
Verdict v(const std::string& status, const std::string& reason = "") {
    Verdict out;
    out.status = status;
    if (!reason.empty()) out.reasons.push_back(reason);
    return out;
}
bool no_call(const LocusFacts& locus) { return locus.call_class != "strict_call" && locus.call_class != "genotype_equivalent_call"; }
bool dosage_ran(const LocusFacts& locus) { return locus.dosage_mode != "off" && locus.dosage_audit_status == "COMPLETE"; }
} // namespace

Verdict genotype_verdict(const LocusFacts& locus, const BlockFacts& block) {
    if (no_call(locus)) return v("NO_CALL", "RUN_NO_CALL:" + locus.no_call_reason);
    std::vector<Verdict> parts;
    const std::string& rep = block.representation_filter;
    if (rep != "PASS" && rep.rfind("LOWQUAL", 0) != 0) parts.push_back(v("NOT_EVALUATED", "REPRESENTATION_NOT_EVALUATED"));
    if (locus.library_sensitivity == "NOT_EVALUATED" || block.library_not_evaluated) {
        parts.push_back(v("NOT_EVALUATED", "LIBRARY_SENSITIVITY_NOT_EVALUATED"));
    }
    if (rep == "LOWQUAL_REPRESENTATION_CONFLICT") parts.push_back(v("LOW_QUALITY", "REPRESENTATION_CONFLICT"));
    if (locus.library_sensitivity == "SENSITIVE" && block.library_sensitive) parts.push_back(v("LOW_QUALITY", "LIBRARY_SENSITIVE"));
    // The reported alleles are the caller's where the dosage stage found them dosage-inconsistent and did not resolve
    // the exact diplotype: a content correction, a triggered block left unrepaired, or a complete repair not applied.
    if (dosage_ran(locus) && (block.dosage_status == "CONTENT_REPAIRED" || block.dosage_status == "TRIGGERED_UNREPAIRED" ||
                              (block.dosage_status == "REPAIRED" && !block.dosage_repair_applied))) {
        parts.push_back(v("LOW_QUALITY", "EXACT_DIPLOTYPE_UNRESOLVED"));
    }
    return combine(parts);
}

Verdict dosage_verdict(const LocusFacts& locus, const BlockFacts& block) {
    if (no_call(locus)) return v("NO_CALL", "RUN_NO_CALL");
    if (locus.dosage_mode == "off") return v("NOT_EVALUATED", "DOSAGE_NOT_EVALUATED:dosage_mode_off");
    if (locus.dosage_audit_status != "COMPLETE") return v("NOT_EVALUATED", "DOSAGE_NOT_EVALUATED:audit_" + locus.dosage_audit_status);
    const std::string& s = block.dosage_status;
    if (s == "CONSISTENT") return v("PASS", "DOSAGE_CONSISTENT");
    if (s == "REPAIRED") {
        const bool complete = block.repair_class == "COMPLETE_DIPLOTYPE_REPAIR";
        if (complete && block.dosage_repair_applied) return v("PASS", "COMPLETE_REPAIR_APPLIED");
        // audit mode: the resolved repair is usable dosage content; only the selected mode kept it from the genotype
        if (complete && locus.dosage_mode == "audit" && block.repair_not_applied_reason == kAuditModeNotApplied) {
            return v("PASS", "COMPLETE_REPAIR_AVAILABLE_NOT_APPLIED");
        }
        // an eligibility, reproduction or safety guard refused it (repair mode, or repair mode's guards in audit mode)
        return v("LOW_QUALITY", "COMPLETE_REPAIR_NOT_APPLIED:" + block.repair_not_applied_reason);
    }
    if (s == "CONTENT_REPAIRED" && block.dosage_content_scope == "COMPLETE") return v("PASS", "CONTENT_CORRECTED_EXACT_DIPLOTYPE_UNRESOLVED");
    // Inconsistent-input guard, not a pipeline state: the dosage stage writes CONTENT_REPAIRED only for the class
    // CONDITIONAL_DOSAGE_CONTENT, which classify_repair_scope assigns only when the dosage content scope is COMPLETE
    // (genotype_dosage_repair.cpp). A CONTENT_REPAIRED row without that scope can come only from a hand-built or damaged
    // table; it fails closed rather than reading as PASS. No producer emits it.
    if (s == "CONTENT_REPAIRED") return v("LOW_QUALITY", "CONTENT_CORRECTION_CONDITIONAL");
    if (s == "SUSPECT") return v("LOW_QUALITY", "DOSAGE_SUSPECT");
    if (s == "TRIGGERED_UNREPAIRED") return v("LOW_QUALITY", "TRIGGERED_UNREPAIRED");
    if (s == "UNDERPOWERED") return v("NOT_EVALUATED", "DOSAGE_UNDERPOWERED");
    if (s == "NOT_SENSITIVE") return v("NOT_APPLICABLE", "NOT_DOSAGE_SENSITIVE");
    return v("NOT_EVALUATED", "DOSAGE_NOT_EVALUATED");
}

Verdict phase_boundary_verdict(const LocusFacts& locus, const BlockFacts& left, const BlockFacts& right) {
    if (no_call(locus)) return v("NO_CALL", "RUN_NO_CALL");
    if (left.zygosity != "HET" || right.zygosity != "HET") return v("NOT_APPLICABLE", "NO_PHASE_INFORMATIVE_BOUNDARY");
    Verdict base;
    const std::string& ps = locus.phase_support;
    if (locus.call_class == "genotype_equivalent_call" || ps == "UNRESOLVED") base = v("NO_CALL", "PHASE_UNRESOLVED");
    else if (ps == "READS_AUDITED" || ps == "READS_WITHIN_CANDIDATES") base = v("PASS");
    else if (ps == "PRIOR_ONLY") base = v("LOW_QUALITY", "MODEL_ONLY_PHASE");
    else base = v("NOT_EVALUATED", "PHASE_ASSESSMENT_INCOMPLETE");
    std::vector<Verdict> parts{base};
    if (left.dosage_repair_applied || right.dosage_repair_applied) parts.push_back(v("NOT_EVALUATED", "PHASE_NOT_ASSESSED_AFTER_DOSAGE_REPAIR"));
    return combine(parts);
}

Aggregate aggregate(const std::vector<std::pair<std::string, Verdict>>& units) {
    Aggregate out;
    for (const char* s : kStatuses) out.counts[s] = 0;
    out.units = units.size();
    bool any_applicable = false;
    for (const auto& unit : units) {
        ++out.counts[unit.second.status];
        if (unit.second.status == "NOT_APPLICABLE") continue;
        ++out.applicable;
        if (!any_applicable || severity(unit.second.status) > severity(out.status)) out.status = unit.second.status;
        any_applicable = true;
    }
    if (!any_applicable) out.status = "NOT_APPLICABLE";
    // REASON@ids in first-occurrence order; a NOT_APPLICABLE unit's reason is listed only when every unit is one, and a
    // PASS unit contributes only a named PASS reason.
    std::vector<std::pair<std::string, std::vector<std::string>>> by;
    for (const auto& unit : units) {
        if (unit.second.status == "NOT_APPLICABLE" && any_applicable) continue;
        for (const std::string& r : unit.second.reasons) {
            auto it = std::find_if(by.begin(), by.end(), [&](const auto& e) { return e.first == r; });
            if (it == by.end()) {
                by.push_back({r, {}});
                it = by.end() - 1;
            }
            it->second.push_back(unit.first);
        }
    }
    std::string text;
    for (const auto& [r, ids] : by) {
        std::string list;
        for (const std::string& id : ids) list += (list.empty() ? "" : ",") + id;
        text += (text.empty() ? "" : ";") + r + "@" + list;
    }
    out.reasons = text.empty() ? "none" : text;
    return out;
}

namespace {
struct Tsv {
    std::vector<std::string> header;
    std::vector<std::vector<std::string>> rows;
    int col(const std::string& name) const {
        const auto it = std::find(header.begin(), header.end(), name);
        return it == header.end() ? -1 : static_cast<int>(it - header.begin());
    }
    std::string cell(std::size_t row, const std::string& name, const std::string& fallback = "NA") const {
        const int c = col(name);
        return c < 0 || static_cast<std::size_t>(c) >= rows[row].size() ? fallback : rows[row][c];
    }
};
std::vector<std::string> split(const std::string& s, char sep) {
    std::vector<std::string> out;
    std::string cur;
    std::stringstream in(s);
    while (std::getline(in, cur, sep)) out.push_back(cur);
    if (!s.empty() && s.back() == sep) out.push_back("");
    return out;
}
bool read_tsv(const std::string& path, Tsv& t) {
    std::ifstream in(path);
    if (!in) return false;
    std::string line;
    bool first = true;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        if (first) { t.header = split(line, '\t'); first = false; }
        else t.rows.push_back(split(line, '\t'));
    }
    return !first;
}
std::map<std::string, std::string> read_status(const std::string& path) {
    std::map<std::string, std::string> out;
    std::ifstream in(path);
    std::string line;
    while (std::getline(in, line)) {
        const std::size_t tab = line.find('\t');
        if (tab != std::string::npos) out.emplace(line.substr(0, tab), line.substr(tab + 1));   // first occurrence wins
    }
    return out;
}
std::set<std::string> id_set(const std::map<std::string, std::string>& status, const std::string& key) {
    std::set<std::string> out;
    const auto found = status.find(key);
    if (found == status.end() || found->second == "none" || found->second.empty()) return out;
    for (const std::string& id : split(found->second, ',')) out.insert(id);
    return out;
}
// The four block columns appended (or replaced) in one block table, rows matched by block_id.
void rewrite_block_table(const std::string& path, const std::map<std::string, std::pair<Verdict, Verdict>>& by_block) {
    Tsv t;
    if (!read_tsv(path, t)) return;
    static const std::vector<std::string> added = {"genotype_filter", "genotype_filter_reasons", "dosage_filter", "dosage_filter_reasons"};
    std::vector<int> keep;
    for (std::size_t c = 0; c < t.header.size(); ++c) {
        if (std::find(added.begin(), added.end(), t.header[c]) == added.end()) keep.push_back(static_cast<int>(c));
    }
    const int id = t.col("block_id");
    if (id < 0) throw std::runtime_error(path + " has no block_id column");
    std::ostringstream out;
    for (std::size_t k = 0; k < keep.size(); ++k) out << (k ? "\t" : "") << t.header[keep[k]];
    for (const std::string& name : added) out << '\t' << name;
    out << '\n';
    for (const auto& row : t.rows) {
        for (std::size_t k = 0; k < keep.size(); ++k) out << (k ? "\t" : "") << (static_cast<std::size_t>(keep[k]) < row.size() ? row[keep[k]] : "");
        const auto found = by_block.find(row[id]);
        if (found == by_block.end()) throw std::runtime_error(path + ": no headline for block " + row[id]);
        out << '\t' << found->second.first.status << '\t' << found->second.first.reason_text() << '\t'
            << found->second.second.status << '\t' << found->second.second.reason_text() << '\n';
    }
    atomic_write(path, out.str());
}
} // namespace

LocusHeadline apply_headline(const std::string& out_prefix, const std::vector<std::string>& block_ids,
                             const std::string& dosage_mode, const std::string& dosage_audit_status,
                             const std::map<std::string, std::string>& repair_not_applied) {
    const std::map<std::string, std::string> status = read_status(out_prefix + ".mosaic_status.tsv");
    LocusFacts locus;
    const auto get = [&](const std::string& key, const std::string& fallback) {
        const auto found = status.find(key);
        return found == status.end() ? fallback : found->second;
    };
    // The call class; a development-view run (no adaptive repair) writes none, so it is derived from the run status:
    // COMPLETE is a strict call, CONDITIONAL_GENOTYPE_EQUIVALENT a genotype-equivalent one, anything else no genotype.
    const std::string run_status = get("status", "NA");
    locus.call_class = get("call_class", run_status == "COMPLETE"                          ? "strict_call"
                                         : run_status == "CONDITIONAL_GENOTYPE_EQUIVALENT" ? "genotype_equivalent_call"
                                                                                           : "no_call");
    if (locus.call_class == "no_call" && get("no_call_reason", "").empty() && run_status != "NA") {
        locus.no_call_reason = "status_" + run_status;
    }
    if (locus.no_call_reason == "NA") locus.no_call_reason = get("no_call_reason", "NA");
    locus.phase_support = get("phase_support", "NA");
    locus.dosage_mode = dosage_mode;
    locus.dosage_audit_status = dosage_audit_status;
    // Library sensitivity: NOT_APPLICABLE for a declared library; otherwise the audit's own status, and an auto run
    // without one never reads as evaluated.
    if (get("library_provenance", "").rfind("declared", 0) == 0) locus.library_sensitivity = "NOT_APPLICABLE";
    else locus.library_sensitivity = get("library_sensitivity_status", get("library_sensitivity", "NOT_EVALUATED"));
    const std::set<std::string> sensitive = id_set(status, "library_sensitive_blocks");
    const std::set<std::string> unaudited = id_set(status, "library_sensitivity_not_evaluated_blocks");

    // Block facts: the reported block table, and the dosage stage's per-block record where the table lacks it (audit).
    std::vector<BlockFacts> blocks;
    Tsv bg, db;
    const bool called = locus.call_class == "strict_call" || locus.call_class == "genotype_equivalent_call";
    if (called && read_tsv(out_prefix + ".block_genotypes.tsv", bg)) {
        std::map<std::string, std::size_t> dosage_row;
        if (read_tsv(out_prefix + ".dosage_blocks.tsv", db)) {
            for (std::size_t r = 0; r < db.rows.size(); ++r) dosage_row[db.cell(r, "block_id")] = r;
        }
        for (std::size_t r = 0; r < bg.rows.size(); ++r) {
            BlockFacts b;
            b.block_id = bg.cell(r, "block_id");
            b.zygosity = bg.cell(r, "zygosity");
            b.representation_filter = bg.cell(r, "filter");
            b.library_sensitive = sensitive.count(b.block_id) != 0;
            b.library_not_evaluated = unaudited.count(b.block_id) != 0;
            const auto d = dosage_row.find(b.block_id);
            const auto pick = [&](const std::string& name) {
                const std::string from_table = bg.cell(r, name, "");
                if (!from_table.empty()) return from_table;
                return d == dosage_row.end() ? std::string("NA") : db.cell(d->second, name);
            };
            b.dosage_status = pick("dosage_status");
            b.repair_class = pick("repair_class");
            b.dosage_content_scope = pick("dosage_content_scope");
            b.dosage_repair_applied = bg.cell(r, "dosage_repair_applied", "NA") == "true";
            const auto why = repair_not_applied.find(b.block_id);
            if (why != repair_not_applied.end()) b.repair_not_applied_reason = why->second;
            else if (bg.cell(r, "dosage_stage_filter", "") == "DOSAGE_REPAIR_NOT_APPLIED") {
                b.repair_not_applied_reason = bg.cell(r, "dosage_stage_filter_reason");
            }
            blocks.push_back(b);
        }
    }
    LocusHeadline out;
    std::vector<std::pair<std::string, Verdict>> g, d, p;
    std::map<std::string, std::pair<Verdict, Verdict>> by_block;
    if (blocks.empty()) {
        // No genotype: every block and boundary is NO_CALL.
        BlockFacts none;
        for (const std::string& id : block_ids) {
            g.push_back({id, genotype_verdict(locus, none)});
            d.push_back({id, dosage_verdict(locus, none)});
        }
        for (std::size_t i = 0; i + 1 < block_ids.size(); ++i) {
            p.push_back({block_ids[i] + ">" + block_ids[i + 1], phase_boundary_verdict(locus, none, none)});
        }
        // A run that stopped before its blocks were built (no model) is still a no-call: one locus-level unit per
        // filter, so the public summary reads NO_CALL with the run's reason, never three empty NOT_APPLICABLE values.
        if (no_call(locus)) {
            if (g.empty()) g.push_back({kLocusUnit, genotype_verdict(locus, none)});
            if (d.empty()) d.push_back({kLocusUnit, dosage_verdict(locus, none)});
            if (p.empty()) p.push_back({kLocusUnit, phase_boundary_verdict(locus, none, none)});
        }
    } else {
        for (const BlockFacts& b : blocks) {
            const Verdict gv = genotype_verdict(locus, b), dv = dosage_verdict(locus, b);
            g.push_back({b.block_id, gv});
            d.push_back({b.block_id, dv});
            by_block[b.block_id] = {gv, dv};
        }
        for (std::size_t i = 0; i + 1 < blocks.size(); ++i) {
            p.push_back({blocks[i].block_id + ">" + blocks[i + 1].block_id, phase_boundary_verdict(locus, blocks[i], blocks[i + 1])});
        }
        rewrite_block_table(out_prefix + ".block_genotypes.tsv", by_block);
        rewrite_block_table(out_prefix + ".genotypes.tsv", by_block);
    }
    out.genotype = aggregate(g);
    out.dosage = aggregate(d);
    out.phase = aggregate(p);
    out.call_summary = "GENOTYPE_" + out.genotype.status + ";DOSAGE_" + out.dosage.status + ";PHASE_" + out.phase.status;
    return out;
}

const std::vector<std::string>& locus_columns() {
    static const std::vector<std::string> columns = [] {
        std::vector<std::string> c = {"headline_rule", "genotype_filter", "genotype_filter_reasons", "dosage_filter",
                                      "dosage_filter_reasons", "phase_filter", "phase_filter_reasons", "call_summary"};
        for (const char* f : {"genotype", "dosage", "phase"}) {
            const std::string p = std::string(f) + "_filter_";
            for (const char* s : {"units", "applicable", "pass", "low_quality", "not_evaluated", "no_call", "not_applicable"}) {
                c.push_back(p + s);
            }
        }
        return c;
    }();
    return columns;
}

std::map<std::string, std::string> locus_values(const LocusHeadline& h) {
    std::map<std::string, std::string> v;
    v["headline_rule"] = kHeadlineRule;
    const std::pair<const char*, const Aggregate*> filters[] = {{"genotype", &h.genotype}, {"dosage", &h.dosage}, {"phase", &h.phase}};
    for (const auto& [name, a] : filters) {
        const std::string p = std::string(name) + "_filter";
        v[p] = a->status;
        v[p + "_reasons"] = a->reasons;
        v[p + "_units"] = std::to_string(a->units);
        v[p + "_applicable"] = std::to_string(a->applicable);
        v[p + "_pass"] = std::to_string(a->counts.at("PASS"));
        v[p + "_low_quality"] = std::to_string(a->counts.at("LOW_QUALITY"));
        v[p + "_not_evaluated"] = std::to_string(a->counts.at("NOT_EVALUATED"));
        v[p + "_no_call"] = std::to_string(a->counts.at("NO_CALL"));
        v[p + "_not_applicable"] = std::to_string(a->counts.at("NOT_APPLICABLE"));
    }
    v["call_summary"] = h.call_summary;
    return v;
}

} // namespace panvar::headline
