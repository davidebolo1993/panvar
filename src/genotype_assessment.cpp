#include "panvar/genotype_assessment.hpp"
#include "panvar/genotype_output_layout.hpp"

#include "panvar/genotype_repair.hpp"   // atomic_write

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace panvar::assessment {

namespace {

using Row = std::map<std::string, std::string>;

const std::vector<std::string> kFindingOrder = {"REPRESENTATION_CONFLICT", "LIBRARY_SENSITIVE", "DOSAGE_DISAGREES", "DOSAGE_SUSPECT",
                                                "PRIOR_DEPENDENT", "LABELS_AMBIGUOUS", "GENOTYPE_UNRESOLVED"};

std::vector<std::string> split(const std::string& s, char sep) {
    std::vector<std::string> out;
    std::string cur;
    for (char c : s) {
        if (c == sep) {
            out.push_back(cur);
            cur.clear();
        } else {
            cur += c;
        }
    }
    out.push_back(cur);
    return out;
}

std::string strip_cr(std::string s) {
    if (!s.empty() && s.back() == '\r') s.pop_back();
    return s;
}

// key<TAB>value, the first occurrence of a key winning (as assessment_spec.kv)
Row read_kv(const std::string& path) {
    Row out;
    std::ifstream in(path);
    for (std::string line; std::getline(in, line);) {
        line = strip_cr(line);
        const std::size_t tab = line.find('\t');
        if (tab == std::string::npos) continue;
        out.emplace(line.substr(0, tab), line.substr(tab + 1));
    }
    return out;
}

struct Table {
    std::vector<std::string> header;
    std::vector<std::vector<std::string>> rows;
};
bool read_table(const std::string& path, Table& t) {
    std::ifstream in(path);
    if (!in) return false;
    bool first = true;
    for (std::string line; std::getline(in, line);) {
        line = strip_cr(line);
        if (first) {
            t.header = split(line, '\t');
            first = false;
            continue;
        }
        if (line.empty()) continue;
        t.rows.push_back(split(line, '\t'));
    }
    return !first;
}
std::vector<Row> read_rows(const std::string& path) {
    Table t;
    std::vector<Row> out;
    if (!read_table(path, t)) return out;
    for (const auto& r : t.rows) {
        Row row;
        for (std::size_t c = 0; c < t.header.size(); ++c) row[t.header[c]] = c < r.size() ? r[c] : "";
        out.push_back(row);
    }
    return out;
}

// dict.get(key, fallback)
std::string get(const Row& m, const std::string& key, const std::string& fallback = "") {
    const auto it = m.find(key);
    return it == m.end() ? fallback : it->second;
}
// `a or b` for strings: the first non-empty
std::string or_(const std::string& a, const std::string& b) { return a.empty() ? b : a; }
bool starts_with(const std::string& s, const std::string& p) { return s.rfind(p, 0) == 0; }
bool contains(const std::string& s, const std::string& p) { return s.find(p) != std::string::npos; }

struct Facts {
    std::string prediction, completion, sequence_status, call_class, no_call_reason, phase_support, library, library_incomplete,
        dosage_mode, audit, audit_note;
    std::set<std::string> sensitive, unaudited;
    std::vector<Row> v2_blocks;
    std::map<std::string, Row> legacy, dosage;
};

std::set<std::string> id_set(const Row& st, const std::string& key) {
    std::set<std::string> out;
    const std::string v = get(st, key, "none");
    if (v == "none" || v.empty()) return out;
    for (const std::string& id : split(v, ',')) out.insert(id);
    return out;
}

Facts read_facts(const std::string& prefix) {
    const Row st = read_kv(prefix + ".mosaic_status.tsv"), sl = read_kv(genotype_output_path(prefix, ".sequence_locus.tsv"));
    Facts f;
    const std::string run_status = get(st, "status", "NA");
    f.call_class = or_(get(st, "call_class"), run_status == "COMPLETE" ? "strict_call"
                                              : run_status == "CONDITIONAL_GENOTYPE_EQUIVALENT" ? "genotype_equivalent_call"
                                                                                                : "no_call");
    if (f.call_class == "no_call" && !get(sl, "legacy_call_class").empty()) f.call_class = get(sl, "legacy_call_class");
    f.no_call_reason = or_(or_(get(st, "no_call_reason"), get(sl, "legacy_no_call_reason")), "status_" + run_status);
    const bool declared = starts_with(get(st, "library_provenance"), "declared");
    f.library = declared ? "NOT_APPLICABLE"
                         : get(st, "library_sensitivity_status", get(st, "library_sensitivity", "NOT_EVALUATED"));
    f.library_incomplete = or_(get(st, "library_sensitivity_incomplete", "not_evaluated"), "not_evaluated");
    f.sensitive = id_set(st, "library_sensitive_blocks");
    f.unaudited = id_set(st, "library_sensitivity_not_evaluated_blocks");
    f.dosage_mode = get(sl, "dosage_mode", "off");
    if (f.dosage_mode.empty() || f.dosage_mode == "NA" || f.dosage_mode == "not_run") f.dosage_mode = "off";
    const Row ds = read_kv(prefix + ".dosage_status.tsv");
    f.audit = f.dosage_mode != "off" ? get(ds, "dosage_audit_status", "NA") : "NA";
    f.audit_note = get(ds, "dosage_audit_note", "none");
    f.prediction = get(sl, "prediction", "NO_PREDICTION");
    f.completion = get(sl, "completion", "NA");
    f.sequence_status = get(sl, "sequence_status", "NA");
    f.phase_support = get(st, "phase_support", get(sl, "legacy_phase_support", "NA"));
    f.v2_blocks = read_rows(genotype_output_path(prefix, ".sequence_blocks.tsv"));
    for (const Row& r : read_rows(prefix + ".block_genotypes.tsv")) f.legacy.emplace(get(r, "block_id"), r);
    for (const Row& r : read_rows(prefix + ".dosage_blocks.tsv")) f.dosage.emplace(get(r, "block_id"), r);
    return f;
}

std::string clean(const std::string& text) {
    std::string out;
    for (char c : text) out += (c == ';' || c == '@' || c == '\t' || c == '\n' || c == '(' || c == ')') ? '_' : c;
    return out;
}

bool parse_int(const std::string& s, long& v) {
    if (s.empty()) return false;
    char* end = nullptr;
    v = std::strtol(s.c_str(), &end, 10);
    return end != nullptr && *end == '\0';
}
std::string sorted_pair(const std::string& a, const std::string& b) {
    long x = 0, y = 0;
    if (!parse_int(a, x) || !parse_int(b, y)) return "NA";
    if (y < x) std::swap(x, y);
    return std::to_string(x) + "/" + std::to_string(y);
}
std::string zygosity(const std::string& g) {
    const std::size_t slash = g.find('/');
    if (g.empty() || g == "NA" || slash == std::string::npos) return "UNKNOWN";
    return g.substr(0, slash) == g.substr(slash + 1) ? "HOM" : "HET";
}

std::string unverified_reason(const Facts& f, const Row& v2, const Row* legacy) {
    if (f.prediction == "PROVISIONAL") return "provisional";
    if (f.prediction != "PREDICTION") return "no_prediction";
    if (get(v2, "block_status") != "IDENTIFIED" || zygosity(get(v2, "genotype", "NA")) == "UNKNOWN") return "genotype_not_identified";
    if (f.call_class != "strict_call" && f.call_class != "genotype_equivalent_call") return "legacy_decision_" + f.no_call_reason;
    if (legacy == nullptr) return "no_legacy_block";
    if (sorted_pair(get(*legacy, "allele1"), get(*legacy, "allele2")) != get(v2, "genotype")) return "legacy_block_differs";
    return "";
}

std::string repair_not_applied(const Row* legacy) {
    const std::string rec = legacy ? get(*legacy, "dosage_filter_reasons") : "";
    if (contains(rec, "COMPLETE_REPAIR_AVAILABLE_NOT_APPLIED")) return "dosage_mode_audit";
    if (starts_with(rec, "COMPLETE_REPAIR_NOT_APPLIED:")) return rec.substr(std::string("COMPLETE_REPAIR_NOT_APPLIED:").size());
    if (legacy && !legacy->empty() && get(*legacy, "dosage_stage_filter") == "DOSAGE_REPAIR_NOT_APPLIED") {
        return get(*legacy, "dosage_stage_filter_reason", "not_recorded");
    }
    return "not_recorded";
}

struct Dosage { std::string check, reason, alternative; };
Dosage block_dosage(const Facts& f, const std::string& why, const Row* legacy, const Row& drow) {
    if (f.dosage_mode == "off") return {"NOT_REQUESTED", "dosage_check_not_requested", "NA"};
    if (f.audit != "COMPLETE") {
        const std::string r = f.audit == "FAILED" ? "AUDIT_FAILED" : f.audit == "NOT_RUN" ? "AUDIT_NOT_RUN"
                              : f.audit == "REFUSED_LIBRARY" ? "AUDIT_UNSUPPORTED_LIBRARY" : "AUDIT_" + f.audit;
        return {"UNAVAILABLE", r, "NA"};
    }
    if (!why.empty()) return {"UNAVAILABLE", "NOT_VERIFIED:" + why, "NA"};
    const auto pick = [&](const std::string& k) {
        if (legacy && !legacy->empty() && !get(*legacy, k).empty()) return get(*legacy, k);
        return get(drow, k, "NA");
    };
    const std::string s = pick("dosage_status");
    const bool applied = legacy != nullptr && get(*legacy, "dosage_repair_applied", "NA") == "true";
    if (s == "NOT_SENSITIVE") return {"NOT_APPLICABLE", "NOT_DOSAGE_SENSITIVE", "NONE"};
    if (s == "CONSISTENT") return {"AGREES", "CONSISTENT", "NONE"};
    if (s == "SUSPECT") return {"SUSPECT", "SUSPECT", "NONE"};
    if (s == "UNDERPOWERED") return {"UNAVAILABLE", "UNDERPOWERED", "NONE"};
    if (s == "REPAIRED" && applied) return {"AGREES", "REPORTED_IS_DOSAGE_REPAIR", "APPLIED"};
    if (s == "REPAIRED") {
        const std::string r = repair_not_applied(legacy);
        return {"DISAGREES", "TRIGGERED", r == "dosage_mode_audit" ? "COMPLETE_REPAIR_AVAILABLE" : "COMPLETE_REPAIR_REFUSED:" + r};
    }
    if (s == "CONTENT_REPAIRED") return {"DISAGREES", "TRIGGERED", "CONTENT_CORRECTION_AVAILABLE"};
    if (s == "TRIGGERED_UNREPAIRED" || s == "AUTO_REPAIR") return {"DISAGREES", "TRIGGERED", "REPAIR_NOT_FOUND:" + or_(pick("dosage_reason"), "none")};
    return {"UNAVAILABLE", "NOT_EVALUATED:" + or_(pick("dosage_reason"), "NA"), "NONE"};
}

using Notes = std::vector<std::pair<std::string, std::vector<std::string>>>;
void note(Notes& list, const std::string& code, const std::string& id) {
    auto it = std::find_if(list.begin(), list.end(), [&](const auto& e) { return e.first == code; });
    if (it == list.end()) {
        list.push_back({code, {}});
        it = list.end() - 1;
    }
    if (std::find(it->second.begin(), it->second.end(), id) == it->second.end()) it->second.push_back(id);
}
std::string join(const std::vector<std::string>& v, const std::string& sep) {
    std::string out;
    for (std::size_t i = 0; i < v.size(); ++i) out += (i ? sep : "") + v[i];
    return out;
}
int count(const std::vector<std::string>& v, std::initializer_list<const char*> values) {
    int n = 0;
    for (const std::string& x : v) {
        for (const char* y : values) n += x == y;
    }
    return n;
}
bool has(const std::vector<std::string>& v, const std::string& x) { return std::find(v.begin(), v.end(), x) != v.end(); }

std::string lower(std::string s) {
    for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

View assess(const Facts& f) {
    const std::string& pred = f.prediction;
    View view;
    Notes findings, missing;
    std::vector<std::string> all_ids;
    for (const Row& r : f.v2_blocks) all_ids.push_back(get(r, "block_id"));
    std::vector<std::string> gvals, dvals, svals, alts, zyg;
    bool verified_all = true;
    if (pred != "NO_PREDICTION") {
        for (const Row& v2 : f.v2_blocks) {
            const std::string bid = get(v2, "block_id");
            const auto lit = f.legacy.find(bid);
            const Row* legacy = lit == f.legacy.end() ? nullptr : &lit->second;
            const std::string why = unverified_reason(f, v2, legacy);
            verified_all = verified_all && why.empty();
            zyg.push_back(zygosity(get(v2, "genotype", "NA")));
            const std::string status = get(v2, "block_status", "NA");
            const std::string res = status == "REPRESENTATION_AMBIGUOUS" ? "LABELS_AMBIGUOUS" : status;
            std::vector<std::string> codes, reasons;
            std::string g;
            if (!why.empty()) {
                g = "NOT_ASSESSED";
                reasons.push_back(why);
                note(missing, "GENOTYPE_CHECK:" + why, bid);
            } else {
                const std::string rep = get(*legacy, "filter", "NA");
                std::vector<std::string> found;
                if (rep == "LOWQUAL_REPRESENTATION_CONFLICT") found.push_back("REPRESENTATION_CONFLICT");
                const bool sens = f.library == "SENSITIVE" && f.sensitive.count(bid) != 0;
                if (sens) found.push_back("LIBRARY_SENSITIVE");
                if (rep != "PASS" && rep != "LOWQUAL_REPRESENTATION_CONFLICT") {
                    reasons.push_back("representation_guard_" + rep);
                    note(missing, "GENOTYPE_CHECK:representation_guard_" + rep, bid);
                }
                if ((f.library == "NOT_EVALUATED" || f.unaudited.count(bid) != 0) && !sens) {
                    note(missing, "LIBRARY_SENSITIVITY:" + f.library_incomplete, bid);
                }
                std::vector<std::string> merged = found;
                merged.insert(merged.end(), reasons.begin(), reasons.end());
                reasons = merged;
                codes.insert(codes.end(), found.begin(), found.end());
                g = !found.empty() ? "CONFLICT"
                    : (rep == "PASS" || rep == "LOWQUAL_REPRESENTATION_CONFLICT") ? "NONE_DETECTED" : "NOT_ASSESSED";
            }
            gvals.push_back(g);
            const auto dit = f.dosage.find(bid);
            const Dosage d = block_dosage(f, why, legacy, dit == f.dosage.end() ? Row{} : dit->second);
            dvals.push_back(d.check);
            if (d.check == "DISAGREES" || d.check == "SUSPECT") codes.push_back("DOSAGE_" + d.check);
            if (d.check == "UNAVAILABLE") {
                if (starts_with(d.reason, "AUDIT_")) {
                    note(missing, "DOSAGE:" + d.reason + (d.reason == "AUDIT_FAILED" ? "(" + clean(f.audit_note) + ")" : ""), "locus");
                } else {
                    note(missing, "DOSAGE:" + d.reason, bid);
                }
            }
            if (d.alternative != "NA" && d.alternative != "NONE") alts.push_back(d.alternative + "@" + bid);
            const std::string support = (legacy && why.empty()) ? get(*legacy, "genotype_support", "NA") : "NA";
            std::string s, sr;
            if (!why.empty()) { s = "NOT_ASSESSED"; sr = why; }
            else if (support == "READS") { s = "READS"; sr = "reads"; }
            else if (support == "PRIOR_ONLY") { s = "PRIOR_DEPENDENT"; sr = "prior_only"; codes.push_back("PRIOR_DEPENDENT"); }
            else { s = "NOT_ASSESSED"; sr = "support_" + support; }
            if (s == "NOT_ASSESSED") note(missing, "READ_SUPPORT:" + sr, bid);
            svals.push_back(s);
            if (pred == "PREDICTION" && res == "LABELS_AMBIGUOUS") codes.push_back("LABELS_AMBIGUOUS");
            if (pred == "PREDICTION" && res == "GENOTYPE_UNRESOLVED") codes.push_back("GENOTYPE_UNRESOLVED");
            for (const std::string& c : kFindingOrder) {
                if (has(codes, c)) note(findings, c, bid);
            }
            view.blocks.push_back({bid, {res, g, reasons.empty() ? "none" : join(reasons, ","), d.check, d.reason, d.alternative, s, sr}});
        }
    }
    Row L;
    const std::size_t n = view.blocks.size();
    if (pred == "NO_PREDICTION") {
        for (const char* k : {"sequence_resolution", "genotype_check", "genotype_check_coverage", "dosage_check", "dosage_check_coverage",
                              "read_support", "read_support_coverage"}) {
            L[k] = "NA";
        }
    } else {
        const int ga = count(gvals, {"CONFLICT", "NONE_DETECTED"});
        L["genotype_check"] = has(gvals, "CONFLICT") ? "CONFLICT" : ga ? "NONE_DETECTED" : "NOT_ASSESSED";
        L["genotype_check_coverage"] = std::to_string(ga) + "/" + std::to_string(n);
        if (f.dosage_mode == "off") {
            L["dosage_check"] = "NOT_REQUESTED";
            L["dosage_check_coverage"] = "NA";
        } else {
            std::vector<std::string> applicable;
            for (const std::string& x : dvals) if (x != "NOT_APPLICABLE") applicable.push_back(x);
            const int da = count(applicable, {"AGREES", "SUSPECT", "DISAGREES"});
            std::string d;
            for (const char* x : {"DISAGREES", "SUSPECT", "AGREES"}) {
                if (d.empty() && has(dvals, x)) d = x;
            }
            if (d.empty()) d = (!dvals.empty() && applicable.empty()) ? "NOT_APPLICABLE" : "UNAVAILABLE";
            L["dosage_check"] = d;
            L["dosage_check_coverage"] = std::to_string(da) + "/" + std::to_string(applicable.size());
        }
        const int sa = count(svals, {"READS", "PRIOR_DEPENDENT"});
        L["read_support"] = has(svals, "PRIOR_DEPENDENT") ? "PRIOR_DEPENDENT" : sa ? "READS" : "NOT_ASSESSED";
        L["read_support_coverage"] = std::to_string(sa) + "/" + std::to_string(n);
        L["sequence_resolution"] = pred == "PREDICTION" ? f.sequence_status : "NOT_DETERMINED";
    }
    const int het = static_cast<int>(std::count(zyg.begin(), zyg.end(), "HET"));
    const int unknown = static_cast<int>(std::count(zyg.begin(), zyg.end(), "UNKNOWN"));
    std::string p, pr;
    if (pred == "NO_PREDICTION") { p = "NA"; pr = "no_prediction"; }
    else if (pred == "PROVISIONAL") { p = "NOT_ASSESSED"; pr = "provisional"; }
    else if (f.sequence_status == "GENOTYPE_RESOLVED_SEQUENCE_SET" || f.sequence_status == "SEQUENCE_TIED") { p = "UNRESOLVED"; pr = "several_sequence_pairs"; }
    else if (unknown == 0 && het < 2) { p = "NOT_APPLICABLE"; pr = "fewer_than_two_heterozygous_blocks"; }
    else if (f.call_class != "strict_call") {
        p = "NOT_ASSESSED";
        pr = "legacy_decision_" + (f.call_class == "genotype_equivalent_call" ? std::string("genotype_equivalent") : f.no_call_reason);
    }
    else if (!verified_all) { p = "NOT_ASSESSED"; pr = "blocks_not_verified"; }
    else if (f.phase_support == "READS_AUDITED" || f.phase_support == "READS_WITHIN_CANDIDATES") { p = "READ_SUPPORTED"; pr = lower(f.phase_support); }
    else if (f.phase_support == "PRIOR_ONLY") { p = "PRIOR_ONLY"; pr = "prior_only"; }
    else { p = "NOT_ASSESSED"; pr = "phase_support_" + f.phase_support; }
    if (p == "PRIOR_ONLY" || p == "UNRESOLVED") note(findings, "PHASE_" + p, "locus");
    if (p == "NOT_ASSESSED") note(missing, "PHASE:" + pr, "locus");
    const auto fmt = [&](const Notes& list) {
        std::vector<std::string> parts;
        for (const auto& [code, ids] : list) {
            const std::string idtext = (all_ids.size() > 1 && ids == all_ids) ? "all" : join(ids, ",");
            parts.push_back(code + "@" + idtext);
        }
        return parts.empty() ? std::string("none") : join(parts, ";");
    };
    const auto cov = [&](const std::string& key) {
        const std::string& v = L[key];
        const std::size_t slash = v.find('/');
        if (v == "NA" || slash == std::string::npos) return std::string();
        const long a = std::strtol(v.substr(0, slash).c_str(), nullptr, 10), b = std::strtol(v.substr(slash + 1).c_str(), nullptr, 10);
        return (0 < a && a < b) ? "[" + v + "]" : std::string();
    };
    L["assessment_rule"] = kRule;
    L["phase"] = p;
    L["phase_reason"] = pr;
    L["findings"] = fmt(findings);
    L["not_assessed"] = fmt(missing);
    L["dosage_alternatives"] = alts.empty() ? "none" : join(alts, ";");
    L["assessment_summary"] = pred + ";SEQUENCE_" + L["sequence_resolution"] + ";GENOTYPE_" + L["genotype_check"] + cov("genotype_check_coverage") +
                              ";DOSAGE_" + L["dosage_check"] + cov("dosage_check_coverage") + ";SUPPORT_" + L["read_support"] +
                              cov("read_support_coverage") + ";PHASE_" + p;
    for (const std::string& k : locus_keys()) view.locus.push_back({k, L[k]});
    return view;
}

// Replaces (or appends) the named columns of a TSV table, matching rows by block_id.
void rewrite_columns(const std::string& path, const std::vector<std::string>& names,
                     const std::map<std::string, std::vector<std::string>>& by_id, const std::string& key_column) {
    Table t;
    if (!read_table(path, t)) return;
    std::vector<std::size_t> keep;
    for (std::size_t c = 0; c < t.header.size(); ++c) {
        if (std::find(names.begin(), names.end(), t.header[c]) == names.end()) keep.push_back(c);
    }
    const auto id_it = std::find(t.header.begin(), t.header.end(), key_column);
    const std::size_t id_col = static_cast<std::size_t>(id_it - t.header.begin());
    std::ostringstream out;
    for (std::size_t k = 0; k < keep.size(); ++k) out << (k ? "\t" : "") << t.header[keep[k]];
    for (const std::string& n : names) out << '\t' << n;
    out << '\n';
    for (const auto& row : t.rows) {
        for (std::size_t k = 0; k < keep.size(); ++k) out << (k ? "\t" : "") << (keep[k] < row.size() ? row[keep[k]] : "");
        const std::string id = (id_it != t.header.end() && id_col < row.size()) ? row[id_col] : "";
        const auto found = by_id.find(id);
        for (std::size_t c = 0; c < names.size(); ++c) out << '\t' << (found == by_id.end() ? "NA" : found->second[c]);
        out << '\n';
    }
    atomic_write(path, out.str());
}

}  // namespace

const std::vector<std::string>& locus_keys() {
    static const std::vector<std::string> keys = {"assessment_rule", "sequence_resolution", "genotype_check", "genotype_check_coverage",
                                                  "dosage_check", "dosage_check_coverage", "read_support", "read_support_coverage",
                                                  "phase", "phase_reason", "findings", "not_assessed", "dosage_alternatives",
                                                  "assessment_summary"};
    return keys;
}

const std::vector<std::string>& block_columns() {
    static const std::vector<std::string> cols = {"block_resolution", "genotype_check", "genotype_check_reasons", "dosage_check",
                                                  "dosage_check_reason", "dosage_alternative", "read_support", "read_support_reason"};
    return cols;
}

View assess_prefix(const std::string& prefix) { return assess(read_facts(prefix)); }

std::string render(const View& view) {
    std::ostringstream out;
    out << "#locus\n";
    for (const auto& [k, v] : view.locus) out << k << '\t' << v << '\n';
    out << "#blocks\nblock_id";
    for (const std::string& c : block_columns()) out << '\t' << c;
    out << '\n';
    for (const auto& [id, cols] : view.blocks) {
        out << id;
        for (const std::string& v : cols) out << '\t' << v;
        out << '\n';
    }
    return out.str();
}

void write_view(const std::string& prefix, const View& view) {
    // the locus rows: earlier assessment rows dropped, the current ones appended (idempotent)
    {
        const std::string path = genotype_output_path(prefix, ".sequence_locus.tsv");
        std::ifstream in(path);
        if (in) {
            std::string text;
            for (std::string line; std::getline(in, line);) {
                const std::string key = line.substr(0, line.find('\t'));
                if (std::find(locus_keys().begin(), locus_keys().end(), key) != locus_keys().end()) continue;
                text += line + "\n";
            }
            in.close();
            for (const auto& [k, v] : view.locus) text += k + "\t" + v + "\n";
            atomic_write(path, text);
        }
    }
    std::map<std::string, std::vector<std::string>> by_id;
    for (const auto& [id, cols] : view.blocks) by_id[id] = cols;
    rewrite_columns(genotype_output_path(prefix, ".sequence_blocks.tsv"), block_columns(), by_id, "block_id");
    rewrite_columns(prefix + ".block_genotypes.tsv", block_columns(), by_id, "block_id");
    rewrite_columns(prefix + ".genotypes.tsv", block_columns(), by_id, "block_id");
    // the locus mirror: one header and one row
    {
        Table t;
        if (read_table(prefix + ".locus_genotype.tsv", t) && !t.rows.empty()) {
            std::vector<std::string> names = locus_keys();
            names.push_back("headline_filters_status");
            std::vector<std::string> values;
            for (const auto& [k, v] : view.locus) values.push_back(v);
            values.push_back(std::string("superseded_by_") + kRule);
            std::map<std::string, std::vector<std::string>> one;
            const std::string first = t.header.empty() ? "" : t.header.front();
            one[t.rows.front().empty() ? "" : t.rows.front().front()] = values;
            rewrite_columns(prefix + ".locus_genotype.tsv", names, one, first);
        }
    }
}

void write_assessment(const std::string& prefix) { write_view(prefix, assess_prefix(prefix)); }

}  // namespace panvar::assessment
