#include "panvar/compound_annotation.hpp"

#include "panvar/genotype_repair.hpp"   // sha256_bytes, sha256_file

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace panvar {

namespace {

std::vector<std::string> split(const std::string& s, char sep) {
    std::vector<std::string> out;
    std::string item;
    std::stringstream in(s);
    while (std::getline(in, item, sep)) out.push_back(item);
    return out;
}

// Neumaier compensated summation, in the caller's (group) order: deterministic.
struct Sum {
    double sum = 0.0, c = 0.0;
    void add(double x) {
        const double t = sum + x;
        c += std::fabs(sum) >= std::fabs(x) ? (sum - t) + x : (x - t) + sum;
        sum = t;
    }
    double value() const { return sum + c; }
};

std::string number(double v) {
    std::ostringstream o;
    o << std::setprecision(17) << v;
    return o.str();
}

}  // namespace

CompoundDeclarations read_compound_declarations(const std::string& path) {
    std::ifstream file(path);
    if (!file) throw std::runtime_error("cannot read compound interval declarations " + path);
    CompoundDeclarations d;
    d.path = path;
    d.sha256 = sha256_file(path);
    bool header = false;
    for (std::string line; std::getline(file, line);) {
        if (line.empty()) continue;
        if (line[0] == '#') {
            if (line.rfind("# compound_interval_declarations_v1", 0) == 0) header = true;
            continue;
        }
        const std::vector<std::string> f = split(line, '\t');
        if (f.size() >= 2 && f[0] == "model_fingerprint") {
            d.model_fingerprint = f[1];
        } else if (f.size() >= 2 && f[0] == "block_ids") {
            d.model_block_ids = split(f[1], ',');
        } else if (f.size() >= 4 && f[0] == "interval") {
            d.intervals.push_back({f[1], split(f[2], ','), f[3]});
        } else {
            throw std::runtime_error("compound interval declarations " + path + ": unrecognised line: " + line);
        }
    }
    if (!header || d.model_fingerprint.empty() || d.model_block_ids.empty()) {
        throw std::runtime_error("compound interval declarations " + path + ": needs the v1 header, model_fingerprint and block_ids");
    }
    return d;
}

std::vector<CompoundSpan> resolve_compound_spans(const CompoundDeclarations* declarations, const std::string& model_fingerprint,
                                                 const std::vector<std::string>& block_ids,
                                                 const std::vector<std::pair<std::size_t, std::size_t>>& representation_intervals,
                                                 std::string* declaration_status) {
    std::vector<CompoundSpan> spans;
    std::string status = "NONE";
    if (declarations != nullptr) {
        if (declarations->model_fingerprint != model_fingerprint) {
            status = "INVALID_FINGERPRINT";
        } else if (declarations->model_block_ids != block_ids) {
            status = "INVALID_BLOCK_IDS";
        } else {
            status = "VALID";
            std::vector<CompoundSpan> declared;
            for (const CompoundDeclaration& d : declarations->intervals) {
                std::vector<std::size_t> idx;
                for (const std::string& id : d.block_ids) {
                    const auto it = std::find(block_ids.begin(), block_ids.end(), id);
                    if (it == block_ids.end()) {
                        status = "INVALID_BLOCK_IDS";
                        break;
                    }
                    idx.push_back(static_cast<std::size_t>(it - block_ids.begin()));
                }
                if (status != "VALID") break;
                if (idx.empty()) {
                    status = "INVALID_BLOCK_IDS";
                    break;
                }
                for (std::size_t k = 1; k < idx.size(); ++k) {
                    if (idx[k] != idx[k - 1] + 1) status = "INVALID_BLOCK_IDS";
                }
                if (status != "VALID") break;
                declared.push_back({idx.front(), idx.back(), {d.id}});
            }
            if (status == "VALID") spans = std::move(declared);
        }
    }
    for (const auto& [a, b] : representation_intervals) spans.push_back({a, b, {"representation_interval"}});
    std::sort(spans.begin(), spans.end(), [](const CompoundSpan& x, const CompoundSpan& y) { return std::tie(x.first, x.last) < std::tie(y.first, y.last); });
    std::vector<CompoundSpan> merged;
    for (CompoundSpan& s : spans) {
        if (!merged.empty() && s.first <= merged.back().last) {   // overlapping (a shared block): coalesce
            merged.back().last = std::max(merged.back().last, s.last);
            for (const std::string& src : s.sources) {
                if (std::find(merged.back().sources.begin(), merged.back().sources.end(), src) == merged.back().sources.end()) {
                    merged.back().sources.push_back(src);
                }
            }
        } else {
            merged.push_back(std::move(s));
        }
    }
    if (declaration_status != nullptr) *declaration_status = status;
    return merged;
}

CompoundAnnotationRow compound_assess(const CompoundAlleleBytes& bytes, std::size_t first, std::size_t last,
                                      const std::vector<CompoundGroup>& groups, const std::vector<std::size_t>& representative_first,
                                      const std::vector<std::size_t>& representative_second) {
    CompoundAnnotationRow row;
    row.first = first;
    row.last = last;
    row.groups = groups.size();
    // exact-sequence interning: label tuple -> sequence id (bytes compared exactly, no hash)
    std::map<std::vector<std::size_t>, std::size_t> tuple_id;
    std::map<std::string, std::size_t> bytes_id;
    std::vector<const std::string*> sequences;
    const auto hap = [&](const std::vector<std::size_t>& labels) {
        std::vector<std::size_t> key(labels.begin() + static_cast<std::ptrdiff_t>(first), labels.begin() + static_cast<std::ptrdiff_t>(last) + 1);
        const auto found = tuple_id.find(key);
        if (found != tuple_id.end()) return found->second;
        std::string s;
        for (std::size_t b = first; b <= last; ++b) s += bytes(b, labels[b]);
        const auto ins = bytes_id.emplace(std::move(s), bytes_id.size());
        if (ins.second) sequences.push_back(&ins.first->first);
        tuple_id.emplace(std::move(key), ins.first->second);
        return ins.first->second;
    };
    const auto pair_of = [&](const std::vector<std::size_t>& a, const std::vector<std::size_t>& b) {
        const std::size_t x = hap(a), y = hap(b);
        return x <= y ? std::make_pair(x, y) : std::make_pair(y, x);   // by value
    };
    const auto target = pair_of(representative_first, representative_second);
    {
        const std::string& a = *sequences[target.first];
        const std::string& b = *sequences[target.second];
        std::string ha = sha256_bytes(a), hb = sha256_bytes(b);
        const bool swap = hb < ha;
        row.target_a_sha256 = swap ? hb : ha;
        row.target_b_sha256 = swap ? ha : hb;
        row.target_a_bp = swap ? b.size() : a.size();
        row.target_b_bp = swap ? a.size() : b.size();
    }
    // decision-set resolution
    std::set<std::pair<std::size_t, std::size_t>> ds;
    for (const CompoundGroup& g : groups) {
        if (!g.decision_set) continue;
        ++row.decision_set_groups;
        if (g.members.empty()) ++row.decision_set_missing_groups;   // M9-J.1: its projection is unknown
        for (const auto& [m1, m2] : g.members) ds.insert(pair_of(*m1, *m2));
    }
    row.decision_set_projections = ds.size();
    if (!ds.empty()) {
        std::vector<std::string> digests;
        for (const auto& [x, y] : ds) {
            std::string a = sha256_bytes(*sequences[x]).substr(0, 16), b = sha256_bytes(*sequences[y]).substr(0, 16);
            if (b < a) std::swap(a, b);
            digests.push_back(a + "|" + b);
        }
        std::sort(digests.begin(), digests.end());
        row.decision_set_projection_digests.clear();
        for (const std::string& d : digests) row.decision_set_projection_digests += (row.decision_set_projection_digests.empty() ? "" : ";") + d;
    }
    // M9-J.1: resolution independent of the mass bounds; an incompletely projected decision set is never RESOLVED
    if (row.decision_set_groups == 0) {
        row.target_resolution = "NO_DECISION_SET";
        row.mass_reference = "representative_projection_only_no_decision_set";
    } else if (row.decision_set_missing_groups > 0) {
        row.target_resolution = "INCOMPLETE_DECISION_SET";
        row.mass_reference = "representative_projection_only_decision_set_incompletely_projected";
    } else if (ds.size() == 1) {
        row.target_resolution = "RESOLVED";
        row.mass_reference = "the_resolved_target";
    } else {
        row.target_resolution = "UNRESOLVED_DECISION_SET";
        row.mass_reference = "representative_projection_only_not_a_resolved_target";
    }
    // masses
    Sum total;
    for (const CompoundGroup& g : groups) {
        if (!std::isfinite(g.mass) || g.mass < 0.0) {
            row.status = "INVALID";
            return row;
        }
        total.add(g.mass);
    }
    row.original_total = total.value();
    if (groups.empty() || row.original_total == 0.0) {
        row.status = "UNASSESSED";
        return row;
    }
    row.normalisation_discrepancy = row.original_total - 1.0;
    Sum sl, su, el, eu;
    std::size_t disagree = 0, disagree_positive = 0;
    for (const CompoundGroup& g : groups) {
        const double q = g.mass / row.original_total;
        row.zero_mass_groups += g.mass == 0.0 ? 1 : 0;
        if (g.members.empty()) {
            ++row.missing_groups;
            su.add(q);
            eu.add(q);
            continue;
        }
        std::set<std::pair<std::size_t, std::size_t>> projections;
        for (const auto& [m1, m2] : g.members) projections.insert(pair_of(*m1, *m2));
        row.ambiguous_groups += projections.size() > 1 ? 1 : 0;
        const bool hit = projections.count(target) != 0;
        const bool miss = projections.size() > (hit ? 1u : 0u);
        if (hit && !miss) sl.add(q);
        if (hit) su.add(q);
        if (!hit) el.add(q);
        if (miss) {
            eu.add(q);
            ++disagree;
            disagree_positive += g.mass > 0.0 ? 1 : 0;
        }
    }
    row.support_lower = sl.value();
    row.support_upper = su.value();
    row.error_lower = el.value();
    row.error_upper = eu.value();
    row.bounds = true;
    row.status = row.missing_groups ? "ASSESSED_WITH_MISSING_PROJECTIONS" : "ASSESSED";
    row.competitor = disagree && disagree_positive ? "EVALUATED_ALTERNATIVE"
                     : disagree                    ? "EVALUATED_ALTERNATIVE_ZERO_MASS"
                     : row.missing_groups          ? "UNKNOWN"
                                                   : "NO_EVALUATED_ALTERNATIVE";
    return row;
}

std::string render_compound_annotation(const CompoundAnnotation& a, const std::vector<std::string>& block_ids) {
    std::ostringstream o;
    o << "# compound interval annotations (M9-I.1 contract v2); declarations " << a.declaration_path << " sha256=" << a.declaration_sha256
      << " status=" << a.declaration_status << "\n"
      << "interval\tfirst_block\tlast_block\tblocks\tsources\ttarget_scope\tstatus\ttarget_resolution\tdecision_set_groups\t"
         "decision_set_missing_groups\tdecision_set_projections\t"
         "decision_set_projection_digests\ttarget_hapA_sha256\ttarget_hapB_sha256\ttarget_hapA_bp\ttarget_hapB_bp\tmass_reference\tcompetitor\t"
         "support_lower\tsupport_upper\terror_lower\terror_upper\tgroups\tmissing_groups\tambiguous_groups\tzero_mass_groups\t"
         "original_total\tnormalisation_discrepancy\tlimitation\n";
    for (const CompoundAnnotationRow& r : a.rows) {
        std::string blocks;
        for (std::size_t b = r.first; b <= r.last && b < block_ids.size(); ++b) blocks += (blocks.empty() ? "" : ",") + block_ids[b];
        std::string sources;
        for (const std::string& s : split(r.sources, ',')) sources += (sources.empty() ? "" : ",") + s;
        const auto num = [&](double v) { return r.bounds ? number(v) : std::string("NA"); };
        o << block_ids[r.first] << ".." << block_ids[r.last] << '\t' << block_ids[r.first] << '\t' << block_ids[r.last] << '\t' << blocks << '\t'
          << (sources.empty() ? "NA" : sources) << '\t'
          << "exact_unordered_pair_of_concatenated_interval_sequences;includes_within_interval_phase;not_structural_equivalence" << '\t'
          << r.status << '\t' << r.target_resolution << '\t' << r.decision_set_groups << '\t' << r.decision_set_missing_groups << '\t'
          << r.decision_set_projections << '\t' << r.decision_set_projection_digests << '\t'
          << r.target_a_sha256 << '\t' << r.target_b_sha256 << '\t' << r.target_a_bp << '\t' << r.target_b_bp << '\t' << r.mass_reference << '\t'
          << r.competitor << '\t' << num(r.support_lower) << '\t' << num(r.support_upper) << '\t' << num(r.error_lower) << '\t'
          << num(r.error_upper) << '\t' << r.groups << '\t' << r.missing_groups << '\t' << r.ambiguous_groups << '\t' << r.zero_mass_groups << '\t'
          << (r.bounds ? number(r.original_total) : std::string("NA")) << '\t' << (r.bounds ? number(r.normalisation_discrepancy) : std::string("NA"))
          << '\t' << "explored_pool_and_recorded_labellings_only;pool_conditional;uncalibrated;no_statement_about_unsearched_alternatives;"
                     "intervals_conservative_not_minimal_or_exhaustive"
          << '\n';
    }
    return o.str();
}

}  // namespace panvar
