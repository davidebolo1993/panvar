// M9-J: compound annotations (compound_annotation.hpp) -- the M9-I.1 contract-v2 fixtures, declaration and coalescing gates,
// determinism, and a replay mode over saved tables (no genotyping) for the C++ / Python reference comparison:
//   compound_annotation_test                      run the gates
//   compound_annotation_test replay <manifest> <out>
// manifest rows: case \t block_alleles.tsv \t run prefix \t rep hap1 \t rep hap2 \t first block id \t last block id \t decision groups
#include "panvar/compound_annotation.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

using namespace panvar;

namespace {

int failures = 0;
void check(bool ok, const std::string& what) {
    std::cout << (ok ? "  PASS " : "  FAIL ") << what << '\n';
    if (!ok) ++failures;
}

struct Pool {
    std::vector<std::vector<std::string>> blocks;
    std::vector<std::vector<std::pair<std::vector<std::size_t>, std::vector<std::size_t>>>> labellings;   // per group
    std::vector<double> mass;
    std::vector<bool> decision;
    CompoundAnnotationRow run(std::size_t first, std::size_t last, const std::vector<std::size_t>& r1, const std::vector<std::size_t>& r2) const {
        std::vector<CompoundGroup> g(mass.size());
        for (std::size_t i = 0; i < mass.size(); ++i) {
            g[i].mass = mass[i];
            g[i].decision_set = decision.empty() ? false : decision[i];
            for (const auto& [a, b] : labellings[i]) g[i].members.push_back({&a, &b});
        }
        const CompoundAlleleBytes bytes = [&](std::size_t b, std::size_t a) -> const std::string& { return blocks[b][a]; };
        return compound_assess(bytes, first, last, g, r1, r2);
    }
};

std::string spell(const std::vector<std::vector<std::string>>& blocks, const std::vector<std::size_t>& l) {
    std::string s;
    for (std::size_t b = 0; b < l.size(); ++b) s += blocks[b][l[b]];
    return s;
}

std::vector<std::string> split(const std::string& s, char sep) {
    std::vector<std::string> out;
    std::string item;
    std::stringstream in(s);
    while (std::getline(in, item, sep)) out.push_back(item);
    return out;
}
std::vector<std::size_t> labels(const std::string& s) {
    std::vector<std::size_t> v;
    for (const std::string& x : split(s, ',')) v.push_back(std::stoul(x));
    return v;
}
std::vector<std::map<std::string, std::string>> table(const std::string& path) {
    std::ifstream in(path);
    std::vector<std::map<std::string, std::string>> rows;
    std::vector<std::string> header;
    for (std::string line; std::getline(in, line);) {
        if (line.empty() || line[0] == '#') continue;
        const std::vector<std::string> f = split(line, '\t');
        if (header.empty()) {
            header = f;
            continue;
        }
        std::map<std::string, std::string> row;
        for (std::size_t i = 0; i < header.size() && i < f.size(); ++i) row[header[i]] = f[i];
        rows.push_back(row);
    }
    return rows;
}

int replay(const std::string& manifest, const std::string& out_path) {
    std::ifstream in(manifest);
    std::ofstream out(out_path);
    out << std::setprecision(17) << "case\tstatus\ttarget_resolution\tdecision_set_projections\tcompetitor\tsupport_lower\tsupport_upper\terror_lower\t"
                                    "error_upper\tgroups\tmissing_groups\tambiguous_groups\tzero_mass_groups\toriginal_total\tnormalisation_discrepancy\n";
    for (std::string line; std::getline(in, line);) {
        if (line.empty()) continue;
        const std::vector<std::string> f = split(line, '\t');
        std::vector<std::vector<std::string>> blocks;
        std::map<std::string, std::size_t> ids;
        for (const auto& r : table(f[1])) {
            const std::size_t b = std::stoul(r.at("block_index"));
            if (blocks.size() <= b) blocks.resize(b + 1);
            blocks[b].push_back(r.at("sequence"));
            ids[r.at("block_id")] = b;
        }
        std::map<std::string, std::size_t> group_index;
        Pool p;
        p.blocks = blocks;
        for (const auto& r : table(f[2] + ".mosaic_scores.tsv")) {
            group_index[r.at("sequence_candidate_id")] = p.mass.size();
            p.mass.push_back(std::stod(r.at("candidate_set_probability")));
        }
        p.labellings.resize(p.mass.size());
        p.decision.assign(p.mass.size(), false);
        for (const auto& r : table(f[2] + ".mosaic_candidates.tsv")) {
            p.labellings[group_index.at(r.at("sequence_candidate_id"))].push_back({labels(r.at("hap1_alleles")), labels(r.at("hap2_alleles"))});
        }
        for (const std::string& g : split(f[7], ',')) {
            if (group_index.count(g)) p.decision[group_index.at(g)] = true;
        }
        const CompoundAnnotationRow row = p.run(ids.at(f[5]), ids.at(f[6]), labels(f[3]), labels(f[4]));
        out << f[0] << '\t' << row.status << '\t' << row.target_resolution << '\t' << row.decision_set_projections << '\t' << row.competitor << '\t'
            << row.support_lower << '\t' << row.support_upper << '\t' << row.error_lower << '\t' << row.error_upper << '\t' << row.groups << '\t'
            << row.missing_groups << '\t' << row.ambiguous_groups << '\t' << row.zero_mass_groups << '\t' << row.original_total << '\t'
            << row.normalisation_discrepancy << '\n';
    }
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc == 4 && std::string(argv[1]) == "replay") return replay(argv[2], argv[3]);
    std::cout << "compound annotations (contract v2):\n";
    const std::vector<std::vector<std::string>> blocks = {{"TT"}, {"A", "AC"}, {"CG", "G"}, {"CCC", "CCA"}};
    const std::vector<std::size_t> x = {0, 0, 0, 0}, y = {0, 1, 1, 0}, z = {0, 1, 0, 0};
    // M1 missing projections: both uppers, neither lower; competitor UNKNOWN
    {
        Pool p{blocks, {{{x, x}}, {}}, {0.6, 0.4}, {}};
        const auto r = p.run(1, 2, x, x);
        check(std::fabs(r.support_lower - 0.6) < 1e-15 && std::fabs(r.support_upper - 1.0) < 1e-15 && r.error_lower == 0.0 &&
                  std::fabs(r.error_upper - 0.4) < 1e-15 && r.missing_groups == 1 && r.competitor == "UNKNOWN" &&
                  r.status == "ASSESSED_WITH_MISSING_PROJECTIONS",
              "M1 missing group: both upper bounds only; counted; UNKNOWN competitor");
    }
    // M2 empty pool / zero total -> UNASSESSED; M3 invalid -> INVALID
    {
        check(Pool{blocks, {}, {}, {}}.run(1, 2, x, x).status == "UNASSESSED" && Pool{blocks, {{{x, x}}}, {0.0}, {}}.run(1, 2, x, x).status == "UNASSESSED",
              "M2 empty or zero-mass pool is UNASSESSED");
        bool all = true;
        for (const double bad : {std::nan(""), static_cast<double>(INFINITY), -1e-12}) all = all && Pool{blocks, {{{x, x}}}, {bad}, {}}.run(1, 2, x, x).status == "INVALID";
        check(all, "M3 NaN, infinite and negative masses are INVALID");
    }
    // M4 rounded totals: normalised, discrepancy recorded
    for (const double tot : {1.0000003, 0.9999997}) {
        Pool p{blocks, {{{x, x}}, {{z, z}}}, {0.7 * tot, 0.3 * tot}, {}};
        const auto r = p.run(1, 2, x, x);
        check(std::fabs(r.support_lower + r.error_upper - 1.0) < 1e-15 && std::fabs(r.normalisation_discrepancy - (tot - 1.0)) < 1e-12,
              "M4 total " + std::to_string(tot) + ": normalised, discrepancy recorded");
    }
    // M5 tiny error retained; M6 zero-mass competitor; M6b wholly supporting pool
    {
        check(Pool{blocks, {{{x, x}}, {{z, z}}}, {1.0, 1e-300}, {}}.run(1, 2, x, x).error_lower == 1e-300, "M5 a 1e-300 error mass is retained exactly");
        const auto r = Pool{blocks, {{{x, x}}, {{z, z}}}, {1.0, 0.0}, {}}.run(1, 2, x, x);
        check(r.competitor == "EVALUATED_ALTERNATIVE_ZERO_MASS" && r.error_upper == 0.0 && r.zero_mass_groups == 1, "M6 zero-mass competitor detected and counted");
        check(Pool{blocks, {{{x, x}}}, {1.0}, {}}.run(1, 2, x, x).competitor == "NO_EVALUATED_ALTERNATIVE", "M6b wholly supporting pool");
    }
    // M7 physically realisable mixed projection ("A"+"CG" == "AC"+"G"), single-block interval; M7b enclosing interval
    {
        const bool same = spell(blocks, x) == spell(blocks, y);
        Pool p{blocks, {{{x, x}, {y, y}}, {{z, z}}}, {0.8, 0.2}, {}};
        const auto r = p.run(1, 1, x, x);
        check(same && r.ambiguous_groups == 1 && r.support_lower == 0.0 && std::fabs(r.support_upper - 0.8) < 1e-15 && std::fabs(r.error_lower - 0.2) < 1e-15 &&
                  std::fabs(r.error_upper - 1.0) < 1e-15,
              "M7 realisable mixed projection: identical whole sequence, interval [0, 0.8] support, mass unsplit");
        const auto r2 = p.run(1, 2, x, x);
        check(r2.ambiguous_groups == 0 && std::fabs(r2.support_lower - 0.8) < 1e-15, "M7b the enclosing interval has one projection");
    }
    // M8 phase counterexample: truth AC/GT, prediction AT/GC
    {
        const std::vector<std::vector<std::string>> pb = {{"A", "G"}, {"C", "T"}};
        Pool p{pb, {{{{0, 0}, {1, 1}}}, {{{0, 1}, {1, 0}}}}, {0.5, 0.5}, {}};
        const auto r = p.run(0, 1, {0, 0}, {1, 1});
        check(std::fabs(r.support_lower - 0.5) < 1e-15 && std::fabs(r.error_lower - 0.5) < 1e-15, "M8 phase: AT/GC is a different compound pair from AC/GT");
    }
    // homologue exchange and duplicate labellings
    {
        Pool a{blocks, {{{x, y}}, {{z, x}}}, {0.6, 0.4}, {}};
        Pool b{blocks, {{{y, x}}, {{x, z}}}, {0.6, 0.4}, {}};
        Pool d{blocks, {{{x, y}, {x, y}, {y, x}}, {{z, x}}}, {0.6, 0.4}, {}};
        const auto ra = a.run(1, 2, x, y), rb = b.run(1, 2, y, x), rd = d.run(1, 2, x, y);
        check(ra.support_lower == rb.support_lower && ra.error_lower == rb.error_lower && ra.target_a_sha256 == rb.target_a_sha256,
              "homologue exchange: identical bounds and target");
        check(rd.support_lower == ra.support_lower && rd.support_upper == ra.support_upper && rd.ambiguous_groups == 0,
              "duplicate labellings create no mass and no ambiguity");
    }
    // decision-set disagreement -> UNRESOLVED_DECISION_SET, representative-only mass reference
    {
        Pool p{blocks, {{{x, x}}, {{z, z}}}, {0.5, 0.5}, {true, true}};
        const auto r = p.run(1, 2, x, x);
        check(r.target_resolution == "UNRESOLVED_DECISION_SET" && r.decision_set_projections == 2 &&
                  r.mass_reference.find("representative_projection_only") == 0,
              "decision-set members with different compound projections: unresolved, representative-only mass");
        Pool q{blocks, {{{x, x}}, {{y, y}}}, {0.5, 0.5}, {true, true}};
        check(q.run(1, 2, x, x).target_resolution == "RESOLVED", "decision-set members with synonymous labellings: resolved");
    }
    // declarations: parse, fingerprint / block-id mismatch, non-consecutive, coalescing overlaps (adjacent stays apart)
    {
        const std::string path = "compound_decl_test.tsv";
        {
            std::ofstream f(path);
            f << "# compound_interval_declarations_v1\nmodel_fingerprint\tabc\nblock_ids\tb0,b1,b2,b3,b4\ninterval\tI1\tb1,b2\tinventory\n"
                 "interval\tI2\tb4\tinventory\n";
        }
        const CompoundDeclarations d = read_compound_declarations(path);
        const std::vector<std::string> ids = {"b0", "b1", "b2", "b3", "b4"};
        std::string st;
        auto spans = resolve_compound_spans(&d, "abc", ids, {{2, 3}}, &st);
        check(st == "VALID" && spans.size() == 2 && spans[0].first == 1 && spans[0].last == 3 && spans[0].sources.size() == 2 && spans[1].first == 4,
              "declarations valid; overlapping representation interval coalesced (b1..b3), adjacent b4 kept apart");
        spans = resolve_compound_spans(&d, "zzz", ids, {{2, 3}}, &st);
        check(st == "INVALID_FINGERPRINT" && spans.size() == 1 && spans[0].sources[0] == "representation_interval",
              "fingerprint mismatch: declarations dropped, representation intervals kept");
        spans = resolve_compound_spans(&d, "abc", {"b0", "b1", "b2", "bX", "b4"}, {}, &st);
        check(st == "INVALID_BLOCK_IDS" && spans.empty(), "block-id mismatch: declarations dropped");
        {
            std::ofstream f(path);
            f << "# compound_interval_declarations_v1\nmodel_fingerprint\tabc\nblock_ids\tb0,b1,b2,b3,b4\ninterval\tI1\tb1,b3\tinventory\n";
        }
        const CompoundDeclarations gap = read_compound_declarations(path);
        resolve_compound_spans(&gap, "abc", ids, {}, &st);
        check(st == "INVALID_BLOCK_IDS", "non-consecutive declared blocks are invalid");
        {
            std::ofstream f(path);
            f << "model_fingerprint\tabc\n";
        }
        bool threw = false;
        try {
            read_compound_declarations(path);
        } catch (const std::runtime_error&) {
            threw = true;
        }
        check(threw, "a malformed declaration file is refused");
        std::remove(path.c_str());
    }
    // determinism: identical rows from concurrent invocations
    {
        Pool p{blocks, {{{x, x}, {y, y}}, {{z, z}}, {{x, z}}}, {0.5, 0.3, 0.2}, {true, false, false}};
        std::vector<CompoundAnnotationRow> rows(4);
        std::vector<std::thread> th;
        for (int i = 0; i < 4; ++i) th.emplace_back([&, i] { rows[static_cast<std::size_t>(i)] = p.run(1, 2, x, x); });
        for (auto& t : th) t.join();
        bool same = true;
        for (const auto& r : rows) {
            same = same && r.support_lower == rows[0].support_lower && r.error_upper == rows[0].error_upper && r.target_a_sha256 == rows[0].target_a_sha256;
        }
        check(same, "identical annotation from 4 concurrent invocations");
    }
    std::cout << (failures == 0 ? "compound_annotation: all passed\n" : "compound_annotation: FAILURES\n");
    return failures == 0 ? 0 : 1;
}
