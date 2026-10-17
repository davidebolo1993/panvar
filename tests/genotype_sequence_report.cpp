// S1-B: the sequence-first reporting layer on hand-built pooled decisions (no reads, no scoring). Toy loci spell
// blocks from short strings; decide_pooled_genotype makes the same decision the caller makes.
//   coupled  -- a unit counted in sv:1 or sv:3 across an UNCHANGED context ctx:2 that shares its period (the CYP2D6
//               sv:3..ctx:4..sv:5 shape): one sequence pair, three explored labellings. The interval must enclose ctx:2,
//               keep ctx:2 IDENTIFIED, and hold exactly the observed tuples (no Cartesian combination).
//   tie      -- two distinct sequence pairs with bit-identical scores that disagree at sv:3: SEQUENCE_TIED, the
//               representative chosen by the sorted SHA-256 pair whatever the input order.
//   ge       -- two distinct sequence pairs (a phase switch) with identical block genotypes: GENOTYPE_RESOLVED_SEQUENCE_SET.
//   bounds   -- a unique winner and a lower group whose two labellings disagree at sv:1: lower/upper error mass count
//               that group once (upper only), never twice.
//   outputs  -- PREDICTION / PROVISIONAL / NO_PREDICTION files, markers and exit codes.
#include "panvar/genotype_sequence_report.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace panvar;

namespace {

void require(bool ok, const std::string& what) {
    if (!ok) {
        std::cerr << "FAILED: " << what << '\n';
        std::exit(1);
    }
}

struct Toy {
    std::vector<std::vector<std::string>> alleles;   // per block
    std::vector<std::string> ids;
    std::string spell(const std::vector<std::size_t>& labels) const {
        std::string out;
        for (std::size_t b = 0; b < labels.size(); ++b) out += alleles[b][labels[b]];
        return out;
    }
};

PoolCandidate candidate(const Toy& toy, std::vector<std::size_t> a, std::vector<std::size_t> b, double score) {
    PoolCandidate c;
    c.first = std::move(a);
    c.second = std::move(b);
    c.first_sequence = sequence_digest(toy.spell(c.first));
    c.second_sequence = sequence_digest(toy.spell(c.second));
    c.score.log_likelihood = score;
    c.score.fragment_term = score;
    c.score.exposure_term = 0.0;
    c.score.absolute_term_sum = std::abs(score);
    c.score.fragments_with_mass = 10;
    c.sequence_prior = 0.0;
    c.path_prior = 0.0;
    return c;
}

std::string slurp(const std::string& path) {
    std::ifstream in(path);
    std::stringstream s;
    s << in.rdbuf();
    return s.str();
}

std::map<std::string, std::string> kv(const std::string& path) {
    std::map<std::string, std::string> out;
    std::istringstream in(slurp(path));
    for (std::string line; std::getline(in, line);) {
        const std::size_t tab = line.find('\t');
        if (tab != std::string::npos) out[line.substr(0, tab)] = line.substr(tab + 1);
    }
    return out;
}

// ctx:0 C, sv:1 {"", A}, ctx:2 {AA}, sv:3 {"", A}, ctx:4 G: (A,-) and (-,A) both spell C AAA G.
Toy coupled_toy() {
    return Toy{{{"C"}, {"", "A"}, {"AA"}, {"", "A"}, {"G"}}, {"ctx:0", "sv:1", "ctx:2", "sv:3", "ctx:4"}};
}

void test_coupled() {
    const Toy toy = coupled_toy();
    const std::vector<std::size_t> x{0, 1, 0, 0, 0}, y{0, 0, 0, 1, 0}, none{0, 0, 0, 0, 0};
    std::vector<PoolCandidate> pool = {candidate(toy, x, x, -10.0), candidate(toy, x, y, -10.0), candidate(toy, y, y, -10.0),
                                       candidate(toy, none, x, -20.0)};
    const PooledDecision d = decide_pooled_genotype(pool, 5, 10, 6);
    require(d.decision == TopSetDecision::DistinctTie, "fixture: the legacy decision is a label-only distinct tie");
    const SequenceReport r = build_sequence_report(d, toy.ids, [&](const auto& l) { return toy.spell(l); });
    require(r.prediction == "PREDICTION" && r.sequence_status == "RESOLVED", "coupled: one sequence pair is RESOLVED");
    require(r.blocks[0].status == "IDENTIFIED" && r.blocks[2].status == "IDENTIFIED" && r.blocks[4].status == "IDENTIFIED",
            "coupled: ctx:0, ctx:2 and ctx:4 stay IDENTIFIED");
    require(r.blocks[1].status == "REPRESENTATION_AMBIGUOUS" && r.blocks[3].status == "REPRESENTATION_AMBIGUOUS",
            "coupled: sv:1 and sv:3 are representation-ambiguous");
    require(r.intervals.size() == 1 && r.intervals[0].first_block == 1 && r.intervals[0].last_block == 3,
            "coupled: ONE interval sv:1..sv:3 enclosing the unchanged context ctx:2 (not two adjacent-block runs)");
    require(r.blocks[2].interval == "sv:1..sv:3", "coupled: ctx:2 is marked inside the interval");
    // observed tuples over (sv:1, ctx:2, sv:3), per haplotype; identical sequences -> each pair sorted
    const std::vector<std::pair<std::vector<std::size_t>, std::vector<std::size_t>>> expect = {
        {{0, 0, 1}, {0, 0, 1}}, {{0, 0, 1}, {1, 0, 0}}, {{1, 0, 0}, {1, 0, 0}}};
    require(r.intervals[0].tuples == expect, "coupled: exactly the three observed tuples");
    for (const auto& t : r.intervals[0].tuples) {
        require(!(t.first == std::vector<std::size_t>{1, 0, 1}) && !(t.second == std::vector<std::size_t>{1, 0, 1}) &&
                    !(t.first == std::vector<std::size_t>{0, 0, 0}),
                "coupled: no Cartesian combination (A at both sv:1 and sv:3, or at neither) is invented");
    }
    require(r.hap1_bases == "CAAAG" && r.hap2_bases == "CAAAG", "coupled: the representative spells the pair");
    require(r.blocks[0].error_mass_lower == 0.0 && r.blocks[0].error_mass_upper == 0.0, "coupled: ctx:0 bounds 0");
    std::cout << "  coupled ambiguity across unchanged context: one enclosing interval sv:1..sv:3, ctx:2 identified inside "
                 "it, three observed tuples, no Cartesian combinations\n";
}

void test_tie_order_independent() {
    const Toy toy{{{"C"}, {"T", "A"}, {"GG"}, {"T", "C", "G"}, {"A"}}, {"ctx:0", "sv:1", "ctx:2", "sv:3", "ctx:4"}};
    const std::vector<std::size_t> h1{0, 0, 0, 0, 0}, h2{0, 0, 0, 1, 0}, h3{0, 0, 0, 2, 0};
    std::vector<PoolCandidate> pool = {candidate(toy, h1, h2, -7.0), candidate(toy, h1, h3, -7.0),
                                       candidate(toy, h2, h3, -30.0)};
    std::string chosen;
    std::mt19937 rng(11);
    for (int trial = 0; trial < 6; ++trial) {
        std::shuffle(pool.begin(), pool.end(), rng);
        const PooledDecision d = decide_pooled_genotype(pool, 5, 10, 6);
        const SequenceReport r = build_sequence_report(d, toy.ids, [&](const auto& l) { return toy.spell(l); });
        require(r.sequence_status == "SEQUENCE_TIED" && r.decision_set_groups == 2, "tie: SEQUENCE_TIED over two groups");
        require(r.exact_score_ties == 2 && r.representative_tie_break == "canonical_sha256", "tie: broken by SHA-256");
        require(r.blocks[3].status == "GENOTYPE_UNRESOLVED" && r.blocks[1].status == "IDENTIFIED", "tie: sv:3 unresolved");
        // the expected representative: the smaller sorted SHA-256 pair
        auto key = [&](const std::vector<std::size_t>& a, const std::vector<std::size_t>& b) {
            std::string x = sha256_bytes(toy.spell(a)), y = sha256_bytes(toy.spell(b));
            return x < y ? x + y : y + x;
        };
        const std::string want = std::min(key(h1, h2), key(h1, h3));
        require(r.hap1_sha256 + r.hap2_sha256 == want, "tie: the representative is the canonical SHA-256 minimum");
        if (trial == 0) chosen = r.hap1_sha256 + r.hap2_sha256;
        require(chosen == r.hap1_sha256 + r.hap2_sha256, "tie: the representative never depends on input order");
    }
    std::cout << "  distinct-sequence exact tie: SEQUENCE_TIED, representative by canonical SHA-256 over 6 input orders\n";
}

void test_genotype_equivalent() {
    const Toy toy{{{"C", "T"}, {"AA"}, {"G", "A"}}, {"sv:0", "ctx:1", "sv:2"}};
    std::vector<PoolCandidate> pool = {candidate(toy, {0, 0, 0}, {1, 0, 1}, -5.0), candidate(toy, {0, 0, 1}, {1, 0, 0}, -5.0)};
    const PooledDecision d = decide_pooled_genotype(pool, 3, 10, 6);
    require(d.decision == TopSetDecision::GenotypeEquivalent, "fixture: legacy genotype-equivalent");
    const SequenceReport r = build_sequence_report(d, toy.ids, [&](const auto& l) { return toy.spell(l); });
    require(r.sequence_status == "GENOTYPE_RESOLVED_SEQUENCE_SET" && r.decision_set_groups == 2,
            "ge: a genotype-resolved sequence set");
    for (const auto& b : r.blocks) require(b.status == "IDENTIFIED", "ge: every block identified");
    require(!r.hap1_bases.empty(), "ge: a representative is retained");
    std::cout << "  genotype-equivalent sequence set: all blocks identified, representative retained\n";
}

void test_bounds_count_groups_once() {
    const Toy toy = coupled_toy();
    const std::vector<std::size_t> none{0, 0, 0, 0, 0}, x{0, 1, 0, 0, 0}, y{0, 0, 0, 1, 0};
    // winner: (none, none); lower group: one extra unit, explored in sv:1 (x) and in sv:3 (y)
    std::vector<PoolCandidate> pool = {candidate(toy, none, none, -10.0), candidate(toy, none, x, -10.0 - std::log(3.0)),
                                       candidate(toy, none, y, -10.0 - std::log(3.0))};
    const PooledDecision d = decide_pooled_genotype(pool, 5, 10, 6);
    const SequenceReport r = build_sequence_report(d, toy.ids, [&](const auto& l) { return toy.spell(l); });
    require(r.sequence_status == "RESOLVED" && d.groups.size() == 2, "bounds: two groups, the winner resolved");
    const double p_lower = 1.0 / (1.0 + 1.0 / 3.0) / 3.0;   // the lower group's candidate-set probability, once
    for (const std::size_t b : {1, 3}) {
        require(r.blocks[b].status == "IDENTIFIED", "bounds: sv:1 and sv:3 identified for the winner");
        require(r.blocks[b].error_mass_lower == 0.0, "bounds: lower 0 (one of the group's labellings agrees)");
        require(std::abs(r.blocks[b].error_mass_upper - p_lower) < 1e-12,
                "bounds: upper = the group's probability counted once, not twice");
    }
    require(std::abs(r.blocks[1].error_mass_upper - r.blocks[3].error_mass_upper) < 1e-15,
            "bounds: the symmetric blocks get symmetric bounds");
    std::cout << "  error-mass bounds: a two-labelling lower group counted once (upper " << p_lower << ", lower 0) at both "
              << "symmetric blocks\n";
}

void test_outputs(const std::string& dir) {
    std::filesystem::create_directories(dir);
    const Toy toy = coupled_toy();
    const std::vector<std::size_t> x{0, 1, 0, 0, 0}, y{0, 0, 0, 1, 0};
    std::vector<PoolCandidate> pool = {candidate(toy, x, x, -10.0), candidate(toy, x, y, -10.0)};
    const PooledDecision d = decide_pooled_genotype(pool, 5, 10, 6);
    const SequenceReport r = build_sequence_report(d, toy.ids, [&](const auto& l) { return toy.spell(l); });
    const std::string p = dir + "/pred";
    write_sequence_report(r, {"RESCORE_INCOMPLETE", "no_call", "distinct_tie"}, p, "unused");
    auto locus = kv(p + ".sequence_locus.tsv");
    require(locus["prediction"] == "PREDICTION" && locus["sequence_status"] == "RESOLVED" && locus["exit_code"] == "0" &&
                locus["structure_confidence"] == "NA" && locus["structure_confidence_reason"] == "UNCALIBRATED" &&
                locus["projection_coverage"] == "EXPLORED" && locus["legacy_no_call_reason"] == "distinct_tie" &&
                locus["representation_intervals"] == "sv:1..sv:3",
            "outputs: the PREDICTION locus record");
    require(slurp(p + ".sequence_haplotypes.fa").find("prediction=PREDICTION sequence_status=RESOLVED") != std::string::npos,
            "outputs: the FASTA is marked");
    require(sequence_report_exit_code(r) == 0, "outputs: PREDICTION exits 0");

    const SequenceReport prov = provisional_sequence_report(x, y, -12.5, "round_cap", "stopped", toy.ids,
                                                            [&](const auto& l) { return toy.spell(l); });
    write_sequence_report(prov, {"REPAIR_INCOMPLETE", "no_call", "round_cap"}, dir + "/prov", "unused");
    locus = kv(dir + "/prov.sequence_locus.tsv");
    require(locus["prediction"] == "PROVISIONAL" && locus["sequence_status"] == "NOT_DETERMINED" &&
                locus["provisional"] == "yes" && locus["completion"] == "round_cap" && locus["exit_code"] == "3" &&
                locus["decision_scope"] == "NOT_DETERMINED" && locus["structure_confidence"] == "NA",
            "outputs: the PROVISIONAL record claims no completed search or decision set");
    require(slurp(dir + "/prov.sequence_blocks.tsv").find("PROVISIONAL_LABEL") != std::string::npos &&
                slurp(dir + "/prov.sequence_blocks.tsv").find("EXPLORED") == std::string::npos,
            "outputs: provisional blocks carry labels only (no projection coverage, no bounds)");
    require(sequence_report_exit_code(prov) == 3, "outputs: PROVISIONAL exits 3");

    write_sequence_report(SequenceReport{}, {"TIME_LIMIT", "no_call", "none"}, dir + "/none", "why");
    locus = kv(dir + "/none.sequence_locus.tsv");
    require(locus["prediction"] == "NO_PREDICTION" && locus["exit_code"] == "1" && locus["completion_note"] == "why" &&
                !std::filesystem::exists(dir + "/none.sequence_haplotypes.fa"),
            "outputs: NO_PREDICTION writes the reason and no FASTA");
    require(sequence_report_exit_code(SequenceReport{}) == 1, "outputs: NO_PREDICTION exits 1");
    std::cout << "  outputs: PREDICTION (0), PROVISIONAL (3, no decision set, no bounds), NO_PREDICTION (1, no FASTA)\n";
}

}  // namespace

int main(int argc, char** argv) {
    const std::string dir = argc > 1 ? argv[1] : "genotype_sequence_report_out";
    test_coupled();
    test_tie_order_independent();
    test_genotype_equivalent();
    test_bounds_count_groups_once();
    test_outputs(dir);
    std::cout << "genotype_sequence_report: all passed\n";
    return 0;
}
