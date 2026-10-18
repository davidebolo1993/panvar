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
#include <functional>
#include <set>
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

// ---- S1-C: labels, bases, hashes and lengths travel together -----------------------------------------------------------
// Sol's reproduction (one block, alleles A / C, one candidate (A, C)) and every two-allele case whose internal-digest
// order disagrees with its SHA-256 order: spelling the reported labels must give the FASTA bases, the bases must hash and
// measure as reported, and swapping the input homologues must not change anything.
void check_aligned(const SequenceReport& r, const std::function<std::string(const std::vector<std::size_t>&)>& spell,
                   const std::string& what) {
    require(spell(r.representative_first) == r.hap1_bases && spell(r.representative_second) == r.hap2_bases,
            what + ": the reported labels spell the reported sequences (hap1 and hap2)");
    require(sha256_bytes(r.hap1_bases) == r.hap1_sha256 && sha256_bytes(r.hap2_bases) == r.hap2_sha256,
            what + ": each sequence hashes as reported");
    require(r.hap1_sha256 <= r.hap2_sha256, what + ": hap1 is the SHA-256-smaller sequence");
}

void test_orientation() {
    {   // the reproduction, exactly
        PoolCandidate c;
        c.first = {0};
        c.second = {1};
        c.first_sequence = sequence_digest("A");
        c.second_sequence = sequence_digest("C");
        c.score.log_likelihood = -1;
        c.score.fragment_term = -1;
        c.score.absolute_term_sum = 1;
        c.score.fragments_with_mass = 10;
        const auto spell = [](const std::vector<std::size_t>& x) { return x[0] == 0 ? std::string("A") : std::string("C"); };
        const SequenceReport r = build_sequence_report(decide_pooled_genotype({c}, 1, 10, 6), {"sv:1"}, spell);
        check_aligned(r, spell, "reproduction");
    }
    // two-allele loci over short strings; keep the cases where digest order and SHA-256 order disagree
    const std::vector<std::string> words = {"A", "C", "G", "T", "AC", "CA", "GT", "TTA", "ACGT", "CCGA", "GATTACA", "TGCA"};
    std::size_t disagreeing = 0, checked = 0;
    for (std::size_t i = 0; i < words.size(); ++i) {
        for (std::size_t j = 0; j < words.size(); ++j) {
            if (i == j) continue;
            const std::string& x = words[i];
            const std::string& y = words[j];
            const bool digest_first = sequence_digest(x) < sequence_digest(y);
            const bool sha_first = sha256_bytes(x) < sha256_bytes(y);
            if (digest_first == sha_first) continue;
            ++disagreeing;
            const Toy toy{{{"G"}, {x, y}, {"T"}}, {"ctx:0", "sv:1", "ctx:2"}};
            const auto spell = [&](const std::vector<std::size_t>& l) { return toy.spell(l); };
            std::vector<std::string> fingerprints;
            for (const bool swap : {false, true}) {
                std::vector<std::size_t> a{0, 0, 0}, b{0, 1, 0};
                if (swap) std::swap(a, b);
                const SequenceReport r =
                    build_sequence_report(decide_pooled_genotype({candidate(toy, a, b, -3.0)}, 3, 10, 6), toy.ids, spell);
                check_aligned(r, spell, "digest/SHA disagreement " + x + "/" + y);
                fingerprints.push_back(r.hap1_sha256 + r.hap2_sha256 + std::to_string(r.hap1_bases.size()) + "/" +
                                       std::to_string(r.hap2_bases.size()));
                ++checked;
            }
            require(fingerprints[0] == fingerprints[1], "an input homologue swap preserves the canonical result");
        }
    }
    require(disagreeing > 0, "fixture: some cases disagree between digest order and SHA-256 order");
    std::cout << "  orientation: the reproduction plus " << disagreeing << " digest/SHA-256-disagreeing pairs (" << checked
              << " reports, both homologue orders): labels spell the FASTA sequences, hashes match, swaps invariant\n";
}

// ---- S1-C: the linear enclosing interval, exhaustively on tiny groups ---------------------------------------------------
// ctx C | sv {"", A} | ctx AA | sv {"", A} | ctx AA | sv {"", A} | ctx G: one extra A has three labellings, two extra
// A have three. Every non-empty subset of a group's labellings (identical-sequence and distinct-sequence groups), in
// two input orders: the reported interval must equal [first, last] of the brute-force union of pairwise differences,
// contain every block where two labellings differ, and hold exactly the observed restricted tuples.
void test_linear_interval() {
    const Toy toy{{{"C"}, {"", "A"}, {"AA"}, {"", "A"}, {"AA"}, {"", "A"}, {"G"}},
                  {"ctx:0", "sv:1", "ctx:2", "sv:3", "ctx:4", "sv:5", "ctx:6"}};
    const auto lab = [](int a, int b, int c) {
        return std::vector<std::size_t>{0, std::size_t(a), 0, std::size_t(b), 0, std::size_t(c), 0};
    };
    const std::vector<std::vector<std::size_t>> one = {lab(1, 0, 0), lab(0, 1, 0), lab(0, 0, 1)};
    const std::vector<std::vector<std::size_t>> two = {lab(1, 1, 0), lab(1, 0, 1), lab(0, 1, 1)};
    std::vector<std::vector<std::pair<std::vector<std::size_t>, std::vector<std::size_t>>>> groups(2);
    for (std::size_t i = 0; i < one.size(); ++i) {
        for (std::size_t j = i; j < one.size(); ++j) groups[0].push_back({one[i], one[j]});   // (S1, S1)
    }
    for (const auto& a : one) {
        for (const auto& b : two) groups[1].push_back({a, b});   // (S1, S2)
    }
    std::size_t cases = 0;
    for (const auto& members : groups) {
        const std::size_t n = members.size();
        for (std::size_t mask = 1; mask < (std::size_t{1} << n); ++mask) {
            std::vector<PoolCandidate> pool;
            for (std::size_t k = 0; k < n; ++k) {
                if (mask >> k & 1) pool.push_back(candidate(toy, members[k].first, members[k].second, -4.0));
            }
            std::string reference;
            for (int order = 0; order < 2; ++order) {
                if (order == 1) std::reverse(pool.begin(), pool.end());
                const PooledDecision d = decide_pooled_genotype(pool, 7, 10, 6);
                const SequenceReport r = build_sequence_report(d, toy.ids, [&](const auto& l) { return toy.spell(l); });
                // brute force over the decision's own members (oriented as the report orients them: by SHA-256)
                std::vector<std::pair<std::vector<std::size_t>, std::vector<std::size_t>>> m;
                for (const std::size_t c : d.groups[r.representative_group]) {
                    auto a = d.candidates[c].first, b = d.candidates[c].second;
                    const bool identical = d.candidates[c].first_sequence == d.candidates[c].second_sequence;
                    if (identical ? b < a : sha256_bytes(toy.spell(b)) < sha256_bytes(toy.spell(a))) std::swap(a, b);
                    m.push_back({a, b});
                }
                std::size_t lo = 7, hi = 0;
                for (std::size_t x = 0; x < m.size(); ++x) {
                    for (std::size_t y = x + 1; y < m.size(); ++y) {
                        for (std::size_t b = 0; b < 7; ++b) {
                            if (m[x].first[b] != m[y].first[b] || m[x].second[b] != m[y].second[b]) {
                                lo = std::min(lo, b);
                                hi = std::max(hi, b);
                            }
                        }
                    }
                }
                std::string got;
                if (lo == 7) {
                    require(r.intervals.empty(), "a group whose labellings never differ has no interval");
                } else {
                    require(r.intervals.size() == 1 && r.intervals[0].first_block == lo && r.intervals[0].last_block == hi,
                            "the linear interval equals the brute-force enclosing range");
                    std::set<std::pair<std::vector<std::size_t>, std::vector<std::size_t>>> want;
                    for (const auto& t : m) {
                        std::vector<std::size_t> a(t.first.begin() + lo, t.first.begin() + hi + 1);
                        std::vector<std::size_t> b(t.second.begin() + lo, t.second.begin() + hi + 1);
                        if (members.front().first.size() && toy.spell(t.first) == toy.spell(t.second) && b < a) std::swap(a, b);
                        want.insert({a, b});
                    }
                    require(std::set<std::pair<std::vector<std::size_t>, std::vector<std::size_t>>>(
                                r.intervals[0].tuples.begin(), r.intervals[0].tuples.end()) == want &&
                                r.intervals[0].tuples.size() == want.size(),
                            "the interval holds exactly the observed tuples (no invention, no loss)");
                    for (std::size_t b = lo; b <= hi; ++b) {
                        require(r.blocks[b].interval != "none", "every block of the range is marked in the interval");
                    }
                    got = std::to_string(lo) + ".." + std::to_string(hi) + ":" + std::to_string(r.intervals[0].tuples.size());
                }
                if (order == 0) reference = got;
                require(got == reference, "the interval never depends on input order");
                ++cases;
            }
        }
    }
    std::cout << "  linear enclosing interval: " << cases << " exhaustive subset/order cases equal the brute-force "
              << "pairwise range, contain every differing block and hold exactly the observed tuples\n";
}

// ---- S1-C: a publication failure leaves no partial report ----------------------------------------------------------------
void test_publication_failure(const std::string& dir) {
    const Toy toy = coupled_toy();
    const std::vector<std::size_t> x{0, 1, 0, 0, 0}, y{0, 0, 0, 1, 0};
    const PooledDecision d = decide_pooled_genotype({candidate(toy, x, x, -10.0), candidate(toy, x, y, -10.0)}, 5, 10, 6);
    SequenceReport r = build_sequence_report(d, toy.ids, [&](const auto& l) { return toy.spell(l); });
    const std::string p = dir + "/fail";
    // a non-empty directory where the alternatives table must be staged: that write fails
    std::filesystem::create_directories(p + ".sequence_alternatives.tsv.writing/blocker");
    const bool ok = publish_sequence_report(r, {"RESCORE_INCOMPLETE", "no_call", "distinct_tie", "off"}, p, "unused");
    require(!ok && r.publication_failed && sequence_report_exit_code(r) == 2, "failure: reported, exit 2");
    for (const char* f : {".sequence_blocks.tsv", ".representation_intervals.tsv", ".sequence_alternatives.tsv",
                          ".sequence_haplotypes.fa"}) {
        require(!std::filesystem::exists(p + f), std::string("failure: no partial ") + f + " is left");
    }
    if (std::filesystem::exists(p + ".sequence_locus.tsv")) {
        const auto locus = kv(p + ".sequence_locus.tsv");
        require(locus.at("prediction") == "NO_PREDICTION" && locus.at("completion") == "report_publication_failed" &&
                    locus.at("exit_code") == "2",
                "failure: a locus record, if any, names the failure and claims no prediction");
    }
    require(sequence_report_status_keys(r).find("sequence_publication\tFAILED") != std::string::npos,
            "failure: the status rows say so");
    std::filesystem::remove_all(p + ".sequence_alternatives.tsv.writing");   // the test's own blocker only
    SequenceReport again = build_sequence_report(d, toy.ids, [&](const auto& l) { return toy.spell(l); });
    require(publish_sequence_report(again, {"RESCORE_INCOMPLETE", "no_call", "distinct_tie", "off"}, p, "unused") &&
                kv(p + ".sequence_locus.tsv").at("prediction") == "PREDICTION" &&
                std::filesystem::exists(p + ".sequence_haplotypes.fa"),
            "after the blocker is gone the same prefix publishes completely");
    std::cout << "  publication failure: no partial report, a NO_PREDICTION record naming the failure, exit 2; a clean "
                 "retry publishes completely (locus record last)\n";
}

}  // namespace

int main(int argc, char** argv) {
    const std::string dir = argc > 1 ? argv[1] : "genotype_sequence_report_out";
    test_coupled();
    test_tie_order_independent();
    test_genotype_equivalent();
    test_bounds_count_groups_once();
    test_outputs(dir);
    test_orientation();
    test_linear_interval();
    test_publication_failure(dir);
    std::cout << "genotype_sequence_report: all passed\n";
    return 0;
}
