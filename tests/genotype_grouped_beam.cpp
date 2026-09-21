// Oracle equivalence for the grouped-prefix proposer.
//
// Two oracles, not one. The first is an independent exhaustive enumeration of every ordered pair of
// template paths, which is what the quantities MEAN; the second is propose_mosaic_diplotypes, the
// existing implementation, which is what production computes today. The grouped proposer must agree
// with both on EVERY final candidate, not merely on the winner, and on every quantity a candidate
// carries: identity, sum-product mass, max-product best history, that history's emission and transition
// components, physical history multiplicity and its saturation flag, and the order candidates come back
// in.
//
// The fixtures are small enough to enumerate exhaustively and are run at r = 0, an intermediate r, and
// r = 1, because the transition degenerates differently at each end: at r = 0 a haplotype cannot leave
// its template and the count recurrence collapses to "the path already here"; at r = 1 every template
// is reachable with equal probability and the same-template and other-template branches coincide.
#include "panvar/mosaic_beam.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr double kNegInf = -std::numeric_limits<double>::infinity();

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void require_close(double observed, double expected, double tolerance, const std::string& message) {
    if (!(std::isfinite(observed) && std::isfinite(expected)) || std::abs(observed - expected) > tolerance) {
        throw std::runtime_error(message + ": observed=" + std::to_string(observed) +
                                 " expected=" + std::to_string(expected));
    }
}

double log_add(double lhs, double rhs) {
    if (lhs == kNegInf) return rhs;
    if (rhs == kNegInf) return lhs;
    const double hi = std::max(lhs, rhs);
    return hi + std::log1p(std::exp(std::min(lhs, rhs) - hi));
}

double transition(std::size_t previous, std::size_t current, double r, std::size_t n) {
    if (previous == current) return std::log((1.0 - r) + r / static_cast<double>(n));
    return r > 0.0 ? std::log(r) - std::log(static_cast<double>(n)) : kNegInf;
}

using OracleKey = std::pair<std::vector<std::size_t>, std::vector<std::size_t>>;

struct OracleValue {
    double logsum = kNegInf;
    double best = kNegInf;
    double emission = 0.0;
    double transition_score = 0.0;
    std::uint64_t histories = 0;
};

// Every ordered pair of template paths, scored history by history. The heterozygous factor of two is
// not written anywhere here: it emerges because both ordered assignments of two distinct allele paths
// are enumerated and land on the same unordered key.
std::map<OracleKey, OracleValue> exhaustive(const panvar::PanelAlleleMatrix& panel,
                                            const panvar::PreparedProposalEmissions& emissions,
                                            const std::vector<double>& recombination) {
    const std::size_t n = panel.template_count();
    const std::size_t blocks = panel.block_count();
    std::map<OracleKey, OracleValue> result;
    std::vector<std::size_t> h1(blocks, 0);
    std::vector<std::size_t> h2(blocks, 0);
    std::function<void(std::size_t, double, double)> extend;
    extend = [&](std::size_t block, double emission_score, double transition_score) {
        if (block == blocks) {
            OracleKey key;
            key.first.resize(blocks);
            key.second.resize(blocks);
            for (std::size_t b = 0; b < blocks; ++b) {
                key.first[b] = panel.template_alleles[h1[b]][b];
                key.second[b] = panel.template_alleles[h2[b]][b];
            }
            if (key.second < key.first) std::swap(key.first, key.second);
            OracleValue& value = result[key];
            const double score = emission_score + transition_score;
            value.logsum = log_add(value.logsum, score);
            ++value.histories;
            if (score > value.best) {
                value.best = score;
                value.emission = emission_score;
                value.transition_score = transition_score;
            }
            return;
        }
        for (std::size_t a = 0; a < n; ++a) {
            h1[block] = a;
            for (std::size_t b = 0; b < n; ++b) {
                h2[block] = b;
                double next_transition = transition_score;
                if (block == 0) {
                    next_transition = -2.0 * std::log(static_cast<double>(n));
                } else {
                    const double t1 = transition(h1[block - 1], a, recombination[block - 1], n);
                    const double t2 = transition(h2[block - 1], b, recombination[block - 1], n);
                    if (t1 == kNegInf || t2 == kNegInf) continue;
                    next_transition += t1 + t2;
                }
                const double next_emission =
                    emission_score + emissions.blocks[block].at(panel.template_alleles[a][block],
                                                                panel.template_alleles[b][block]);
                extend(block + 1, next_emission, next_transition);
            }
        }
    };
    extend(0, 0.0, 0.0);
    return result;
}

panvar::PreparedProposalEmissions make_emissions(const std::vector<std::vector<double>>& tables,
                                                 std::size_t alleles) {
    panvar::PreparedProposalEmissions emissions;
    emissions.total_alleles = alleles * tables.size();
    emissions.scored_alleles = emissions.total_alleles;
    emissions.blocks.resize(tables.size());
    for (std::size_t b = 0; b < tables.size(); ++b) {
        emissions.blocks[b].allele_count = alleles;
        emissions.blocks[b].scored_alleles = alleles;
        emissions.blocks[b].log_score = tables[b];
    }
    return emissions;
}

panvar::MosaicBeamParameters unbounded(std::size_t target) {
    panvar::MosaicBeamParameters parameters;
    parameters.unique_candidate_target = target;
    parameters.state_limit = 1u << 16;
    return parameters;
}

// Every candidate, against both oracles.
void compare(const panvar::PanelAlleleMatrix& panel, const panvar::PreparedProposalEmissions& emissions,
             const std::vector<double>& recombination, const std::string& label) {
    const auto oracle = exhaustive(panel, emissions, recombination);
    const auto grouped =
        panvar::propose_mosaic_diplotypes_grouped(panel, emissions, recombination, unbounded(oracle.size() + 5));
    const auto existing =
        panvar::propose_mosaic_diplotypes(panel, emissions, recombination, unbounded(oracle.size() + 5));

    require(grouped.search_exhaustive && grouped.states_pruned == 0,
            label + ": an unbounded grouped search must not prune");
    require(grouped.candidates.size() == oracle.size(),
            label + ": the grouped beam must return every exhaustive diplotype, got " +
                std::to_string(grouped.candidates.size()) + " against " + std::to_string(oracle.size()));
    require(grouped.candidates.size() == existing.candidates.size(),
            label + ": the two implementations must return the same number of candidates");

    // Candidates are matched by IDENTITY, not by position. The two implementations sum the same
    // histories in different association orders, so a mass can differ by one unit in the last place --
    // measured, on this very fixture -- and two candidates that tie bitwise in one implementation are
    // separated by that ULP in the other, where the documented allele-path tie-break then never
    // engages. Position is therefore not comparable bitwise; the order gate below asserts only what the
    // data actually determines.
    std::map<OracleKey, std::size_t> existing_at;
    for (std::size_t i = 0; i < existing.candidates.size(); ++i) {
        existing_at[{existing.candidates[i].diplotype.first.alleles,
                     existing.candidates[i].diplotype.second.alleles}] = i;
    }
    require(existing_at.size() == existing.candidates.size(),
            label + ": fixture: the existing implementation must return distinct candidates");

    for (std::size_t i = 0; i < grouped.candidates.size(); ++i) {
        const panvar::MosaicCandidate& got = grouped.candidates[i];
        const OracleKey key{got.diplotype.first.alleles, got.diplotype.second.alleles};
        const std::string at = label + " [" + std::to_string(i) + "]";

        const auto where = existing_at.find(key);
        require(where != existing_at.end(),
                at + ": every grouped candidate must be one the existing implementation returned");
        const panvar::MosaicCandidate& want = existing.candidates[where->second];
        require(got.candidate_id == i && got.unique_rank == i + 1, at + ": ranks must be dense and in order");

        const auto found = oracle.find(key);
        require(found != oracle.end(), at + ": candidate must exist in the exhaustive enumeration");
        const OracleValue& truth = found->second;

        require_close(got.retained_logsum_proposal_score, truth.logsum, 1e-9,
                      at + ": sum-product mass must equal the exhaustive log-sum over histories");
        require_close(got.retained_logsum_proposal_score, want.retained_logsum_proposal_score, 1e-9,
                      at + ": sum-product mass must equal the existing implementation");
        require_close(got.best_proposal_score, truth.best, 1e-9,
                      at + ": max-product score must equal the exhaustive best single history");
        require_close(got.best_proposal_score, want.best_proposal_score, 1e-9,
                      at + ": max-product score must equal the existing implementation");
        require_close(got.proposal_emission_component, truth.emission, 1e-9,
                      at + ": the best history's emission component must match");
        require_close(got.proposal_transition_component, truth.transition_score, 1e-9,
                      at + ": the best history's transition component must match");
        require(got.equivalent_retained_histories == truth.histories,
                at + ": history multiplicity must be exact, got " +
                    std::to_string(got.equivalent_retained_histories) + " against " +
                    std::to_string(truth.histories));
        require(got.equivalent_retained_histories == want.equivalent_retained_histories,
                at + ": history multiplicity must equal the existing implementation");
        require(got.history_count_saturated == want.history_count_saturated,
                at + ": the saturation flag must equal the existing implementation");
        require(!got.history_count_saturated, at + ": fixture: nothing here may saturate a uint64 count");
        require(got.every_allele_scored == want.every_allele_scored, at + ": scored-allele flag must match");

        // The ordering rule itself, not just agreement with the other implementation.
        if (i > 0) {
            const panvar::MosaicCandidate& previous = grouped.candidates[i - 1];
            const bool ordered =
                previous.retained_logsum_proposal_score > got.retained_logsum_proposal_score ||
                (previous.retained_logsum_proposal_score == got.retained_logsum_proposal_score &&
                 (previous.best_proposal_score > got.best_proposal_score ||
                  (previous.best_proposal_score == got.best_proposal_score &&
                   std::tie(previous.diplotype.first.alleles, previous.diplotype.second.alleles) <
                       std::tie(got.diplotype.first.alleles, got.diplotype.second.alleles))));
            require(ordered, at + ": candidates must follow mass, then witness, then allele-path order");
        }
    }

    // Cross-implementation order, asserted only between candidates the data separates. Where two masses
    // agree to within the comparison tolerance their relative order is not determined by the model, and
    // an ULP of association order decides it; where they genuinely differ, both implementations must
    // agree on which comes first.
    for (std::size_t i = 0; i < grouped.candidates.size(); ++i) {
        for (std::size_t j = i + 1; j < grouped.candidates.size(); ++j) {
            const double gap = grouped.candidates[i].retained_logsum_proposal_score -
                               grouped.candidates[j].retained_logsum_proposal_score;
            if (!(gap > 1e-9)) continue;
            const OracleKey earlier{grouped.candidates[i].diplotype.first.alleles,
                                    grouped.candidates[i].diplotype.second.alleles};
            const OracleKey later{grouped.candidates[j].diplotype.first.alleles,
                                  grouped.candidates[j].diplotype.second.alleles};
            require(existing_at.at(earlier) < existing_at.at(later),
                    label + ": candidates the data separates must be ordered the same way by both "
                            "implementations");
        }
    }
}

panvar::PanelAlleleMatrix two_allele_panel() {
    panvar::PanelAlleleMatrix panel;
    panel.template_names = {"h0", "h1", "h2", "h3"};
    panel.template_alleles = {{0, 0, 0}, {0, 1, 1}, {1, 0, 1}, {1, 1, 0}};
    return panel;
}

// Three alleles, classes of different sizes, and templates that change class between blocks -- the case
// a single template scan with an allele-class scatter has to get right.
panvar::PanelAlleleMatrix three_allele_panel() {
    panvar::PanelAlleleMatrix panel;
    panel.template_names = {"t0", "t1", "t2", "t3", "t4"};
    panel.template_alleles = {{0, 0, 0}, {1, 1, 1}, {2, 2, 2}, {0, 1, 2}, {2, 1, 0}};
    return panel;
}

void test_matches_the_oracles() {
    const panvar::PreparedProposalEmissions two = make_emissions(
        {{0.0, -0.2, -0.2, -1.0}, {-0.3, 0.0, 0.0, -0.7}, {-0.6, 0.1, 0.1, -0.4}}, 2);
    const panvar::PanelAlleleMatrix panel = two_allele_panel();
    compare(panel, two, {0.0, 0.0}, "two alleles, r=0");
    compare(panel, two, {0.25, 0.6}, "two alleles, intermediate r");
    compare(panel, two, {1.0, 1.0}, "two alleles, r=1");

    const panvar::PreparedProposalEmissions three = make_emissions(
        {{0.0, -0.2, -0.5, -0.2, -0.9, -0.3, -0.5, -0.3, -1.1},
         {-0.3, 0.0, -0.4, 0.0, -0.7, -0.15, -0.4, -0.15, -0.8},
         {-0.6, 0.1, -0.25, 0.1, -0.4, 0.05, -0.25, 0.05, -0.55}},
        3);
    const panvar::PanelAlleleMatrix wide = three_allele_panel();
    compare(wide, three, {0.0, 0.0}, "three alleles, r=0");
    compare(wide, three, {0.25, 0.6}, "three alleles, intermediate r");
    compare(wide, three, {1.0, 1.0}, "three alleles, r=1");
    std::cout << "  grouped proposer matches both oracles at r=0, intermediate r and r=1\n";
}

// Mutation gates. Each builds what the comparison WOULD have expected under a specific wrong
// implementation, and requires that it disagrees with the real result -- so the assertions above are
// shown to discriminate rather than merely to pass.
void test_mutations_are_caught() {
    const panvar::PreparedProposalEmissions emissions = make_emissions(
        {{0.0, -0.2, -0.2, -1.0}, {-0.3, 0.0, 0.0, -0.7}, {-0.6, 0.1, 0.1, -0.4}}, 2);
    const panvar::PanelAlleleMatrix panel = two_allele_panel();
    const std::vector<double> recombination = {0.25, 0.6};
    const auto oracle = exhaustive(panel, emissions, recombination);
    const auto actual =
        panvar::propose_mosaic_diplotypes_grouped(panel, emissions, recombination, unbounded(oracle.size() + 5));

    // Mutation 1: the heterozygous factor of two omitted. An unordered pair of DISTINCT prefixes stands
    // for two ordered histories; dropping that halves its mass and its multiplicity.
    bool factor_of_two_bites = false;
    for (const panvar::MosaicCandidate& candidate : actual.candidates) {
        const auto found = oracle.find({candidate.diplotype.first.alleles, candidate.diplotype.second.alleles});
        require(found != oracle.end(), "mutation fixture: candidate must exist in the oracle");
        if (candidate.diplotype.first.alleles == candidate.diplotype.second.alleles) continue;
        const double without = found->second.logsum - std::log(2.0);
        const std::uint64_t halved = found->second.histories / 2;
        if (std::abs(candidate.retained_logsum_proposal_score - without) > 1e-9 ||
            candidate.equivalent_retained_histories != halved) {
            factor_of_two_bites = true;
        }
    }
    require(factor_of_two_bites,
            "the gates must fail an implementation omitting the heterozygous factor of two");

    // Mutation 2: only the best template attribution carried, so the group's mass would be its best
    // single history rather than the sum over every attribution.
    bool best_only_bites = false;
    for (const panvar::MosaicCandidate& candidate : actual.candidates) {
        const auto found = oracle.find({candidate.diplotype.first.alleles, candidate.diplotype.second.alleles});
        if (std::abs(candidate.retained_logsum_proposal_score - found->second.best) > 1e-9) {
            best_only_bites = true;
        }
    }
    require(best_only_bites,
            "the gates must fail an implementation carrying only the best template attribution");
    std::cout << "  mutation gates bite: heterozygous factor of two, and best-attribution-only\n";
}

// Pre-allocation is transactional: the budgets are computed from the same loops the search runs and are
// checked before anything is allocated, so a refused search refuses outright.
void test_preallocation_refuses_before_allocating() {
    const panvar::PreparedProposalEmissions emissions = make_emissions(
        {{0.0, -0.2, -0.2, -1.0}, {-0.3, 0.0, 0.0, -0.7}, {-0.6, 0.1, 0.1, -0.4}}, 2);
    const panvar::PanelAlleleMatrix panel = two_allele_panel();
    const std::vector<std::size_t> alleles = {2, 2, 2};

    const std::uint64_t work = panvar::grouped_extension_upper_bound(panel.template_count(), 64, alleles);
    const std::uint64_t cells = panvar::grouped_prefix_vector_cells(panel.template_count(), 64);
    require(work > 0 && cells > 0, "both budgets must be positive for a real search");

    panvar::MosaicBeamParameters admitted;
    admitted.unique_candidate_target = 8;
    admitted.state_limit = 64;
    admitted.max_extension_evaluations = work;
    admitted.max_prefix_vector_cells = cells;
    const auto ran = panvar::propose_mosaic_diplotypes_grouped(panel, emissions, {0.25, 0.6}, admitted);
    require(!ran.candidates.empty(), "a search admitted at exactly its budget must run");
    require(ran.extensions_evaluated <= work,
            "the work actually done must not exceed the budget computed from the same loops");

    for (const bool starve_work : {true, false}) {
        panvar::MosaicBeamParameters refused = admitted;
        if (starve_work) {
            refused.max_extension_evaluations = work - 1;
        } else {
            refused.max_prefix_vector_cells = cells - 1;
        }
        bool threw = false;
        try {
            panvar::propose_mosaic_diplotypes_grouped(panel, emissions, {0.25, 0.6}, refused);
        } catch (const std::length_error&) {
            threw = true;
        }
        require(threw, std::string("a budget below the declared ") + (starve_work ? "work" : "vector cells") +
                           " must refuse before searching");
    }
    std::cout << "  pre-allocation: admitted at budget, refused one below, on both counts\n";
}

}   // namespace

int main() {
    try {
        test_matches_the_oracles();
        test_mutations_are_caught();
        test_preallocation_refuses_before_allocating();
        std::cout << "genotype grouped beam: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "genotype grouped beam: FAIL: " << error.what() << '\n';
        return 1;
    }
}
