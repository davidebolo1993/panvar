#include "panvar/genotype_model.hpp"
#include "panvar/marker_model.hpp"
#include "panvar/mosaic_beam.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <set>
#include <cstring>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {

constexpr double kNegInf = -std::numeric_limits<double>::infinity();

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void require_close(double observed, double expected, double tolerance, const std::string& message) {
    if (!std::isfinite(observed) || !std::isfinite(expected) ||
        std::abs(observed - expected) > tolerance) {
        throw std::runtime_error(message + ": observed=" + std::to_string(observed) +
                                 " expected=" + std::to_string(expected));
    }
}

void require_throws(const std::function<void()>& operation, const std::string& message) {
    try {
        operation();
    } catch (const std::invalid_argument&) {
        return;
    }
    throw std::runtime_error(message);
}

std::size_t marker_slot(const panvar::PreparedMarkerPanel& panel, std::uint64_t code) {
    const auto found = std::lower_bound(panel.codes.begin(), panel.codes.end(), code);
    if (found == panel.codes.end() || *found != code) {
        throw std::runtime_error("fixture marker is absent from prepared universe");
    }
    return static_cast<std::size_t>(found - panel.codes.begin());
}

double poisson_kernel(std::uint64_t observed, double mean) {
    return static_cast<double>(observed) * std::log(mean) - mean;
}

panvar::StructuralLocusInput marker_fixture(bool two_varying_blocks = false) {
    panvar::StructuralLocusInput input;
    input.locus_id = "marker-fixture";
    input.blocks = {
        {"b0", {{"a-rich", "AAAAA"}, {"mixed", "AAACC"}, {"a-rich-alias", "AAAAA"}}},
    };
    input.invariant_segments = {"CGT", ""};
    input.templates = {
        {"h0", panvar::TemplateFrameStatus::Complete, {"a-rich"}},
        {"h1", panvar::TemplateFrameStatus::Complete, {"mixed"}},
        {"h2", panvar::TemplateFrameStatus::Complete, {"a-rich-alias"}},
    };
    if (two_varying_blocks) {
        input.blocks.push_back(
            {"b1", {{"b1-a-rich", "AAAAA"}, {"b1-mixed", "AAACC"}}});
        input.invariant_segments = {"CGT", "", ""};
        input.templates[0].allele_source_ids.push_back("b1-a-rich");
        input.templates[1].allele_source_ids.push_back("b1-mixed");
        input.templates[2].allele_source_ids.push_back("b1-a-rich");
    }
    return input;
}

panvar::SequencingRead read(const std::string& name, const std::string& bases) {
    return {name, bases, std::string(bases.size(), 'I')};
}

void test_marker_model() {
    const panvar::PreparedGenotypeModel model =
        panvar::prepare_genotype_model(marker_fixture());
    const panvar::PreparedMarkerPanel panel = panvar::build_marker_panel(model, 3, 3);
    require(panel.total_alleles == 2 && panel.scored_alleles == 2,
            "every canonical block allele must be retained for proposal scoring");
    require(panel.blocks.size() == 1 && panel.blocks[0].marker_slots.size() == 3,
            "hand fixture must retain its three varying canonical 3-mers");
    require(panel.anchor_slots.size() == 1 && panel.anchor_haploid_multiplicity[0] == 1,
            "constant context marker must remain a separate depth anchor");

    // Independent two-bit expectations: AAA=0, AAC=1, ACC=5, canonical(CGT)=ACG=6.
    const std::size_t aaa = marker_slot(panel, 0);
    const std::size_t aac = marker_slot(panel, 1);
    const std::size_t acc = marker_slot(panel, 5);
    const std::size_t acg = marker_slot(panel, 6);
    require(panel.anchor_slots[0] == acg, "hand-coded invariant marker key must match panel path");

    const std::vector<panvar::ReadPair> fragments = {
        {read("f/1", "AAAAA"), read("f/2", "AAACC")},
    };
    const panvar::MarkerObservations observations =
        panvar::count_marker_observations(panel, fragments);
    require(observations.reads == 2 && observations.bases == 10 &&
                observations.syncmer_occurrences == 6 && observations.matched_occurrences == 6,
            "read marker accounting must count every observed occurrence");
    require(observations.counts[aaa] == 4 && observations.counts[aac] == 1 &&
                observations.counts[acc] == 1 && observations.counts[acg] == 0,
            "hand-counted marker multiplicities must match read projection exactly");

    panvar::MarkerModelParameters parameters;
    parameters.kmer_size = 3;
    parameters.syncmer_s = 3;
    parameters.expected_observations_per_copy = 0.5;
    parameters.background_observations = 0.5;
    parameters.provenance = "fixture: explicit";
    const panvar::PreparedProposalEmissions emissions =
        panvar::build_proposal_emissions(panel, observations, parameters);
    require(emissions.every_allele_scored() && emissions.blocks[0].scored_alleles == 2,
            "proposal table must score every allele without a hidden cap");
    for (std::size_t first = 0; first < 2; ++first) {
        for (std::size_t second = 0; second < 2; ++second) {
            require(emissions.blocks[0].at(first, second) ==
                        emissions.blocks[0].at(second, first),
                    "diploid marker emission must be exactly symmetric under homologue swap");
        }
    }

    const double base = 0.5;
    const double expected_aa =
        poisson_kernel(4, base + 0.5 * 6.0) - poisson_kernel(4, base);
    const double expected_am =
        poisson_kernel(4, base + 0.5 * 4.0) - poisson_kernel(4, base) +
        poisson_kernel(1, base + 0.5) - poisson_kernel(1, base) +
        poisson_kernel(1, base + 0.5) - poisson_kernel(1, base);
    const std::size_t a_rich = model.blocks()[0].allele_id_for_source("a-rich");
    const std::size_t mixed = model.blocks()[0].allele_id_for_source("mixed");
    require_close(emissions.blocks[0].at(a_rich, a_rich), expected_aa, 1e-15,
                  "homozygous multiplicity score must equal hand Poisson calculation");
    require_close(emissions.blocks[0].at(a_rich, mixed), expected_am, 1e-15,
                  "heterozygous multiplicity score must equal hand Poisson calculation");

    const double constant = emissions.blocks[0].candidate_independent_log_constant;
    const bool original_order = emissions.blocks[0].at(a_rich, a_rich) <
                                emissions.blocks[0].at(a_rich, mixed);
    const bool shifted_order = emissions.blocks[0].at(a_rich, a_rich) + constant <
                               emissions.blocks[0].at(a_rich, mixed) + constant;
    require(original_order == shifted_order,
            "candidate-independent marker constant must not change proposal ranking");

    const panvar::PreparedGenotypeModel multiblock_model =
        panvar::prepare_genotype_model(marker_fixture(true));
    const panvar::PreparedMarkerPanel multiblock =
        panvar::build_marker_panel(multiblock_model, 3, 3);
    require(multiblock.dropped_multiblock_markers == 3 &&
                multiblock.blocks[0].marker_slots.empty() &&
                multiblock.blocks[1].marker_slots.empty(),
            "a marker varying in two blocks must not be counted in two proposal emissions");

    require_throws([&] { panvar::build_marker_panel(model, 3, 4); },
                   "syncmer s greater than k must refuse");
    panvar::MarkerModelParameters bad = parameters;
    bad.syncmer_s = 0;
    require_throws([&] { panvar::validate_marker_parameters(bad); },
                   "zero syncmer s must refuse instead of producing an empty model");
}

double log_add_exp(double lhs, double rhs) {
    if (lhs == kNegInf) return rhs;
    if (rhs == kNegInf) return lhs;
    const double hi = std::max(lhs, rhs);
    return hi + std::log1p(std::exp(std::min(lhs, rhs) - hi));
}

double transition(std::size_t from, std::size_t to, double r, std::size_t n) {
    const double p = (from == to ? 1.0 - r : 0.0) + r / static_cast<double>(n);
    return p > 0.0 ? std::log(p) : kNegInf;
}

struct OracleKey {
    std::vector<std::size_t> first;
    std::vector<std::size_t> second;
    bool operator<(const OracleKey& other) const noexcept {
        return std::tie(first, second) < std::tie(other.first, other.second);
    }
};

struct OracleValue {
    double best = kNegInf;
    double logsum = kNegInf;
    double emission = 0.0;
    double transition = 0.0;
    std::size_t histories = 0;
};

std::map<OracleKey, OracleValue> exhaustive_candidates(
    const panvar::PanelAlleleMatrix& panel,
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
            value.logsum = log_add_exp(value.logsum, score);
            ++value.histories;
            if (score > value.best) {
                value.best = score;
                value.emission = emission_score;
                value.transition = transition_score;
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
                const double next_emission = emission_score + emissions.blocks[block].at(
                    panel.template_alleles[a][block], panel.template_alleles[b][block]);
                extend(block + 1, next_emission, next_transition);
            }
        }
    };
    extend(0, 0.0, 0.0);
    return result;
}

panvar::PreparedProposalEmissions beam_emissions() {
    panvar::PreparedProposalEmissions emissions;
    emissions.total_alleles = 6;
    emissions.scored_alleles = 6;
    emissions.blocks.resize(3);
    const std::vector<std::vector<double>> values = {
        {0.0, -0.2, -0.2, -1.0},
        {-0.3, 0.0, 0.0, -0.7},
        {-0.6, 0.1, 0.1, -0.4},
    };
    for (std::size_t block = 0; block < 3; ++block) {
        emissions.blocks[block].allele_count = 2;
        emissions.blocks[block].scored_alleles = 2;
        emissions.blocks[block].log_score = values[block];
    }
    return emissions;
}

void compare_to_oracle(
    const panvar::PanelAlleleMatrix& panel,
    const panvar::PreparedProposalEmissions& emissions,
    const std::vector<double>& recombination,
    const std::string& label) {
    const auto oracle = exhaustive_candidates(panel, emissions, recombination);
    panvar::MosaicBeamParameters parameters;
    parameters.unique_candidate_target = oracle.size() + 5;
    parameters.state_limit = 1u << 20;
    const panvar::MosaicBeamResult actual =
        panvar::propose_mosaic_diplotypes(panel, emissions, recombination, parameters);
    require(actual.search_exhaustive && actual.states_pruned == 0,
            label + ": an unbounded merged beam must not prune");
    require(actual.candidates.size() == oracle.size() && actual.unique_candidates_found == oracle.size(),
            label + ": the merged beam must recover every exhaustive unique mosaic");
    require(!actual.unique_target_reached && !actual.candidate_set_conditional,
            label + ": an exhaustive search with fewer than K candidates is complete, not conditional");
    for (const panvar::MosaicCandidate& candidate : actual.candidates) {
        const OracleKey key{candidate.diplotype.first.alleles, candidate.diplotype.second.alleles};
        const auto found = oracle.find(key);
        require(found != oracle.end(), label + ": beam candidate must exist in exhaustive enumeration");
        // Merging sums histories in a different order from the oracle, so the log-sum is compared with a
        // tolerance for that reordering; it is a sum, never a maximum.
        require_close(candidate.retained_logsum_proposal_score, found->second.logsum, 1e-12,
                      label + ": merged mass must equal the exhaustive log-sum over histories");
        require_close(candidate.best_proposal_score, found->second.best, 1e-12,
                      label + ": the witness must be the exhaustive best history");
        require_close(candidate.proposal_emission_component, found->second.emission, 1e-12,
                      label + ": the witness emission component must equal the exhaustive value");
        require_close(candidate.proposal_transition_component, found->second.transition, 1e-12,
                      label + ": the witness transition component must equal the exhaustive value");
        require(candidate.equivalent_retained_histories == found->second.histories &&
                    !candidate.history_count_saturated,
                label + ": merged physical history multiplicity must equal the exhaustive count");
        require(!(candidate.diplotype.second.alleles < candidate.diplotype.first.alleles),
                label + ": the global homologue swap must be canonicalised in the final candidate");
        require(candidate.every_allele_scored, label + ": the manifest must record complete allele scoring");
    }
    bool summed_histories_are_live = false;
    for (const panvar::MosaicCandidate& candidate : actual.candidates) {
        if (candidate.equivalent_retained_histories > 1 &&
            candidate.retained_logsum_proposal_score > candidate.best_proposal_score + 1e-6) {
            summed_histories_are_live = true;
        }
    }
    require(summed_histories_are_live,
            label + ": the fixture must contain merged histories whose sum exceeds their maximum");
}

void test_mosaic_beam() {
    const panvar::PreparedProposalEmissions emissions = beam_emissions();

    panvar::PanelAlleleMatrix three;
    three.template_names = {"h0", "h1", "h2"};
    three.template_alleles = {{0, 0, 0}, {0, 1, 1}, {1, 0, 1}};
    compare_to_oracle(three, emissions, {0.0, 0.0}, "three templates, r=0");
    compare_to_oracle(three, emissions, {0.25, 0.6}, "three templates, intermediate r");
    compare_to_oracle(three, emissions, {1.0, 1.0}, "three templates, r=1");

    // Four templates in which the pairing of a prefix with its own template decides the future: h0 and h3
    // spell different prefixes at block 0 and different alleles at block 1. A state key that
    // canonicalised the prefix pair and the template pair separately would exchange them, and the
    // exhaustive comparison below would fail.
    panvar::PanelAlleleMatrix four;
    four.template_names = {"h0", "h1", "h2", "h3"};
    four.template_alleles = {{0, 0, 0}, {0, 1, 1}, {1, 0, 1}, {1, 1, 0}};
    bool association_matters = false;
    for (std::size_t t = 0; t < 4; ++t) {
        for (std::size_t u = 0; u < 4; ++u) {
            if (four.template_alleles[t][0] != four.template_alleles[u][0] &&
                four.template_alleles[t][1] != four.template_alleles[u][1]) {
                association_matters = true;
            }
        }
    }
    require(association_matters, "fixture: template-prefix association must change some future");
    compare_to_oracle(four, emissions, {0.0, 0.0}, "four templates, r=0");
    compare_to_oracle(four, emissions, {0.25, 0.6}, "four templates, intermediate r");
    compare_to_oracle(four, emissions, {1.0, 1.0}, "four templates, r=1");

    // The transition term is live: replacing the intermediate arm with r=1 changes its oracle.
    const auto intermediate = exhaustive_candidates(three, emissions, {0.25, 0.6});
    const auto mutated_transition = exhaustive_candidates(three, emissions, {1.0, 1.0});
    bool transition_mutation_bites = intermediate.size() != mutated_transition.size();
    for (const auto& item : intermediate) {
        const auto found = mutated_transition.find(item.first);
        if (found == mutated_transition.end() || std::abs(item.second.logsum - found->second.logsum) > 1e-8) {
            transition_mutation_bites = true;
            break;
        }
    }
    require(transition_mutation_bites, "fixture must fail a mutation replacing the Li-Stephens transition");

    // Search status. An exhaustive search that finds more than K returns a conditional set; a pruned
    // search that reaches K returns a conditional set; a pruned search that misses K has not reached its
    // target; only an exhaustive search with at most K candidates is unconditional.
    const std::size_t total = exhaustive_candidates(four, emissions, {0.25, 0.6}).size();
    require(total > 3, "fixture: the four-template panel must have more than three unique mosaics");

    const auto over = panvar::propose_mosaic_diplotypes(four, emissions, {0.25, 0.6}, {2, 1u << 20});
    require(over.search_exhaustive && over.unique_target_reached && over.candidate_set_conditional &&
                over.candidates.size() == 2 && over.unique_candidates_found == total,
            "an exhaustive search that finds more than K must return K and mark the set conditional");

    const auto pruned_hit = panvar::propose_mosaic_diplotypes(four, emissions, {0.25, 0.6}, {2, 6});
    require(!pruned_hit.search_exhaustive && pruned_hit.states_pruned > 0 && pruned_hit.unique_target_reached &&
                pruned_hit.candidate_set_conditional,
            "a pruned search that reaches K must mark its candidate set conditional");

    const auto starved = panvar::propose_mosaic_diplotypes(four, emissions, {0.25, 0.6}, {50, 1});
    require(!starved.search_exhaustive && starved.states_pruned > 0 && !starved.unique_target_reached,
            "a pruned search that misses K must report neither exhaustiveness nor its target");

    // Tracing is a diagnostic and must not change the search: with a target declared, a pruned search
    // returns the same candidates, masses and work. The trace must account for every pruned state, keep a
    // surviving diplotype's states in every layer, and show the layer where a lost diplotype was dropped.
    {
        const auto all = panvar::propose_mosaic_diplotypes(four, emissions, {0.25, 0.6}, {1000, 1u << 20});
        panvar::MosaicBeamParameters narrow{1000, 6};
        const auto untraced = panvar::propose_mosaic_diplotypes(four, emissions, {0.25, 0.6}, narrow);
        require(!untraced.search_exhaustive && untraced.trace.empty(),
                "fixture: the narrow search must prune, and an untraced search must carry no trace");
        const auto same_diplotype = [](const panvar::MosaicDiplotype& a, const panvar::MosaicDiplotype& b) {
            return (a.first == b.first && a.second == b.second) || (a.first == b.second && a.second == b.first);
        };
        const auto proposed = [&](const panvar::MosaicBeamResult& result, const panvar::MosaicDiplotype& d) {
            for (const auto& candidate : result.candidates) {
                if (same_diplotype(candidate.diplotype, d)) return true;
            }
            return false;
        };
        std::optional<panvar::MosaicDiplotype> kept;
        std::optional<panvar::MosaicDiplotype> lost;
        for (const auto& candidate : all.candidates) {
            if (proposed(untraced, candidate.diplotype)) {
                if (!kept) kept = candidate.diplotype;
            } else if (!lost) {
                lost = candidate.diplotype;
            }
        }
        require(kept.has_value() && lost.has_value(),
                "fixture: the narrow search must keep one exhaustive diplotype and prune another");
        for (const auto& item : {std::make_pair(*kept, true), std::make_pair(*lost, false)}) {
            narrow.trace_target = item.first;
            const auto traced = panvar::propose_mosaic_diplotypes(four, emissions, {0.25, 0.6}, narrow);
            bool identical = traced.candidates.size() == untraced.candidates.size() &&
                             traced.states_pruned == untraced.states_pruned &&
                             traced.extensions_evaluated == untraced.extensions_evaluated;
            for (std::size_t i = 0; identical && i < traced.candidates.size(); ++i) {
                identical = traced.candidates[i].diplotype.first == untraced.candidates[i].diplotype.first &&
                            traced.candidates[i].diplotype.second == untraced.candidates[i].diplotype.second &&
                            traced.candidates[i].retained_logsum_proposal_score ==
                                untraced.candidates[i].retained_logsum_proposal_score &&
                            traced.candidates[i].best_proposal_score == untraced.candidates[i].best_proposal_score;
            }
            require(identical, "a traced search must return exactly the untraced search");
            require(traced.trace.size() == four.block_count(), "a trace must record one layer per block");
            std::uint64_t pruned = 0;
            bool always_retained = true;
            bool dropped_where_offered = false;
            for (const auto& layer : traced.trace) {
                pruned += layer.states_offered - layer.states_retained;
                always_retained = always_retained && layer.target_states_retained > 0;
                if (layer.target_states_offered > 0 && layer.target_states_retained == 0) dropped_where_offered = true;
            }
            require(pruned == traced.states_pruned, "the trace must account for every pruned state");
            if (item.second) {
                require(always_retained, "a proposed diplotype must be retained in every traced layer");
            } else {
                require(!always_retained && dropped_where_offered,
                        "a pruned diplotype's trace must show the layer where it was offered and dropped");
            }
        }
    }

    // Ties: with all emissions equal, the order is fully determined by the documented key, and two runs
    // are identical.
    panvar::PreparedProposalEmissions tied = emissions;
    for (auto& block : tied.blocks) std::fill(block.log_score.begin(), block.log_score.end(), 0.0);
    const auto first = panvar::propose_mosaic_diplotypes(four, tied, {1.0, 1.0}, {5, 7});
    const auto second = panvar::propose_mosaic_diplotypes(four, tied, {1.0, 1.0}, {5, 7});
    require(first.candidates.size() == second.candidates.size(), "tie fixture must be deterministic");
    for (std::size_t i = 0; i < first.candidates.size(); ++i) {
        require(first.candidates[i].diplotype.first == second.candidates[i].diplotype.first &&
                    first.candidates[i].diplotype.second == second.candidates[i].diplotype.second,
                "repeated runs must return identical candidates in identical order");
        if (i > 0) {
            const auto& a = first.candidates[i - 1];
            const auto& b = first.candidates[i];
            const bool ordered =
                a.retained_logsum_proposal_score > b.retained_logsum_proposal_score ||
                (a.retained_logsum_proposal_score == b.retained_logsum_proposal_score &&
                 (a.best_proposal_score > b.best_proposal_score ||
                  (a.best_proposal_score == b.best_proposal_score &&
                   std::tie(a.diplotype.first.alleles, a.diplotype.second.alleles) <
                       std::tie(b.diplotype.first.alleles, b.diplotype.second.alleles))));
            require(ordered, "candidates must follow mass, then witness, then allele-path order");
        }
    }

    // All templates spell the same alleles: merging must collapse every history into one mosaic, found
    // exhaustively within a tiny state budget.
    panvar::PanelAlleleMatrix collapsed = four;
    collapsed.template_alleles = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    const auto one = panvar::propose_mosaic_diplotypes(collapsed, tied, {0.3, 0.3}, {4, 16});
    require(one.candidates.size() == 1 && one.search_exhaustive && !one.unique_target_reached &&
                !one.candidate_set_conditional && one.candidates[0].equivalent_retained_histories == 16 * 16 * 16,
            "identical templates must merge into one mosaic carrying every physical history");

    require_throws([&] { panvar::propose_mosaic_diplotypes(three, emissions, {0.2}, {2, 10}); },
                   "a missing boundary transition must refuse");
    require_throws([&] { panvar::propose_mosaic_diplotypes(three, emissions, {0.2, 0.2}, {2, 0}); },
                   "a zero state limit must refuse");
    panvar::MosaicBeamParameters starved_budget{2, 10, 1};
    bool budget_refused = false;
    try {
        panvar::propose_mosaic_diplotypes(three, emissions, {0.2, 0.2}, starved_budget);
    } catch (const std::length_error&) {
        budget_refused = true;
    }
    require(budget_refused, "an extension budget below the declared work must refuse before searching");
}

// Marker depth comes from the invariant anchors, whose diploid multiplicity is known whatever the
// sample's genotype. The previous estimate divided every matched occurrence by twice the size of the
// marker universe, which ignored anchor multiplicity, variable-marker genotype and copy number.
void test_marker_depth_estimation() {
    const auto make = [](const std::vector<std::uint32_t>& multiplicity,
                         const std::vector<std::uint64_t>& counts) {
        panvar::PreparedMarkerPanel panel;
        panvar::MarkerObservations observations;
        for (std::size_t k = 0; k < multiplicity.size(); ++k) {
            panel.codes.push_back(k);
            panel.anchor_slots.push_back(static_cast<std::uint32_t>(k));
            panel.anchor_haploid_multiplicity.push_back(multiplicity[k]);
            observations.counts.push_back(counts[k]);
        }
        return std::make_pair(panel, observations);
    };

    // A declared background below the model floor is refused, not held at the floor. A declared zero
    // used to take a closed form here and then be replaced by 1e-6 before proposal scoring.
    {
        const auto fixture = make({1, 1, 2, 3}, {10, 14, 40, 57});
        bool refused = false;
        try {
            panvar::estimate_marker_depth(fixture.first, fixture.second, 0.0);
        } catch (const std::invalid_argument&) {
            refused = true;
        }
        require(refused, "a declared background below the model floor must refuse");
        const panvar::MarkerDepthEstimate estimate =
            panvar::estimate_marker_depth(fixture.first, fixture.second, panvar::kMarkerBackgroundFloor);
        require(estimate.converged && estimate.background == panvar::kMarkerBackgroundFloor &&
                    !estimate.background_estimated,
                "a declared background at the floor must be kept exactly and not re-estimated");
    }

    // One-dimensional root with a declared positive background.
    {
        const auto fixture = make({1, 2, 2, 4}, {9, 17, 15, 33});
        const double background = 1.5;
        const panvar::MarkerDepthEstimate estimate =
            panvar::estimate_marker_depth(fixture.first, fixture.second, background);
        double score = 0.0;
        const std::vector<double> c = {9, 17, 15, 33};
        const std::vector<double> d = {2, 4, 4, 8};
        for (std::size_t k = 0; k < 4; ++k) {
            score += c[k] * d[k] / (background + estimate.rate_per_copy * d[k]) - d[k];
        }
        require(std::abs(score) < 1e-8, "the declared-background rate must be the likelihood root");
    }

    // Joint recovery from simulated anchors with known generating values.
    {
        std::mt19937_64 rng(1234);
        const double true_background = 3.0;
        const double true_rate = 7.5;
        std::vector<std::uint32_t> multiplicity;
        std::vector<std::uint64_t> counts;
        for (int k = 0; k < 6000; ++k) {
            const std::uint32_t m = 1 + static_cast<std::uint32_t>(k % 3);
            std::poisson_distribution<long> draw(true_background + true_rate * 2.0 * m);
            multiplicity.push_back(m);
            counts.push_back(static_cast<std::uint64_t>(draw(rng)));
        }
        const auto fixture = make(multiplicity, counts);
        const panvar::MarkerDepthEstimate estimate =
            panvar::estimate_marker_depth(fixture.first, fixture.second, std::nullopt);
        require(estimate.background_estimated && estimate.converged && estimate.termination == "interior_stationary",
                "with several anchor multiplicities both parameters must be fitted and converge in the interior");
        require(std::abs(estimate.rate_per_copy - true_rate) < 0.05 * true_rate,
                "the per-copy rate must be recovered from anchors: got " +
                    std::to_string(estimate.rate_per_copy));
        require(std::abs(estimate.background - true_background) < 0.25 * true_background,
                "the background must be recovered from anchors: got " +
                    std::to_string(estimate.background));
        // The old estimator, fed the same anchors, is off by the multiplicity structure.
        double total = 0.0;
        for (const auto c : counts) total += static_cast<double>(c);
        const double old_rate = total / (2.0 * static_cast<double>(counts.size()));
        require(std::abs(old_rate - true_rate) > 0.2 * true_rate,
                "fixture: the multiplicity-blind estimator must be visibly wrong here");
    }

    // The boundary optimum -- the NA19240 case. Counts that fit c = 22m - 4 exactly imply a negative
    // intercept, so the maximum is on the face background = floor: the rate must be the root there, the
    // background gradient must be negative, and the fit must say so rather than hit its iteration cap.
    {
        std::vector<std::uint32_t> multiplicity;
        std::vector<std::uint64_t> counts;
        for (int k = 0; k < 3000; ++k) {
            const std::uint32_t m = 1 + static_cast<std::uint32_t>(k % 3);
            multiplicity.push_back(m);
            counts.push_back(22 * m - 4);
        }
        const auto fixture = make(multiplicity, counts);
        const panvar::MarkerDepthEstimate estimate =
            panvar::estimate_marker_depth(fixture.first, fixture.second, std::nullopt);
        require(estimate.converged && estimate.background_at_floor && estimate.termination == "boundary_kkt",
                "a boundary optimum must converge on the boundary face, got " + estimate.termination);
        require(estimate.background == panvar::kMarkerBackgroundFloor,
                "the fitted boundary background must be the model floor");
        double copy_sum = 0.0;
        for (const auto m : multiplicity) copy_sum += 2.0 * m;
        require(estimate.gradient_background < 0.0 && std::abs(estimate.gradient_rate) <= 1e-8 * copy_sum,
                "the boundary estimate must satisfy KKT: negative background gradient, stationary rate");
        const auto rate_score = [&](double r) {
            double value = -copy_sum;
            for (std::size_t k = 0; k < counts.size(); ++k) {
                const double d = 2.0 * multiplicity[k];
                value += static_cast<double>(counts[k]) * d / (estimate.background + r * d);
            }
            return value;
        };
        require(std::abs(rate_score(estimate.rate_per_copy)) <= 1e-8 * copy_sum,
                "the boundary rate must be the likelihood root with the background on its floor");
    }

    // Refusals: an unidentifiable background, and anchors with no observation at all.
    {
        const auto same = make({2, 2, 2}, {20, 22, 19});
        bool refused = false;
        try {
            panvar::estimate_marker_depth(same.first, same.second, std::nullopt);
        } catch (const std::invalid_argument&) {
            refused = true;
        }
        require(refused, "a single anchor multiplicity must refuse to estimate the background");
        const auto empty = make({1, 2}, {0, 0});
        refused = false;
        try {
            panvar::estimate_marker_depth(empty.first, empty.second, panvar::kMarkerBackgroundFloor);
        } catch (const std::invalid_argument&) {
            refused = true;
        }
        require(refused, "anchors with no observation must refuse rather than report zero depth");
    }
}

// Real C4 invariant-anchor histograms, {haploid multiplicity, observed count, anchors}, dumped by
// genotype_probe anchors: NA20355, whose depth fit was falsely refused (no_ascent) in the round-5 cohort, and
// NA19240, whose optimum is on the background floor.
struct AnchorRow {
    std::uint32_t multiplicity;
    std::uint64_t count;
    std::uint64_t anchors;
};
const std::vector<AnchorRow> kNa20355Anchors = {
#include "genotype_depth_anchors_na20355.inc"
};
const std::vector<AnchorRow> kNa19240Anchors = {
#include "genotype_depth_anchors_na19240.inc"
};

struct ProfileFit {
    double background = 0.0;
    double rate = 0.0;
    long double log_likelihood = 0.0L;
    bool at_floor = false;
};

// An oracle for the joint depth fit that shares nothing with its Newton iteration. At a fixed background the
// rate is the root of the strictly decreasing rate score. The profile log-likelihood's slope in the background
// is the background score at that rate (envelope theorem), and it decreases because the joint log-likelihood is
// concave: the maximum is on the floor when that slope is already non-positive there, and otherwise at its root.
// Both roots are bisected on the histogram in long double (which is double on some platforms).
ProfileFit profile_oracle(const std::vector<AnchorRow>& rows) {
    using ld = long double;
    ld n = 0.0L;
    ld copy_sum = 0.0L;
    for (const AnchorRow& row : rows) {
        n += static_cast<ld>(row.anchors);
        copy_sum += static_cast<ld>(row.anchors) * 2.0L * static_cast<ld>(row.multiplicity);
    }
    const auto rate_at = [&](ld bg) {
        const auto score = [&](ld r) {
            ld value = -copy_sum;
            for (const AnchorRow& row : rows) {
                const ld d = 2.0L * static_cast<ld>(row.multiplicity);
                value += static_cast<ld>(row.anchors) * static_cast<ld>(row.count) * d / (bg + r * d);
            }
            return value;
        };
        ld lo = 0.0L;
        ld hi = 1.0L;
        while (score(hi) > 0.0L) hi *= 2.0L;
        for (int i = 0; i < 400; ++i) {
            const ld mid = 0.5L * (lo + hi);
            (score(mid) > 0.0L ? lo : hi) = mid;
        }
        return 0.5L * (lo + hi);
    };
    const auto slope = [&](ld bg) {
        const ld r = rate_at(bg);
        ld value = -n;
        for (const AnchorRow& row : rows) {
            const ld d = 2.0L * static_cast<ld>(row.multiplicity);
            value += static_cast<ld>(row.anchors) * static_cast<ld>(row.count) / (bg + r * d);
        }
        return value;
    };
    const ld floor = static_cast<ld>(panvar::kMarkerBackgroundFloor);
    ProfileFit fit;
    ld bg = floor;
    fit.at_floor = slope(floor) <= 0.0L;
    if (!fit.at_floor) {
        ld lo = floor;
        ld hi = 1.0L;
        while (slope(hi) > 0.0L) hi *= 2.0L;
        for (int i = 0; i < 200; ++i) {
            const ld mid = 0.5L * (lo + hi);
            (slope(mid) > 0.0L ? lo : hi) = mid;
        }
        bg = 0.5L * (lo + hi);
    }
    const ld r = rate_at(bg);
    ld log_likelihood = -n * bg - r * copy_sum;
    for (const AnchorRow& row : rows) {
        if (row.count == 0) continue;
        const ld d = 2.0L * static_cast<ld>(row.multiplicity);
        log_likelihood += static_cast<ld>(row.anchors) * static_cast<ld>(row.count) * std::log(bg + r * d);
    }
    fit.background = static_cast<double>(bg);
    fit.rate = static_cast<double>(r);
    fit.log_likelihood = log_likelihood;
    return fit;
}

// The production fit must stop by its own rule where the independent oracle puts the maximum: on the same side
// of the floor, at the same rate and background, with the same log-likelihood within its resolution.
void test_depth_fit_matches_profile_oracle() {
    for (const auto& item : {std::make_pair(std::string("NA20355"), &kNa20355Anchors),
                             std::make_pair(std::string("NA19240"), &kNa19240Anchors)}) {
        panvar::PreparedMarkerPanel panel;
        panvar::MarkerObservations observations;
        std::uint32_t slot = 0;
        for (const AnchorRow& row : *item.second) {
            for (std::uint64_t i = 0; i < row.anchors; ++i, ++slot) {
                panel.codes.push_back(slot);
                panel.anchor_slots.push_back(slot);
                panel.anchor_haploid_multiplicity.push_back(row.multiplicity);
                observations.counts.push_back(row.count);
            }
        }
        require(panel.anchor_slots.size() == 25167, "fixture: " + item.first + " must carry every C4 anchor");
        const panvar::MarkerDepthEstimate estimate =
            panvar::estimate_marker_depth(panel, observations, std::nullopt);
        const ProfileFit oracle = profile_oracle(*item.second);
        require(estimate.converged && estimate.newton_decrement <= estimate.likelihood_resolution,
                item.first + ": the depth fit must converge by its decrement rule, stopped by " + estimate.termination);
        require(estimate.background_at_floor == oracle.at_floor &&
                    estimate.termination == (oracle.at_floor ? "boundary_kkt" : "interior_stationary"),
                item.first + ": the fit must stop on the oracle's side of the background floor");
        require(std::abs(estimate.rate_per_copy - oracle.rate) <= 1e-7 * oracle.rate,
                item.first + ": the rate must match the profile oracle");
        require(oracle.at_floor ? estimate.background == panvar::kMarkerBackgroundFloor
                                : std::abs(estimate.background - oracle.background) <= 1e-5 * oracle.background,
                item.first + ": the background must match the profile oracle");
        const double likelihood_gap = std::abs(estimate.log_likelihood - static_cast<double>(oracle.log_likelihood));
        require(likelihood_gap <= std::max(4.0 * estimate.likelihood_resolution, 1e-7),
                item.first + ": the log-likelihood must match the profile oracle within its resolution, off by " +
                    std::to_string(likelihood_gap));
        std::cout << "  depth fit on " << item.first << " anchors: rate " << estimate.rate_per_copy << " (oracle "
                  << oracle.rate << "), background " << estimate.background << " (oracle " << oracle.background
                  << "), " << estimate.termination << " in " << estimate.iterations << " iterations, decrement "
                  << estimate.newton_decrement << " <= resolution " << estimate.likelihood_resolution << '\n';
    }
}

// One background from fit to emission. The NA19240-like boundary fit ends on the model floor, and
// proposal scoring must use exactly that value. The reads observe AAC and ACC once each; the a-rich
// allele carries neither, so a-rich/a-rich gives both markers zero copies and each observation
// contributes log(background). The command used to raise a fitted 1e-9 to 1e-6 before scoring, which
// moves this pair's score by log(1000) per observation; the gate shows the table sees that change.
void test_marker_background_reaches_proposal() {
    panvar::PreparedMarkerPanel anchors;
    panvar::MarkerObservations anchor_counts;
    for (std::uint32_t k = 0; k < 3000; ++k) {
        const std::uint32_t m = 1 + k % 3;
        anchors.codes.push_back(k);
        anchors.anchor_slots.push_back(k);
        anchors.anchor_haploid_multiplicity.push_back(m);
        anchor_counts.counts.push_back(22 * m - 4);
    }
    const panvar::MarkerDepthEstimate depth =
        panvar::estimate_marker_depth(anchors, anchor_counts, std::nullopt);
    require(depth.converged && depth.background_at_floor && depth.background == panvar::kMarkerBackgroundFloor,
            "the boundary depth fit must end on the model floor");
    const panvar::MarkerModelParameters parameters = panvar::proposal_marker_parameters(depth, 3, 3);
    require(parameters.background_observations == depth.background &&
                parameters.expected_observations_per_copy == depth.rate_per_copy,
            "proposal scoring must use the fitted background and rate unchanged");

    const panvar::PreparedGenotypeModel model = panvar::prepare_genotype_model(marker_fixture());
    const panvar::PreparedMarkerPanel panel = panvar::build_marker_panel(model, 3, 3);
    const std::vector<panvar::ReadPair> fragments = {
        {read("f/1", "AAAAA"), read("f/2", "AAACC")},
    };
    const panvar::MarkerObservations observations =
        panvar::count_marker_observations(panel, fragments);
    require(observations.counts[marker_slot(panel, 1)] == 1 && observations.counts[marker_slot(panel, 5)] == 1,
            "fixture: AAC and ACC must each be observed once");
    const std::size_t a_rich = model.blocks()[0].allele_id_for_source("a-rich");
    const std::size_t mixed = model.blocks()[0].allele_id_for_source("mixed");
    const auto contrast = [&](const panvar::MarkerModelParameters& p) {
        const panvar::PreparedProposalEmissions e = panvar::build_proposal_emissions(panel, observations, p);
        return e.blocks[0].at(mixed, mixed) - e.blocks[0].at(a_rich, a_rich);
    };
    const double bg = depth.background;
    const double r = depth.rate_per_copy;
    // mixed/mixed: AAA, AAC, ACC two copies each; a-rich/a-rich: AAA six copies, AAC and ACC none.
    const double expected = poisson_kernel(4, bg + 2.0 * r) - poisson_kernel(4, bg + 6.0 * r) +
                            2.0 * (poisson_kernel(1, bg + 2.0 * r) - poisson_kernel(1, bg));
    require_close(contrast(parameters), expected, 1e-9,
                  "the zero-copy pair contrast must be scored at the fitted background");
    panvar::MarkerModelParameters raised = parameters;
    raised.background_observations = 1000.0 * bg;
    const double shift = contrast(parameters) - contrast(raised);
    require(std::abs(shift - 2.0 * std::log(1000.0)) < 1e-3,
            "a 1000-fold background change must move the zero-copy contrast by log(1000) per observation, got " +
                std::to_string(shift));
    raised.background_observations = 0.5 * panvar::kMarkerBackgroundFloor;
    bool refused = false;
    try {
        panvar::build_proposal_emissions(panel, observations, raised);
    } catch (const std::invalid_argument&) {
        refused = true;
    }
    require(refused, "proposal emissions must refuse a background below the model floor");
}

// The protected panel baseline. A zero-switch panel pair that looks bad at the first block and is the best
// diplotype overall is dropped there by a narrow beam -- the held-out C4 failure -- and can never be
// recovered by a larger K. The panel baseline contains it by construction, and the off-panel mosaics the
// beam adds never duplicate a panel pair.
void test_panel_baseline_retains_a_pruned_pair() {
    panvar::PanelAlleleMatrix panel;
    panel.template_names = {"t0", "t1", "t2", "t3"};
    panel.template_alleles = {{0, 0, 0}, {1, 1, 1}, {2, 2, 2}, {3, 3, 3}};
    panvar::PreparedProposalEmissions emissions;
    emissions.total_alleles = 12;
    emissions.scored_alleles = 12;
    for (std::size_t b = 0; b < 3; ++b) {
        panvar::BlockProposalEmissions block;
        block.allele_count = 4;
        block.scored_alleles = 4;
        block.log_score.assign(16, 0.0);
        // Locally bad at block 0 and best overall. With switches this rare (about -15.2 nats each), a mosaic
        // that reaches the true alleles after block 0 cannot recover the two switches it needs.
        const double truth = b == 0 ? -10.0 : 8.0;
        block.log_score[0 * 4 + 1] = truth;
        block.log_score[1 * 4 + 0] = truth;
        emissions.blocks.push_back(block);
    }
    const panvar::MosaicDiplotype truth{panvar::MosaicPath{{0, 0, 0}}, panvar::MosaicPath{{1, 1, 1}}};
    const auto is_truth = [&](const panvar::MosaicDiplotype& d) {
        return (d.first == truth.first && d.second == truth.second) || (d.first == truth.second && d.second == truth.first);
    };

    const std::vector<double> rare_switches = {1e-6, 1e-6};
    const auto wide = panvar::propose_mosaic_diplotypes(panel, emissions, rare_switches, {64, 1u << 20});
    require(wide.search_exhaustive && is_truth(wide.candidates.front().diplotype),
            "fixture: exhaustively, the true panel pair must be the best diplotype");

    panvar::MosaicBeamParameters narrow{64, 2};
    narrow.trace_target = truth;
    const auto beam = panvar::propose_mosaic_diplotypes(panel, emissions, rare_switches, narrow);
    bool proposed = false;
    for (const auto& candidate : beam.candidates) proposed = proposed || is_truth(candidate.diplotype);
    require(!beam.search_exhaustive && !proposed, "fixture: the narrow beam must not propose the true panel pair");
    require(beam.trace.front().target_states_offered > 0 && beam.trace.front().target_states_retained == 0,
            "fixture: the trace must show the true pair dropped at block 0");

    const auto baseline = panvar::panel_diplotypes(panel);
    require(baseline.size() == 10, "four distinct templates give ten unordered panel diplotypes");
    bool retained = false;
    for (std::size_t i = 0; i < baseline.size(); ++i) {
        retained = retained || is_truth(baseline[i].diplotype);
        require(baseline[i].candidate_id == i &&
                    !(baseline[i].diplotype.second.alleles < baseline[i].diplotype.first.alleles),
                "panel diplotypes must be canonical and numbered in order");
    }
    require(retained, "the panel baseline must contain the pruned true pair");

    const auto mosaics = panvar::off_panel_mosaics(beam, baseline);
    std::size_t beam_panel_pairs = 0;
    for (const auto& candidate : beam.candidates) {
        bool in_panel = false;
        for (const auto& pair : baseline) {
            in_panel = in_panel || (pair.diplotype.first == candidate.diplotype.first &&
                                    pair.diplotype.second == candidate.diplotype.second);
        }
        beam_panel_pairs += in_panel ? 1 : 0;
    }
    require(mosaics.candidates.size() + beam_panel_pairs == beam.candidates.size() &&
                mosaics.unique_candidates_found == mosaics.candidates.size(),
            "off-panel mosaics must be exactly the beam's diplotypes that are not panel pairs");
    for (std::size_t rank = 0; rank < mosaics.candidates.size(); ++rank) {
        const auto& candidate = mosaics.candidates[rank];
        require(candidate.unique_rank == rank + 1 && candidate.candidate_id == baseline.size() + rank,
                "off-panel mosaics must be ranked in proposal order after the panel ids");
        for (const auto& pair : baseline) {
            require(!(pair.diplotype.first == candidate.diplotype.first && pair.diplotype.second == candidate.diplotype.second),
                    "an off-panel mosaic must not duplicate a panel diplotype");
        }
    }

    // Identical template rows give one panel diplotype, not two.
    panvar::PanelAlleleMatrix duplicated = panel;
    duplicated.template_names.push_back("t0-copy");
    duplicated.template_alleles.push_back({0, 0, 0});
    require(panvar::panel_diplotypes(duplicated).size() == 10,
            "templates with identical allele rows must not add panel diplotypes");
}

} // namespace


// ---- I2-E: exact reuse of marker-proposal contributions ----------------------------------------------------------------
// Synthetic panels exercising repeated markers, varying copy multiplicities, fixed copies outside the block, empty
// inventories, zero observations and markers absent from one allele. Reuse on must equal reuse off bit for bit (every
// emission entry, every constant, every count), under successive calls with different observations and parameters, at
// every budget (tiny budgets fall back exactly), and downstream proposals must be identical.
panvar::PreparedMarkerPanel synthetic_panel(std::uint64_t seed, std::size_t blocks) {
    std::mt19937_64 rng(seed);
    panvar::PreparedMarkerPanel panel;
    panel.kmer_size = 31;
    panel.syncmer_s = 15;
    std::uint32_t next_slot = 0;
    for (std::size_t b = 0; b < blocks; ++b) {
        panvar::BlockMarkerInventory block;
        block.block_id = "b" + std::to_string(b);
        const std::size_t markers = b == 1 ? 0 : 5 + rng() % 40;              // block 1: an empty marker inventory
        block.allele_count = b == 2 ? 1 : 2 + rng() % 9;                       // block 2: a single allele
        for (std::size_t m = 0; m < markers; ++m) {
            block.marker_slots.push_back(next_slot++);
            block.fixed_haploid_multiplicity.push_back(static_cast<std::uint32_t>(rng() % 4));   // fixed copies elsewhere
        }
        block.allele_markers.resize(block.allele_count);
        for (std::size_t a = 0; a < block.allele_count; ++a) {
            if (a == 0 && block.allele_count > 2) continue;                    // an allele carrying no marker
            for (std::uint32_t m = 0; m < markers; ++m) {
                if (rng() % 3 == 0) continue;                                  // absent from this allele
                const std::uint32_t copies = 1 + static_cast<std::uint32_t>(rng() % 5 == 0 ? rng() % 7 : rng() % 2);
                block.allele_markers[a].emplace_back(m, copies);
            }
        }
        panel.total_alleles += block.allele_count;
        panel.scored_alleles += block.allele_count;
        panel.blocks.push_back(std::move(block));
    }
    panel.codes.resize(next_slot);
    for (std::uint32_t i = 0; i < next_slot; ++i) panel.codes[i] = i;
    return panel;
}

panvar::MarkerObservations synthetic_observations(const panvar::PreparedMarkerPanel& panel, std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    panvar::MarkerObservations o;
    o.counts.resize(panel.codes.size());
    for (auto& c : o.counts) c = rng() % 4 == 0 ? 0 : rng() % 30;   // zero observations included
    return o;
}

void require_bitwise_emissions(const panvar::PreparedProposalEmissions& a, const panvar::PreparedProposalEmissions& b,
                               const std::string& label) {
    require(a.total_alleles == b.total_alleles && a.scored_alleles == b.scored_alleles && a.blocks.size() == b.blocks.size(),
            label + ": emission metadata differs");
    for (std::size_t k = 0; k < a.blocks.size(); ++k) {
        const auto& x = a.blocks[k];
        const auto& y = b.blocks[k];
        require(x.allele_count == y.allele_count && x.scored_alleles == y.scored_alleles && x.log_score.size() == y.log_score.size(),
                label + ": block metadata differs");
        require(std::memcmp(&x.candidate_independent_log_constant, &y.candidate_independent_log_constant, sizeof(double)) == 0,
                label + ": candidate-independent constant differs");
        require(x.log_score.empty() ||
                    std::memcmp(x.log_score.data(), y.log_score.data(), x.log_score.size() * sizeof(double)) == 0,
                label + ": an emission entry differs in block " + std::to_string(k));
    }
}

std::uint64_t distinct_terms(const panvar::PreparedMarkerPanel& panel) {
    std::uint64_t total = 0;
    for (const auto& block : panel.blocks) {
        std::vector<std::set<std::uint64_t>> single(block.marker_slots.size()), need(block.marker_slots.size());
        for (const auto& allele : block.allele_markers) for (const auto& [m, c] : allele) single[m].insert(c);
        // copies a pair can present at a marker: one allele's alone (the other lacks it) or both alleles' summed, the
        // homozygous pair included -- only for allele pairs that actually occur
        for (std::size_t x = 0; x < block.allele_count; ++x) {
            for (std::size_t y = x; y < block.allele_count; ++y) {
                std::map<std::uint32_t, std::uint64_t> sum;
                for (const auto& [m, c] : block.allele_markers[x]) sum[m] += c;
                for (const auto& [m, c] : block.allele_markers[y]) sum[m] += c;
                for (const auto& [m, c] : sum) need[m].insert(c);
            }
        }
        for (const auto& n : need) total += n.size();
    }
    return total;
}

void test_contribution_reuse() {
    for (const std::uint64_t seed : {11u, 12u, 13u}) {
        const panvar::PreparedMarkerPanel panel = synthetic_panel(seed, 6);
        const std::vector<std::pair<panvar::MarkerObservations, panvar::MarkerModelParameters>> calls = [&] {
            std::vector<std::pair<panvar::MarkerObservations, panvar::MarkerModelParameters>> out;
            for (int k = 0; k < 3; ++k) {
                panvar::MarkerModelParameters p;
                p.kmer_size = 31;
                p.syncmer_s = 15;
                p.expected_observations_per_copy = k == 1 ? 3.7 : 1.25;   // different fitted parameters in successive calls
                p.background_observations = k == 2 ? 0.02 : 1e-6;
                p.provenance = "fixture: synthetic call " + std::to_string(k);
                out.push_back({synthetic_observations(panel, seed * 10 + static_cast<std::uint64_t>(k % 2)), p});
            }
            return out;
        }();
        for (std::size_t c = 0; c < calls.size(); ++c) {
            const auto& [obs, params] = calls[c];
            const std::string label = "seed " + std::to_string(seed) + ", call " + std::to_string(c);
            panvar::ProposalReuseReport off_report, on_report, tiny_report, partial_report;
            const auto direct = panvar::build_proposal_emissions(panel, obs, params, panvar::ProposalReuseOptions{}, &off_report);
            require_bitwise_emissions(direct, panvar::build_proposal_emissions(panel, obs, params), label + ": overloads");
            const auto reused = panvar::build_proposal_emissions(panel, obs, params, {true, std::size_t{1} << 20}, &on_report);
            require_bitwise_emissions(direct, reused, label + ": reuse on");
            require(on_report.executed + on_report.reused == on_report.visits && on_report.visits == off_report.visits &&
                        off_report.executed == off_report.visits && on_report.fallback == 0 && on_report.reused > 0,
                    label + ": accounting");
            require(on_report.executed == distinct_terms(panel),
                    label + ": each required term must be evaluated exactly once");
            // tiny budget: nothing tabled, every term direct, identical
            const auto tiny = panvar::build_proposal_emissions(panel, obs, params, {true, 16}, &tiny_report);
            require_bitwise_emissions(direct, tiny, label + ": tiny budget");
            require(tiny_report.markers_tabled == 0 && tiny_report.fallback == tiny_report.visits && tiny_report.reused == 0 &&
                        tiny_report.peak_table_bytes <= 16,
                    label + ": a tiny budget must fall back for every term and retain nothing beyond it");
            // I2-F correction: the planning arrays are temporary, allocated even when nothing is admitted, and reported
            // apart from the budgeted retained storage
            std::size_t largest_markers = 0;
            for (const auto& block : panel.blocks) largest_markers = std::max(largest_markers, block.marker_slots.size());
            require(tiny_report.peak_planning_bytes == largest_markers * (sizeof(std::uint32_t) + sizeof(std::uint64_t)) &&
                        tiny_report.peak_planning_bytes > 16,
                    label + ": the tiny budget's planning scratch must be reported, outside the budget");
            require(off_report.peak_planning_bytes == 0 && off_report.peak_table_bytes == 0, label + ": reuse off allocates nothing");
            require(on_report.peak_table_bytes <= (std::size_t{1} << 20), label + ": the table must stay within its budget");
            // a partial budget: some markers tabled, the rest direct, identical
            const std::size_t partial = static_cast<std::size_t>(on_report.peak_table_bytes / 2);
            const auto half = panvar::build_proposal_emissions(panel, obs, params, {true, partial}, &partial_report);
            require_bitwise_emissions(direct, half, label + ": partial budget");
            require(partial_report.markers_tabled > 0 && partial_report.markers_not_tabled > 0 && partial_report.fallback > 0 &&
                        partial_report.peak_table_bytes <= partial,
                    label + ": a partial budget must table some markers, fall back for the rest and stay within it");
        }
        // downstream proposals: identical candidates and scores from either emission table
        panvar::PanelAlleleMatrix templates;
        std::mt19937_64 rng(seed + 99);
        for (std::size_t t = 0; t < 5; ++t) {
            templates.template_names.push_back("t" + std::to_string(t));
            std::vector<std::size_t> path;
            for (const auto& block : panel.blocks) path.push_back(rng() % block.allele_count);
            templates.template_alleles.push_back(path);
        }
        const auto& [obs, params] = calls[0];
        const auto a = panvar::propose_mosaic_diplotypes(templates, panvar::build_proposal_emissions(panel, obs, params),
                                                         {0.2, 0.2, 0.2, 0.2, 0.2}, {12, 4096});
        const auto b = panvar::propose_mosaic_diplotypes(
            templates, panvar::build_proposal_emissions(panel, obs, params, {true, std::size_t{1} << 20}), {0.2, 0.2, 0.2, 0.2, 0.2},
            {12, 4096});
        require(a.candidates.size() == b.candidates.size() && a.unique_candidates_found == b.unique_candidates_found &&
                    a.states_pruned == b.states_pruned,
                "downstream proposal shapes must be identical");
        for (std::size_t i = 0; i < a.candidates.size(); ++i) {
            const auto& x = a.candidates[i];
            const auto& y = b.candidates[i];
            require(x.diplotype.first.alleles == y.diplotype.first.alleles && x.diplotype.second.alleles == y.diplotype.second.alleles &&
                        std::memcmp(&x.best_proposal_score, &y.best_proposal_score, sizeof(double)) == 0 &&
                        std::memcmp(&x.retained_logsum_proposal_score, &y.retained_logsum_proposal_score, sizeof(double)) == 0 &&
                        x.equivalent_retained_histories == y.equivalent_retained_histories,
                    "downstream proposals must be identical, rank by rank");
        }
    }
    // the hand fixture's real model, and error handling identical in both modes
    {
        const panvar::PreparedGenotypeModel model = panvar::prepare_genotype_model(marker_fixture(true));
        const panvar::PreparedMarkerPanel panel = panvar::build_marker_panel(model, 3, 3);
        const std::vector<panvar::ReadPair> fragments = {{read("f/1", "AAAAACGTAC"), read("f/2", "AAACCGGTTT")}};
        const panvar::MarkerObservations obs = panvar::count_marker_observations(panel, fragments);
        panvar::MarkerModelParameters p;
        p.kmer_size = 3;
        p.syncmer_s = 3;
        p.expected_observations_per_copy = 0.5;
        p.background_observations = 0.5;
        p.provenance = "fixture: explicit";
        require_bitwise_emissions(panvar::build_proposal_emissions(panel, obs, p),
                                  panvar::build_proposal_emissions(panel, obs, p, {true, 4096}), "hand fixture");
        panvar::PreparedMarkerPanel broken = synthetic_panel(21, 3);
        broken.blocks[0].allele_markers[1].emplace_back(static_cast<std::uint32_t>(broken.blocks[0].marker_slots.size() + 5), 1);
        const panvar::MarkerObservations bobs = synthetic_observations(broken, 5);
        panvar::MarkerModelParameters bp = p;
        bp.kmer_size = 31;
        bp.syncmer_s = 15;
        std::string off_error, on_error;
        try { (void)panvar::build_proposal_emissions(broken, bobs, bp); } catch (const std::exception& e) { off_error = e.what(); }
        try { (void)panvar::build_proposal_emissions(broken, bobs, bp, {true, 1u << 20}); } catch (const std::exception& e) { on_error = e.what(); }
        require(!off_error.empty() && off_error == on_error, "an out-of-block marker must be refused identically in both modes");
        panvar::MarkerModelParameters invalid = bp;
        invalid.expected_observations_per_copy = -1.0;
        off_error.clear();
        on_error.clear();
        try { (void)panvar::build_proposal_emissions(broken, bobs, invalid); } catch (const std::exception& e) { off_error = e.what(); }
        try { (void)panvar::build_proposal_emissions(broken, bobs, invalid, {true, 1u << 20}); } catch (const std::exception& e) { on_error = e.what(); }
        require(!off_error.empty() && off_error == on_error, "invalid parameters must be refused identically in both modes");
    }
    std::cout << "  contribution reuse: bitwise-identical emissions and constants over 3 panels x 3 successive calls (different "
                 "observations and parameters), each required term evaluated once; tiny and partial budgets fall back exactly; "
                 "downstream proposals identical; errors identical\n";
}

int main() {
    try {
        test_marker_model();
        test_mosaic_beam();
        test_marker_depth_estimation();
        test_marker_background_reaches_proposal();
        test_depth_fit_matches_profile_oracle();
        test_panel_baseline_retains_a_pruned_pair();
        test_contribution_reuse();
        std::cout << "genotype mosaic proposal: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "genotype mosaic proposal: FAIL: " << error.what() << '\n';
        return 1;
    }
}
