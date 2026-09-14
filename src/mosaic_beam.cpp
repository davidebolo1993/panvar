#include "panvar/mosaic_beam.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace panvar {
namespace {

constexpr double kNegInf = -std::numeric_limits<double>::infinity();
constexpr std::uint32_t kNoPrevious = std::numeric_limits<std::uint32_t>::max();

struct BeamNode {
    double score = kNegInf;
    double emission_score = 0.0;
    double transition_score = 0.0;
    std::uint32_t previous = kNoPrevious;
    std::uint32_t template1 = 0;
    std::uint32_t template2 = 0;
};

struct CandidateKey {
    std::vector<std::size_t> first;
    std::vector<std::size_t> second;

    bool operator<(const CandidateKey& other) const noexcept {
        return std::tie(first, second) < std::tie(other.first, other.second);
    }
};

struct CandidateAccumulator {
    double best_score = kNegInf;
    double logsum_score = kNegInf;
    double best_emission = 0.0;
    double best_transition = 0.0;
    std::size_t histories = 0;
};

struct BeamRun {
    std::vector<std::vector<BeamNode>> layers;
    std::uint64_t extensions = 0;
};

double log_add_exp(double lhs, double rhs) {
    if (lhs == kNegInf) return rhs;
    if (rhs == kNegInf) return lhs;
    const double hi = std::max(lhs, rhs);
    const double lo = std::min(lhs, rhs);
    return hi + std::log1p(std::exp(lo - hi));
}

bool better_node(const BeamNode& lhs, const BeamNode& rhs) {
    if (lhs.score != rhs.score) return lhs.score > rhs.score;
    if (lhs.previous != rhs.previous) return lhs.previous < rhs.previous;
    if (lhs.template1 != rhs.template1) return lhs.template1 < rhs.template1;
    return lhs.template2 < rhs.template2;
}

void keep_best(std::vector<BeamNode>& nodes, std::size_t width) {
    if (nodes.size() > width) {
        std::nth_element(nodes.begin(), nodes.begin() + static_cast<std::ptrdiff_t>(width),
                         nodes.end(), better_node);
        nodes.resize(width);
    }
    std::sort(nodes.begin(), nodes.end(), better_node);
}

void add_extension_count(std::uint64_t& total, std::uint64_t amount) {
    if (amount > std::numeric_limits<std::uint64_t>::max() - total) {
        throw std::length_error("mosaic beam extension count overflows uint64");
    }
    total += amount;
}

BeamRun run_beam(
    const PanelAlleleMatrix& panel,
    const PreparedProposalEmissions& emissions,
    const std::vector<std::vector<std::size_t>>& allele_by_block,
    const std::vector<std::vector<double>>& log_transition_by_boundary,
    std::size_t width,
    std::size_t chunk_size) {
    const std::size_t templates = panel.template_count();
    const std::size_t blocks = panel.block_count();
    BeamRun result;
    result.layers.resize(blocks);

    std::vector<BeamNode>& initial = result.layers[0];
    if (templates > std::numeric_limits<std::uint32_t>::max()) {
        throw std::length_error("mosaic beam template index exceeds uint32");
    }
    const double initial_prior = -2.0 * std::log(static_cast<double>(templates));
    const BlockProposalEmissions& initial_emission = emissions.blocks[0];
    for (std::size_t h1 = 0; h1 < templates; ++h1) {
        for (std::size_t h2 = 0; h2 < templates; ++h2) {
            add_extension_count(result.extensions, 1);
            const double emission = initial_emission.log_score[
                allele_by_block[0][h1] * initial_emission.allele_count +
                allele_by_block[0][h2]];
            const double score = initial_prior + emission;
            if (!std::isfinite(score)) continue;
            initial.push_back({score, emission, initial_prior, kNoPrevious,
                               static_cast<std::uint32_t>(h1),
                               static_cast<std::uint32_t>(h2)});
        }
    }
    keep_best(initial, width);

    for (std::size_t block = 1; block < blocks; ++block) {
        const std::vector<BeamNode>& previous = result.layers[block - 1];
        std::vector<BeamNode> best;
        std::vector<BeamNode> pool;
        pool.reserve(width + chunk_size);
        const BlockProposalEmissions& block_emission = emissions.blocks[block];
        const std::vector<double>& log_transition = log_transition_by_boundary[block - 1];

        const auto flush = [&] {
            if (pool.empty()) return;
            if (!best.empty()) {
                pool.insert(pool.end(), best.begin(), best.end());
            }
            keep_best(pool, width);
            best.swap(pool);
            pool.clear();
            pool.reserve(width + chunk_size);
        };

        for (std::size_t previous_id = 0; previous_id < previous.size(); ++previous_id) {
            if (previous_id > std::numeric_limits<std::uint32_t>::max()) {
                throw std::length_error("mosaic beam backpointer exceeds uint32");
            }
            for (std::size_t h1 = 0; h1 < templates; ++h1) {
                const double transition1 = log_transition[
                    static_cast<std::size_t>(previous[previous_id].template1) * templates + h1];
                if (transition1 == kNegInf) {
                    add_extension_count(result.extensions, templates);
                    continue;
                }
                for (std::size_t h2 = 0; h2 < templates; ++h2) {
                    add_extension_count(result.extensions, 1);
                    const double transition2 = log_transition[
                        static_cast<std::size_t>(previous[previous_id].template2) * templates + h2];
                    if (transition2 == kNegInf) continue;
                    const double emission = block_emission.log_score[
                        allele_by_block[block][h1] * block_emission.allele_count +
                        allele_by_block[block][h2]];
                    const double transition = transition1 + transition2;
                    const double score = previous[previous_id].score + transition + emission;
                    if (!std::isfinite(score)) continue;
                    pool.push_back({score,
                                    previous[previous_id].emission_score + emission,
                                    previous[previous_id].transition_score + transition,
                                    static_cast<std::uint32_t>(previous_id),
                                    static_cast<std::uint32_t>(h1),
                                    static_cast<std::uint32_t>(h2)});
                    if (pool.size() >= chunk_size) flush();
                }
            }
        }
        flush();
        // `best` carries the expansion pool's buffer after the swap in flush(), whose capacity is
        // the chunk size. Layers are retained for every block, so moving that buffer in would keep
        // one chunk-sized allocation per block alive: at C4's complete decomposition that is
        // gigabytes of untouched capacity. Store exactly the retained nodes instead.
        result.layers[block].assign(best.begin(), best.end());
        best.clear();
        best.shrink_to_fit();
        if (result.layers[block].empty()) break;
    }
    return result;
}

std::uint64_t extension_upper_bound(
    std::size_t templates,
    std::size_t blocks,
    std::size_t width) {
    if (templates > std::numeric_limits<std::uint64_t>::max() / templates) {
        throw std::length_error("mosaic beam template square overflows uint64");
    }
    const std::uint64_t template_pairs =
        static_cast<std::uint64_t>(templates) * static_cast<std::uint64_t>(templates);
    std::uint64_t result = template_pairs;
    if (blocks <= 1) return result;
    if (width > std::numeric_limits<std::uint64_t>::max() / template_pairs) {
        throw std::length_error("mosaic beam per-block extension bound overflows uint64");
    }
    const std::uint64_t per_block = static_cast<std::uint64_t>(width) * template_pairs;
    if (blocks - 1 > std::numeric_limits<std::uint64_t>::max() / per_block) {
        throw std::length_error("mosaic beam total extension bound overflows uint64");
    }
    const std::uint64_t tail = static_cast<std::uint64_t>(blocks - 1) * per_block;
    if (tail > std::numeric_limits<std::uint64_t>::max() - result) {
        throw std::length_error("mosaic beam total extension bound overflows uint64");
    }
    return result + tail;
}

std::map<CandidateKey, CandidateAccumulator> collapse_final_histories(
    const PanelAlleleMatrix& panel,
    const BeamRun& run) {
    std::map<CandidateKey, CandidateAccumulator> candidates;
    if (run.layers.empty() || run.layers.back().empty()) return candidates;
    const std::size_t blocks = panel.block_count();
    const std::vector<BeamNode>& final = run.layers.back();
    for (std::size_t final_id = 0; final_id < final.size(); ++final_id) {
        CandidateKey key;
        key.first.resize(blocks);
        key.second.resize(blocks);
        std::size_t node_id = final_id;
        for (std::size_t reverse_block = blocks; reverse_block > 0; --reverse_block) {
            const std::size_t block = reverse_block - 1;
            if (node_id >= run.layers[block].size()) {
                throw std::logic_error("mosaic beam backpointer is outside its layer");
            }
            const BeamNode& node = run.layers[block][node_id];
            key.first[block] = panel.template_alleles[node.template1][block];
            key.second[block] = panel.template_alleles[node.template2][block];
            if (block > 0) node_id = node.previous;
        }
        if (key.second < key.first) std::swap(key.first, key.second);

        CandidateAccumulator& accumulator = candidates[key];
        const BeamNode& best = final[final_id];
        accumulator.logsum_score = log_add_exp(accumulator.logsum_score, best.score);
        ++accumulator.histories;
        if (best.score > accumulator.best_score) {
            accumulator.best_score = best.score;
            accumulator.best_emission = best.emission_score;
            accumulator.best_transition = best.transition_score;
        }
    }
    return candidates;
}

std::vector<MosaicCandidate> materialize_candidates(
    const std::map<CandidateKey, CandidateAccumulator>& collapsed,
    bool every_allele_scored,
    std::size_t target) {
    std::vector<MosaicCandidate> result;
    result.reserve(collapsed.size());
    for (const auto& item : collapsed) {
        MosaicCandidate candidate;
        candidate.diplotype.first.alleles = item.first.first;
        candidate.diplotype.second.alleles = item.first.second;
        candidate.best_proposal_score = item.second.best_score;
        candidate.retained_logsum_proposal_score = item.second.logsum_score;
        candidate.proposal_emission_component = item.second.best_emission;
        candidate.proposal_transition_component = item.second.best_transition;
        candidate.equivalent_retained_histories = item.second.histories;
        candidate.every_allele_scored = every_allele_scored;
        result.push_back(std::move(candidate));
    }
    std::sort(result.begin(), result.end(), [](const MosaicCandidate& lhs,
                                               const MosaicCandidate& rhs) {
        if (lhs.retained_logsum_proposal_score != rhs.retained_logsum_proposal_score) {
            return lhs.retained_logsum_proposal_score > rhs.retained_logsum_proposal_score;
        }
        if (lhs.best_proposal_score != rhs.best_proposal_score) {
            return lhs.best_proposal_score > rhs.best_proposal_score;
        }
        if (lhs.diplotype.first.alleles != rhs.diplotype.first.alleles) {
            return lhs.diplotype.first.alleles < rhs.diplotype.first.alleles;
        }
        return lhs.diplotype.second.alleles < rhs.diplotype.second.alleles;
    });
    if (result.size() > target) result.resize(target);
    for (std::size_t rank = 0; rank < result.size(); ++rank) {
        result[rank].candidate_id = rank;
        result[rank].unique_rank = rank + 1;
    }
    return result;
}

void validate_beam_inputs(
    const PanelAlleleMatrix& panel,
    const PreparedProposalEmissions& emissions,
    const std::vector<double>& switch_probability,
    const MosaicBeamParameters& parameters) {
    validate_panel(panel);
    if (emissions.blocks.size() != panel.block_count()) {
        throw std::invalid_argument("proposal emissions and panel have different block counts");
    }
    if (switch_probability.size() + 1 != panel.block_count()) {
        throw std::invalid_argument("mosaic beam needs one switch probability per boundary");
    }
    for (double r : switch_probability) {
        if (!std::isfinite(r) || r < 0.0 || r > 1.0) {
            throw std::invalid_argument("mosaic beam switch probabilities must be in [0,1]");
        }
    }
    if (parameters.unique_candidate_target == 0 || parameters.raw_history_limit == 0) {
        throw std::invalid_argument("mosaic beam candidate target and raw limit must be positive");
    }
    if (parameters.expansion_chunk_size == 0) {
        throw std::invalid_argument("mosaic beam expansion chunk size must be positive");
    }
    if (parameters.raw_history_limit > std::numeric_limits<std::uint32_t>::max()) {
        throw std::invalid_argument("mosaic beam raw history limit exceeds uint32 backpointers");
    }
    for (std::size_t block = 0; block < panel.block_count(); ++block) {
        const BlockProposalEmissions& emission = emissions.blocks[block];
        if (emission.allele_count == 0 ||
            emission.allele_count > std::numeric_limits<std::size_t>::max() /
                                        emission.allele_count ||
            emission.log_score.size() != emission.allele_count * emission.allele_count) {
            throw std::invalid_argument("proposal emission table is malformed");
        }
        for (double value : emission.log_score) {
            if (!std::isfinite(value)) {
                throw std::invalid_argument("proposal emission table contains a non-finite value");
            }
        }
        for (const auto& row : panel.template_alleles) {
            if (row[block] >= emission.allele_count) {
                throw std::invalid_argument("panel allele is outside proposal emission table");
            }
        }
    }
}

} // namespace

MosaicBeamResult propose_mosaic_diplotypes(
    const PanelAlleleMatrix& panel,
    const PreparedProposalEmissions& emissions,
    const std::vector<double>& switch_probability,
    const MosaicBeamParameters& parameters) {
    validate_beam_inputs(panel, emissions, switch_probability, parameters);

    MosaicBeamResult result;
    const std::size_t templates = panel.template_count();
    const std::size_t blocks = panel.block_count();
    if (templates > std::numeric_limits<std::size_t>::max() / templates) {
        throw std::length_error("mosaic beam transition table size overflows size_t");
    }

    // These two block-major tables keep vector-of-vector panel lookups and logarithms out of the
    // tens-of-millions extension loop. They are shared across adaptive-width reruns.
    std::vector<std::vector<std::size_t>> allele_by_block(
        blocks, std::vector<std::size_t>(templates, 0));
    for (std::size_t block = 0; block < blocks; ++block) {
        for (std::size_t h = 0; h < templates; ++h) {
            allele_by_block[block][h] = panel.template_alleles[h][block];
        }
    }
    std::vector<std::vector<double>> log_transition_by_boundary;
    log_transition_by_boundary.reserve(switch_probability.size());
    for (double r : switch_probability) {
        std::vector<double> transition_table(templates * templates, kNegInf);
        for (std::size_t from = 0; from < templates; ++from) {
            for (std::size_t to = 0; to < templates; ++to) {
                const double probability = (from == to ? 1.0 - r : 0.0) +
                                           r / static_cast<double>(templates);
                if (probability > 0.0) {
                    transition_table[from * templates + to] = std::log(probability);
                }
            }
        }
        log_transition_by_boundary.push_back(std::move(transition_table));
    }

    std::size_t width = std::min(
        parameters.raw_history_limit,
        std::max<std::size_t>(1, parameters.unique_candidate_target));
    while (true) {
        const std::uint64_t predicted = extension_upper_bound(templates, blocks, width);
        if (result.extension_upper_bound > parameters.max_extension_evaluations ||
            predicted > parameters.max_extension_evaluations - result.extension_upper_bound) {
            throw std::length_error(
                "mosaic beam extension work exceeds the declared pre-allocation budget");
        }
        result.extension_upper_bound += predicted;
        BeamRun run = run_beam(
            panel, emissions, allele_by_block, log_transition_by_boundary,
            width, parameters.expansion_chunk_size);
        add_extension_count(result.extensions_evaluated, run.extensions);
        const auto collapsed = collapse_final_histories(panel, run);
        result.beam_width_used = width;
        result.final_raw_histories = run.layers.empty() ? 0 : run.layers.back().size();
        result.candidates = materialize_candidates(
            collapsed, emissions.every_allele_scored(), parameters.unique_candidate_target);
        result.unique_target_reached =
            collapsed.size() >= parameters.unique_candidate_target;
        if (result.unique_target_reached || width == parameters.raw_history_limit) {
            result.raw_history_limit_reached =
                width == parameters.raw_history_limit && !result.unique_target_reached;
            return result;
        }
        const std::size_t doubled = width > parameters.raw_history_limit / 2
            ? parameters.raw_history_limit
            : width * 2;
        width = std::max(width + 1, doubled);
    }
}

} // namespace panvar
