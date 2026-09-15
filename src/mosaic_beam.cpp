#include "panvar/mosaic_beam.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <utility>

namespace panvar {
namespace {

constexpr double kNegInf = -std::numeric_limits<double>::infinity();
constexpr std::uint32_t kRoot = std::numeric_limits<std::uint32_t>::max();
constexpr std::uint64_t kCountCap = std::numeric_limits<std::uint64_t>::max();

// prefix_id = intern(previous_prefix_id, allele_id). Equal ids are equal allele prefixes, and a prefix is
// reconstructed by walking parents, so no allele vector is ever copied during the search.
class PrefixInterner {
public:
    std::uint32_t intern(std::uint32_t parent, std::size_t allele) {
        const std::uint64_t key =
            (static_cast<std::uint64_t>(parent) << 32) | static_cast<std::uint64_t>(allele);
        const auto found = index_.find(key);
        if (found != index_.end()) return found->second;
        if (parent_.size() >= kRoot) throw std::length_error("mosaic prefix interner exhausted uint32 ids");
        const std::uint32_t id = static_cast<std::uint32_t>(parent_.size());
        parent_.push_back(parent);
        allele_.push_back(static_cast<std::uint32_t>(allele));
        index_.emplace(key, id);
        return id;
    }

    std::vector<std::size_t> path(std::uint32_t id, std::size_t blocks) const {
        std::vector<std::size_t> out(blocks, 0);
        for (std::size_t b = blocks; b > 0; --b) {
            if (id == kRoot) throw std::logic_error("mosaic prefix is shorter than the block chain");
            out[b - 1] = allele_[id];
            id = parent_[id];
        }
        if (id != kRoot) throw std::logic_error("mosaic prefix is longer than the block chain");
        return out;
    }

private:
    std::vector<std::uint32_t> parent_;
    std::vector<std::uint32_t> allele_;
    std::unordered_map<std::uint64_t, std::uint32_t> index_;
};

struct State {
    // Canonical: (prefix1, template1) <= (prefix2, template2). The prefix travels with its own template.
    std::uint32_t prefix1 = 0;
    std::uint32_t template1 = 0;
    std::uint32_t prefix2 = 0;
    std::uint32_t template2 = 0;
    double log_mass = kNegInf;
    double best = kNegInf;
    double best_emission = 0.0;
    double best_transition = 0.0;
    std::uint64_t histories = 0;
    bool saturated = false;
};

double log_add_exp_local(double lhs, double rhs) {
    if (lhs == kNegInf) return rhs;
    if (rhs == kNegInf) return lhs;
    const double hi = std::max(lhs, rhs);
    return hi + std::log1p(std::exp(std::min(lhs, rhs) - hi));
}

void add_count(std::uint64_t& total, bool& saturated, std::uint64_t amount, bool amount_saturated) {
    saturated = saturated || amount_saturated;
    if (amount > kCountCap - total) {
        total = kCountCap;
        saturated = true;
    } else {
        total += amount;
    }
}

std::tuple<std::uint32_t, std::uint32_t, std::uint32_t, std::uint32_t> key_of(const State& state) {
    return std::make_tuple(state.prefix1, state.template1, state.prefix2, state.template2);
}

void canonicalise(State& state) {
    if (std::make_pair(state.prefix2, state.template2) < std::make_pair(state.prefix1, state.template1)) {
        std::swap(state.prefix1, state.prefix2);
        std::swap(state.template1, state.template2);
    }
}

// Total order used for pruning: more mass first, then a better witness, then the key.
bool better_state(const State& lhs, const State& rhs) {
    if (lhs.log_mass != rhs.log_mass) return lhs.log_mass > rhs.log_mass;
    if (lhs.best != rhs.best) return lhs.best > rhs.best;
    return key_of(lhs) < key_of(rhs);
}

void merge_into(State& into, const State& from) {
    into.log_mass = log_add_exp_local(into.log_mass, from.log_mass);
    if (from.best > into.best) {
        into.best = from.best;
        into.best_emission = from.best_emission;
        into.best_transition = from.best_transition;
    }
    add_count(into.histories, into.saturated, from.histories, from.saturated);
}

// Keeps the state_limit best states of a layer after every state has been fully merged. Children of
// different parent prefix pairs never share a key, so a group's children are final when offered here.
class LayerSelector {
public:
    explicit LayerSelector(std::size_t limit) : limit_(limit), heap_(worse_on_top) {}

    void offer(const State& state) {
        ++offered_;
        if (heap_.size() < limit_) {
            heap_.push(state);
        } else if (better_state(state, heap_.top())) {
            heap_.pop();
            heap_.push(state);
        }
    }

    std::vector<State> take(std::uint64_t& pruned) {
        pruned += offered_ - heap_.size();
        std::vector<State> out;
        out.reserve(heap_.size());
        while (!heap_.empty()) {
            out.push_back(heap_.top());
            heap_.pop();
        }
        std::sort(out.begin(), out.end(),
                  [](const State& lhs, const State& rhs) { return key_of(lhs) < key_of(rhs); });
        return out;
    }

private:
    static bool worse_on_top(const State& lhs, const State& rhs) { return better_state(lhs, rhs); }
    std::size_t limit_;
    std::uint64_t offered_ = 0;
    std::priority_queue<State, std::vector<State>, bool (*)(const State&, const State&)> heap_;
};

void add_extensions(std::uint64_t& total, std::uint64_t amount) {
    if (amount > std::numeric_limits<std::uint64_t>::max() - total) {
        throw std::length_error("mosaic beam extension count overflows uint64");
    }
    total += amount;
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
    if (parameters.unique_candidate_target == 0 || parameters.state_limit == 0) {
        throw std::invalid_argument("mosaic beam candidate target and state limit must be positive");
    }
    if (panel.template_count() >= kRoot) {
        throw std::invalid_argument("mosaic beam template index exceeds uint32");
    }
    for (std::size_t block = 0; block < panel.block_count(); ++block) {
        const BlockProposalEmissions& emission = emissions.blocks[block];
        if (emission.allele_count == 0 ||
            emission.allele_count > std::numeric_limits<std::size_t>::max() / emission.allele_count ||
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

std::uint64_t extension_upper_bound(std::size_t templates, std::size_t blocks, std::size_t limit) {
    const std::uint64_t pairs = static_cast<std::uint64_t>(templates) * templates;
    if (blocks <= 1) return pairs;
    if (limit > std::numeric_limits<std::uint64_t>::max() / pairs) {
        throw std::length_error("mosaic beam per-block extension bound overflows uint64");
    }
    const std::uint64_t per_block = static_cast<std::uint64_t>(limit) * pairs;
    if (blocks - 1 > (std::numeric_limits<std::uint64_t>::max() - pairs) / per_block) {
        throw std::length_error("mosaic beam total extension bound overflows uint64");
    }
    return pairs + static_cast<std::uint64_t>(blocks - 1) * per_block;
}

} // namespace

MosaicBeamResult propose_mosaic_diplotypes(
    const PanelAlleleMatrix& panel,
    const PreparedProposalEmissions& emissions,
    const std::vector<double>& switch_probability,
    const MosaicBeamParameters& parameters) {
    validate_beam_inputs(panel, emissions, switch_probability, parameters);
    const std::size_t templates = panel.template_count();
    const std::size_t blocks = panel.block_count();

    MosaicBeamResult result;
    result.state_limit_used = parameters.state_limit;
    result.extension_upper_bound = extension_upper_bound(templates, blocks, parameters.state_limit);
    if (result.extension_upper_bound > parameters.max_extension_evaluations) {
        throw std::length_error("mosaic beam extension work exceeds the declared pre-allocation budget");
    }

    std::vector<std::vector<std::size_t>> allele_by_block(blocks, std::vector<std::size_t>(templates, 0));
    for (std::size_t block = 0; block < blocks; ++block) {
        for (std::size_t h = 0; h < templates; ++h) allele_by_block[block][h] = panel.template_alleles[h][block];
    }

    PrefixInterner prefixes;
    const double log_templates = std::log(static_cast<double>(templates));

    // ---- block 0: every ordered template pair, merged by key -------------------------------------------
    std::vector<State> layer;
    {
        std::vector<State> initial;
        initial.reserve(templates * templates);
        const BlockProposalEmissions& emission = emissions.blocks[0];
        std::vector<std::uint32_t> prefix_of(templates);
        for (std::size_t h = 0; h < templates; ++h) prefix_of[h] = prefixes.intern(kRoot, allele_by_block[0][h]);
        for (std::size_t h1 = 0; h1 < templates; ++h1) {
            for (std::size_t h2 = 0; h2 < templates; ++h2) {
                add_extensions(result.extensions_evaluated, 1);
                const double e = emission.log_score[allele_by_block[0][h1] * emission.allele_count +
                                                    allele_by_block[0][h2]];
                State state;
                state.prefix1 = prefix_of[h1];
                state.template1 = static_cast<std::uint32_t>(h1);
                state.prefix2 = prefix_of[h2];
                state.template2 = static_cast<std::uint32_t>(h2);
                state.log_mass = -2.0 * log_templates + e;
                state.best = state.log_mass;
                state.best_emission = e;
                state.best_transition = -2.0 * log_templates;
                state.histories = 1;
                canonicalise(state);
                initial.push_back(state);
            }
        }
        std::sort(initial.begin(), initial.end(),
                  [](const State& lhs, const State& rhs) { return key_of(lhs) < key_of(rhs); });
        LayerSelector selector(parameters.state_limit);
        for (std::size_t i = 0; i < initial.size();) {
            State merged = initial[i];
            std::size_t j = i + 1;
            for (; j < initial.size() && key_of(initial[j]) == key_of(merged); ++j) merge_into(merged, initial[j]);
            selector.offer(merged);
            i = j;
        }
        layer = selector.take(result.states_pruned);
    }

    // ---- blocks 1..B-1: extend each parent prefix-pair group over every ordered template pair ------------
    std::vector<double> cell_mass(templates * templates, kNegInf);
    std::vector<double> cell_best(templates * templates, kNegInf);
    std::vector<double> cell_best_emission(templates * templates, 0.0);
    std::vector<double> cell_best_transition(templates * templates, 0.0);
    std::vector<std::uint64_t> cell_histories(templates * templates, 0);
    std::vector<char> cell_saturated(templates * templates, 0);
    std::vector<std::size_t> touched;

    for (std::size_t block = 1; block < blocks && !layer.empty(); ++block) {
        const BlockProposalEmissions& emission = emissions.blocks[block];
        const double r = switch_probability[block - 1];
        const double log_same = std::log((1.0 - r) + r / static_cast<double>(templates));
        const double log_other = r > 0.0 ? std::log(r) - log_templates : kNegInf;

        // Group parents by their prefix pair. States are sorted by key, so prefix1 then template1 order;
        // re-sort by (prefix1, prefix2) to make each group contiguous.
        std::sort(layer.begin(), layer.end(), [](const State& lhs, const State& rhs) {
            return std::tie(lhs.prefix1, lhs.prefix2, lhs.template1, lhs.template2) <
                   std::tie(rhs.prefix1, rhs.prefix2, rhs.template1, rhs.template2);
        });

        LayerSelector selector(parameters.state_limit);
        std::vector<std::uint32_t> child_prefix1(emission.allele_count);
        std::vector<std::uint32_t> child_prefix2(emission.allele_count);
        for (std::size_t g = 0; g < layer.size();) {
            std::size_t end = g + 1;
            while (end < layer.size() && layer[end].prefix1 == layer[g].prefix1 &&
                   layer[end].prefix2 == layer[g].prefix2) {
                ++end;
            }
            const std::uint32_t parent1 = layer[g].prefix1;
            const std::uint32_t parent2 = layer[g].prefix2;
            const bool symmetric = parent1 == parent2;
            for (std::size_t a = 0; a < emission.allele_count; ++a) {
                child_prefix1[a] = prefixes.intern(parent1, a);
                child_prefix2[a] = symmetric ? child_prefix1[a] : prefixes.intern(parent2, a);
            }

            // Extending one representative of each merged parent over every ORDERED new template pair is
            // exact: the swapped ordered history, extended by the swapped pair, lands on the same unordered
            // child with the same transition product.
            for (std::size_t p = g; p < end; ++p) {
                const State& parent = layer[p];
                for (std::size_t x = 0; x < templates; ++x) {
                    const double t1 = parent.template1 == x ? log_same : log_other;
                    if (t1 == kNegInf) {
                        continue;
                    }
                    for (std::size_t y = 0; y < templates; ++y) {
                        const double t2 = parent.template2 == y ? log_same : log_other;
                        if (t2 == kNegInf) continue;
                        add_extensions(result.extensions_evaluated, 1);
                        const double e = emission.log_score[allele_by_block[block][x] * emission.allele_count +
                                                            allele_by_block[block][y]];
                        // A symmetric parent's children (x,y) and (y,x) are one unordered child.
                        const std::size_t cell = symmetric ? std::min(x, y) * templates + std::max(x, y)
                                                           : x * templates + y;
                        if (cell_mass[cell] == kNegInf && cell_histories[cell] == 0) touched.push_back(cell);
                        const double transition = t1 + t2;
                        cell_mass[cell] = log_add_exp_local(cell_mass[cell], parent.log_mass + transition + e);
                        const double candidate_best = parent.best + transition + e;
                        if (candidate_best > cell_best[cell]) {
                            cell_best[cell] = candidate_best;
                            cell_best_emission[cell] = parent.best_emission + e;
                            cell_best_transition[cell] = parent.best_transition + transition;
                        }
                        bool saturated = cell_saturated[cell] != 0;
                        add_count(cell_histories[cell], saturated, parent.histories, parent.saturated);
                        cell_saturated[cell] = saturated ? 1 : 0;
                    }
                }
            }

            std::sort(touched.begin(), touched.end());
            for (const std::size_t cell : touched) {
                const std::size_t x = cell / templates;
                const std::size_t y = cell % templates;
                State child;
                child.prefix1 = child_prefix1[allele_by_block[block][x]];
                child.template1 = static_cast<std::uint32_t>(x);
                child.prefix2 = child_prefix2[allele_by_block[block][y]];
                child.template2 = static_cast<std::uint32_t>(y);
                child.log_mass = cell_mass[cell];
                child.best = cell_best[cell];
                child.best_emission = cell_best_emission[cell];
                child.best_transition = cell_best_transition[cell];
                child.histories = cell_histories[cell];
                child.saturated = cell_saturated[cell] != 0;
                canonicalise(child);
                selector.offer(child);
                cell_mass[cell] = kNegInf;
                cell_best[cell] = kNegInf;
                cell_histories[cell] = 0;
                cell_saturated[cell] = 0;
            }
            touched.clear();
            g = end;
        }
        layer = selector.take(result.states_pruned);
    }

    // ---- final collapse: drop templates, merge by unordered pair of complete allele paths ----------------
    result.final_states = layer.size();
    std::sort(layer.begin(), layer.end(), [](const State& lhs, const State& rhs) {
        return std::tie(lhs.prefix1, lhs.prefix2) < std::tie(rhs.prefix1, rhs.prefix2);
    });
    std::vector<MosaicCandidate> candidates;
    for (std::size_t i = 0; i < layer.size();) {
        State merged = layer[i];
        std::size_t j = i + 1;
        for (; j < layer.size() && layer[j].prefix1 == merged.prefix1 && layer[j].prefix2 == merged.prefix2; ++j) {
            merge_into(merged, layer[j]);
        }
        MosaicCandidate candidate;
        candidate.diplotype.first.alleles = prefixes.path(merged.prefix1, blocks);
        candidate.diplotype.second.alleles = prefixes.path(merged.prefix2, blocks);
        if (candidate.diplotype.second.alleles < candidate.diplotype.first.alleles) {
            std::swap(candidate.diplotype.first, candidate.diplotype.second);
        }
        candidate.best_proposal_score = merged.best;
        candidate.proposal_emission_component = merged.best_emission;
        candidate.proposal_transition_component = merged.best_transition;
        candidate.retained_logsum_proposal_score = merged.log_mass;
        candidate.equivalent_retained_histories = merged.histories;
        candidate.history_count_saturated = merged.saturated;
        candidate.every_allele_scored = emissions.every_allele_scored();
        candidates.push_back(std::move(candidate));
        i = j;
    }
    std::sort(candidates.begin(), candidates.end(), [](const MosaicCandidate& lhs, const MosaicCandidate& rhs) {
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

    result.unique_candidates_found = candidates.size();
    result.search_exhaustive = result.states_pruned == 0;
    result.unique_target_reached = candidates.size() >= parameters.unique_candidate_target;
    result.candidate_set_conditional =
        !result.search_exhaustive || candidates.size() > parameters.unique_candidate_target;
    if (candidates.size() > parameters.unique_candidate_target) {
        candidates.resize(parameters.unique_candidate_target);
    }
    for (std::size_t rank = 0; rank < candidates.size(); ++rank) {
        candidates[rank].candidate_id = rank;
        candidates[rank].unique_rank = rank + 1;
    }
    result.candidates = std::move(candidates);
    return result;
}

} // namespace panvar
