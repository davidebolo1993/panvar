#include "panvar/mosaic_beam.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <queue>
#include <set>
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

    // Lookup without interning, for tracing: never changes an id.
    std::optional<std::uint32_t> find(std::uint32_t parent, std::size_t allele) const {
        const std::uint64_t key =
            (static_cast<std::uint64_t>(parent) << 32) | static_cast<std::uint64_t>(allele);
        const auto found = index_.find(key);
        if (found == index_.end()) return std::nullopt;
        return found->second;
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

// Keeps at most state_limit states of a layer after every state has been fully merged, PREFIX FIRST.
//
// A layer's width is useful only to the extent it holds distinct diplotype prefixes. States sharing a
// prefix pair differ solely in which templates currently carry that prefix; they are the same history
// counted many ways, and they rise and fall together when the next block's emission moves. Ranking
// states alone therefore lets one prefix's template variants take the whole width -- measured at 1,018
// of 1,024 retained states on C4 ctx:4 -- so a rival prefix is dropped while the leader is held over a
// thousand times.
//
// The rule, with no parameter and no fixed per-prefix cap: group the layer's states by unordered prefix
// pair, rank groups by aggregate (log-sum-exp) mass, keep the best state of each group in rank order
// until the limit is reached, then spend whatever capacity is left on the best remaining states under
// the existing state order. Where capacity is not scarce this retains exactly what state-only selection
// retained; under pressure every prefix gets one slot before any prefix gets two.
//
// Memory stays O(limit) because groups close early. A child prefix is interned from exactly one parent
// prefix, so every state of a child group is offered inside one parent's extension and the group is
// final the moment that parent group closes -- which is why close_group() exists rather than a map over
// the whole layer, which would have to hold every one of the tens of millions of states offered.
//
// Determinism, on which gate 2 rests: a group's states are sorted by key before its mass is summed, so
// the aggregate is bitwise independent of the order the states were offered in; groups are ranked by
// mass, then by their best state under the existing total order, then -- through that order's own key
// tie-break -- by prefix and template ids. Distinct groups have distinct best-state keys, so the group
// order is total and no tie is left to arrival order.
class LayerSelector {
public:
    explicit LayerSelector(std::size_t limit)
        : limit_(limit), groups_(worse_group_on_top), remainder_(worse_state_on_top) {}

    void offer(const State& state) {
        ++offered_;
        pending_.push_back(state);
    }

    // Every state of the parent prefix-pair group just extended has now been offered, so each child
    // group among them is complete and can be reduced to (aggregate mass, best state).
    void close_group() {
        if (pending_.empty()) return;
        std::sort(pending_.begin(), pending_.end(), [](const State& lhs, const State& rhs) {
            return std::tie(lhs.prefix1, lhs.prefix2, lhs.template1, lhs.template2) <
                   std::tie(rhs.prefix1, rhs.prefix2, rhs.template1, rhs.template2);
        });
        for (std::size_t i = 0; i < pending_.size();) {
            std::size_t j = i + 1;
            while (j < pending_.size() && pending_[j].prefix1 == pending_[i].prefix1 &&
                   pending_[j].prefix2 == pending_[i].prefix2) {
                ++j;
            }
            Group group;
            std::size_t best = i;
            for (std::size_t k = i; k < j; ++k) {
                group.mass = log_add_exp_local(group.mass, pending_[k].log_mass);
                if (k != i && better_state(pending_[k], pending_[best])) best = k;
            }
            group.best = pending_[best];
            ++groups_offered_;
            // A group evicted here can never be wanted again: eviction only happens once more groups
            // than the limit exist, and that is exactly the case where the top groups fill the width
            // and no remaining state is ever drawn.
            if (groups_.size() < limit_) {
                groups_.push(group);
            } else if (!groups_.empty() && better_group(group, groups_.top())) {
                groups_.pop();
                groups_.push(group);
            }
            for (std::size_t k = i; k < j; ++k) {
                if (k == best) continue;
                if (remainder_.size() < limit_) {
                    remainder_.push(pending_[k]);
                } else if (!remainder_.empty() && better_state(pending_[k], remainder_.top())) {
                    remainder_.pop();
                    remainder_.push(pending_[k]);
                }
            }
            i = j;
        }
        pending_.clear();
    }

    std::vector<State> take(std::uint64_t& pruned) {
        close_group();   // the final parent group of the layer, and block 0's single batch
        std::vector<State> out;
        out.reserve(std::min<std::size_t>(limit_, static_cast<std::size_t>(offered_)));
        while (!groups_.empty()) {
            out.push_back(groups_.top().best);
            groups_.pop();
        }
        // Unused capacity goes to the best states left over, which is what makes the rule collapse to
        // the old one when the width is not scarce. The heap yields its worst first, so it is drained
        // and ordered before the best are taken.
        if (out.size() < limit_) {
            std::vector<State> rest;
            rest.reserve(remainder_.size());
            while (!remainder_.empty()) {
                rest.push_back(remainder_.top());
                remainder_.pop();
            }
            std::sort(rest.begin(), rest.end(), better_state);
            for (std::size_t i = 0; i < rest.size() && out.size() < limit_; ++i) out.push_back(rest[i]);
        }
        pruned += offered_ - out.size();
        std::sort(out.begin(), out.end(),
                  [](const State& lhs, const State& rhs) { return key_of(lhs) < key_of(rhs); });
        return out;
    }

    // Distinct prefix pairs OFFERED to this layer, counted exactly and for free: every group is closed
    // here, so the count needs no set over the offered states.
    std::size_t groups_offered() const { return groups_offered_; }

private:
    struct Group {
        double mass = kNegInf;
        State best;
    };
    static bool better_group(const Group& lhs, const Group& rhs) {
        if (lhs.mass != rhs.mass) return lhs.mass > rhs.mass;
        return better_state(lhs.best, rhs.best);
    }
    static bool worse_group_on_top(const Group& lhs, const Group& rhs) { return better_group(lhs, rhs); }
    static bool worse_state_on_top(const State& lhs, const State& rhs) { return better_state(lhs, rhs); }

    std::size_t limit_;
    std::uint64_t offered_ = 0;
    std::size_t groups_offered_ = 0;
    std::vector<State> pending_;   // one parent group's children, never the whole layer
    std::priority_queue<Group, std::vector<Group>, bool (*)(const Group&, const Group&)> groups_;
    std::priority_queue<State, std::vector<State>, bool (*)(const State&, const State&)> remainder_;
};

// Follows a declared diplotype through the layers. It only reads states and the interner.
class LayerTracer {
public:
    LayerTracer(const MosaicBeamParameters& parameters, const PrefixInterner& prefixes)
        : parameters_(parameters), prefixes_(prefixes) {}

    bool on() const { return parameters_.trace_target.has_value(); }

    void begin(std::size_t block) {
        if (!on()) return;
        current_ = MosaicBeamTraceLayer{};
        current_.block = block;
        masses_.clear();
        block_ = block;
        resolved_ = false;
    }

    void offer(const State& state) {
        if (!on()) return;
        ++current_.states_offered;
        masses_.push_back(state.log_mass);
        if (!resolve()) return;
        if (!matches(state)) return;
        ++current_.target_states_offered;
        current_.target_logsum_log_mass = log_add_exp_local(current_.target_logsum_log_mass, state.log_mass);
        if (state.log_mass > current_.target_best_log_mass) {
            current_.target_best_log_mass = state.log_mass;
            current_.target_best_template1 = state.template1;
            current_.target_best_template2 = state.template2;
        }
        if (declared(state)) {
            current_.declared_templates_offered = true;
            current_.declared_templates_log_mass = state.log_mass;
        }
    }

    void finish(const std::vector<State>& retained, std::size_t groups_offered,
                std::vector<MosaicBeamTraceLayer>& out) {
        if (!on()) return;
        current_.prefix_pairs_offered = groups_offered;
        current_.states_retained = retained.size();
        current_.target_prefixes_present = resolve();
        {
            std::vector<std::pair<std::uint32_t, std::uint32_t>> pairs;
            pairs.reserve(retained.size());
            for (const State& state : retained) pairs.emplace_back(state.prefix1, state.prefix2);
            std::sort(pairs.begin(), pairs.end());
            // One pass for both: the number of distinct prefix pairs, and the size of the largest
            // group sharing one. Run lengths over the sorted pairs give each prefix's share of the
            // width, so the crowding is measured rather than inferred from the count.
            std::size_t distinct = 0;
            std::size_t largest = 0;
            for (std::size_t i = 0; i < pairs.size();) {
                std::size_t j = i;
                while (j < pairs.size() && pairs[j] == pairs[i]) ++j;
                ++distinct;
                largest = std::max(largest, j - i);
                i = j;
            }
            current_.prefix_pairs_retained = distinct;
            current_.max_states_per_prefix = largest;
        }
        current_.states_pruned = current_.states_offered - static_cast<std::uint64_t>(retained.size());
        for (const State& state : retained) {
            current_.retained_min_log_mass = current_.retained_min_log_mass == kNegInf
                                                 ? state.log_mass
                                                 : std::min(current_.retained_min_log_mass, state.log_mass);
            if (current_.target_prefixes_present && matches(state)) {
                ++current_.target_states_retained;
                if (declared(state)) current_.declared_templates_retained = true;
            }
        }
        if (current_.target_states_offered > 0) {
            for (const double mass : masses_) {
                if (mass > current_.target_best_log_mass) ++current_.offered_above_target_best;
            }
        }
        previous1_ = current_.target_prefixes_present ? std::optional<std::uint32_t>(target1_) : std::nullopt;
        previous2_ = current_.target_prefixes_present ? std::optional<std::uint32_t>(target2_) : std::nullopt;
        out.push_back(current_);
    }

private:
    // The target's prefix ids at this block exist once their parents' groups have been extended.
    bool resolve() {
        if (resolved_) return true;
        const MosaicDiplotype& target = *parameters_.trace_target;
        std::optional<std::uint32_t> one;
        std::optional<std::uint32_t> two;
        if (block_ == 0) {
            one = prefixes_.find(kRoot, target.first.alleles[0]);
            two = prefixes_.find(kRoot, target.second.alleles[0]);
        } else if (previous1_ && previous2_) {
            one = prefixes_.find(*previous1_, target.first.alleles[block_]);
            two = prefixes_.find(*previous2_, target.second.alleles[block_]);
        }
        if (!one || !two) return false;
        target1_ = *one;
        target2_ = *two;
        resolved_ = true;
        return true;
    }

    bool matches(const State& state) const {
        return (state.prefix1 == target1_ && state.prefix2 == target2_) ||
               (state.prefix1 == target2_ && state.prefix2 == target1_);
    }

    bool declared(const State& state) const {
        if (!parameters_.trace_templates) return false;
        const std::uint32_t d1 = static_cast<std::uint32_t>(parameters_.trace_templates->first);
        const std::uint32_t d2 = static_cast<std::uint32_t>(parameters_.trace_templates->second);
        return (state.prefix1 == target1_ && state.template1 == d1 && state.prefix2 == target2_ &&
                state.template2 == d2) ||
               (state.prefix1 == target2_ && state.template1 == d2 && state.prefix2 == target1_ &&
                state.template2 == d1);
    }

    const MosaicBeamParameters& parameters_;
    const PrefixInterner& prefixes_;
    MosaicBeamTraceLayer current_;
    std::vector<double> masses_;
    std::size_t block_ = 0;
    bool resolved_ = false;
    std::uint32_t target1_ = 0;
    std::uint32_t target2_ = 0;
    std::optional<std::uint32_t> previous1_;
    std::optional<std::uint32_t> previous2_;
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
    if (parameters.trace_target) {
        if (parameters.trace_target->first.alleles.size() != panel.block_count() ||
            parameters.trace_target->second.alleles.size() != panel.block_count()) {
            throw std::invalid_argument("mosaic beam trace target must have one allele per block");
        }
        if (parameters.trace_templates && (parameters.trace_templates->first >= panel.template_count() ||
                                           parameters.trace_templates->second >= panel.template_count())) {
            throw std::invalid_argument("mosaic beam trace templates are outside the panel");
        }
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

} // namespace

// The one planner. Production checks the beam's pre-allocation budget against this; anything that
// predicts the same work -- the preflight -- must call it rather than restate it. A second formula in
// floating point is not the same function near the cap, which is exactly where the answer matters.
std::uint64_t extension_upper_bound(std::size_t templates, std::size_t blocks, std::size_t limit) {
    const std::uint64_t pairs = static_cast<std::uint64_t>(templates) * templates;
    if (pairs == 0) return 0;
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

std::vector<MosaicCandidate> panel_diplotypes(const PanelAlleleMatrix& panel) {
    validate_panel(panel);
    std::vector<MosaicCandidate> out;
    std::set<std::pair<std::vector<std::size_t>, std::vector<std::size_t>>> seen;
    for (std::size_t i = 0; i < panel.template_count(); ++i) {
        for (std::size_t j = i; j < panel.template_count(); ++j) {
            MosaicCandidate candidate;
            candidate.diplotype.first.alleles.assign(panel.template_alleles[i].begin(), panel.template_alleles[i].end());
            candidate.diplotype.second.alleles.assign(panel.template_alleles[j].begin(), panel.template_alleles[j].end());
            if (candidate.diplotype.second.alleles < candidate.diplotype.first.alleles) {
                std::swap(candidate.diplotype.first, candidate.diplotype.second);
            }
            if (!seen.emplace(candidate.diplotype.first.alleles, candidate.diplotype.second.alleles).second) continue;
            candidate.candidate_id = out.size();
            out.push_back(std::move(candidate));
        }
    }
    return out;
}

MosaicBeamResult off_panel_mosaics(const MosaicBeamResult& beam, const std::vector<MosaicCandidate>& panel) {
    std::set<std::pair<std::vector<std::size_t>, std::vector<std::size_t>>> keys;
    for (const MosaicCandidate& candidate : panel) {
        keys.emplace(candidate.diplotype.first.alleles, candidate.diplotype.second.alleles);
    }
    MosaicBeamResult out = beam;
    out.candidates.clear();
    for (const MosaicCandidate& candidate : beam.candidates) {
        if (keys.count({candidate.diplotype.first.alleles, candidate.diplotype.second.alleles}) != 0) continue;
        out.candidates.push_back(candidate);
    }
    for (std::size_t rank = 0; rank < out.candidates.size(); ++rank) {
        out.candidates[rank].candidate_id = panel.size() + rank;
        out.candidates[rank].unique_rank = rank + 1;
    }
    out.unique_candidates_found = out.candidates.size();
    return out;
}

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
    LayerTracer tracer(parameters, prefixes);
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
        tracer.begin(0);
        for (std::size_t i = 0; i < initial.size();) {
            State merged = initial[i];
            std::size_t j = i + 1;
            for (; j < initial.size() && key_of(initial[j]) == key_of(merged); ++j) merge_into(merged, initial[j]);
            selector.offer(merged);
            tracer.offer(merged);
            i = j;
        }
        selector.close_group();
        layer = selector.take(result.states_pruned);
        tracer.finish(layer, selector.groups_offered(), result.trace);
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
        tracer.begin(block);
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
                tracer.offer(child);
                cell_mass[cell] = kNegInf;
                cell_best[cell] = kNegInf;
                cell_histories[cell] = 0;
                cell_saturated[cell] = 0;
            }
            touched.clear();
            // This parent's children are complete, and no other parent can add to any of their prefix
            // groups, so the groups can be reduced now and the states released.
            selector.close_group();
            g = end;
        }
        layer = selector.take(result.states_pruned);
        tracer.finish(layer, selector.groups_offered(), result.trace);
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
