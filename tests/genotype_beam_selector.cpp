// Gate 2 for prefix-first selection: what a layer RETAINS must be a function of the states offered to
// it, never of the order they were offered in.
//
// This cannot be tested through propose_mosaic_diplotypes. The only input that reorders the offered
// states is the panel's template order, and relabelling templates does not merely reorder -- template
// ids are the final tie-break of the beam's state order, so a relabelled panel legitimately keeps a
// different one of two equal-mass states. That sensitivity predates prefix-first selection and is a
// property of the tie-break, not of the selection rule.
//
// So this test drives the selector directly, which is why it compiles mosaic_beam.cpp into its own
// translation unit: LayerSelector and State live in an unnamed namespace and have no other reader. The
// states offered are held fixed and only their order varies.
#include "../src/mosaic_beam.cpp"

#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

panvar::State make_state(std::uint32_t prefix1, std::uint32_t template1, std::uint32_t prefix2,
                         std::uint32_t template2, double log_mass) {
    panvar::State state;
    state.prefix1 = prefix1;
    state.template1 = template1;
    state.prefix2 = prefix2;
    state.template2 = template2;
    state.log_mass = log_mass;
    state.best = log_mass;
    state.best_emission = log_mass;
    state.histories = 1;
    panvar::canonicalise(state);
    return state;
}

// One batch is one parent prefix-pair group: the unit the beam closes groups on. Every state of a child
// prefix group is offered inside a single batch, which is what lets the selector close groups early.
using Batch = std::vector<panvar::State>;

std::vector<panvar::State> select(const std::vector<Batch>& batches, std::size_t limit,
                                  std::uint64_t& pruned) {
    panvar::LayerSelector selector(limit);
    for (const Batch& batch : batches) {
        for (const panvar::State& state : batch) selector.offer(state);
        selector.close_group();
    }
    return selector.take(pruned);
}

bool identical(const std::vector<panvar::State>& lhs, const std::vector<panvar::State>& rhs) {
    if (lhs.size() != rhs.size()) return false;
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        // Bitwise on the mass: an order-dependent log-sum-exp would show up here and nowhere else.
        if (panvar::key_of(lhs[i]) != panvar::key_of(rhs[i]) || lhs[i].log_mass != rhs[i].log_mass ||
            lhs[i].best != rhs[i].best || lhs[i].histories != rhs[i].histories) {
            return false;
        }
    }
    return true;
}

// Four prefix groups of uneven size and interleaved mass, so that ranking groups and ranking states
// disagree: group 10/11 holds five states including the single best, while groups 20/21 and 30/31 hold
// one apiece whose mass sits below several members of the crowded group.
std::vector<Batch> fixture() {
    Batch crowded;
    crowded.push_back(make_state(10, 0, 11, 1, -1.0));
    crowded.push_back(make_state(10, 2, 11, 3, -1.5));
    crowded.push_back(make_state(10, 4, 11, 5, -1.7));
    crowded.push_back(make_state(10, 6, 11, 7, -1.9));
    crowded.push_back(make_state(10, 8, 11, 9, -2.1));
    Batch others;
    others.push_back(make_state(20, 0, 21, 1, -2.6));
    others.push_back(make_state(20, 2, 21, 3, -3.4));
    others.push_back(make_state(30, 1, 31, 2, -2.9));
    Batch singleton;
    singleton.push_back(make_state(40, 0, 41, 0, -4.2));
    return {crowded, others, singleton};
}

void test_offer_order_does_not_change_retention() {
    const std::vector<Batch> batches = fixture();
    std::size_t offered = 0;
    for (const Batch& batch : batches) offered += batch.size();

    for (const std::size_t limit : {std::size_t(1), std::size_t(3), std::size_t(4), std::size_t(6),
                                    std::size_t(9), std::size_t(64)}) {
        std::uint64_t pruned = 0;
        const std::vector<panvar::State> reference = select(batches, limit, pruned);
        require(reference.size() <= limit, "a layer must never retain more than the limit");
        require(pruned == offered - reference.size(), "every state not retained must be counted pruned");

        std::mt19937 generator(20260921u + static_cast<unsigned>(limit));
        for (int trial = 0; trial < 40; ++trial) {
            // Shuffle the states inside each batch, and the batches among themselves. Both are order
            // changes the selector must be blind to; neither moves a state into another group.
            std::vector<Batch> shuffled = batches;
            for (Batch& batch : shuffled) std::shuffle(batch.begin(), batch.end(), generator);
            std::shuffle(shuffled.begin(), shuffled.end(), generator);
            std::uint64_t again = 0;
            const std::vector<panvar::State> retained = select(shuffled, limit, again);
            require(again == pruned && identical(retained, reference),
                    "reordering the offered states must not change what the layer retains, at limit " +
                        std::to_string(limit));
        }
    }
}

// Gate 1 at the selector: with capacity to spare nothing is pruned and every offered state survives,
// so a layer that never hits the limit is exactly the exhaustive layer.
void test_spare_capacity_retains_everything() {
    const std::vector<Batch> batches = fixture();
    std::size_t offered = 0;
    for (const Batch& batch : batches) offered += batch.size();
    std::uint64_t pruned = 0;
    const std::vector<panvar::State> retained = select(batches, offered, pruned);
    require(retained.size() == offered && pruned == 0,
            "a width at least as large as the offering must retain all of it and prune nothing");
}

// The rule itself: under pressure every prefix gets a slot before any prefix gets a second one, and
// spare capacity then goes to the best remaining states under the existing order.
void test_prefix_first_then_best_remaining() {
    const std::vector<Batch> batches = fixture();
    const auto prefixes_in = [](const std::vector<panvar::State>& states) {
        std::set<std::pair<std::uint32_t, std::uint32_t>> seen;
        for (const panvar::State& state : states) seen.emplace(state.prefix1, state.prefix2);
        return seen;
    };

    std::uint64_t pruned = 0;
    // Four groups, width four: one slot each, and the crowded group cannot take a second.
    const std::vector<panvar::State> tight = select(batches, 4, pruned);
    require(tight.size() == 4 && prefixes_in(tight).size() == 4,
            "with the width equal to the number of prefixes, every prefix must be represented");

    // Width six: four group bests, then the two best remaining states, both from the crowded group.
    pruned = 0;
    const std::vector<panvar::State> loose = select(batches, 6, pruned);
    require(loose.size() == 6 && prefixes_in(loose).size() == 4,
            "spare capacity must go to the best remaining states, not to new prefixes that do not exist");
    std::size_t crowded_states = 0;
    for (const panvar::State& state : loose) {
        if (state.prefix1 == 10 && state.prefix2 == 11) ++crowded_states;
    }
    require(crowded_states == 3,
            "the crowded prefix must take exactly the spare slots: one as a group best, two as fill");

    // Width one: the single best group, which is the group holding the best state.
    pruned = 0;
    const std::vector<panvar::State> one = select(batches, 1, pruned);
    require(one.size() == 1 && one.front().prefix1 == 10 && one.front().prefix2 == 11,
            "a width of one must go to the highest-ranked prefix group");
}

}   // namespace

int main() {
    try {
        test_offer_order_does_not_change_retention();
        test_spare_capacity_retains_everything();
        test_prefix_first_then_best_remaining();
        std::cout << "genotype beam selector: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "genotype beam selector: FAIL: " << error.what() << '\n';
        return 1;
    }
}
