// Branch-local lifetime of tie-branching state in one shared PreparedScoringContext.
//
// A scripted branching search -- a root tie seeding six branches, each running rounds of neighbours and a promoted
// incumbent, one of them splitting into two nested child seeds -- is replayed under three lifetimes:
//   census_free      production: a waiting branch is a descriptor (its allele choice) and holds nothing in the
//                    context; its start is placed in full, censused, just before it runs, and at every settle point
//                    (the root recorded, each branch processed) every repair slot is released;
//   branch_local     the eager frontier: every seed is a census-holding slot built by block delta at its tie, and
//                    at every settle point every repair slot is released except the seeds still waiting;
//   retain_parents   the pre-fix lifetime: eager seeds, every parent incumbent stays pinned and nothing is released
//                    at a branch's end.
// Required: every seed, neighbour and incumbent is bit-identical (mass column and exposure) to a fresh full
// placement of its sequence, in every lifetime, and bit-identical between the census-free and the eager frontier;
// a seed stays usable after its parent is released, and so do nested child seeds after theirs; released slots are
// reused; under branch_local the censuses held at a settle point are exactly the waiting seeds; and under
// census_free no census, sequence, mass column or live slot is held at any settle point however large the frontier,
// and at most the running branch's start holds a seed census.
#include "panvar/fragment_likelihood.hpp"
#include "panvar/genotype_model.hpp"

#include <algorithm>
#include <cstring>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

std::string revcomp(const std::string& bases) {
    std::string out(bases.size(), 'N');
    for (std::size_t i = 0; i < bases.size(); ++i) {
        const char base = bases[bases.size() - 1 - i];
        out[i] = base == 'A' ? 'T' : base == 'C' ? 'G' : base == 'G' ? 'C' : 'A';
    }
    return out;
}

std::string random_bases(std::mt19937_64& rng, std::size_t length) {
    std::string out(length, 'A');
    for (char& base : out) base = "ACGT"[rng() % 4];
    return out;
}

panvar::LibraryModel library() {
    panvar::LibraryModel model;
    model.insert_mean_bp = 260.0;
    model.insert_sd_bp = 20.0;
    model.support_sigmas = 3.0;
    model.substitution_error = 0.02;
    model.outlier_weight = 0.005;
    model.fragment_start_intensity = 0.01;
    model.provenance = "fixture: declared";
    return model;
}

panvar::FragmentScoringOptions options() {
    panvar::FragmentScoringOptions out;
    out.mismatch_band = 6;
    out.minimum_piece_bp = 10;
    out.threads = 1;
    out.piece_hit_publish_batch = 1;
    out.emission_cache_bytes = 1u << 20;
    return out;
}

// The locus: a 6 kb haplotype cut into 10 blocks of 600 bp; each block has three alternative alleles (a
// substitution, an insertion, a deletion).
struct Locus {
    std::string base;
    std::vector<std::pair<std::size_t, std::size_t>> blocks;   // [begin, end) on the base
    std::vector<std::vector<std::string>> alleles;              // per block: alternatives to the base's bytes
    std::vector<panvar::ReadPair> fragments;
};

Locus locus() {
    std::mt19937_64 rng(20260925);
    Locus l;
    l.base = random_bases(rng, 6000);
    for (std::size_t b = 0; b < 10; ++b) {
        const std::size_t begin = b * 600, end = begin + 600;
        l.blocks.push_back({begin, end});
        std::string sub = l.base.substr(begin, 600);
        sub[300] = sub[300] == 'A' ? 'C' : 'A';
        std::string ins = l.base.substr(begin, 300) + random_bases(rng, 9) + l.base.substr(begin + 300, 300);
        std::string del = l.base.substr(begin, 280) + l.base.substr(begin + 292, 308);
        l.alleles.push_back({sub, ins, del});
    }
    // Reads from the base and from single-block variants, so every allele has support somewhere.
    std::vector<std::string> sources = {l.base};
    for (std::size_t b = 0; b < 10; ++b) {
        sources.push_back(l.base.substr(0, l.blocks[b].first) + l.alleles[b][b % 3] + l.base.substr(l.blocks[b].second));
    }
    std::size_t n = 0;
    for (const std::string& source : sources) {
        for (std::size_t start = 0; start + 260 <= source.size(); start += 37) {
            std::string forward = source.substr(start, 100);
            std::string reverse = revcomp(source.substr(start + 160, 100));
            if (rng() % 4 == 0) forward[rng() % 100] = "ACGT"[rng() % 4];
            l.fragments.push_back({{"f" + std::to_string(n) + "/1", forward, std::string(100, 'I')},
                                   {"f" + std::to_string(n) + "/2", reverse, std::string(100, 'I')}});
            ++n;
        }
    }
    return l;
}

// A haplotype as the base with some blocks replaced: block -> alternative index.
using Choice = std::map<std::size_t, std::size_t>;

std::string spell(const Locus& l, const Choice& choice) {
    std::string out;
    std::size_t at = 0;
    for (std::size_t b = 0; b < l.blocks.size(); ++b) {
        const auto found = choice.find(b);
        out += found == choice.end() ? l.base.substr(l.blocks[b].first, 600) : l.alleles[b][found->second];
        at = l.blocks[b].second;
    }
    (void)at;
    return out;
}

// Where block b sits in the spelled haplotype of `choice`.
std::pair<std::size_t, std::size_t> span_of(const Locus& l, const Choice& choice, std::size_t block) {
    std::size_t at = 0;
    for (std::size_t b = 0; b < block; ++b) {
        const auto found = choice.find(b);
        at += found == choice.end() ? 600 : l.alleles[b][found->second].size();
    }
    const auto found = choice.find(block);
    return {at, at + (found == choice.end() ? 600 : l.alleles[block][found->second].size())};
}

// Only one PreparedScoringContext may hold the process's emission cache, so the checker records what each slot held
// (its sequence, mass column and exposure) while a run's context is alive, and compares all of them with a fresh
// reference context -- a full placement of every recorded sequence -- once that context is gone.
struct Checker {
    const Locus& l;
    struct Recorded {
        std::string sequence, label;
        std::vector<double> mass;
        double exposure;
    };
    std::vector<Recorded> recorded;
    std::vector<Recorded> last;   // the previous run's records, kept for a cross-lifetime comparison
    std::size_t compared = 0;
    explicit Checker(const Locus& locus) : l(locus) {}
    void same(const panvar::PreparedScoringContext& context, std::size_t slot, const std::string& sequence,
              const std::string& label) {
        Recorded r{sequence, label, std::vector<double>(l.fragments.size()), context.exposure_of(slot)};
        for (std::size_t f = 0; f < l.fragments.size(); ++f) r.mass[f] = context.mass_of(f, slot);
        recorded.push_back(std::move(r));
    }
    void verify() {
        panvar::PreparedScoringContext reference(l.fragments, library(), options());
        for (const Recorded& r : recorded) {
            const std::size_t fresh = reference.add_haplotypes({r.sequence}).front();
            require(r.exposure == reference.exposure_of(fresh), r.label + ": exposure differs from a fresh placement");
            for (std::size_t f = 0; f < l.fragments.size(); ++f) {
                const double b = reference.mass_of(f, fresh);
                require(std::memcmp(&r.mass[f], &b, sizeof(double)) == 0,
                        r.label + ": mass differs from a fresh placement at fragment " + std::to_string(f));
            }
            ++compared;
        }
        last = std::move(recorded);
        recorded.clear();
    }
};

struct Trace {
    std::vector<std::size_t> census_at_settle, pending_at_settle, live_at_settle, total_at_settle;
    std::vector<std::uint64_t> census_bytes_at_settle, sequence_bytes_at_settle, mass_bytes_at_settle;
    std::size_t reused_slots = 0;
    std::size_t seeds_materialised = 0;
    std::size_t peak_live_seed_censuses = 0;
    panvar::FineWindowReport fine;   // I2-B: the fine run's window discovery
};

enum class Lifetime { CensusFree, BranchLocal, RetainParents };

// One scripted branching search in its own context, under one lifetime.
// `fine` (I2-B): each round's neighbours are scored as one add_repaired_haplotypes batch with round-local fine window
// discovery, as the production search issues them.
// `fused_root` (I2-D): the root is placed and censused in one pass (add_haplotypes with retain_censuses).
Trace run(const Locus& l, Lifetime lifetime, Checker& check, bool fine = false, bool fused_root = false) {
    const bool branch_local = lifetime != Lifetime::RetainParents;
    const bool census_free = lifetime == Lifetime::CensusFree;
    panvar::PreparedScoringContext context(l.fragments, library(), options());
    if (fine) context.enable_fine_window_discovery(std::size_t{1} << 26, std::uint64_t{1} << 40);
    Trace trace;
    constexpr std::size_t kNoSlot = static_cast<std::size_t>(-1);
    context.begin_transient_round();
    const Choice root_choice = {};
    const std::size_t root = context.add_haplotypes({spell(l, root_choice)}, nullptr, false, fused_root).front();
    const std::uint64_t placements_before_retain = context.full_placement_count();
    context.retain_census(root);
    require(!fused_root || context.full_placement_count() == placements_before_retain,
            "I2-D: retain_census after a fused root must place nothing");
    std::set<std::size_t> legacy_pins = {root};
    std::set<std::size_t> seen_slots = {root};
    // The root tie: six seeds, each the root with one block replaced (a promoted repair: census retained).
    struct Seed {
        Choice choice;
        std::size_t slot;
        bool split;
    };
    std::vector<Seed> worklist;
    for (std::size_t b = 0; b < 6; ++b) {
        Choice choice = root_choice;
        choice[b] = b % 3;
        const auto [begin, end] = span_of(l, root_choice, b);
        const std::size_t slot = census_free ? kNoSlot
                                             : context.add_repaired_haplotype(root, begin, end, l.alleles[b][b % 3],
                                                                              nullptr, true);
        worklist.push_back({choice, slot, b == 1});   // the second branch will split
    }
    const auto pending_slots = [&](std::size_t next) {
        std::vector<std::size_t> pins;
        for (std::size_t i = next; i < worklist.size(); ++i) {
            if (worklist[i].slot != kNoSlot) pins.push_back(worklist[i].slot);
        }
        return pins;
    };
    const auto settle = [&](std::size_t next) {
        std::vector<std::size_t> keep = pending_slots(next);
        if (!branch_local) keep.insert(keep.end(), legacy_pins.begin(), legacy_pins.end());
        if (branch_local) context.release_transient(keep);
        const panvar::ContextMemoryAccount account = context.memory_account();
        trace.census_at_settle.push_back(account.census_count);
        trace.census_bytes_at_settle.push_back(account.census_bytes);
        trace.sequence_bytes_at_settle.push_back(account.sequence_bytes);
        trace.mass_bytes_at_settle.push_back(account.mass_column_bytes);
        trace.pending_at_settle.push_back(worklist.size() - next);
        trace.live_at_settle.push_back(account.slots_live);
        trace.total_at_settle.push_back(account.slots_total);
        for (std::size_t i = next; i < worklist.size(); ++i) {
            require(census_free ? worklist[i].slot == kNoSlot : context.holds_census(worklist[i].slot),
                    census_free ? "a waiting census-free seed must hold no slot" : "a waiting seed must keep its census");
        }
    };
    settle(0);
    if (branch_local) require(!context.holds_census(root), "the root parent is released once its seeds exist");
    for (std::size_t next = 0; next < worklist.size(); ++next) {
        Seed seed = worklist[next];
        if (census_free) {
            // Just in time: the start spelled from its descriptor, placed in full and censused.
            require(context.memory_account().census_count == 0, "no census may be held before a branch materialises");
            seed.slot = context.add_haplotypes({spell(l, seed.choice)}, nullptr, true).front();
            context.retain_census(seed.slot);
            if (!seen_slots.insert(seed.slot).second) ++trace.reused_slots;
            ++trace.seeds_materialised;
        }
        trace.peak_live_seed_censuses =
            std::max(trace.peak_live_seed_censuses, pending_slots(next + 1).size() + 1);
        check.same(context, seed.slot, spell(l, seed.choice), "seed " + std::to_string(next));
        // Two rounds: eight neighbours (no census), then a promoted incumbent (census), then the round's release.
        std::size_t incumbent = seed.slot;
        Choice current = seed.choice;
        for (std::size_t round = 0; round < 2; ++round) {
            std::vector<panvar::RepairRequest> requests;
            std::vector<std::pair<Choice, bool>> requested;
            for (std::size_t j = 0; j < 8; ++j) {
                const std::size_t block = (next + round + j) % l.blocks.size();
                Choice choice = current;
                choice[block] = (j + round) % 3;
                if (choice == current) continue;
                const auto [begin, end] = span_of(l, current, block);
                if (fine) {
                    requests.push_back({incumbent, begin, end, l.alleles[block][choice[block]]});
                    requested.push_back({choice, j % 3 == 0});
                    continue;
                }
                const std::size_t slot =
                    context.add_repaired_haplotype(incumbent, begin, end, l.alleles[block][choice[block]], nullptr, false);
                if (!seen_slots.insert(slot).second) ++trace.reused_slots;
                if (j % 3 == 0) check.same(context, slot, spell(l, choice), "neighbour");
            }
            if (fine) {
                const std::vector<std::size_t> slots = context.add_repaired_haplotypes(requests);
                for (std::size_t k = 0; k < slots.size(); ++k) {
                    if (!seen_slots.insert(slots[k]).second) ++trace.reused_slots;
                    if (requested[k].second) check.same(context, slots[k], spell(l, requested[k].first), "neighbour");
                }
            }
            const std::size_t block = (next + 5 + round) % l.blocks.size();
            Choice won = current;
            won[block] = (round + 1) % 3;
            const auto [begin, end] = span_of(l, current, block);
            const std::size_t winner = context.add_repaired_haplotype(incumbent, begin, end, l.alleles[block][won[block]],
                                                                      nullptr, true);
            if (!seen_slots.insert(winner).second) ++trace.reused_slots;
            check.same(context, winner, spell(l, won), "promoted incumbent");
            // The round's own release: the new incumbent and every pin survive (as CoordinateRepairSearch does).
            std::vector<std::size_t> keep = {winner};
            const std::vector<std::size_t> pins = pending_slots(next + 1);
            keep.insert(keep.end(), pins.begin(), pins.end());
            keep.push_back(seed.slot);   // a running branch's own seed stays pinned until it has run
            if (!branch_local) keep.insert(keep.end(), legacy_pins.begin(), legacy_pins.end());
            context.release_transient(keep);
            incumbent = winner;
            current = won;
        }
        if (seed.split) {
            // A split: its terminal incumbent is the parent of two nested child seeds.
            legacy_pins.insert(incumbent);
            for (std::size_t c = 0; c < 2; ++c) {
                const std::size_t block = (next + 3 + c) % l.blocks.size();
                Choice child = current;
                child[block] = (c + 2) % 3;
                const auto [begin, end] = span_of(l, current, block);
                const std::size_t slot =
                    census_free ? kNoSlot
                                : context.add_repaired_haplotype(incumbent, begin, end, l.alleles[block][child[block]],
                                                                 nullptr, true);
                worklist.push_back({child, slot, false});
            }
        }
        settle(next + 1);
        if (branch_local && !census_free && seed.split) {
            require(!context.holds_census(incumbent), "a split parent is released once its children are seeded");
            for (std::size_t i = worklist.size() - 2; i < worklist.size(); ++i) {
                const std::size_t block = (i + 7) % l.blocks.size();
                Choice from_child = worklist[i].choice;
                from_child[block] = 1;
                const auto [begin, end] = span_of(l, worklist[i].choice, block);
                const std::size_t probe =
                    context.add_repaired_haplotype(worklist[i].slot, begin, end, l.alleles[block][1], nullptr, false);
                check.same(context, probe, spell(l, from_child), "a delta from a nested seed after its parent's release");
            }
        }
    }
    trace.fine = context.fine_window_report();
    return trace;
}

}  // namespace


// I4-C: the bounded seed cache's admission and release, on the context primitives the production driver uses. Seeds are
// built by exact block delta (census retained) from a census-holding root; a budget admits the first ones whose exact
// payload fits; every rejected seed gives back EXACTLY its own new slots at once (never a pre-existing one); a later
// just-in-time placement then reuses a recycled slot. Required: admitted and recycled columns bit-identical to fresh
// placements; release_slots touches nothing else; a seed whose sequence the context already holds (an identical
// sequence under another label) resolves to that slot and is never released by a rejection; reversed homologue order
// scores bit-identically; release_slots refuses a non-live slot. At 1 and 4 threads.
void test_bounded_seed_release(const Locus& l, Checker& check, std::size_t threads) {
    panvar::FragmentScoringOptions o = options();
    o.threads = threads;
    panvar::PreparedScoringContext context(l.fragments, library(), o);
    context.begin_transient_round();
    const std::size_t root = context.add_haplotypes({spell(l, {})}, nullptr, false).front();
    context.retain_census(root);
    std::vector<std::size_t> admitted;
    std::vector<Choice> rejected;
    std::uint64_t cache_bytes = 0, budget = 0;
    std::size_t released_total = 0;
    for (std::size_t b = 0; b < 6; ++b) {
        Choice choice;
        choice[b] = b % 3;
        const auto [begin, end] = span_of(l, {}, b);
        const std::vector<std::size_t> before = context.live_transient_slots();
        const std::uint64_t held_before = context.transient_payload_bytes();
        const std::size_t slot = context.add_repaired_haplotype(root, begin, end, l.alleles[b][b % 3], nullptr, true);
        const std::uint64_t payload = context.slot_payload_bytes(slot);
        if (b == 0) budget = 2 * payload + payload / 2;   // admits two seeds of this size, rejects the rest
        if (cache_bytes + payload <= budget) {
            cache_bytes += payload;
            admitted.push_back(slot);
            check.same(context, slot, spell(l, choice), "I4-C admitted seed " + std::to_string(b));
            continue;
        }
        // rejected: exactly the build's new slots go back, and the context returns to what it held before
        std::vector<std::size_t> fresh;
        for (const std::size_t x : context.live_transient_slots())
            if (!std::binary_search(before.begin(), before.end(), x)) fresh.push_back(x);
        require(fresh == std::vector<std::size_t>{slot}, "I4-C: a rejected seed's build creates exactly its own slot");
        context.release_slots(fresh);
        require(context.slot_payload_bytes(slot) == 0 && !context.holds_census(slot) &&
                    context.live_transient_slots() == before && context.transient_payload_bytes() == held_before,
                "I4-C: releasing a rejected seed gives back exactly its slot and nothing else");
        ++released_total;
        rejected.push_back(choice);
    }
    require(admitted.size() == 2 && rejected.size() == 4, "I4-C fixture: two seeds admitted, four rejected");
    for (const std::size_t slot : admitted) require(context.holds_census(slot), "I4-C: admitted seeds keep their census");
    // an identical sequence under another label: the same replacement again resolves to the admitted seed's slot, so a
    // "rejection" of it finds no new slot to release
    {
        const auto [begin, end] = span_of(l, {}, 0);
        const std::vector<std::size_t> before = context.live_transient_slots();
        const std::size_t again = context.add_repaired_haplotype(root, begin, end, l.alleles[0][0], nullptr, true);
        require(again == admitted[0] && context.live_transient_slots() == before,
                "I4-C: an identical sequence resolves to the slot already held and creates nothing to release");
    }
    // just in time: each rejected start placed in full reuses a recycled slot and is bit-identical to a fresh placement
    std::size_t recycled = 0;
    for (const Choice& choice : rejected) {
        const std::size_t total_before = context.memory_account().slots_total;
        const std::size_t slot = context.add_haplotypes({spell(l, choice)}, nullptr, true).front();
        context.retain_census(slot);
        if (context.memory_account().slots_total == total_before) ++recycled;
        check.same(context, slot, spell(l, choice), "I4-C recycled just-in-time seed");
    }
    require(recycled > 0, "I4-C: the released slots must be recycled");
    // reversed homologue order
    const panvar::FragmentScoringResult ab = context.score_candidates({{admitted[0], admitted[1]}});
    const panvar::FragmentScoringResult ba = context.score_candidates({{admitted[1], admitted[0]}});
    const auto bits = [](double x, double y) { return std::memcmp(&x, &y, sizeof(double)) == 0; };
    require(bits(ab.scores[0].log_likelihood, ba.scores[0].log_likelihood) && bits(ab.scores[0].fragment_term, ba.scores[0].fragment_term) &&
                bits(ab.scores[0].exposure_term, ba.scores[0].exposure_term) && ab.scores[0].fragments_with_mass == ba.scores[0].fragments_with_mass,
            "I4-C: reversed homologue order must score bit-identically");
    // refusals
    bool refused = false;
    try {
        context.release_slots({context.memory_account().slots_total + 5});
    } catch (const std::invalid_argument&) {
        refused = true;
    }
    require(refused, "I4-C: release_slots must refuse a slot that is not a live transient haplotype");
    std::cout << "  ok: I4-C bounded seed admission at " << threads << " thread(s): 2 admitted, " << released_total
              << " rejected seeds released at once (exactly their own slots), " << recycled << " recycled just-in-time "
                 "starts; identical sequence shares its slot; reversed homologue order bit-identical\n";
}

int main() {
    try {
        const Locus l = locus();
        Checker check(l);
        const Trace local = run(l, Lifetime::BranchLocal, check);
        check.verify();
        std::vector<Checker::Recorded> eager_records;
        for (const Checker::Recorded& r : check.last) {
            if (r.label.rfind("a delta from a nested seed", 0) != 0) eager_records.push_back(r);   // eager-only probes
        }
        const Trace legacy = run(l, Lifetime::RetainParents, check);
        check.verify();
        const Trace jit = run(l, Lifetime::CensusFree, check);
        check.verify();
        // The census-free frontier against the eager one: the same seeds, neighbours and incumbents, bit for bit.
        require(check.last.size() == eager_records.size(), "census-free and eager runs recorded different slot sets");
        for (std::size_t i = 0; i < eager_records.size(); ++i) {
            const Checker::Recorded& a = eager_records[i];
            const Checker::Recorded& b = check.last[i];
            require(a.label == b.label && a.sequence == b.sequence, "census-free and eager runs diverge at " + a.label);
            require(std::memcmp(&a.exposure, &b.exposure, sizeof(double)) == 0, a.label + ": exposure differs eager vs census-free");
            require(std::memcmp(a.mass.data(), b.mass.data(), a.mass.size() * sizeof(double)) == 0,
                    a.label + ": mass column differs eager vs census-free");
        }
        require(jit.pending_at_settle == local.pending_at_settle, "the frontier differs between census-free and eager");
        // I2-B: the census-free branching search again, every round a batch with fine window discovery, each window's
        // rows checked against the reference scans; then every record against fresh placements and against the
        // census-free run's own, bit for bit, with the same frontier, slot reuse and seed materialisation.
        {
            const std::vector<Checker::Recorded> jit_records = check.last;
            panvar::set_fine_window_row_check(true);
            const Trace fine = run(l, Lifetime::CensusFree, check, true, true);   // I2-D: with a fused root too
            panvar::set_fine_window_row_check(false);
            check.verify();
            require(check.last.size() == jit_records.size(), "I2-B: the fine branching run recorded a different slot set");
            for (std::size_t i = 0; i < jit_records.size(); ++i) {
                const Checker::Recorded& a = jit_records[i];
                const Checker::Recorded& b = check.last[i];
                require(a.label == b.label && a.sequence == b.sequence, "I2-B: the fine branching run diverges at " + a.label);
                require(std::memcmp(&a.exposure, &b.exposure, sizeof(double)) == 0 &&
                            std::memcmp(a.mass.data(), b.mass.data(), a.mass.size() * sizeof(double)) == 0,
                        "I2-B: " + a.label + " differs between the fine and the reference branching run");
            }
            require(fine.pending_at_settle == jit.pending_at_settle && fine.seeds_materialised == jit.seeds_materialised &&
                        fine.reused_slots > 0 && fine.census_at_settle == jit.census_at_settle,
                    "I2-B: the frontier, materialisation, slot reuse and censuses must be the census-free run's");
            require(fine.fine.rounds > 0 && fine.fine.units_served > 0 && fine.fine.fallback_windows == 0,
                    "I2-B: the branching run must serve shared fine units");
            std::cout << "  ok: I2-B/I2-D tie-branching replay with fine window discovery and a fused root seed -- " << fine.fine.rounds << " rounds, "
                      << fine.fine.worker_windows << " windows (row-checked), " << fine.fine.units_served
                      << " units served; every record equals the census-free reference run and fresh placements\n";
        }
        require(jit.seeds_materialised == local.pending_at_settle.size() - 1, "every branch materialises exactly one seed");
        std::size_t peak_frontier = 0;
        for (std::size_t i = 0; i < jit.census_at_settle.size(); ++i) {
            require(jit.census_at_settle[i] == 0 && jit.census_bytes_at_settle[i] == 0,
                    "census_free: a census is held at settle point " + std::to_string(i));
            require(jit.live_at_settle[i] == 0, "census_free: a slot is live at settle point " + std::to_string(i));
            require(jit.sequence_bytes_at_settle[i] == jit.sequence_bytes_at_settle[0] &&
                        jit.mass_bytes_at_settle[i] == jit.mass_bytes_at_settle[0],
                    "census_free: sequence or mass-column bytes change with the frontier at settle point " + std::to_string(i));
            peak_frontier = std::max(peak_frontier, jit.pending_at_settle[i]);
        }
        require(peak_frontier > 1, "the fixture's frontier must grow for the memory gate to mean anything");
        require(jit.peak_live_seed_censuses == 1, "census_free: only the running branch's start may hold a seed census");
        require(jit.reused_slots > 0, "census_free: released slots must be reused");
        require(jit.total_at_settle.back() <= local.total_at_settle.back(),
                "census_free: the slot table outgrows the eager frontier's");
        for (std::size_t i = 0; i < local.census_at_settle.size(); ++i) {
            require(local.census_at_settle[i] == local.pending_at_settle[i],
                    "branch_local: the censuses held at settle point " + std::to_string(i) + " are not exactly the waiting seeds (" +
                        std::to_string(local.census_at_settle[i]) + " vs " + std::to_string(local.pending_at_settle[i]) + ")");
            require(local.live_at_settle[i] == local.pending_at_settle[i],
                    "branch_local: the live slots at a settle point are not exactly the waiting seeds");
        }
        // The excess over the frontier: censuses held beyond the waiting seeds. Zero throughout under branch_local;
        // under retain_parents it grows with the work already done (pinned parents, the last finished incumbent).
        const auto excess = [](const Trace& t, std::size_t i) { return t.census_at_settle[i] - t.pending_at_settle[i]; };
        require(excess(legacy, legacy.census_at_settle.size() - 1) > excess(legacy, 0) &&
                    excess(legacy, legacy.census_at_settle.size() - 1) > 0,
                "the retain-parents lifetime must show the accumulation the branch-local lifetime removes");
        require(local.reused_slots > 0, "released slots must be reused");
        require(local.total_at_settle.back() <= legacy.total_at_settle.back(), "the branch-local slot table does not outgrow the old one");
        std::cout << "  ok: " << check.compared << " seeds, neighbours, promoted incumbents and nested-seed deltas bit-identical "
                     "to fresh placements, under all three lifetimes\n";
        std::cout << "  ok: " << eager_records.size() << " seeds, neighbours and incumbents bit-identical between the "
                     "census-free and the eager frontier\n";
        std::cout << "  ok: census_free at each settle point -- frontier:";
        for (const std::size_t n : jit.pending_at_settle) std::cout << ' ' << n;
        std::cout << "; censuses:";
        for (const std::size_t n : jit.census_at_settle) std::cout << ' ' << n;
        std::cout << "; live slots:";
        for (const std::size_t n : jit.live_at_settle) std::cout << ' ' << n;
        std::cout << "; sequence bytes " << jit.sequence_bytes_at_settle.back() << ", mass-column bytes "
                  << jit.mass_bytes_at_settle.back() << " throughout\n";
        std::cout << "  ok: seeds materialised " << jit.seeds_materialised << "; peak live seed censuses census_free "
                  << jit.peak_live_seed_censuses << ", eager " << local.peak_live_seed_censuses << "; census_free slot table "
                  << jit.total_at_settle.back() << ", slot reuses " << jit.reused_slots << "\n";
        std::cout << "  ok: censuses beyond the waiting seeds at each settle point -- branch_local:";
        for (std::size_t i = 0; i < local.census_at_settle.size(); ++i) std::cout << ' ' << excess(local, i);
        std::cout << "; retain_parents:";
        for (std::size_t i = 0; i < legacy.census_at_settle.size(); ++i) std::cout << ' ' << excess(legacy, i);
        std::cout << "\n  ok: branch_local censuses at each settle point = waiting seeds:";
        for (std::size_t i = 0; i < local.census_at_settle.size(); ++i) std::cout << ' ' << local.census_at_settle[i];
        std::cout << "; retain_parents:";
        for (std::size_t i = 0; i < legacy.census_at_settle.size(); ++i) std::cout << ' ' << legacy.census_at_settle[i];
        std::cout << "\n  ok: released slots reused " << local.reused_slots << " times; slot table " << local.total_at_settle.back()
                  << " (retain_parents " << legacy.total_at_settle.back() << ")\n";
        for (const std::size_t threads : {std::size_t{1}, std::size_t{4}}) {
            test_bounded_seed_release(l, check, threads);
            check.verify();
        }
        std::cout << "branch lifetime gates passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "branch lifetime gate FAILED: " << error.what() << '\n';
        return 1;
    }
}
