// Branch-local lifetime of tie-branching state in one shared PreparedScoringContext.
//
// A scripted branching search -- a root tie seeding six branches, each running rounds of neighbours and a promoted
// incumbent, one of them splitting into two nested child seeds -- is replayed under two lifetimes:
//   branch_local     at every settle point (the root seeded, each branch processed) every repair slot is released
//                    except the seeds still waiting;
//   retain_parents   the pre-fix lifetime: every parent incumbent stays pinned and nothing is released at a
//                    branch's end.
// Required: every seed, neighbour and incumbent is bit-identical (mass column and exposure) to a fresh full
// placement of its sequence, in both lifetimes; a seed stays usable after its parent is released, and so do nested
// child seeds after theirs; released slots are reused; and under branch_local the censuses held at a settle point
// are exactly the waiting seeds -- they do not grow with the number of branches already completed.
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
        recorded.clear();
    }
};

struct Trace {
    std::vector<std::size_t> census_at_settle, pending_at_settle, live_at_settle, total_at_settle;
    std::vector<std::uint64_t> census_bytes_at_settle;
    std::size_t reused_slots = 0;
};

// One scripted branching search in `context`. branch_local selects the lifetime.
Trace run(const Locus& l, bool branch_local, Checker& check) {
    panvar::PreparedScoringContext context(l.fragments, library(), options());
    Trace trace;
    context.begin_transient_round();
    const Choice root_choice = {};
    const std::size_t root = context.add_haplotypes({spell(l, root_choice)}).front();
    context.retain_census(root);
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
        const std::size_t slot = context.add_repaired_haplotype(root, begin, end, l.alleles[b][b % 3], nullptr, true);
        worklist.push_back({choice, slot, b == 1});   // the second branch will split
    }
    const auto pending_slots = [&](std::size_t next) {
        std::vector<std::size_t> pins;
        for (std::size_t i = next; i < worklist.size(); ++i) pins.push_back(worklist[i].slot);
        return pins;
    };
    const auto settle = [&](std::size_t next) {
        std::vector<std::size_t> keep = pending_slots(next);
        if (!branch_local) keep.insert(keep.end(), legacy_pins.begin(), legacy_pins.end());
        if (branch_local) context.release_transient(keep);
        const panvar::ContextMemoryAccount account = context.memory_account();
        trace.census_at_settle.push_back(account.census_count);
        trace.census_bytes_at_settle.push_back(account.census_bytes);
        trace.pending_at_settle.push_back(worklist.size() - next);
        trace.live_at_settle.push_back(account.slots_live);
        trace.total_at_settle.push_back(account.slots_total);
        for (std::size_t i = next; i < worklist.size(); ++i) {
            require(context.holds_census(worklist[i].slot), "a waiting seed must keep its census");
        }
    };
    settle(0);
    if (branch_local) require(!context.holds_census(root), "the root parent is released once its seeds exist");
    for (std::size_t next = 0; next < worklist.size(); ++next) {
        const Seed seed = worklist[next];
        check.same(context, seed.slot, spell(l, seed.choice), "seed " + std::to_string(next));
        // Two rounds: eight neighbours (no census), then a promoted incumbent (census), then the round's release.
        std::size_t incumbent = seed.slot;
        Choice current = seed.choice;
        for (std::size_t round = 0; round < 2; ++round) {
            for (std::size_t j = 0; j < 8; ++j) {
                const std::size_t block = (next + round + j) % l.blocks.size();
                Choice choice = current;
                choice[block] = (j + round) % 3;
                if (choice == current) continue;
                const auto [begin, end] = span_of(l, current, block);
                const std::size_t slot =
                    context.add_repaired_haplotype(incumbent, begin, end, l.alleles[block][choice[block]], nullptr, false);
                if (!seen_slots.insert(slot).second) ++trace.reused_slots;
                if (j % 3 == 0) check.same(context, slot, spell(l, choice), "neighbour");
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
                    context.add_repaired_haplotype(incumbent, begin, end, l.alleles[block][child[block]], nullptr, true);
                worklist.push_back({child, slot, false});
            }
        }
        settle(next + 1);
        if (branch_local && seed.split) {
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
    return trace;
}

}  // namespace

int main() {
    try {
        const Locus l = locus();
        Checker check(l);
        const Trace local = run(l, true, check);
        check.verify();
        const Trace legacy = run(l, false, check);
        check.verify();
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
                     "to fresh placements, under both lifetimes\n";
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
        std::cout << "branch lifetime gates passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "branch lifetime gate FAILED: " << error.what() << '\n';
        return 1;
    }
}
