// Exactness of the union seeding path's optimisations against the paths they replaced:
//   1  the early-refusing feasibility prefilter decides exactly as the full-scan kernel, prepared and unprepared,
//      over adversarial edit scripts (substitutions, insertions, deletions, mixed, on piece boundaries, at and
//      around the true start), repeats, unrelated sequence, windows cut short by the haplotype end, reads longer
//      than the bit-parallel word, bands 0..8 and gaps 0..4;
//   2  the counting-by-form start set is exactly std::sort + std::unique's list on randomised proposals with heavy
//      duplication, empty and single-form lists, forms above the insertion-sort threshold, and sparse forms;
//   3  the ACGT refusal still fires, with the same message, for a non-ACGT base in the window or the read;
//   4  the emission key resumed from a prepared form's carried hash state equals the full hash in every field,
//      and scoring with it gives identical mass cells and identical emission-cache hits, misses, executed DPs,
//      entries, retention refusals and budget bytes as hashing the whole read -- under a full, an undersized
//      and a zero-room budget, and with every key forced into one hash.
#include "panvar/fragment_likelihood.hpp"
#include "panvar/genotype_model.hpp"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

panvar::SequencingRead read_of(const std::string& name, const std::string& bases) {
    return {name, bases, std::string(bases.size(), 'I')};
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

// A read copied from hap[start, start + length) with the given edits applied: substitutions, insertions of a
// random base, deletions.
std::string edited_read(const std::string& hap, std::size_t start, std::size_t length, std::size_t edits, int mode,
                        std::size_t pieces, std::mt19937_64& rng) {
    std::string read = hap.substr(start, length + 8);
    for (std::size_t e = 0; e < edits; ++e) {
        const char kind = mode == 1 ? 's' : mode == 2 ? 'i' : mode == 3 ? 'd' : "sid"[rng() % 3];
        std::size_t at = rng() % length;
        if (rng() % 2 == 0) at = std::min<std::size_t>(length - 1, (rng() % pieces) * length / pieces);   // boundary
        if (at >= read.size()) at = read.size() - 1;
        if (kind == 's') {
            const char was = read[at];
            do read[at] = "ACGT"[rng() % 4]; while (read[at] == was);
        } else if (kind == 'i') {
            read.insert(read.begin() + static_cast<std::ptrdiff_t>(at), "ACGT"[rng() % 4]);
        } else if (read.size() > 1) {
            read.erase(read.begin() + static_cast<std::ptrdiff_t>(at));
        }
    }
    read.resize(std::min(read.size(), length));
    return read;
}

}  // namespace

int main() {
    try {
        std::mt19937_64 rng(20260925);
        // ---- 1. the prefilter --------------------------------------------------------------------------
        std::uint64_t checked = 0, admitted = 0, refused = 0;
        const auto check = [&](const std::string& hap, std::size_t start, const std::string& read, std::size_t band,
                               std::size_t gap, const std::string& label) {
            if (read.empty() || start > hap.size()) return;
            const panvar::PrefilterDecisions d = panvar::feasibility_prefilter_decisions(hap, start, read, band, gap);
            require(d.production_prepared == d.reference_prepared && d.production_unprepared == d.reference_unprepared &&
                        d.reference_prepared == d.reference_unprepared,
                    label + ": prefilter decisions differ (start " + std::to_string(start) + ", band " +
                        std::to_string(band) + ", gap " + std::to_string(gap) + ")");
            ++checked;
            (d.reference_prepared ? admitted : refused) += 1;
        };
        const std::string unique = random_bases(rng, 1500);
        std::string repeat;
        {
            const std::string unit = random_bases(rng, 7);
            while (repeat.size() < 1500) repeat += unit;
            for (int v = 0; v < 12; ++v) repeat[rng() % repeat.size()] = "ACGT"[rng() % 4];
        }
        std::string period41;
        {
            const std::string unit = random_bases(rng, 41);
            while (period41.size() < 1500) period41 += unit;
        }
        const std::string* haps[3] = {&unique, &repeat, &period41};
        for (int trial = 0; trial < 6000; ++trial) {
            const std::string& hap = *haps[trial % 3];
            const std::size_t band = rng() % 9;
            const std::size_t gap = rng() % 5;
            const std::size_t length = trial % 11 == 0 ? 60 + rng() % 200 : 100;   // some above the 128-bit word
            const std::size_t start = 100 + rng() % (hap.size() - length - 120);
            const std::size_t edits = rng() % (band + 5);
            const int mode = static_cast<int>(rng() % 4);
            const std::string read = edited_read(hap, start, length, edits, mode, band + 1, rng);
            for (int shift = -4; shift <= 4; ++shift) {
                check(hap, static_cast<std::size_t>(static_cast<long>(start) + shift), read, band, gap,
                      "edit trial " + std::to_string(trial));
            }
            // Windows cut short by the haplotype end: the read placed flush against it and past it.
            const std::size_t tail = hap.size() - read.size();
            for (std::size_t back = 0; back <= gap + 2 && back <= tail; ++back) {
                check(hap, tail - back + (back % 2), read, band, gap, "end trial " + std::to_string(trial));
            }
        }
        for (int trial = 0; trial < 2000; ++trial) {   // unrelated sequence: refusals, decided early or late
            const std::string read = random_bases(rng, trial % 7 == 0 ? 150 : 100);
            check(unique, rng() % 1200, read, rng() % 9, rng() % 5, "unrelated " + std::to_string(trial));
        }
        require(admitted > 1000 && refused > 1000, "the prefilter fixture must exercise both decisions");
        std::cout << "  ok: prefilter -- " << checked << " decisions (" << admitted << " admitted, " << refused
                  << " refused), early-refusing = full-scan kernel, prepared and unprepared\n";

        // ---- 2. the start set -------------------------------------------------------------------------
        std::uint64_t lists = 0;
        for (int trial = 0; trial < 3000; ++trial) {
            const std::size_t forms = 1 + rng() % (trial % 5 == 0 ? 3 : 5000);
            const std::size_t count = trial % 97 == 0 ? 0 : rng() % 20000;
            std::vector<std::pair<std::uint32_t, std::uint32_t>> proposals;
            for (std::size_t i = 0; i < count; ++i) {
                const std::uint32_t form = static_cast<std::uint32_t>(rng() % forms);
                const std::uint32_t start = static_cast<std::uint32_t>(rng() % (trial % 3 == 0 ? 50 : 300000));
                const std::size_t copies = 1 + rng() % 4;   // duplicate proposals, as several seeds give
                for (std::size_t c = 0; c < copies; ++c) proposals.push_back({form, start});
                for (std::uint32_t s = start > 3 ? start - 3 : 0; s <= start + 3 && rng() % 2 == 0; ++s) {
                    proposals.push_back({form, s});   // a +/- max_gap neighbourhood
                }
            }
            std::shuffle(proposals.begin(), proposals.end(), rng);
            const auto reference = panvar::sort_unique_proposals_for_test(proposals, forms, false);
            const auto bucketed = panvar::sort_unique_proposals_for_test(proposals, forms, true);
            require(reference == bucketed, "start set " + std::to_string(trial) + ": the bucketed list differs");
            ++lists;
        }
        std::cout << "  ok: start set -- " << lists << " randomised proposal lists, bucketed = std::sort + std::unique\n";

        // ---- 3. the ACGT refusal ----------------------------------------------------------------------
        panvar::LibraryModel library;
        library.insert_mean_bp = 260.0;
        library.insert_sd_bp = 20.0;
        library.substitution_error = 0.01;
        library.fragment_start_intensity = 0.01;
        library.provenance = "fixture: declared";
        const std::string hap = unique.substr(0, 400);
        std::string hap_n = hap;
        hap_n[250] = 'N';
        const std::string read = hap.substr(200, 100);
        std::string read_n = read;
        read_n[10] = 'N';
        const auto refuses = [&](const std::string& h, const std::string& r) {
            try {
                (void)panvar::gapped_emission(h, r, 200, library, 6);
            } catch (const std::invalid_argument& error) {
                return std::string(error.what()).find("ACGT") != std::string::npos;
            }
            return false;
        };
        require(refuses(hap_n, read) && refuses(hap, read_n), "a non-ACGT base must still refuse");
        require(!refuses(hap, read), "an all-ACGT window must not refuse");
        std::cout << "  ok: a non-ACGT base in the window or the read still refuses\n";
        // ---- 4. the emission key ----------------------------------------------------------------------
        std::uint64_t keys = 0;
        for (int trial = 0; trial < 20000; ++trial) {
            const std::string read = random_bases(rng, 1 + rng() % 180);
            const std::string window = random_bases(rng, rng() % 200);
            require(panvar::emission_key_prefix_matches(read, window), "a resumed emission key differs from the full hash");
            ++keys;
        }
        panvar::set_emission_cache_hash_collapse(true);
        require(panvar::emission_key_prefix_matches(unique.substr(0, 100), unique.substr(50, 103)),
                "under forced collisions the resumed key must equal the full one too");
        panvar::set_emission_cache_hash_collapse(false);
        {
            // A small locus: three haplotypes sharing stretches, fragments drawn from them with errors.
            std::vector<std::string> haplotypes = {unique.substr(0, 900), unique.substr(300, 600) + repeat.substr(0, 300),
                                                   period41.substr(0, 500) + unique.substr(0, 400)};
            std::vector<panvar::ReadPair> fragments;
            for (std::size_t h = 0; h < haplotypes.size(); ++h) {
                for (std::size_t start = 0; start + 260 <= haplotypes[h].size(); start += 17) {
                    std::string forward = haplotypes[h].substr(start, 100);
                    std::string reverse = revcomp(haplotypes[h].substr(start + 160, 100));
                    if (rng() % 3 == 0) forward[rng() % 100] = "ACGT"[rng() % 4];
                    if (rng() % 3 == 0) reverse[rng() % 100] = "ACGT"[rng() % 4];
                    fragments.push_back({read_of("f/1", forward), read_of("f/2", reverse)});
                }
            }
            std::vector<panvar::MosaicPath> paths;
            for (std::size_t h = 0; h < haplotypes.size(); ++h) paths.push_back(panvar::MosaicPath{{h}});
            const panvar::InternedHaplotypes interned = panvar::intern_haplotypes(paths, haplotypes);
            std::vector<panvar::HaplotypeCandidate> candidates;
            for (std::size_t x = 0; x < interned.sequences.size(); ++x) {
                for (std::size_t y = x; y < interned.sequences.size(); ++y) candidates.push_back({x, y});
            }
            panvar::LibraryModel pairs = library;
            pairs.support_sigmas = 3.0;
            pairs.outlier_weight = 0.005;
            const auto run = [&](bool prefix, std::size_t budget, bool collapse) {
                panvar::FragmentScoringOptions o;
                o.minimum_piece_bp = 10;
                o.threads = 1;
                o.emission_cache_bytes = budget;
                panvar::set_emission_key_prefix(prefix);
                panvar::set_emission_cache_hash_collapse(collapse);
                panvar::FragmentScoringResult result = panvar::score_diplotypes(interned, candidates, fragments, pairs, o);
                panvar::set_emission_cache_hash_collapse(false);
                panvar::set_emission_key_prefix(true);
                return result;
            };
            std::uint64_t compared = 0;
            for (const std::size_t budget : {std::size_t{1} << 24, std::size_t{20000}, std::size_t{1}}) {
                for (const bool collapse : {false, true}) {
                    const panvar::FragmentScoringResult reference = run(false, budget, collapse);
                    const panvar::FragmentScoringResult resumed = run(true, budget, collapse);
                    require(reference.log_mass.size() == resumed.log_mass.size() &&
                                std::memcmp(reference.log_mass.data(), resumed.log_mass.data(),
                                            reference.log_mass.size() * sizeof(double)) == 0,
                            "the resumed key changed a mass cell");
                    const panvar::EmissionCacheCounters& a = reference.emission_cache;
                    const panvar::EmissionCacheCounters& b = resumed.emission_cache;
                    require(a.executed_prefilters == b.executed_prefilters && a.dp_cache_requests == b.dp_cache_requests &&
                                a.dp_cache_hits == b.dp_cache_hits && a.dp_cache_misses == b.dp_cache_misses &&
                                a.executed_dps == b.executed_dps && a.entries == b.entries &&
                                a.entries_not_retained == b.entries_not_retained && a.payload_bytes == b.payload_bytes &&
                                a.budgeted_bytes_used == b.budgeted_bytes_used,
                            "the resumed key changed a hit, miss, executed DP, retention or budget decision (budget " +
                                std::to_string(budget) + (collapse ? ", collapsed)" : ")"));
                    require(a.dp_cache_requests > 0 && (budget == 1 || a.dp_cache_hits > 0),
                            "the key fixture must exercise the cache (requests " + std::to_string(a.dp_cache_requests) +
                                ", hits " + std::to_string(a.dp_cache_hits) + ", entries " + std::to_string(a.entries) + ", misses " + std::to_string(a.dp_cache_misses) +
                                ", executed " + std::to_string(a.executed_dps) + ", not retained " +
                                std::to_string(a.entries_not_retained) + ", budgeted " +
                                std::to_string(a.budgeted_bytes_used) + ")");
                    if (budget == 1) require(a.entries_not_retained > 0, "a one-byte budget must decline retention");
                    ++compared;
                }
            }
            std::cout << "  ok: emission key -- " << keys << " resumed keys equal to the full hash; " << compared
                      << " scorings (full, undersized and one-byte budgets, with and without forced collisions): mass "
                         "cells and every cache counter identical to hashing the whole read\n";
        }
        std::cout << "seeding exactness gates passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "seeding exactness gate FAILED: " << error.what() << "\n";
        return 1;
    }
}
