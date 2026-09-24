// Exact union placement reuse: every haploid mass cell, exposure, candidate component and ranking must equal
// the reference scorer's (reuse off) bit for bit, whatever the segmentation, budget, thread count or order.
#include "panvar/fragment_likelihood.hpp"
#include "panvar/genotype_model.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <map>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

bool same_double(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }

std::string revcomp(const std::string& bases) {
    std::string out(bases.size(), 'N');
    for (std::size_t i = 0; i < bases.size(); ++i) {
        const char base = bases[bases.size() - 1 - i];
        out[i] = base == 'A' ? 'T' : base == 'C' ? 'G' : base == 'G' ? 'C' : 'A';
    }
    return out;
}

panvar::SequencingRead read_of(const std::string& name, const std::string& bases) {
    return {name, bases, std::string(bases.size(), 'I')};
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

panvar::FragmentScoringOptions options(std::size_t threads) {
    panvar::FragmentScoringOptions out;
    out.mismatch_band = 6;
    out.minimum_piece_bp = 10;
    out.threads = threads;
    out.piece_hit_publish_batch = 1;
    out.emission_cache_bytes = 1u << 20;
    return out;
}

std::string random_bases(std::mt19937_64& rng, std::size_t length) {
    static const char* alphabet = "ACGT";
    std::string out(length, 'A');
    for (char& base : out) base = alphabet[rng() % 4];
    return out;
}

panvar::ReadPair fragment_at(const std::string& source, std::size_t start, std::size_t insert, std::size_t read_length,
                             const std::string& name, std::mt19937_64& rng) {
    std::string forward = source.substr(start, read_length);
    std::string reverse = revcomp(source.substr(start + insert - read_length, read_length));
    // A sprinkling of sequencing errors, so verification and the gapped emission really work.
    static const char* alphabet = "ACGT";
    if (rng() % 3 == 0) forward[rng() % forward.size()] = alphabet[rng() % 4];
    if (rng() % 3 == 0) reverse[rng() % reverse.size()] = alphabet[rng() % 4];
    return {read_of(name + "/1", forward), read_of(name + "/2", reverse)};
}

// A locus of five haplotypes built from segments: shared long segments at different coordinates and twice
// within one haplotype; substitution, insertion, deletion and empty alleles; short segments a read crosses
// several of; a tandem repeat; shared segments at both locus ends.
struct Fixture {
    std::vector<std::string> haplotypes;
    std::vector<std::vector<std::uint32_t>> segments;
    std::vector<panvar::ReadPair> fragments;
};

Fixture fixture() {
    std::mt19937_64 rng(20260925);
    const std::string A = random_bases(rng, 900);
    const std::string B = random_bases(rng, 700);
    std::string C = random_bases(rng, 650);
    std::string C_sub = C;
    C_sub[300] = C_sub[300] == 'A' ? 'C' : 'A';
    std::string C_ins = C.substr(0, 320) + random_bases(rng, 37) + C.substr(320);
    std::string C_del = C.substr(0, 280) + C.substr(330);
    const std::string S1 = random_bases(rng, 40);
    const std::string S2 = random_bases(rng, 25);
    const std::string unit = random_bases(rng, 41);
    std::string T;
    for (int i = 0; i < 10; ++i) T += unit;
    const std::string empty;
    const std::vector<std::vector<const std::string*>> layout = {
        {&A, &S1, &C, &T, &B},
        {&A, &S2, &C_sub, &T, &B},
        {&B, &S1, &C_ins, &A},
        {&A, &C_del, &empty, &A, &B},
        {&T, &S1, &S2, &C, &A},
    };
    Fixture out;
    for (const auto& parts : layout) {
        std::string haplotype;
        std::vector<std::uint32_t> lengths;
        for (const std::string* part : parts) {
            haplotype += *part;
            lengths.push_back(static_cast<std::uint32_t>(part->size()));
        }
        out.haplotypes.push_back(haplotype);
        out.segments.push_back(lengths);
    }
    std::size_t n = 0;
    for (std::size_t h = 0; h < out.haplotypes.size(); ++h) {
        const std::string& source = out.haplotypes[h];
        for (std::size_t start = 0; start + 260 <= source.size(); start += 29) {
            out.fragments.push_back(fragment_at(source, start, 260, 100, "f" + std::to_string(n++), rng));
        }
        // Flush against both locus ends.
        out.fragments.push_back(fragment_at(source, 0, 260, 100, "first" + std::to_string(h), rng));
        out.fragments.push_back(fragment_at(source, source.size() - 260, 260, 100, "last" + std::to_string(h), rng));
    }
    return out;
}

struct Run {
    panvar::InternedHaplotypes interned;
    std::vector<panvar::HaplotypeCandidate> candidates;
    panvar::FragmentScoringResult result;
    panvar::CertifiedRanking ranking;
};

Run score(const Fixture& f, const std::vector<std::size_t>& order, std::size_t threads,
          const std::vector<std::vector<std::uint32_t>>* segments, std::size_t budget, std::uint64_t cap) {
    Run run;
    std::vector<panvar::MosaicPath> paths;
    std::vector<std::string> sequences;
    for (std::size_t i = 0; i < order.size(); ++i) {
        paths.push_back(panvar::MosaicPath{{i}});
        sequences.push_back(f.haplotypes[order[i]]);
    }
    run.interned = panvar::intern_haplotypes(paths, sequences);
    for (std::size_t x = 0; x < run.interned.sequences.size(); ++x) {
        for (std::size_t y = x; y < run.interned.sequences.size(); ++y) run.candidates.push_back({x, y});
    }
    std::vector<std::vector<std::uint32_t>> ordered;
    if (segments != nullptr) {
        ordered.resize(run.interned.sequences.size());
        for (std::size_t i = 0; i < order.size(); ++i) ordered[run.interned.index_of[i]] = (*segments)[order[i]];
    }
    panvar::FragmentScoringOptions o = options(threads);
    if (segments != nullptr) {
        o.union_segment_lengths = &ordered;
        o.union_placement_cache_bytes = budget;
        o.max_union_cached_placements_consumed = cap;
    }
    run.result = panvar::score_diplotypes(run.interned, run.candidates, f.fragments, library(), o);
    run.ranking = panvar::certify_scored_ranking(run.result.scores, f.fragments.size(), o.mismatch_band,
                                                 run.result.counters, panvar::CertificationOptions{});
    return run;
}

// Compares two runs by haplotype SEQUENCE and candidate sequence pair, so permuted orders can be compared.
void require_identical(const Run& reference, const Run& other, const std::string& label) {
    std::map<std::string, std::size_t> slot_of;
    for (std::size_t h = 0; h < other.interned.sequences.size(); ++h) slot_of[other.interned.sequences[h]] = h;
    const std::size_t fragments = reference.result.log_mass.size() / reference.interned.sequences.size();
    for (std::size_t h = 0; h < reference.interned.sequences.size(); ++h) {
        const std::size_t k = slot_of.at(reference.interned.sequences[h]);
        require(same_double(reference.result.exposure[h], other.result.exposure[k]), label + ": exposure differs");
        for (std::size_t f = 0; f < fragments; ++f) {
            require(same_double(reference.result.mass_of(f, h), other.result.mass_of(f, k)),
                    label + ": mass cell differs at fragment " + std::to_string(f) + " haplotype " + std::to_string(h));
        }
    }
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> candidate_of;
    for (std::size_t c = 0; c < other.candidates.size(); ++c) {
        candidate_of[std::minmax(other.candidates[c].first, other.candidates[c].second)] = c;
    }
    std::vector<std::size_t> reference_to_other(reference.candidates.size());
    for (std::size_t c = 0; c < reference.candidates.size(); ++c) {
        const std::size_t x = slot_of.at(reference.interned.sequences[reference.candidates[c].first]);
        const std::size_t y = slot_of.at(reference.interned.sequences[reference.candidates[c].second]);
        const std::size_t d = candidate_of.at(std::minmax(x, y));
        reference_to_other[c] = d;
        const panvar::DiplotypeScore& a = reference.result.scores[c];
        const panvar::DiplotypeScore& b = other.result.scores[d];
        require(same_double(a.log_likelihood, b.log_likelihood) && same_double(a.fragment_term, b.fragment_term) &&
                    same_double(a.exposure_term, b.exposure_term) && same_double(a.absolute_term_sum, b.absolute_term_sum) &&
                    a.fragments_with_mass == b.fragments_with_mass,
                label + ": candidate components differ");
    }
    require(reference.ranking.status == other.ranking.status, label + ": ranking status differs");
    require(reference.ranking.candidates.size() == other.ranking.candidates.size(), label + ": ranking size differs");
    for (std::size_t i = 0; i < reference.ranking.candidates.size(); ++i) {
        require(reference_to_other[reference.ranking.candidates[i].candidate_id] == other.ranking.candidates[i].candidate_id &&
                    same_double(reference.ranking.candidates[i].score, other.ranking.candidates[i].score),
                label + ": certified ranking differs at rank " + std::to_string(i));
    }
    std::vector<std::size_t> ws_ref, ws_other;
    for (const std::size_t w : reference.ranking.winner_set) ws_ref.push_back(reference_to_other[w]);
    ws_other = other.ranking.winner_set;
    std::sort(ws_ref.begin(), ws_ref.end());
    std::sort(ws_other.begin(), ws_other.end());
    require(ws_ref == ws_other, label + ": winner set differs");
    // The logical placement and origin work is the reference's; only executed seeding may be lower.
    require(reference.result.counters.accepted_mate_placements == other.result.counters.accepted_mate_placements &&
                reference.result.counters.physical_origins == other.result.counters.physical_origins &&
                reference.result.counters.statistical_state_classes == other.result.counters.statistical_state_classes,
            label + ": accepted placements, physical origins or state classes differ");
}

}  // namespace

int main() {
    try {
        const Fixture f = fixture();
        const std::vector<std::size_t> natural = {0, 1, 2, 3, 4};
        const Run reference = score(f, natural, 1, nullptr, 0, 0);
        std::cout << "  reference: " << f.fragments.size() << " fragments, " << reference.interned.sequences.size()
                  << " haplotypes, " << reference.result.counters.accepted_mate_placements << " accepted placements, "
                  << reference.result.counters.physical_origins << " physical origins\n";

        const std::size_t full = std::size_t{1} << 30;
        const Run retained = score(f, natural, 1, &f.segments, full, 0);
        const panvar::UnionPlacementCacheReport& report = retained.result.union_placement_cache;
        require(report.active && report.admitted_segments > 0 && report.hits > 0 && report.placements_replayed > 0,
                "the fixture must actually reuse placements");
        require(retained.result.counters.full_read_verifications < reference.result.counters.full_read_verifications,
                "reuse must execute fewer verifications");
        require_identical(reference, retained, "fully retained");
        std::cout << "  ok: fully retained -- " << report.admitted_segments << " segments, " << report.hits << " hits, "
                  << report.placements_replayed << " placements replayed, verifications "
                  << retained.result.counters.full_read_verifications << " of "
                  << reference.result.counters.full_read_verifications << "\n";

        for (const std::size_t threads : {std::size_t{2}, std::size_t{4}, std::size_t{8}}) {
            require_identical(reference, score(f, natural, threads, &f.segments, full, 0), std::to_string(threads) + " threads");
        }
        std::cout << "  ok: 1, 2, 4 and 8 threads\n";

        const Run disabled = score(f, natural, 1, &f.segments, 0, 0);
        require(!disabled.result.union_placement_cache.active &&
                    disabled.result.counters.full_read_verifications == reference.result.counters.full_read_verifications,
                "a zero budget must be the reference path");
        require_identical(reference, disabled, "disabled");
        std::cout << "  ok: disabled (budget 0) is the reference path\n";

        // Byte budget exact and one below: exact retains everything; one below declines at least one segment.
        const std::uint64_t exact = report.payload_bytes;
        const Run at_exact = score(f, natural, 1, &f.segments, exact, 0);
        require(at_exact.result.union_placement_cache.declined_segments == 0 &&
                    at_exact.result.union_placement_cache.payload_bytes == exact,
                "exactly the needed bytes must retain every eligible segment");
        require_identical(reference, at_exact, "budget exact");
        const Run below = score(f, natural, 1, &f.segments, exact - 1, 0);
        require(below.result.union_placement_cache.declined_segments >= 1 &&
                    below.result.union_placement_cache.payload_bytes <= exact - 1,
                "one byte below must decline a segment and hold the budget");
        require_identical(reference, below, "budget one below");
        const Run undersized = score(f, natural, 1, &f.segments, exact / 3, 0);
        require(undersized.result.union_placement_cache.payload_bytes <= exact / 3, "an undersized cache holds its budget");
        require_identical(reference, undersized, "undersized");
        std::cout << "  ok: byte budget exact (" << exact << ") and one below, undersized (" << exact / 3 << ")\n";

        // Consumption cap exact and one below.
        const std::uint64_t consumed = report.placements_replayed;
        require_identical(reference, score(f, natural, 1, &f.segments, full, consumed), "consumption cap exact");
        bool refused = false;
        try {
            (void)score(f, natural, 1, &f.segments, full, consumed - 1);
        } catch (const panvar::ResourceExhausted&) {
            refused = true;
        }
        require(refused, "one cached placement below the consumption cap must refuse with no ranking");
        std::cout << "  ok: cached-placement consumption cap exact (" << consumed << ") accepts, one below refuses\n";

        panvar::set_union_placement_cache_hash_collapse(true);
        const Run collided = score(f, natural, 1, &f.segments, full, 0);
        panvar::set_union_placement_cache_hash_collapse(false);
        require(collided.result.union_placement_cache.distinct_segments == report.distinct_segments &&
                    collided.result.union_placement_cache.hits == report.hits,
                "forced hash collisions must separate segments by bytes alone");
        require_identical(reference, collided, "forced hash collisions");
        std::cout << "  ok: forced hash collisions\n";

        const std::vector<std::size_t> permuted = {3, 0, 4, 2, 1};
        require_identical(reference, score(f, permuted, 4, &f.segments, full, 0), "permuted order");
        std::cout << "  ok: haplotype and candidate order permuted\n";

        // The segmentation is a hint: a wrong one (arbitrary cut points) stays exact, it only reuses less.
        std::vector<std::vector<std::uint32_t>> wrong;
        for (const std::string& h : f.haplotypes) {
            const std::uint32_t L = static_cast<std::uint32_t>(h.size());
            wrong.push_back({L / 3, L / 3, L - 2 * (L / 3)});
        }
        require_identical(reference, score(f, natural, 2, &wrong, full, 0), "wrong segmentation");
        std::vector<std::vector<std::uint32_t>> whole;
        for (const std::string& h : f.haplotypes) whole.push_back({static_cast<std::uint32_t>(h.size())});
        require_identical(reference, score(f, natural, 2, &whole, full, 0), "one segment per haplotype");
        std::cout << "  ok: a wrong segmentation and whole-haplotype segments stay exact\n";
        std::cout << "union placement cache gates passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "union placement cache gate FAILED: " << error.what() << "\n";
        return 1;
    }
}
