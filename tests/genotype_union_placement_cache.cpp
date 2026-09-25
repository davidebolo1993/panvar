// Exact union placement reuse (phase 1) and origin reuse (phase 2): every haploid mass cell, exposure, candidate
// component and ranking must equal the reference scorer's (reuse off) bit for bit, whatever the segmentation,
// budget, cap, thread count or order, and the logical placement, origin and state totals must reconcile with it.
#include "panvar/fragment_likelihood.hpp"
#include "panvar/genotype_model.hpp"

#include <algorithm>
#include <cmath>
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
// several of; a tandem repeat; shared segments at both locus ends; and a segment holding two copies of one
// stretch a base apart, so one fragment has origins of equal insert and orientation but different latent mass.
struct Fixture {
    std::vector<std::string> haplotypes;
    std::vector<std::vector<std::uint32_t>> segments;
    std::vector<panvar::ReadPair> fragments;
    std::vector<std::pair<std::size_t, std::size_t>> near_duplicate;   // (haplotype, offset) of each D
    std::size_t near_duplicate_bp = 0;
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
    const std::string X = random_bases(rng, 330);
    std::string X_mut = X;
    X_mut[170] = X_mut[170] == 'G' ? 'T' : 'G';
    const std::string D = X + random_bases(rng, 60) + X_mut;
    const std::vector<std::vector<const std::string*>> layout = {
        {&A, &S1, &C, &T, &B},
        {&A, &S2, &C_sub, &T, &D, &B},
        {&B, &S1, &C_ins, &D, &A},
        {&A, &C_del, &empty, &A, &B},
        {&T, &S1, &S2, &C, &D, &A},
    };
    Fixture out;
    out.near_duplicate_bp = D.size();
    for (const auto& parts : layout) {
        std::string haplotype;
        std::vector<std::uint32_t> lengths;
        for (const std::string* part : parts) {
            if (part == &D) out.near_duplicate.push_back({out.haplotypes.size(), haplotype.size()});
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

struct Config {
    std::vector<std::size_t> order = {0, 1, 2, 3, 4};
    std::size_t threads = 1;
    const std::vector<std::vector<std::uint32_t>>* segments = nullptr;
    std::size_t placement_budget = 0;
    std::uint64_t placement_cap = 0;
    std::size_t origin_budget = 0;
    std::uint64_t origins_cap = 0;
    std::uint64_t states_cap = 0;
    bool reverse_candidates = false;
    std::size_t prefilter_budget = 0;
    std::uint64_t prefilter_cap = 0;
    std::uint64_t max_physical_origins = 0;   // 0: the default
    std::uint64_t max_state_classes = 0;
    bool record_all_origins = false;
};

Run score(const Fixture& f, const Config& config) {
    Run run;
    std::vector<panvar::MosaicPath> paths;
    std::vector<std::string> sequences;
    for (std::size_t i = 0; i < config.order.size(); ++i) {
        paths.push_back(panvar::MosaicPath{{i}});
        sequences.push_back(f.haplotypes[config.order[i]]);
    }
    run.interned = panvar::intern_haplotypes(paths, sequences);
    for (std::size_t x = 0; x < run.interned.sequences.size(); ++x) {
        for (std::size_t y = x; y < run.interned.sequences.size(); ++y) run.candidates.push_back({x, y});
    }
    if (config.reverse_candidates) std::reverse(run.candidates.begin(), run.candidates.end());
    std::vector<std::vector<std::uint32_t>> ordered;
    if (config.segments != nullptr) {
        ordered.resize(run.interned.sequences.size());
        for (std::size_t i = 0; i < config.order.size(); ++i) {
            ordered[run.interned.index_of[i]] = (*config.segments)[config.order[i]];
        }
    }
    panvar::FragmentScoringOptions o = options(config.threads);
    if (config.segments != nullptr) {
        o.union_segment_lengths = &ordered;
        o.union_placement_cache_bytes = config.placement_budget;
        o.max_union_cached_placements_consumed = config.placement_cap;
        o.union_origin_cache_bytes = config.origin_budget;
        o.max_union_cached_origins_consumed = config.origins_cap;
        o.max_union_cached_states_consumed = config.states_cap;
    }
    o.union_prefilter_cache_bytes = config.prefilter_budget;
    o.max_union_cached_prefilters_consumed = config.prefilter_cap;
    if (config.max_physical_origins != 0) o.max_physical_origins = config.max_physical_origins;
    if (config.max_state_classes != 0) o.max_state_classes = config.max_state_classes;
    if (config.record_all_origins) {
        for (std::size_t i = 0; i < f.fragments.size(); ++i) o.record_origins_for_fragments.push_back(i);
    }
    run.result = panvar::score_diplotypes(run.interned, run.candidates, f.fragments, library(), o);
    run.ranking = panvar::certify_scored_ranking(run.result.scores, f.fragments.size(), o.mismatch_band,
                                                 run.result.counters, panvar::CertificationOptions{});
    return run;
}

// The reference's totals recovered from a run: executed, minus what phase 2 constructed once, plus what it served.
struct Logical {
    std::uint64_t accepted = 0, origins = 0, pairs = 0, valid = 0, classes = 0;
    bool operator==(const Logical& o) const {
        return accepted == o.accepted && origins == o.origins && pairs == o.pairs && valid == o.valid &&
               classes == o.classes;
    }
};
Logical logical(const Run& run) {
    const panvar::FragmentScoringCounters& c = run.result.counters;
    const panvar::UnionOriginCacheReport& o = run.result.union_origin_cache;
    return {c.accepted_mate_placements, c.physical_origins - o.origins_constructed + c.reused_physical_origins,
            c.latent_alignment_pairs - o.latent_pairs_constructed + c.reused_latent_alignment_pairs,
            c.valid_fr_states - o.valid_states_constructed + c.reused_valid_fr_states, c.merged_state_classes};
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
    // The logical placement, origin, pair and state work is the reference's; only executed work may be lower.
    require(logical(reference) == logical(other),
            label + ": logical accepted placements, physical origins, latent pairs, valid states or state classes differ");
}

// Prefilter reuse against the same configuration without it: identical scores (require_identical) and identical
// logical verification counters; only kernels executed may fall.
void require_prefilter_exact(const Run& without, const Run& with, const std::string& label) {
    require_identical(without, with, label);
    const panvar::FragmentScoringCounters& a = without.result.counters;
    const panvar::FragmentScoringCounters& b = with.result.counters;
    require(a.full_read_verifications == b.full_read_verifications && a.prefilter_rejections == b.prefilter_rejections &&
                a.dp_evaluations == b.dp_evaluations && a.proposed_starts == b.proposed_starts &&
                a.unique_starts == b.unique_starts,
            label + ": a logical verification counter moved");
    const panvar::UnionPrefilterCacheReport& p = with.result.union_prefilter_cache;
    require(with.result.emission_cache.executed_prefilters + p.hits == without.result.emission_cache.executed_prefilters,
            label + ": executed prefilters + decisions served must equal the kernels run without reuse");
}

bool refuses(const Fixture& f, const Config& config) {
    try {
        (void)score(f, config);
    } catch (const panvar::ResourceExhausted&) {
        return true;
    }
    return false;
}

}  // namespace

int main() {
    try {
        const Fixture f = fixture();
        const Config off;
        const Run reference = score(f, off);
        require(reference.result.counters.merged_state_classes == reference.result.counters.statistical_state_classes &&
                    reference.result.counters.reused_physical_origins == 0,
                "the reference path's merged classes are its classes, and it serves nothing");
        std::cout << "  reference: " << f.fragments.size() << " fragments, " << reference.interned.sequences.size()
                  << " haplotypes, " << reference.result.counters.accepted_mate_placements << " accepted placements, "
                  << reference.result.counters.physical_origins << " physical origins, "
                  << reference.result.counters.statistical_state_classes << " state classes\n";

        // ---- phase 1: placement reuse -------------------------------------------------------------------
        const std::size_t full = std::size_t{1} << 30;
        Config p1;
        p1.segments = &f.segments;
        p1.placement_budget = full;
        const Run retained = score(f, p1);
        const panvar::UnionPlacementCacheReport& report = retained.result.union_placement_cache;
        require(report.active && report.admitted_segments > 0 && report.hits > 0 && report.placements_replayed > 0,
                "the fixture must actually reuse placements");
        require(retained.result.counters.full_read_verifications < reference.result.counters.full_read_verifications,
                "reuse must execute fewer verifications");
        require(!retained.result.union_origin_cache.active, "phase 2 is off without its budget");
        require_identical(reference, retained, "fully retained");
        std::cout << "  ok: fully retained -- " << report.admitted_segments << " segments, " << report.hits << " hits, "
                  << report.placements_replayed << " placements replayed, verifications "
                  << retained.result.counters.full_read_verifications << " of "
                  << reference.result.counters.full_read_verifications << "\n";

        for (const std::size_t threads : {std::size_t{2}, std::size_t{4}, std::size_t{8}}) {
            Config c = p1;
            c.threads = threads;
            const Run threaded = score(f, c);
            require_identical(reference, threaded, std::to_string(threads) + " threads");
            // The fill's contents, payload and admission do not depend on who filled what.
            const panvar::UnionPlacementCacheReport& t = threaded.result.union_placement_cache;
            require(t.admitted_segments == report.admitted_segments && t.declined_segments == report.declined_segments &&
                        t.placements_enumerated == report.placements_enumerated &&
                        t.placements_cached == report.placements_cached && t.payload_bytes == report.payload_bytes &&
                        t.hits == report.hits && t.placements_replayed == report.placements_replayed,
                    std::to_string(threads) + " threads: the placement cache's contents or payload differ");
        }
        std::cout << "  ok: 1, 2, 4 and 8 threads (cache contents, payload and admission invariant)\n";

        Config zero = p1;
        zero.placement_budget = 0;
        const Run disabled = score(f, zero);
        require(!disabled.result.union_placement_cache.active &&
                    disabled.result.counters.full_read_verifications == reference.result.counters.full_read_verifications,
                "a zero budget must be the reference path");
        require_identical(reference, disabled, "disabled");
        std::cout << "  ok: disabled (budget 0) is the reference path\n";

        // Byte budget exact and one below: exact retains everything; one below declines at least one segment.
        const std::uint64_t exact = report.payload_bytes;
        Config c_exact = p1;
        c_exact.placement_budget = exact;
        const Run at_exact = score(f, c_exact);
        require(at_exact.result.union_placement_cache.declined_segments == 0 &&
                    at_exact.result.union_placement_cache.payload_bytes == exact,
                "exactly the needed bytes must retain every eligible segment");
        require_identical(reference, at_exact, "budget exact");
        Config c_below = p1;
        c_below.placement_budget = exact - 1;
        const Run below = score(f, c_below);
        require(below.result.union_placement_cache.declined_segments >= 1 &&
                    below.result.union_placement_cache.payload_bytes <= exact - 1,
                "one byte below must decline a segment and hold the budget");
        require_identical(reference, below, "budget one below");
        Config c_under = p1;
        c_under.placement_budget = exact / 3;
        const Run undersized = score(f, c_under);
        require(undersized.result.union_placement_cache.payload_bytes <= exact / 3, "an undersized cache holds its budget");
        require_identical(reference, undersized, "undersized");
        std::cout << "  ok: byte budget exact (" << exact << ") and one below, undersized (" << exact / 3 << ")\n";

        // Consumption cap exact and one below.
        const std::uint64_t consumed = report.placements_replayed;
        Config c_cap = p1;
        c_cap.placement_cap = consumed;
        require_identical(reference, score(f, c_cap), "consumption cap exact");
        c_cap.placement_cap = consumed - 1;
        require(refuses(f, c_cap), "one cached placement below the consumption cap must refuse with no ranking");
        std::cout << "  ok: cached-placement consumption cap exact (" << consumed << ") accepts, one below refuses\n";

        panvar::set_union_placement_cache_hash_collapse(true);
        const Run collided = score(f, p1);
        panvar::set_union_placement_cache_hash_collapse(false);
        require(collided.result.union_placement_cache.distinct_segments == report.distinct_segments &&
                    collided.result.union_placement_cache.hits == report.hits,
                "forced hash collisions must separate segments by bytes alone");
        require_identical(reference, collided, "forced hash collisions");
        std::cout << "  ok: forced hash collisions\n";

        Config c_perm = p1;
        c_perm.order = {3, 0, 4, 2, 1};
        c_perm.threads = 4;
        require_identical(reference, score(f, c_perm), "permuted order");
        std::cout << "  ok: haplotype and candidate order permuted\n";

        // The segmentation is a hint: a wrong one (arbitrary cut points) stays exact, it only reuses less.
        std::vector<std::vector<std::uint32_t>> wrong;
        for (const std::string& h : f.haplotypes) {
            const std::uint32_t L = static_cast<std::uint32_t>(h.size());
            wrong.push_back({L / 3, L / 3, L - 2 * (L / 3)});
        }
        Config c_wrong = p1;
        c_wrong.segments = &wrong;
        c_wrong.threads = 2;
        require_identical(reference, score(f, c_wrong), "wrong segmentation");
        std::vector<std::vector<std::uint32_t>> whole;
        for (const std::string& h : f.haplotypes) whole.push_back({static_cast<std::uint32_t>(h.size())});
        Config c_whole = c_wrong;
        c_whole.segments = &whole;
        require_identical(reference, score(f, c_whole), "one segment per haplotype");
        std::cout << "  ok: a wrong segmentation and whole-haplotype segments stay exact\n";

        // ---- phase 2: origin reuse ----------------------------------------------------------------------
        const std::uint64_t unbounded = std::uint64_t{1} << 62;
        Config p2 = p1;
        p2.origin_budget = full;
        p2.origins_cap = unbounded;
        p2.states_cap = unbounded;
        const Run served = score(f, p2);
        const panvar::UnionOriginCacheReport& oreport = served.result.union_origin_cache;
        const panvar::FragmentScoringCounters& sc = served.result.counters;
        require(oreport.active && oreport.admitted_segments > 0 && oreport.hits > 0 && sc.reused_physical_origins > 0 &&
                    sc.reused_state_aggregates > 0,
                "the fixture must actually serve origins");
        require(oreport.hits > oreport.admitted_segments, "every admitted segment is served at several occurrences");
        require(sc.physical_origins < reference.result.counters.physical_origins &&
                    sc.statistical_state_classes < reference.result.counters.statistical_state_classes,
                "origin reuse must execute fewer origins and state classes");
        require(oreport.declined_segments == 0 && oreport.evictions == 0 &&
                    oreport.misses == oreport.segment_occurrences - oreport.hits,
                "a full budget declines and evicts nothing");
        require_identical(reference, served, "origin reuse");
        std::cout << "  ok: origin reuse -- " << oreport.admitted_segments << " of " << oreport.eligible_segments
                  << " segments, " << oreport.hits << " hits, " << oreport.misses << " misses; origins executed "
                  << sc.physical_origins << " (" << oreport.origins_constructed << " constructed once) + "
                  << sc.reused_physical_origins << " served = reference " << reference.result.counters.physical_origins
                  << "; state classes executed " << sc.statistical_state_classes << ", aggregates served "
                  << sc.reused_state_aggregates << ", merged " << sc.merged_state_classes << "\n";

        // A fragment with two origins of equal insert and orientation but different latent mass, both inside the
        // served interior of one near-duplicate segment D (copies a base apart): the reference's census, recorded
        // origin by origin, must contain one, so the merge is shown to keep such classes apart.
        {
            Config recorded;
            recorded.record_all_origins = true;
            const Run named = score(f, recorded);
            const std::size_t R = 100;
            const std::size_t g = library().max_gap_bp;
            std::size_t found = 0;
            for (const auto& [key, origins] : named.result.recorded_origins) {
                const std::string& sequence = named.interned.sequences[key.second];
                for (const auto& [h, offset] : f.near_duplicate) {
                    if (f.haplotypes[h] != sequence) continue;
                    std::map<std::pair<std::size_t, int>, std::set<double>> masses;
                    for (const panvar::FragmentPlacement& origin : origins) {
                        const std::size_t start = origin.forward_start;
                        const std::size_t end = start + origin.insert_bp;
                        if (start < offset + g || start > offset + f.near_duplicate_bp - R - g) continue;
                        if (end < offset + R + g || end > offset + f.near_duplicate_bp - g) continue;
                        masses[{origin.insert_bp, static_cast<int>(origin.orientation)}].insert(origin.log_latent_mass);
                    }
                    for (const auto& item : masses) found += item.second.size() > 1 ? 1 : 0;
                }
            }
            require(found > 0, "the near-duplicate segment must give some fragment origins of equal insert and "
                               "orientation but different latent mass inside a served interior");
            std::cout << "  ok: " << found << " (fragment, occurrence, insert, orientation) classes inside served "
                      << "near-duplicate interiors carry more than one latent mass\n";
        }

        for (const std::size_t threads : {std::size_t{2}, std::size_t{4}, std::size_t{8}}) {
            Config c = p2;
            c.threads = threads;
            require_identical(reference, score(f, c), "origin reuse, " + std::to_string(threads) + " threads");
        }
        std::cout << "  ok: origin reuse at 1, 2, 4 and 8 threads\n";

        Config p2_zero = p2;
        p2_zero.origin_budget = 0;
        const Run p2_disabled = score(f, p2_zero);
        require(!p2_disabled.result.union_origin_cache.active && p2_disabled.result.counters.reused_physical_origins == 0,
                "a zero origin budget is phase 1 alone");
        require_identical(reference, p2_disabled, "origin reuse disabled");

        const std::uint64_t oexact = oreport.payload_bytes;
        Config p2_exact = p2;
        p2_exact.origin_budget = oexact;
        const Run o_exact = score(f, p2_exact);
        require(o_exact.result.union_origin_cache.declined_segments == 0 &&
                    o_exact.result.union_origin_cache.payload_bytes == oexact,
                "exactly the needed origin bytes retain every eligible segment");
        require_identical(reference, o_exact, "origin budget exact");
        Config p2_below = p2;
        p2_below.origin_budget = oexact - 1;
        const Run o_below = score(f, p2_below);
        require(o_below.result.union_origin_cache.declined_segments >= 1 &&
                    o_below.result.union_origin_cache.payload_bytes <= oexact - 1,
                "one origin byte below declines a segment and holds the budget");
        require_identical(reference, o_below, "origin budget one below");
        Config p2_under = p2;
        p2_under.origin_budget = oexact / 4;
        const Run o_under = score(f, p2_under);
        require(o_under.result.union_origin_cache.payload_bytes <= oexact / 4, "an undersized origin cache holds its budget");
        require_identical(reference, o_under, "origin cache undersized");
        // Phase 2 needs both strands of a segment from phase 1: an undersized placement cache serves fewer origins.
        Config p2_p1under = p2;
        p2_p1under.placement_budget = exact / 3;
        require_identical(reference, score(f, p2_p1under), "origin reuse over an undersized placement cache");
        std::cout << "  ok: origin cache disabled, exact (" << oexact << ") and one below, undersized, and over an "
                  << "undersized placement cache\n";

        Config p2_cap = p2;
        p2_cap.origins_cap = sc.reused_physical_origins;
        p2_cap.states_cap = sc.reused_state_aggregates;
        require_identical(reference, score(f, p2_cap), "origin and state consumption caps exact");
        Config p2_origins_below = p2_cap;
        p2_origins_below.origins_cap = sc.reused_physical_origins - 1;
        require(refuses(f, p2_origins_below), "one served origin below its cap must refuse with no ranking");
        Config p2_states_below = p2_cap;
        p2_states_below.states_cap = sc.reused_state_aggregates - 1;
        require(refuses(f, p2_states_below), "one served state aggregate below its cap must refuse with no ranking");
        for (const std::size_t threads : {std::size_t{4}}) {
            Config c = p2_origins_below;
            c.threads = threads;
            require(refuses(f, c), "the served-origin cap refuses at 4 threads too");
        }
        Config p2_nocap = p2;
        p2_nocap.states_cap = 0;
        bool invalid = false;
        try {
            (void)score(f, p2_nocap);
        } catch (const std::invalid_argument&) {
            invalid = true;
        }
        require(invalid, "origin reuse without a positive consumption cap is refused as a configuration");
        std::cout << "  ok: served-origin cap exact (" << sc.reused_physical_origins << ") and served-state cap exact ("
                  << sc.reused_state_aggregates << ") accept, one below each refuses; a zero cap is invalid\n";

        // The executed caps bind executed work only: the executed origins and classes are exactly enough under
        // reuse, one fewer refuses, and the reference refuses at the same caps.
        Config p2_exec = p2;
        p2_exec.max_physical_origins = sc.physical_origins;
        p2_exec.max_state_classes = sc.statistical_state_classes;
        require_identical(reference, score(f, p2_exec), "executed caps exact under origin reuse");
        Config p2_exec_origins = p2_exec;
        p2_exec_origins.max_physical_origins = sc.physical_origins - 1;
        require(refuses(f, p2_exec_origins), "one executed origin below max_physical_origins refuses");
        Config p2_exec_states = p2_exec;
        p2_exec_states.max_state_classes = sc.statistical_state_classes - 1;
        require(refuses(f, p2_exec_states), "one executed class below max_state_classes refuses");
        Config ref_at_exec = off;
        ref_at_exec.max_physical_origins = sc.physical_origins;
        require(refuses(f, ref_at_exec), "the reference, which executes every origin, refuses at the reused run's cap");
        std::cout << "  ok: max_physical_origins (" << sc.physical_origins << ") and max_state_classes ("
                  << sc.statistical_state_classes << ") bind executed work: exact accepts, one below refuses\n";

        panvar::set_union_placement_cache_hash_collapse(true);
        const Run o_collided = score(f, p2);
        panvar::set_union_placement_cache_hash_collapse(false);
        require(o_collided.result.union_origin_cache.hits == oreport.hits, "forced collisions serve the same occurrences");
        require_identical(reference, o_collided, "origin reuse, forced hash collisions");

        Config p2_perm = p2;
        p2_perm.order = {3, 0, 4, 2, 1};
        p2_perm.threads = 4;
        require_identical(reference, score(f, p2_perm), "origin reuse, permuted haplotypes");
        Config p2_req = p2;
        p2_req.reverse_candidates = true;
        p2_req.threads = 2;
        require_identical(reference, score(f, p2_req), "origin reuse, reversed candidate requests");
        Config p2_wrong = p2;
        p2_wrong.segments = &wrong;
        require_identical(reference, score(f, p2_wrong), "origin reuse, wrong segmentation");
        Config p2_whole = p2;
        p2_whole.segments = &whole;
        require_identical(reference, score(f, p2_whole), "origin reuse, whole-haplotype segments");
        std::cout << "  ok: origin reuse under forced collisions, permuted haplotypes, reversed requests, wrong and "
                  << "whole-haplotype segmentations\n";

        // ---- prefilter reuse --------------------------------------------------------------------------
        const std::uint64_t many = std::uint64_t{1} << 40;
        for (const Config* base : std::initializer_list<const Config*>{&off, &p1, &p2}) {
            Config with = *base;
            with.prefilter_budget = std::size_t{1} << 24;
            with.prefilter_cap = many;
            const Run plain = score(f, *base);
            const Run served_run = score(f, with);
            const panvar::UnionPrefilterCacheReport& p = served_run.result.union_prefilter_cache;
            require(p.active && p.hits > 0 && p.window_contents > 0 && p.decisions_held > 0,
                    "the fixture must actually serve prefilter decisions");
            require_prefilter_exact(plain, served_run, "prefilter reuse");
            require_identical(reference, served_run, "prefilter reuse against the reference");
            if (base == &p2) {
                std::cout << "  ok: prefilter reuse over phases 1+2 -- " << p.window_contents << " window contents, "
                          << p.decisions_held << " decisions held, " << p.hits << " served, " << p.misses
                          << " kernels run (without reuse "
                          << plain.result.emission_cache.executed_prefilters << ")\n";
                for (const std::size_t threads : {std::size_t{2}, std::size_t{4}, std::size_t{8}}) {
                    Config c = with;
                    c.threads = threads;
                    require_prefilter_exact(plain, score(f, c), "prefilter reuse, " + std::to_string(threads) + " threads");
                }
                Config permuted = with;
                permuted.order = {3, 0, 4, 2, 1};
                permuted.threads = 4;
                require_identical(reference, score(f, permuted), "prefilter reuse, permuted haplotypes");
                Config requests = with;
                requests.reverse_candidates = true;
                require_identical(reference, score(f, requests), "prefilter reuse, reversed requests");
                panvar::set_union_prefilter_cache_hash_collapse(true);
                const Run collided_run = score(f, with);
                panvar::set_union_prefilter_cache_hash_collapse(false);
                require(collided_run.result.union_prefilter_cache.window_contents == p.window_contents &&
                            collided_run.result.union_prefilter_cache.hits == p.hits,
                        "forced window-hash collisions must intern the same contents and serve the same decisions");
                require_prefilter_exact(plain, collided_run, "prefilter reuse, forced hash collisions");
                Config tiny = with;
                tiny.prefilter_budget = 1;   // the smallest table: admission stops at once, the kernel runs instead
                const Run tiny_run = score(f, tiny);
                require(tiny_run.result.union_prefilter_cache.decisions_not_retained > 0,
                        "a tiny decision table must decline to retain");
                require_prefilter_exact(plain, tiny_run, "prefilter reuse, undersized table");
                Config at_cap = with;
                at_cap.prefilter_cap = p.hits;
                require_prefilter_exact(plain, score(f, at_cap), "prefilter consumption cap exact");
                Config below = with;
                below.prefilter_cap = p.hits - 1;
                require(refuses(f, below), "one served decision below the cap must refuse with no ranking");
                below.threads = 4;
                const Run four = score(f, [&] { Config c = with; c.threads = 4; return c; }());
                below.prefilter_cap = four.result.union_prefilter_cache.hits - 1;
                require(refuses(f, below), "the prefilter cap refuses at 4 threads too");
                Config zero_cap = with;
                zero_cap.prefilter_cap = 0;
                bool invalid = false;
                try {
                    (void)score(f, zero_cap);
                } catch (const std::invalid_argument&) {
                    invalid = true;
                }
                require(invalid, "prefilter reuse without a positive cap is refused as a configuration");
            }
        }
        std::cout << "  ok: prefilter reuse alone, over phase 1 and over phases 1+2; 1/2/4/8 threads, permuted haplotypes, "
                     "reversed requests, forced window-hash collisions, an undersized table, the consumption cap exact "
                     "and one below (1 and 4 threads), a zero cap invalid\n";

        panvar::set_origin_cache_self_check(true);
        for (const std::size_t threads : {std::size_t{1}, std::size_t{4}}) {
            Config c = p2;
            c.threads = threads;
            require_identical(reference, score(f, c), "origin reuse self-check");
        }
        panvar::set_origin_cache_self_check(false);
        std::cout << "  ok: every served fragment census equals a fresh enumeration (self-check, 1 and 4 threads)\n";
        std::cout << "union placement cache gates passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "union placement cache gate FAILED: " << error.what() << "\n";
        return 1;
    }
}
