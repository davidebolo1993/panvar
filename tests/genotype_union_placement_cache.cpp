// Exact union placement reuse (phase 1) and origin reuse (phase 2): every haploid mass cell, exposure, candidate
// component and ranking must equal the reference scorer's (reuse off) bit for bit, whatever the segmentation,
// budget, cap, thread count or order, and the logical placement, origin and state totals must reconcile with it.
#include "panvar/fragment_likelihood.hpp"
#include "panvar/genotype_model.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <functional>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
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
    std::size_t sort_scratch = 0;
    std::uint64_t max_physical_origins = 0;   // 0: the default
    std::uint64_t max_state_classes = 0;
    bool record_all_origins = false;
    // I2-H
    panvar::GapRecursionBackend backend = panvar::GapRecursionBackend::scaled_reference;
    std::size_t band = 0;                                           // 0: the default band
    std::size_t emission_cache_bytes = std::size_t(-1);             // -1: the default
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
    o.start_sort_scratch_bytes = config.sort_scratch;
    o.max_union_cached_prefilters_consumed = config.prefilter_cap;
    if (config.max_physical_origins != 0) o.max_physical_origins = config.max_physical_origins;
    if (config.max_state_classes != 0) o.max_state_classes = config.max_state_classes;
    o.gap_recursion_backend = config.backend;
    if (config.band != 0) o.mismatch_band = config.band;
    if (config.emission_cache_bytes != std::size_t(-1)) o.emission_cache_bytes = config.emission_cache_bytes;
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
void require_identical(const Run& reference, const Run& other, const std::string& label, bool logical_counters = true) {
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
    // (The content index keeps no placement counters; its origin and state totals are compared separately.)
    require(!logical_counters || logical(reference) == logical(other),
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


// ---- EXPERIMENTAL (I1): shared-content placement index --------------------------------------------------------------
// Every case compares every haploid mass, exposure, candidate component, the certified ranking and the winner set with
// the reference scorer bit for bit (require_identical), plus the physical-origin and state-class totals.

// One index run: registers the run's unique sequences (in interned order, with the given segmentation), discovers,
// scores the same candidates. Ranking certified with the reference's counters (certification reads only the scores
// and the band here; the index does not keep the scorer's placement counters).
// I1.1: `prefilter_bytes` > 0 gives the batch the exact prefilter-decision cache (the reference's configuration).
panvar::FragmentScoringOptions index_options(std::size_t threads, std::size_t prefilter_bytes, std::size_t origin_bytes = 0) {
    panvar::FragmentScoringOptions o = options(threads);
    o.union_origin_cache_bytes = origin_bytes;   // I1.2: the batch's interior-origin aggregate budget (0: off)
    o.union_prefilter_cache_bytes = prefilter_bytes;
    o.max_union_cached_prefilters_consumed = prefilter_bytes > 0 ? (std::uint64_t{1} << 40) : 0;
    return o;
}

Run index_run(const Run& reference, const std::vector<panvar::ReadPair>& fragments,
              const std::vector<std::vector<std::uint32_t>>& segments_by_slot, std::size_t threads,
              panvar::ContentBatchReport* report_out = nullptr, std::size_t prefilter_bytes = 0, std::size_t origin_bytes = 0) {
    panvar::ContentPlacementIndex index(100, library().max_gap_bp);
    for (std::size_t h = 0; h < reference.interned.sequences.size(); ++h) {
        const std::size_t slot = index.add_haplotype(reference.interned.sequences[h], segments_by_slot[h]);
        require(slot == h, "index slots follow registration order");
    }
    panvar::ContentIndexBatch batch(index, fragments, library(), index_options(threads, prefilter_bytes, origin_bytes));
    batch.discover();
    Run run;
    run.interned = reference.interned;
    run.candidates = reference.candidates;
    run.result = batch.score(run.candidates);
    run.ranking = panvar::certify_scored_ranking(run.result.scores, fragments.size(), options(threads).mismatch_band,
                                                 reference.result.counters, panvar::CertificationOptions{});
    if (report_out != nullptr) *report_out = batch.report();
    return run;
}

std::vector<std::vector<std::uint32_t>> segments_in_interned_order(const Fixture& f, const Run& run,
                                                                   const std::vector<std::size_t>& order) {
    std::vector<std::vector<std::uint32_t>> out(run.interned.sequences.size());
    for (std::size_t i = 0; i < order.size(); ++i) out[run.interned.index_of[i]] = f.segments[order[i]];
    return out;
}

void require_origin_totals(const Run& reference, const panvar::ContentBatchReport& r, const std::string& label) {
    require(r.physical_origins == reference.result.counters.physical_origins &&
                r.state_classes == reference.result.counters.statistical_state_classes,
            label + ": physical origins " + std::to_string(r.physical_origins) + " vs " +
                std::to_string(reference.result.counters.physical_origins) + ", state classes " +
                std::to_string(r.state_classes) + " vs " + std::to_string(reference.result.counters.statistical_state_classes));
}

// Tiny and empty blocks: A | t1 t2 t3 t4 | B where every t is 0..9 bp, so one read crosses all four. The panel rows
// are the first four; the mosaic (added after the index is prepared) takes a combination no row has.
Fixture tiny_fixture(std::vector<std::string>& mosaic, std::vector<std::uint32_t>& mosaic_segments) {
    std::mt19937_64 rng(20261016);
    const std::string A = random_bases(rng, 420), B = random_bases(rng, 430);
    const std::vector<std::vector<std::string>> tiny = {{"", random_bases(rng, 5)},
                                                       {random_bases(rng, 7), random_bases(rng, 3)},
                                                       {"", random_bases(rng, 9)},
                                                       {random_bases(rng, 4), ""}};
    const auto spell = [&](const std::vector<int>& alleles, std::vector<std::uint32_t>& lengths) {
        std::string h = A;
        lengths = {static_cast<std::uint32_t>(A.size())};
        for (std::size_t b = 0; b < tiny.size(); ++b) {
            h += tiny[b][alleles[b]];
            lengths.push_back(static_cast<std::uint32_t>(tiny[b][alleles[b]].size()));
        }
        h += B;
        lengths.push_back(static_cast<std::uint32_t>(B.size()));
        return h;
    };
    Fixture out;
    for (const std::vector<int>& row : {std::vector<int>{0, 0, 0, 0}, {1, 1, 1, 1}, {0, 1, 0, 1}, {1, 0, 1, 0}}) {
        std::vector<std::uint32_t> lengths;
        out.haplotypes.push_back(spell(row, lengths));
        out.segments.push_back(lengths);
    }
    mosaic = {spell({1, 1, 0, 0}, mosaic_segments)};
    std::mt19937_64 reads(11);
    std::size_t n = 0;
    // reads from the mosaic (its junctions are on no panel row) and from two panel rows, both orientations
    for (const std::string* source : {&mosaic[0], &out.haplotypes[1], &out.haplotypes[2]}) {
        for (std::size_t start = 0; start + 260 <= source->size(); start += 11) {
            panvar::ReadPair fr = fragment_at(*source, start, 260, 100, "t" + std::to_string(n), reads);
            if (n++ % 2 == 1) std::swap(fr.first, fr.second);
            out.fragments.push_back(fr);
        }
    }
    return out;
}

void content_index_gates(const Fixture& f, const Run& reference) {
    const std::vector<std::size_t> order = {0, 1, 2, 3, 4};
    // 1. The base fixture (shared segments at shifted coordinates, near-duplicate copies, a tandem repeat, indels,
    // empty alleles, both locus ends), at 1 and 4 threads.
    for (const std::size_t threads : {std::size_t{1}, std::size_t{4}}) {
        panvar::ContentBatchReport r;
        require_identical(reference, index_run(reference, f.fragments, segments_in_interned_order(f, reference, order), threads, &r),
                          "content index, " + std::to_string(threads) + " thread(s)", false);
        require_origin_totals(reference, r, "content index");
        if (threads == 1) {
            std::cout << "  ok: content index, base fixture -- " << r.units_evaluated << " units, " << r.unit_bases_scanned
                      << " unit bases scanned, " << r.unique_starts << " starts verified (reference "
                      << reference.result.counters.unique_starts << "), " << r.expanded_placements << " placements expanded, "
                      << r.duplicate_placements_removed << " duplicates removed, origins " << r.physical_origins << "\n";
        }
    }
    // 2. Both mate orientations: every other fragment's mates swapped.
    {
        Fixture g = f;
        for (std::size_t i = 1; i < g.fragments.size(); i += 2) std::swap(g.fragments[i].first, g.fragments[i].second);
        const Run ref = score(g, Config{});
        panvar::ContentBatchReport r;
        require_identical(ref, index_run(ref, g.fragments, segments_in_interned_order(g, ref, order), 2, &r), "swapped mates", false);
        require_origin_totals(ref, r, "swapped mates");
    }
    // 3. Any segmentation is exact: arbitrary cut points, whole haplotypes, and a sequence registered twice under
    // different segmentations (two labellings of one sequence) -- both slots must carry the reference's masses.
    {
        std::vector<std::vector<std::uint32_t>> wrong, whole;
        for (const std::string& h : reference.interned.sequences) {
            const std::uint32_t L = static_cast<std::uint32_t>(h.size());
            wrong.push_back({L / 3, 0, L / 3, L - 2 * (L / 3)});
            whole.push_back({L});
        }
        require_identical(reference, index_run(reference, f.fragments, wrong, 2), "arbitrary segmentation", false);
        require_identical(reference, index_run(reference, f.fragments, whole, 2), "whole-haplotype segments", false);
        panvar::ContentPlacementIndex index(100, library().max_gap_bp);
        const auto by_order = segments_in_interned_order(f, reference, order);
        for (std::size_t h = 0; h < reference.interned.sequences.size(); ++h) index.add_haplotype(reference.interned.sequences[h], by_order[h]);
        const std::size_t twin = index.add_haplotype(reference.interned.sequences[0], wrong[0]);
        panvar::ContentIndexBatch batch(index, f.fragments, library(), options(2));
        batch.discover();
        const panvar::FragmentScoringResult got = batch.score({{0, twin}});
        for (std::size_t frag = 0; frag < f.fragments.size(); ++frag) {
            require(same_double(got.mass_of(frag, twin), reference.result.mass_of(frag, 0)) &&
                        same_double(got.mass_of(frag, 0), reference.result.mass_of(frag, 0)),
                    "a sequence under two segmentations must carry one set of masses");
        }
        std::cout << "  ok: swapped mates; arbitrary, whole-haplotype and twin segmentations of one sequence\n";
    }
    // 4. Forced hash collisions: every unit key hashes alike; byte comparison alone separates them.
    {
        panvar::set_content_index_hash_collapse(true);
        panvar::ContentBatchReport r;
        const Run collided = index_run(reference, f.fragments, segments_in_interned_order(f, reference, order), 2, &r);
        panvar::set_content_index_hash_collapse(false);
        require_identical(reference, collided, "content index, forced hash collisions", false);
        std::cout << "  ok: forced unit-hash collisions\n";
    }
    // 5. Tiny and empty blocks, and a novel mosaic added after the index was prepared (incremental discovery).
    {
        std::vector<std::string> mosaic;
        std::vector<std::uint32_t> mosaic_segments;
        Fixture t = tiny_fixture(mosaic, mosaic_segments);
        Fixture full = t;
        full.haplotypes.push_back(mosaic[0]);
        full.segments.push_back(mosaic_segments);
        Config five;
        five.order = {0, 1, 2, 3, 4};
        const Run ref = score(full, five);
        // the index: the four panel rows first, discovered; then the mosaic, discovered incrementally
        panvar::ContentPlacementIndex index(100, library().max_gap_bp);
        std::vector<std::size_t> slot_of_input(5);
        for (std::size_t i = 0; i < 4; ++i) slot_of_input[i] = index.add_haplotype(full.haplotypes[i], full.segments[i]);
        panvar::ContentIndexBatch batch(index, full.fragments, library(), options(2));
        batch.discover();
        const std::uint64_t before_units = batch.report().units_evaluated;
        slot_of_input[4] = index.add_haplotype(full.haplotypes[4], full.segments[4]);
        batch.discover();
        const std::uint64_t new_units = batch.report().units_evaluated - before_units;
        require(new_units > 0 && batch.report().discovery_calls == 2, "the mosaic must bring new units, discovered incrementally");
        // candidates in the reference's unique-sequence space, mapped to index slots
        std::vector<panvar::HaplotypeCandidate> cands;
        std::vector<std::size_t> slot_of_unique(ref.interned.sequences.size());
        for (std::size_t i = 0; i < 5; ++i) slot_of_unique[ref.interned.index_of[i]] = slot_of_input[i];
        for (const auto& c : ref.candidates) cands.push_back({slot_of_unique[c.first], slot_of_unique[c.second]});
        const panvar::FragmentScoringResult got = batch.score(cands);
        for (std::size_t u = 0; u < ref.interned.sequences.size(); ++u) {
            require(same_double(got.exposure[slot_of_unique[u]], ref.result.exposure[u]), "tiny blocks: exposure differs");
            for (std::size_t frag = 0; frag < full.fragments.size(); ++frag) {
                require(same_double(got.mass_of(frag, slot_of_unique[u]), ref.result.mass_of(frag, u)),
                        "tiny blocks / incremental mosaic: mass differs at fragment " + std::to_string(frag));
            }
        }
        for (std::size_t c = 0; c < cands.size(); ++c) {
            require(same_double(got.scores[c].log_likelihood, ref.result.scores[c].log_likelihood),
                    "tiny blocks / incremental mosaic: candidate score differs");
        }
        require(batch.report().physical_origins == ref.result.counters.physical_origins, "tiny blocks: origin total differs");
        // the mosaic's junction reads really cross several tiny blocks and have mass only on haplotypes carrying them
        std::size_t mosaic_only = 0;
        const std::size_t m = ref.interned.index_of[4];
        for (std::size_t frag = 0; frag < full.fragments.size(); ++frag) {
            bool elsewhere = false;
            for (std::size_t u = 0; u < ref.interned.sequences.size(); ++u) {
                if (u != m && ref.result.mass_of(frag, u) >= ref.result.mass_of(frag, m)) elsewhere = true;
            }
            mosaic_only += !elsewhere && !std::isinf(ref.result.mass_of(frag, m));
        }
        require(mosaic_only > 0, "control: some fragments must fit the mosaic strictly best");
        std::cout << "  ok: tiny/empty blocks crossed by one read, novel mosaic added after preparation (" << new_units
                  << " new units discovered incrementally; " << mosaic_only << " fragments fit the mosaic strictly best)\n";
    }
    // 6. One prepared index, two read batches (even and odd fragments), each against its own reference.
    {
        panvar::ContentPlacementIndex index(100, library().max_gap_bp);
        const auto by_order = segments_in_interned_order(f, reference, order);
        for (std::size_t h = 0; h < reference.interned.sequences.size(); ++h) index.add_haplotype(reference.interned.sequences[h], by_order[h]);
        const std::size_t units = index.unit_count();
        for (int parity = 0; parity < 2; ++parity) {
            Fixture half = f;
            half.fragments.clear();
            for (std::size_t i = parity; i < f.fragments.size(); i += 2) half.fragments.push_back(f.fragments[i]);
            const Run ref = score(half, Config{});
            panvar::ContentIndexBatch batch(index, half.fragments, library(), options(2));
            batch.discover();
            Run run;
            run.interned = ref.interned;
            run.candidates = ref.candidates;
            run.result = batch.score(run.candidates);
            run.ranking = panvar::certify_scored_ranking(run.result.scores, half.fragments.size(), options(2).mismatch_band,
                                                         ref.result.counters, panvar::CertificationOptions{});
            require_identical(ref, run, "read batch " + std::to_string(parity), false);
        }
        require(index.unit_count() == units, "scoring a batch must not change the index");
        std::cout << "  ok: one prepared index reused for two read batches\n";
    }
}


// ---- I1.1: indexed discovery with the exact prefilter-decision cache ------------------------------------------------
// Each case: the cached reference, I1 with the cache off and I1 with it on, all bit-identical (masses, exposure,
// candidate components, certified ranking, winner set; origin and state totals). Executed-work counters are reported
// apart: the cache may only lower executed prefilters.
void content_index_prefilter_gates(const Fixture& f, const Run& reference) {
    const std::vector<std::size_t> order = {0, 1, 2, 3, 4};
    const auto segs = segments_in_interned_order(f, reference, order);
    const std::size_t on = std::size_t{1} << 27;
    for (const std::size_t threads : {std::size_t{1}, std::size_t{4}}) {
        panvar::ContentBatchReport off_r, on_r;
        require_identical(reference, index_run(reference, f.fragments, segs, threads, &off_r, 0), "I1 cache off", false);
        require_identical(reference, index_run(reference, f.fragments, segs, threads, &on_r, on), "I1.1 cache on", false);
        require_origin_totals(reference, on_r, "I1.1 cache on");
        require(off_r.prefilter_status == "off" && on_r.prefilter_status == "on" && on_r.prefilter.hits > 0,
                "the cache must be off, then on and serving");
        require(on_r.executed_prefilters + on_r.prefilter.hits == off_r.executed_prefilters &&
                    on_r.unique_starts == off_r.unique_starts && on_r.dp_evaluations == off_r.dp_evaluations,
                "executed prefilters + decisions served must equal the kernels run without the cache; logical counters fixed");
        if (threads == 1) {
            std::cout << "  ok: I1.1 base fixture -- prefilters executed " << on_r.executed_prefilters << " (+" << on_r.prefilter.hits
                      << " served) vs " << off_r.executed_prefilters << " without; " << on_r.prefilter.window_contents
                      << " window contents indexed over " << on_r.prefilter.positions_indexed << " positions\n";
        }
    }
    // swapped mates and an arbitrary segmentation, cache on
    {
        Fixture g = f;
        for (std::size_t i = 1; i < g.fragments.size(); i += 2) std::swap(g.fragments[i].first, g.fragments[i].second);
        const Run ref = score(g, Config{});
        require_identical(ref, index_run(ref, g.fragments, segments_in_interned_order(g, ref, order), 2, nullptr, on), "I1.1 swapped mates", false);
        std::vector<std::vector<std::uint32_t>> wrong;
        for (const std::string& h : reference.interned.sequences) {
            const std::uint32_t L = static_cast<std::uint32_t>(h.size());
            wrong.push_back({L / 3, 0, L / 3, L - 2 * (L / 3)});
        }
        require_identical(reference, index_run(reference, f.fragments, wrong, 2, nullptr, on), "I1.1 arbitrary segmentation", false);
    }
    // forced collisions in both the window index and the unit index; a tiny decision table
    {
        panvar::set_union_prefilter_cache_hash_collapse(true);
        panvar::set_content_index_hash_collapse(true);
        const Run collided = index_run(reference, f.fragments, segs, 2, nullptr, on);
        panvar::set_content_index_hash_collapse(false);
        panvar::set_union_prefilter_cache_hash_collapse(false);
        require_identical(reference, collided, "I1.1 forced window and unit hash collisions", false);
        panvar::ContentBatchReport tiny;
        require_identical(reference, index_run(reference, f.fragments, segs, 2, &tiny, 1), "I1.1 one-byte cache budget", false);
        require(tiny.prefilter.decisions_not_retained > 0, "a one-byte budget must decline decisions");
        std::cout << "  ok: I1.1 swapped mates, arbitrary segmentation, forced window/unit collisions, one-byte budget ("
                  << tiny.prefilter.decisions_not_retained << " decisions not retained)\n";
    }
    // tiny/empty blocks and a mosaic added after discovery: one batch, one cache across both discover() calls
    {
        std::vector<std::string> mosaic;
        std::vector<std::uint32_t> mosaic_segments;
        Fixture t = tiny_fixture(mosaic, mosaic_segments);
        t.haplotypes.push_back(mosaic[0]);
        t.segments.push_back(mosaic_segments);
        Config five;
        five.order = {0, 1, 2, 3, 4};
        const Run ref = score(t, five);
        panvar::ContentPlacementIndex index(100, library().max_gap_bp);
        std::vector<std::size_t> slot_of_input(5);
        for (std::size_t i = 0; i < 4; ++i) slot_of_input[i] = index.add_haplotype(t.haplotypes[i], t.segments[i]);
        panvar::ContentIndexBatch batch(index, t.fragments, library(), index_options(2, on));
        batch.discover();
        const std::uint64_t hits_before = batch.report().prefilter.hits;
        slot_of_input[4] = index.add_haplotype(t.haplotypes[4], t.segments[4]);
        batch.discover();
        std::vector<std::size_t> slot_of_unique(ref.interned.sequences.size());
        for (std::size_t i = 0; i < 5; ++i) slot_of_unique[ref.interned.index_of[i]] = slot_of_input[i];
        std::vector<panvar::HaplotypeCandidate> cands;
        for (const auto& c : ref.candidates) cands.push_back({slot_of_unique[c.first], slot_of_unique[c.second]});
        const panvar::FragmentScoringResult got = batch.score(cands);
        for (std::size_t u = 0; u < ref.interned.sequences.size(); ++u) {
            require(same_double(got.exposure[slot_of_unique[u]], ref.result.exposure[u]), "I1.1 tiny blocks: exposure");
            for (std::size_t frag = 0; frag < t.fragments.size(); ++frag)
                require(same_double(got.mass_of(frag, slot_of_unique[u]), ref.result.mass_of(frag, u)), "I1.1 tiny blocks / incremental: mass");
        }
        for (std::size_t c = 0; c < cands.size(); ++c)
            require(same_double(got.scores[c].log_likelihood, ref.result.scores[c].log_likelihood), "I1.1 tiny blocks: score");
        std::cout << "  ok: I1.1 tiny/empty blocks and an incremental mosaic under one batch cache (decisions served: "
                  << hits_before << " before the addition, " << batch.report().prefilter.hits << " after)\n";
    }
    // two read batches over one index, each with its own cache
    {
        panvar::ContentPlacementIndex index(100, library().max_gap_bp);
        for (std::size_t h = 0; h < reference.interned.sequences.size(); ++h) index.add_haplotype(reference.interned.sequences[h], segs[h]);
        for (int parity = 0; parity < 2; ++parity) {
            Fixture half = f;
            half.fragments.clear();
            for (std::size_t i = parity; i < f.fragments.size(); i += 2) half.fragments.push_back(f.fragments[i]);
            const Run ref = score(half, Config{});
            panvar::ContentIndexBatch batch(index, half.fragments, library(), index_options(2, on));
            batch.discover();
            Run run;
            run.interned = ref.interned;
            run.candidates = ref.candidates;
            run.result = batch.score(run.candidates);
            run.ranking = panvar::certify_scored_ranking(run.result.scores, half.fragments.size(), options(2).mismatch_band,
                                                         ref.result.counters, panvar::CertificationOptions{});
            require_identical(ref, run, "I1.1 read batch " + std::to_string(parity), false);
        }
        std::cout << "  ok: I1.1 two read batches, each with its own cache\n";
    }
    // unsupported configuration: mixed read lengths -> explicit fallback to uncached discovery, still exact
    {
        Fixture m = f;
        m.fragments[3].first.bases.pop_back();
        m.fragments[3].first.qualities.pop_back();
        const Run ref = score(m, Config{});
        panvar::ContentBatchReport r;
        require_identical(ref, index_run(ref, m.fragments, segments_in_interned_order(m, ref, order), 2, &r, on), "I1.1 mixed read lengths", false);
        require(r.prefilter_status == "fallback_uncached_mixed_read_lengths" && r.prefilter.hits == 0,
                "mixed read lengths must fall back to uncached discovery, and say so");
        std::cout << "  ok: I1.1 mixed read lengths fall back to uncached discovery (" << r.prefilter_status << ")\n";
    }
}


// ---- I1.2: interior-origin reuse on the indexed route ---------------------------------------------------------------
// Each case: the reference, the index with reuse off and with reuse on (prefilter cache on in both), bit-identical
// in masses, exposure, candidate components, certified ranking and winner set; physical origins (direct + served) and
// state classes equal to the reference's.

// One segment S (160 bp) twice with a 20 bp spacer, so pairs join the two copies (forward in one, reverse in the
// other): those must stay direct; within-copy pairs are served.
Fixture cross_copy_fixture() {
    std::mt19937_64 rng(20261017);
    const std::string A = random_bases(rng, 500), S = random_bases(rng, 160), sp = random_bases(rng, 20), B = random_bases(rng, 500);
    std::string S2 = S;
    S2[80] = S2[80] == 'A' ? 'C' : 'A';
    Fixture out;
    const std::vector<std::vector<const std::string*>> layout = {{&A, &S, &sp, &S, &B}, {&A, &S, &sp, &S2, &B}, {&B, &S, &sp, &S, &A}};
    for (const auto& parts : layout) {
        std::string h;
        std::vector<std::uint32_t> lengths;
        for (const std::string* part : parts) { h += *part; lengths.push_back(static_cast<std::uint32_t>(part->size())); }
        out.haplotypes.push_back(h);
        out.segments.push_back(lengths);
    }
    std::mt19937_64 reads(13);
    std::size_t n = 0;
    for (const std::string& h : out.haplotypes) {
        for (std::size_t start = 0; start + 260 <= h.size(); start += 7) {
            panvar::ReadPair fr = fragment_at(h, start, 230 + (start % 41), 100, "c" + std::to_string(n), reads);
            if (n++ % 2 == 1) std::swap(fr.first, fr.second);
            out.fragments.push_back(fr);
        }
    }
    return out;
}

void content_index_origin_gates(const Fixture& f, const Run& reference) {
    const std::vector<std::size_t> order = {0, 1, 2, 3, 4};
    const auto segs = segments_in_interned_order(f, reference, order);
    const std::size_t pre = std::size_t{1} << 27, budget = std::size_t{256} << 20;
    for (const std::size_t threads : {std::size_t{1}, std::size_t{4}}) {
        panvar::ContentBatchReport off, on;
        require_identical(reference, index_run(reference, f.fragments, segs, threads, &off, pre, 0), "I1.2 reuse off", false);
        require_identical(reference, index_run(reference, f.fragments, segs, threads, &on, pre, budget), "I1.2 reuse on", false);
        require_origin_totals(reference, on, "I1.2 reuse on");
        require(on.origin_reuse_status == "on" && on.origin_entities_admitted > 0 && on.served_physical_origins > 0 &&
                    on.direct_physical_origins < off.physical_origins,
                "reuse must admit aggregates, serve origins and pair fewer directly");
        if (threads == 1) {
            std::cout << "  ok: I1.2 base fixture -- " << on.origin_entities_admitted << " of " << on.origin_entities_eligible
                      << " eligible entities admitted (" << on.origin_held_bytes << " bytes), origins direct "
                      << on.direct_physical_origins << " + served " << on.served_physical_origins << " = reference "
                      << reference.result.counters.physical_origins << " (reuse off: all " << off.physical_origins << " direct)\n";
        }
    }
    // cross-copy pairs: one segment twice within insert range
    {
        const Fixture c = cross_copy_fixture();
        Config three;
        three.order = {0, 1, 2};
        const Run ref = score(c, three);
        const auto cs = segments_in_interned_order(c, ref, three.order);
        panvar::ContentBatchReport on;
        require_identical(ref, index_run(ref, c.fragments, cs, 2, nullptr, pre, 0), "I1.2 cross-copy, reuse off", false);
        require_identical(ref, index_run(ref, c.fragments, cs, 2, &on, pre, budget), "I1.2 cross-copy, reuse on", false);
        require_origin_totals(ref, on, "I1.2 cross-copy");
        require(on.served_physical_origins > 0 && on.direct_physical_origins > 0, "cross-copy: both served and direct origins");
        std::cout << "  ok: I1.2 repeated copies within insert range (direct " << on.direct_physical_origins << ", served "
                  << on.served_physical_origins << ")\n";
    }
    // alternative segmentation; forced collisions; budgets that decline everything and that decline some
    {
        std::vector<std::vector<std::uint32_t>> wrong;
        for (const std::string& h : reference.interned.sequences) {
            const std::uint32_t L = static_cast<std::uint32_t>(h.size());
            wrong.push_back({L / 3, 0, L / 3, L - 2 * (L / 3)});
        }
        require_identical(reference, index_run(reference, f.fragments, wrong, 2, nullptr, pre, budget), "I1.2 arbitrary segmentation", false);
        panvar::set_union_prefilter_cache_hash_collapse(true);
        panvar::set_content_index_hash_collapse(true);
        const Run collided = index_run(reference, f.fragments, segs, 2, nullptr, pre, budget);
        panvar::set_content_index_hash_collapse(false);
        panvar::set_union_prefilter_cache_hash_collapse(false);
        require_identical(reference, collided, "I1.2 forced collisions", false);
        panvar::ContentBatchReport none, some, full;
        require_identical(reference, index_run(reference, f.fragments, segs, 2, &full, pre, budget), "I1.2 full budget", false);
        require_identical(reference, index_run(reference, f.fragments, segs, 2, &none, pre, 1), "I1.2 one-byte budget", false);
        require(none.origin_entities_admitted == 0 && none.origin_entities_declined > 0 && none.served_physical_origins == 0,
                "a one-byte budget must decline every aggregate and pair directly");
        const std::size_t partial = full.origin_held_bytes / 2;
        require_identical(reference, index_run(reference, f.fragments, segs, 2, &some, pre, partial), "I1.2 half budget", false);
        require_origin_totals(reference, some, "I1.2 half budget");
        require(some.origin_entities_admitted > 0 && some.origin_entities_declined > 0 && some.origin_held_bytes <= partial,
                "half the budget must admit some, decline some and hold its bytes");
        std::cout << "  ok: I1.2 arbitrary segmentation, forced collisions, one-byte budget (all " << none.origin_entities_declined
                  << " declined, direct), half budget (" << some.origin_entities_admitted << " admitted, "
                  << some.origin_entities_declined << " declined)\n";
    }
    // tiny/empty blocks and a mosaic added after a first score: aggregates built again for entities it makes eligible
    {
        std::vector<std::string> mosaic;
        std::vector<std::uint32_t> mosaic_segments;
        Fixture t = tiny_fixture(mosaic, mosaic_segments);
        t.haplotypes.push_back(mosaic[0]);
        t.segments.push_back(mosaic_segments);
        Config five;
        five.order = {0, 1, 2, 3, 4};
        const Run ref = score(t, five);
        panvar::ContentPlacementIndex index(100, library().max_gap_bp);
        std::vector<std::size_t> slot_of_input(5);
        for (std::size_t i = 0; i < 4; ++i) slot_of_input[i] = index.add_haplotype(t.haplotypes[i], t.segments[i]);
        panvar::ContentIndexBatch batch(index, t.fragments, library(), index_options(2, pre, budget));
        batch.discover();
        (void)batch.score({{0, 0}});   // aggregates built before the addition
        const std::uint64_t built_before = batch.report().origin_entities_built;
        slot_of_input[4] = index.add_haplotype(t.haplotypes[4], t.segments[4]);
        batch.discover();
        std::vector<std::size_t> slot_of_unique(ref.interned.sequences.size());
        for (std::size_t i = 0; i < 5; ++i) slot_of_unique[ref.interned.index_of[i]] = slot_of_input[i];
        std::vector<panvar::HaplotypeCandidate> cands;
        for (const auto& c : ref.candidates) cands.push_back({slot_of_unique[c.first], slot_of_unique[c.second]});
        const panvar::FragmentScoringResult got = batch.score(cands);
        for (std::size_t u = 0; u < ref.interned.sequences.size(); ++u)
            for (std::size_t frag = 0; frag < t.fragments.size(); ++frag)
                require(same_double(got.mass_of(frag, slot_of_unique[u]), ref.result.mass_of(frag, u)), "I1.2 tiny blocks / incremental: mass");
        for (std::size_t c = 0; c < cands.size(); ++c)
            require(same_double(got.scores[c].log_likelihood, ref.result.scores[c].log_likelihood), "I1.2 tiny blocks: score");
        require(batch.report().physical_origins == ref.result.counters.physical_origins, "I1.2 tiny blocks: origin total");
        std::cout << "  ok: I1.2 tiny/empty blocks and an incremental mosaic (" << built_before << " aggregates before the addition, "
                  << batch.report().origin_entities_built << " after; " << batch.report().served_physical_origins << " origins served)\n";
    }
    // two read batches, each with its own aggregates
    {
        panvar::ContentPlacementIndex index(100, library().max_gap_bp);
        for (std::size_t h = 0; h < reference.interned.sequences.size(); ++h) index.add_haplotype(reference.interned.sequences[h], segs[h]);
        for (int parity = 0; parity < 2; ++parity) {
            Fixture half = f;
            half.fragments.clear();
            for (std::size_t i = parity; i < f.fragments.size(); i += 2) half.fragments.push_back(f.fragments[i]);
            const Run ref = score(half, Config{});
            panvar::ContentIndexBatch batch(index, half.fragments, library(), index_options(2, pre, budget));
            batch.discover();
            Run run;
            run.interned = ref.interned;
            run.candidates = ref.candidates;
            run.result = batch.score(run.candidates);
            run.ranking = panvar::certify_scored_ranking(run.result.scores, half.fragments.size(), options(2).mismatch_band,
                                                         ref.result.counters, panvar::CertificationOptions{});
            require_identical(ref, run, "I1.2 read batch " + std::to_string(parity), false);
            require(batch.report().physical_origins == ref.result.counters.physical_origins, "I1.2 read batch: origin total");
        }
        std::cout << "  ok: I1.2 two read batches, each with its own aggregates\n";
    }
}


// ---- I1.3: retained (incremental) indexed scoring ----------------------------------------------------------------------
// Every stage is compared with a fresh reference: masses and exposure of every registered slot, every likelihood
// component of every requested pair, bit for bit; counters show completed columns and pairs are not recomputed while
// genuinely new ones are.

const panvar::DiplotypeScore& reference_pair(const Run& ref, std::size_t a, std::size_t b) {
    for (std::size_t c = 0; c < ref.candidates.size(); ++c) {
        if (std::minmax(ref.candidates[c].first, ref.candidates[c].second) == std::minmax(a, b)) return ref.result.scores[c];
    }
    throw std::runtime_error("reference has no such pair");
}

void require_same_score(const panvar::DiplotypeScore& x, const panvar::DiplotypeScore& y, const std::string& label) {
    require(same_double(x.log_likelihood, y.log_likelihood) && same_double(x.fragment_term, y.fragment_term) &&
                same_double(x.exposure_term, y.exposure_term) && same_double(x.absolute_term_sum, y.absolute_term_sum) &&
                x.fragments_with_mass == y.fragments_with_mass,
            label + ": pair likelihood components differ");
}

// `unique_of[slot]` names the reference's unique sequence each slot carries; requests are in slot space.
void require_stage(const panvar::FragmentScoringResult& got, const Run& ref, const std::vector<std::size_t>& unique_of,
                   const std::vector<panvar::HaplotypeCandidate>& requests, std::size_t fragments, const std::string& label) {
    for (std::size_t slot = 0; slot < unique_of.size(); ++slot) {
        require(same_double(got.exposure[slot], ref.result.exposure[unique_of[slot]]), label + ": exposure differs");
        for (std::size_t frag = 0; frag < fragments; ++frag) {
            require(same_double(got.mass_of(frag, slot), ref.result.mass_of(frag, unique_of[slot])),
                    label + ": mass differs at slot " + std::to_string(slot) + " fragment " + std::to_string(frag));
        }
    }
    for (std::size_t c = 0; c < requests.size(); ++c) {
        require_same_score(got.scores[c], reference_pair(ref, unique_of[requests[c].first], unique_of[requests[c].second]), label);
    }
}

void content_index_incremental_gates(const Fixture& f, const Run& reference) {
    const std::vector<std::size_t> order = {0, 1, 2, 3, 4};
    const auto segs = segments_in_interned_order(f, reference, order);
    const std::size_t pre = std::size_t{1} << 27, budget = std::size_t{256} << 20;
    const std::size_t nf = f.fragments.size();
    for (const std::size_t threads : {std::size_t{1}, std::size_t{4}}) {
        const std::string T = " (" + std::to_string(threads) + " thread(s))";
        panvar::ContentPlacementIndex index(100, library().max_gap_bp);
        panvar::ContentIndexBatch batch(index, f.fragments, library(), index_options(threads, pre, budget));
        std::vector<std::size_t> unique_of;
        const auto add = [&](std::size_t u, const std::vector<std::uint32_t>& seg) {
            const std::size_t slot = index.add_haplotype(reference.interned.sequences[u], seg);
            require(slot == unique_of.size(), "slots are appended in order");
            unique_of.push_back(u);
            return slot;
        };
        // 1. three slots, three requested pairs: all computed
        for (std::size_t u = 0; u < 3; ++u) add(u, segs[u]);
        batch.discover();
        const std::vector<panvar::HaplotypeCandidate> first = {{0, 0}, {0, 1}, {1, 2}};
        require_stage(batch.score(first), reference, unique_of, first, nf, "I1.3 first request" + T);
        require(batch.report().columns_computed_last == 3 && batch.report().pairs_computed_last == 3, "first request computes 3 columns, 3 pairs");
        // 2. the identical request: nothing recomputed
        require_stage(batch.score(first), reference, unique_of, first, nf, "I1.3 repeated request" + T);
        require(batch.report().columns_computed_last == 0 && batch.report().pairs_computed_last == 0 &&
                    batch.report().columns_retained_last == 3,
                "a repeated request recomputes no column and no pair");
        // 3. reordered subset, reversed orientation: served, identical
        const std::vector<panvar::HaplotypeCandidate> reordered = {{2, 1}, {0, 0}};
        require_stage(batch.score(reordered), reference, unique_of, reordered, nf, "I1.3 reordered subset" + T);
        require(batch.report().columns_computed_last == 0 && batch.report().pairs_computed_last == 0, "a reordered subset recomputes nothing");
        // 4. new pairs of old slots: combination only
        const std::vector<panvar::HaplotypeCandidate> new_pairs = {{0, 2}, {2, 2}, {1, 0}};
        require_stage(batch.score(new_pairs), reference, unique_of, new_pairs, nf, "I1.3 new pairs of old slots" + T);
        require(batch.report().columns_computed_last == 0 && batch.report().pairs_computed_last == 2,
                "new pairs of old slots: no column, exactly the two new pairs");
        // 5. an appended slot with no new content: the same sequence again, same segmentation (an equivalent label)
        const std::uint64_t units_before = batch.report().units_evaluated;
        const std::size_t twin = add(1, segs[1]);
        batch.discover();
        require(batch.report().units_evaluated == units_before, "a repeated sequence brings no new unit");
        const std::vector<panvar::HaplotypeCandidate> twin_pairs = {{twin, twin}, {1, twin}, {0, 1}};
        require_stage(batch.score(twin_pairs), reference, unique_of, twin_pairs, nf, "I1.3 appended slot, no new content" + T);
        require(batch.report().columns_computed_last == 1 && batch.report().pairs_computed_last == 2,
                "the appended slot is computed once; (0,1) is served");
        // 6. appended slots with genuine new content, then every pair: only the new columns and pairs computed
        add(3, segs[3]);
        add(4, segs[4]);
        batch.discover();
        require(batch.report().units_evaluated > units_before, "new sequences bring new units");
        std::vector<std::size_t> slot_of_unique(5);
        for (std::size_t u = 0; u < 5; ++u) slot_of_unique[u] = u < 3 ? u : u + 1;   // slot 3 is the twin
        std::vector<panvar::HaplotypeCandidate> all;
        for (const auto& c : reference.candidates) all.push_back({slot_of_unique[c.first], slot_of_unique[c.second]});
        const panvar::FragmentScoringResult full = batch.score(all);
        require_stage(full, reference, unique_of, all, nf, "I1.3 appended new content, every pair" + T);
        require(batch.report().columns_computed_last == 2 && batch.report().pairs_computed_last == 10 &&
                    batch.report().columns_computed_total == 6,
                "only the 2 new columns and the 10 pairs not yet seen are computed (6 columns ever)");
        // the final incremental state equals bulk scoring of the same slots from scratch
        {
            panvar::ContentPlacementIndex bulk_index(100, library().max_gap_bp);
            for (std::size_t slot = 0; slot < unique_of.size(); ++slot) bulk_index.add_haplotype(reference.interned.sequences[unique_of[slot]], segs[unique_of[slot]]);
            panvar::ContentIndexBatch bulk(bulk_index, f.fragments, library(), index_options(threads, pre, budget));
            bulk.discover();
            const panvar::FragmentScoringResult b = bulk.score(all);
            require(b.log_mass.size() == full.log_mass.size(), "bulk dimensions");
            for (std::size_t i = 0; i < b.log_mass.size(); ++i) require(same_double(b.log_mass[i], full.log_mass[i]), "I1.3 incremental != bulk mass");
            for (std::size_t c = 0; c < all.size(); ++c) require_same_score(b.scores[c], full.scores[c], "I1.3 incremental != bulk");
            require(bulk.report().physical_origins == batch.report().physical_origins, "I1.3 incremental origin total != bulk");
        }
        // 7. priors stay outside the retained likelihood: two prior vectors, each reproducing the reference's winner
        for (int which = 0; which < 2; ++which) {
            std::vector<double> prior(all.size());
            for (std::size_t c = 0; c < all.size(); ++c) prior[c] = which == 0 ? 0.0 : -37.0 * static_cast<double>(c % 4);
            std::size_t best_got = 0, best_ref = 0;
            for (std::size_t c = 0; c < all.size(); ++c) {
                if (full.scores[c].log_likelihood + prior[c] > full.scores[best_got].log_likelihood + prior[best_got]) best_got = c;
                if (reference.result.scores[c].log_likelihood + prior[c] > reference.result.scores[best_ref].log_likelihood + prior[best_ref]) best_ref = c;
            }
            require(best_got == best_ref, "I1.3 winner under altered priors");
        }
        // 8. separate contexts: another library, another read set -- their own results; this batch's are untouched
        {
            panvar::LibraryModel other = library();
            other.substitution_error = 0.03;
            const panvar::FragmentScoringResult ref_other = panvar::score_diplotypes(reference.interned, reference.candidates, f.fragments, other, options(threads));
            panvar::ContentIndexBatch b2(index, f.fragments, other, index_options(threads, pre, budget));
            b2.discover();
            const panvar::FragmentScoringResult got_other = b2.score(all);
            bool differs = false;
            for (std::size_t slot = 0; slot < unique_of.size(); ++slot) {
                for (std::size_t frag = 0; frag < nf; ++frag) {
                    require(same_double(got_other.mass_of(frag, slot), ref_other.mass_of(frag, unique_of[slot])), "I1.3 other library: mass");
                    differs |= !same_double(got_other.mass_of(frag, slot), full.mass_of(frag, slot));
                }
            }
            require(differs, "control: another library must give other masses");
            Fixture half = f;
            half.fragments.clear();
            for (std::size_t i = 0; i < nf; i += 2) half.fragments.push_back(f.fragments[i]);
            const Run ref_half = score(half, Config{});
            panvar::ContentIndexBatch b3(index, half.fragments, library(), index_options(threads, pre, budget));
            b3.discover();
            require_stage(b3.score(all), ref_half, unique_of, all, half.fragments.size(), "I1.3 other read set" + T);
            require_stage(batch.score(all), reference, unique_of, all, nf, "I1.3 after other contexts" + T);
            require(batch.report().columns_computed_last == 0 && batch.report().pairs_computed_last == 0,
                    "other contexts must not disturb this one's retained results");
        }
        // 9. tiny/empty blocks: a junction mosaic appended after its panel rows were scored
        {
            std::vector<std::string> mosaic;
            std::vector<std::uint32_t> mosaic_segments;
            Fixture t = tiny_fixture(mosaic, mosaic_segments);
            t.haplotypes.push_back(mosaic[0]);
            t.segments.push_back(mosaic_segments);
            Config five;
            five.order = {0, 1, 2, 3, 4};
            const Run ref = score(t, five);
            panvar::ContentPlacementIndex ti(100, library().max_gap_bp);
            panvar::ContentIndexBatch tb(ti, t.fragments, library(), index_options(threads, pre, budget));
            std::vector<std::size_t> tu;
            for (std::size_t i = 0; i < 4; ++i) { ti.add_haplotype(t.haplotypes[i], t.segments[i]); tu.push_back(ref.interned.index_of[i]); }
            tb.discover();
            std::vector<panvar::HaplotypeCandidate> old_pairs;
            for (std::size_t a = 0; a < 4; ++a) for (std::size_t b = a; b < 4; ++b) old_pairs.push_back({a, b});
            require_stage(tb.score(old_pairs), ref, tu, old_pairs, t.fragments.size(), "I1.3 tiny blocks, panel rows" + T);
            const std::uint64_t units = tb.report().units_evaluated;
            ti.add_haplotype(t.haplotypes[4], t.segments[4]);
            tu.push_back(ref.interned.index_of[4]);
            tb.discover();
            require(tb.report().units_evaluated > units, "the junction mosaic brings new units");
            std::vector<panvar::HaplotypeCandidate> every;
            for (std::size_t a = 0; a < 5; ++a) for (std::size_t b = a; b < 5; ++b) every.push_back({a, b});
            require_stage(tb.score(every), ref, tu, every, t.fragments.size(), "I1.3 tiny blocks, mosaic appended" + T);
            require(tb.report().columns_computed_last == 1 && tb.report().pairs_computed_last == 5,
                    "the mosaic: one new column and its five new pairs");
        }
        if (threads == 1) {
            std::cout << "  ok: I1.3 retained scoring -- repeated, reordered and subset requests recompute nothing; new pairs of "
                         "old slots combine only; an appended repeat brings no unit and one column; new content brings "
                         "units and only its columns and pairs; incremental = bulk; priors outside; other library and "
                         "read contexts independent; a tiny-block junction mosaic appended (1 column, 5 pairs)\n";
        }
    }
    std::cout << "  ok: I1.3 at 1 and 4 threads\n";
}


// ---- I1.4: fine (content-defined) discovery units, coarse rows reconstructed ---------------------------------------
// Coarse and fine discovery must give the same coarse interior placement rows (every tuple, emission bits included),
// and the fine route must equal the reference (masses, exposure, candidate components, ranking, winner set, origin
// totals). Origin reuse, prefilter cache and retention stay on throughout.

void require_same_rows(const panvar::ContentIndexBatch& a, const panvar::ContentIndexBatch& b, const std::string& label) {
    const auto ra = a.interior_placement_rows(), rb = b.interior_placement_rows();
    require(ra.size() == rb.size(), label + ": coarse interior unit sets differ");
    for (const auto& [key, rows] : ra) {
        const auto found = rb.find(key);
        require(found != rb.end() && found->second.size() == rows.size(), label + ": reconstructed row count differs");
        for (std::size_t i = 0; i < rows.size(); ++i) {
            const auto& x = rows[i];
            const auto& y = found->second[i];
            require(std::get<0>(x) == std::get<0>(y) && std::get<1>(x) == std::get<1>(y) && std::get<2>(x) == std::get<2>(y) &&
                        std::get<3>(x) == std::get<3>(y) && same_double(std::get<4>(x), std::get<4>(y)),
                    label + ": reconstructed placement differs");
        }
    }
}

// Coarse and fine batches over the same sequences; returns the fine run (for require_identical) and checks the rows.
Run fine_and_coarse(const Run& ref, const std::vector<panvar::ReadPair>& fragments, const std::vector<std::vector<std::uint32_t>>& segs,
                    std::size_t threads, const std::string& label, panvar::ContentIndexReport* index_report = nullptr,
                    panvar::ContentBatchReport* batch_report = nullptr) {
    const std::size_t pre = std::size_t{1} << 27, budget = std::size_t{256} << 20;
    panvar::ContentPlacementIndex coarse(100, library().max_gap_bp, false), fine(100, library().max_gap_bp, true);
    for (std::size_t h = 0; h < ref.interned.sequences.size(); ++h) {
        coarse.add_haplotype(ref.interned.sequences[h], segs[h]);
        fine.add_haplotype(ref.interned.sequences[h], segs[h]);
    }
    panvar::ContentIndexBatch bc(coarse, fragments, library(), index_options(threads, pre, budget));
    panvar::ContentIndexBatch bf(fine, fragments, library(), index_options(threads, pre, budget));
    bc.discover();
    bf.discover();
    require_same_rows(bc, bf, label);
    Run run;
    run.interned = ref.interned;
    run.candidates = ref.candidates;
    run.result = bf.score(run.candidates);
    run.ranking = panvar::certify_scored_ranking(run.result.scores, fragments.size(), options(threads).mismatch_band,
                                                 ref.result.counters, panvar::CertificationOptions{});
    require(bf.report().physical_origins == ref.result.counters.physical_origins, label + ": origin total differs");
    if (index_report != nullptr) *index_report = fine.report();
    if (batch_report != nullptr) *batch_report = bf.report();
    return run;
}

void content_index_fine_gates(const Fixture& f, const Run& reference) {
    const std::vector<std::size_t> order = {0, 1, 2, 3, 4};
    const auto segs = segments_in_interned_order(f, reference, order);
    // base fixture: near-identical alleles (a substitution, a 37 bp insertion, a 50 bp deletion), copies, both ends
    for (const std::size_t threads : {std::size_t{1}, std::size_t{4}}) {
        panvar::ContentIndexReport ir;
        panvar::ContentBatchReport br;
        require_identical(reference, fine_and_coarse(reference, f.fragments, segs, threads, "I1.4 base", &ir, &br), "I1.4 base fixture", false);
        require(ir.fine_occurrences > ir.interior_units && ir.fine_occurrences > ir.fine_units,
                "the fixture must cut segments into several cores and share fine units across alleles");
        require(br.owned_starts_skipped_before_verification > 0 && br.coarse_units_reconstructed > 0,
                "ownership must drop starts before verification; coarse rows must be reconstructed");
        if (threads == 1) {
            std::cout << "  ok: I1.4 base fixture -- " << ir.interior_units << " coarse interiors (" << ir.coarse_interior_bytes
                      << " bytes) as " << ir.fine_occurrences << " fine pieces over " << ir.fine_units << " distinct fine units ("
                      << ir.fine_unit_bytes << " text bytes, " << ir.fine_core_bytes << " owned); " << br.owned_starts_skipped_before_verification
                      << " starts dropped before verification\n";
        }
    }
    // both orientations; mixed read lengths (which also exercises the prefilter fallback); arbitrary segmentation
    {
        Fixture g = f;
        for (std::size_t i = 1; i < g.fragments.size(); i += 2) std::swap(g.fragments[i].first, g.fragments[i].second);
        for (std::size_t i = 0; i < g.fragments.size(); i += 7) {
            g.fragments[i].second.bases.resize(93);
            g.fragments[i].second.qualities.resize(93);
        }
        const Run ref = score(g, Config{});
        require_identical(ref, fine_and_coarse(ref, g.fragments, segments_in_interned_order(g, ref, order), 2, "I1.4 swapped, mixed lengths"),
                          "I1.4 swapped mates and mixed read lengths", false);
        std::vector<std::vector<std::uint32_t>> wrong;
        for (const std::string& h : reference.interned.sequences) {
            const std::uint32_t L = static_cast<std::uint32_t>(h.size());
            wrong.push_back({L / 3, 0, L / 3, L - 2 * (L / 3)});
        }
        require_identical(reference, fine_and_coarse(reference, f.fragments, wrong, 2, "I1.4 arbitrary segmentation"), "I1.4 arbitrary segmentation", false);
        panvar::set_union_prefilter_cache_hash_collapse(true);
        panvar::set_content_index_hash_collapse(true);
        const Run collided = fine_and_coarse(reference, f.fragments, segs, 2, "I1.4 collisions");
        panvar::set_content_index_hash_collapse(false);
        panvar::set_union_prefilter_cache_hash_collapse(false);
        require_identical(reference, collided, "I1.4 forced collisions", false);
        std::cout << "  ok: I1.4 swapped mates with mixed read lengths, arbitrary segmentation, forced collisions\n";
    }
    // tiny/empty blocks and a junction mosaic appended, with repeated requests (retained scoring) under fine discovery
    {
        std::vector<std::string> mosaic;
        std::vector<std::uint32_t> mosaic_segments;
        Fixture t = tiny_fixture(mosaic, mosaic_segments);
        t.haplotypes.push_back(mosaic[0]);
        t.segments.push_back(mosaic_segments);
        Config five;
        five.order = {0, 1, 2, 3, 4};
        const Run ref = score(t, five);
        panvar::ContentPlacementIndex ti(100, library().max_gap_bp, true);
        panvar::ContentIndexBatch tb(ti, t.fragments, library(), index_options(4, std::size_t{1} << 27, std::size_t{256} << 20));
        std::vector<std::size_t> tu;
        for (std::size_t i = 0; i < 4; ++i) { ti.add_haplotype(t.haplotypes[i], t.segments[i]); tu.push_back(ref.interned.index_of[i]); }
        tb.discover();
        std::vector<panvar::HaplotypeCandidate> old_pairs;
        for (std::size_t a = 0; a < 4; ++a) for (std::size_t b = a; b < 4; ++b) old_pairs.push_back({a, b});
        require_stage(tb.score(old_pairs), ref, tu, old_pairs, t.fragments.size(), "I1.4 tiny blocks");
        require_stage(tb.score(old_pairs), ref, tu, old_pairs, t.fragments.size(), "I1.4 tiny blocks, repeated");
        require(tb.report().columns_computed_last == 0 && tb.report().pairs_computed_last == 0, "I1.4: a repeat recomputes nothing");
        ti.add_haplotype(t.haplotypes[4], t.segments[4]);
        tu.push_back(ref.interned.index_of[4]);
        tb.discover();
        std::vector<panvar::HaplotypeCandidate> every;
        for (std::size_t a = 0; a < 5; ++a) for (std::size_t b = a; b < 5; ++b) every.push_back({b, a});   // reversed orientation
        require_stage(tb.score(every), ref, tu, every, t.fragments.size(), "I1.4 tiny blocks, mosaic appended");
        require(tb.report().columns_computed_last == 1 && tb.report().pairs_computed_last == 5, "I1.4: the mosaic adds 1 column and 5 pairs");
        std::cout << "  ok: I1.4 tiny/empty blocks, an appended junction mosaic, repeated and reversed requests\n";
    }
    // retained scoring with nonempty audit-library components: repeated and reordered requests against the reference
    {
        panvar::FragmentScoringOptions audited = index_options(2, std::size_t{1} << 27, std::size_t{256} << 20);
        panvar::LibraryModel a1 = library(), a2 = library();
        a1.fragment_start_intensity *= 2.0;
        a2.outlier_weight = 0.02;
        audited.audit_libraries = {a1, a2};
        panvar::FragmentScoringOptions plain = options(2);
        plain.audit_libraries = audited.audit_libraries;
        const panvar::FragmentScoringResult ref = panvar::score_diplotypes(reference.interned, reference.candidates, f.fragments, library(), plain);
        require(ref.audit_scores.size() == 2, "the reference must carry two audit libraries");
        panvar::ContentPlacementIndex ai(100, library().max_gap_bp, true);
        for (std::size_t h = 0; h < reference.interned.sequences.size(); ++h) ai.add_haplotype(reference.interned.sequences[h], segs[h]);
        panvar::ContentIndexBatch ab(ai, f.fragments, library(), audited);
        ab.discover();
        const auto check = [&](const std::vector<std::size_t>& rows, const std::string& label) {
            std::vector<panvar::HaplotypeCandidate> req;
            for (const std::size_t c : rows) req.push_back({reference.candidates[c].second, reference.candidates[c].first});
            const panvar::FragmentScoringResult got = ab.score(req);
            require(got.audit_scores.size() == 2, label + ": audit components missing");
            for (std::size_t i = 0; i < rows.size(); ++i) {
                require_same_score(got.scores[i], ref.scores[rows[i]], label);
                for (std::size_t k = 0; k < 2; ++k) require_same_score(got.audit_scores[k][i], ref.audit_scores[k][rows[i]], label + " audit");
            }
        };
        check({0, 3, 7}, "I1.4 audit, first");
        check({7, 0}, "I1.4 audit, reordered repeat");
        require(ab.report().pairs_computed_last == 0, "I1.4 audit: a reordered repeat recomputes no pair");
        std::vector<std::size_t> all_rows(reference.candidates.size());
        for (std::size_t c = 0; c < all_rows.size(); ++c) all_rows[c] = c;
        check(all_rows, "I1.4 audit, every pair");
        require(ab.report().pairs_computed_last == all_rows.size() - 3, "I1.4 audit: only unseen pairs computed");
        std::cout << "  ok: I1.4 retained scoring with two audit libraries -- repeated, reordered and full requests\n";
    }
}

}  // namespace


// ---- I2-A: the caller adapter (score_diplotypes_indexed_union) -----------------------------------------------------
// The one-shot indexed route behind score_diplotypes' contract: identical masses, exposures, components (audit views
// included), ranking and winners; deterministic executed work at 1 and 4 threads; every applicable guard fires exactly
// above its exact executed total and accepts at it; refusals throw with no result and report what they did.
panvar::FragmentScoringOptions adapter_options(std::size_t threads) {
    panvar::FragmentScoringOptions o = options(threads);
    o.union_origin_cache_bytes = std::size_t{256} << 20;
    o.max_union_cached_origins_consumed = std::uint64_t{1} << 40;
    o.max_union_cached_states_consumed = std::uint64_t{1} << 40;
    o.max_union_cached_placements_consumed = std::uint64_t{1} << 40;
    o.union_prefilter_cache_bytes = std::size_t{1} << 27;
    o.max_union_cached_prefilters_consumed = std::uint64_t{1} << 40;
    return o;
}
constexpr std::size_t kAdapterStore = std::size_t{1} << 30;

Run adapter_run(const Run& ref, const std::vector<panvar::ReadPair>& fragments, const std::vector<std::vector<std::uint32_t>>& segs,
                const panvar::FragmentScoringOptions& o, panvar::IndexedUnionReport* report = nullptr,
                std::size_t store = kAdapterStore, panvar::FragmentScoringObservation* observed = nullptr) {
    Run run;
    run.interned = ref.interned;
    run.candidates = ref.candidates;
    run.result = panvar::score_diplotypes_indexed_union(run.interned, run.candidates, fragments, library(), o, segs, store,
                                                        report, observed);
    run.ranking = panvar::certify_scored_ranking(run.result.scores, fragments.size(), o.mismatch_band, run.result.counters,
                                                 panvar::CertificationOptions{});
    return run;
}

template <typename F>
bool adapter_refuses(F&& call, panvar::FragmentScoringObservation* observed = nullptr) {
    try {
        call(observed);
    } catch (const panvar::ResourceExhausted&) {
        return true;
    }
    return false;
}

template <typename F>
bool adapter_invalid(F&& call) {
    try {
        call();
    } catch (const std::invalid_argument&) {
        return true;
    }
    return false;
}

void require_same_executed(const panvar::FragmentScoringCounters& a, const panvar::FragmentScoringCounters& b,
                           const std::string& label) {
    require(a.piece_hits == b.piece_hits && a.proposed_starts == b.proposed_starts && a.unique_starts == b.unique_starts &&
                a.full_read_verifications == b.full_read_verifications && a.dp_evaluations == b.dp_evaluations &&
                a.prefilter_rejections == b.prefilter_rejections && a.accepted_mate_placements == b.accepted_mate_placements &&
                a.physical_origins == b.physical_origins && a.statistical_state_classes == b.statistical_state_classes &&
                a.latent_alignment_pairs == b.latent_alignment_pairs && a.valid_fr_states == b.valid_fr_states &&
                a.reused_physical_origins == b.reused_physical_origins &&
                a.reused_state_aggregates == b.reused_state_aggregates && a.merged_state_classes == b.merged_state_classes,
            label + ": executed work differs");
}

void content_index_adapter_gates(const Fixture& f, const Run& reference) {
    const std::vector<std::size_t> order = {0, 1, 2, 3, 4};
    const auto segs = segments_in_interned_order(f, reference, order);
    // agreement and determinism: 1 and 4 threads, bit for bit, with identical executed work
    panvar::IndexedUnionReport r1, r4;
    const Run one = adapter_run(reference, f.fragments, segs, adapter_options(1), &r1);
    const Run four = adapter_run(reference, f.fragments, segs, adapter_options(4), &r4);
    require_identical(reference, one, "I2-A adapter, 1 thread", false);
    require_identical(reference, four, "I2-A adapter, 4 threads", false);
    require_identical(one, four, "I2-A adapter, 1 against 4 threads", false);
    require_same_executed(one.result.counters, four.result.counters, "I2-A adapter, 1 against 4 threads");
    const panvar::FragmentScoringCounters& c = one.result.counters;
    require(r1.batch.origin_entities_admitted > 0 && c.reused_physical_origins > 0 && r1.batch.store_charged_bytes > 0 &&
                r1.batch.placements_consumed_last > 0 && r1.batch.prefilter_status == "on" && r1.index.fine_units > 0,
            "the fixture must exercise fine discovery, origin reuse, the store and the prefilter cache");
    // the reference's logical census, recovered from executed + served (the reference's own phase-2 identity)
    require(c.physical_origins - r1.batch.origins_constructed + c.reused_physical_origins ==
                    reference.result.counters.physical_origins &&
                c.merged_state_classes == reference.result.counters.statistical_state_classes,
            "executed origins minus constructed plus served must be the reference's physical origins");
    require(one.result.threads_used_placement >= 1 && four.result.threads_used_placement == 4 &&
                four.result.threads_used_combination >= 1 && one.result.union_origin_cache.active &&
                one.result.union_prefilter_cache.active && !one.result.union_placement_cache.active,
            "the adapter must report its workers and which caches it used");
    std::cout << "  ok: I2-A adapter = reference at 1 and 4 threads; executed origins " << c.physical_origins << " + served "
              << c.reused_physical_origins << " (constructed " << r1.batch.origins_constructed << "); store "
              << r1.batch.store_charged_bytes << " bytes charged\n";

    // nonempty audit-library components
    {
        panvar::FragmentScoringOptions plain = options(2), audited = adapter_options(2);
        panvar::LibraryModel a1 = library(), a2 = library();
        a1.fragment_start_intensity *= 2.0;
        a2.outlier_weight = 0.02;
        plain.audit_libraries = audited.audit_libraries = {a1, a2};
        const panvar::FragmentScoringResult ref =
            panvar::score_diplotypes(reference.interned, reference.candidates, f.fragments, library(), plain);
        const Run got = adapter_run(reference, f.fragments, segs, audited);
        require(ref.audit_scores.size() == 2 && got.result.audit_scores.size() == 2, "both routes must carry two audit views");
        for (std::size_t k = 0; k < 2; ++k)
            for (std::size_t i = 0; i < ref.scores.size(); ++i) {
                require_same_score(got.result.scores[i], ref.scores[i], "I2-A audit run");
                require_same_score(got.result.audit_scores[k][i], ref.audit_scores[k][i], "I2-A audit view");
            }
        std::cout << "  ok: I2-A adapter with two audit libraries\n";
    }
    // both orientations and mixed read lengths: the reference's refusal under the prefilter cache, agreement without it
    {
        Fixture g = f;
        for (std::size_t i = 1; i < g.fragments.size(); i += 2) std::swap(g.fragments[i].first, g.fragments[i].second);
        for (std::size_t i = 0; i < g.fragments.size(); i += 7) {
            g.fragments[i].second.bases.resize(93);
            g.fragments[i].second.qualities.resize(93);
        }
        const Run ref = score(g, Config{});
        const auto gsegs = segments_in_interned_order(g, ref, order);
        Config prefiltered;
        prefiltered.prefilter_budget = std::size_t{1} << 27;
        prefiltered.prefilter_cap = std::uint64_t{1} << 40;
        bool reference_refused = false;
        try {
            (void)score(g, prefiltered);
        } catch (const std::invalid_argument&) {
            reference_refused = true;
        }
        require(reference_refused, "fixture: the reference refuses mixed read lengths under the prefilter cache");
        require(adapter_invalid([&] { (void)adapter_run(ref, g.fragments, gsegs, adapter_options(2)); }),
                "the adapter must refuse mixed read lengths under the prefilter cache, as the reference does");
        panvar::FragmentScoringOptions no_prefilter = adapter_options(2);
        no_prefilter.union_prefilter_cache_bytes = 0;
        no_prefilter.max_union_cached_prefilters_consumed = 0;
        require_identical(ref, adapter_run(ref, g.fragments, gsegs, no_prefilter), "I2-A swapped mates, mixed lengths", false);
        std::cout << "  ok: I2-A both orientations and mixed read lengths (refused under the prefilter cache like the reference)\n";
    }
    // tiny/empty blocks and a novel junction mosaic, scored in one shot with the panel rows
    {
        std::vector<std::string> mosaic;
        std::vector<std::uint32_t> mosaic_segments;
        Fixture t = tiny_fixture(mosaic, mosaic_segments);
        t.haplotypes.push_back(mosaic[0]);
        t.segments.push_back(mosaic_segments);
        Config five;
        five.order = {0, 1, 2, 3, 4};
        const Run ref = score(t, five);
        require_identical(ref, adapter_run(ref, t.fragments, segments_in_interned_order(t, ref, five.order), adapter_options(4)),
                          "I2-A tiny blocks and a junction mosaic", false);
        std::cout << "  ok: I2-A tiny/empty blocks and a novel junction\n";
    }
    // guards: accepted at the exact executed total, refused one above it (single worker, so every total is exact)
    {
        const panvar::FragmentScoringOptions base = adapter_options(1);
        const panvar::ContentBatchReport& b = r1.batch;
        struct Guard {
            const char* name;
            std::uint64_t total;
            std::function<void(panvar::FragmentScoringOptions&, std::uint64_t)> set;
        };
        const std::uint64_t cells = static_cast<std::uint64_t>(reference.interned.sequences.size()) * f.fragments.size();
        const std::vector<Guard> guards = {
            {"piece hits", c.piece_hits, [](panvar::FragmentScoringOptions& o, std::uint64_t v) { o.max_piece_hits = v; }},
            {"physical origins", c.physical_origins, [](panvar::FragmentScoringOptions& o, std::uint64_t v) { o.max_physical_origins = v; }},
            {"state classes", c.statistical_state_classes, [](panvar::FragmentScoringOptions& o, std::uint64_t v) { o.max_state_classes = v; }},
            {"latent pairs", c.latent_alignment_pairs, [](panvar::FragmentScoringOptions& o, std::uint64_t v) { o.max_latent_alignment_pairs = v; }},
            {"mass cells", cells, [](panvar::FragmentScoringOptions& o, std::uint64_t v) { o.max_mass_cells = v; }},
            {"cached placements", b.placements_consumed_last,
             [](panvar::FragmentScoringOptions& o, std::uint64_t v) { o.max_union_cached_placements_consumed = v; }},
            {"cached origins", b.origins_consumed_last,
             [](panvar::FragmentScoringOptions& o, std::uint64_t v) { o.max_union_cached_origins_consumed = v; }},
            {"cached states", b.states_consumed_last,
             [](panvar::FragmentScoringOptions& o, std::uint64_t v) { o.max_union_cached_states_consumed = v; }},
            {"cached prefilters", b.prefilter.hits,
             [](panvar::FragmentScoringOptions& o, std::uint64_t v) { o.max_union_cached_prefilters_consumed = v; }},
        };
        for (const Guard& g : guards) {
            require(g.total > 1, std::string("fixture: the ") + g.name + " total must be exercised");
            panvar::FragmentScoringOptions at = base, below = base;
            g.set(at, g.total);
            g.set(below, g.total - 1);
            require(!adapter_refuses([&](panvar::FragmentScoringObservation*) { (void)adapter_run(reference, f.fragments, segs, at); }),
                    std::string("I2-A: a ") + g.name + " cap equal to the executed total must accept");
            panvar::FragmentScoringObservation seen;
            panvar::IndexedUnionReport refused;
            require(adapter_refuses([&](panvar::FragmentScoringObservation* o) {
                        (void)adapter_run(reference, f.fragments, segs, below, &refused, kAdapterStore, o);
                    }, &seen),
                    std::string("I2-A: a ") + g.name + " cap one below the executed total must refuse");
            require(!refused.refusal.empty() && seen.threads_used_placement >= 1 && seen.threads_used_combination == 0 &&
                        seen.seconds_diploid_combination == 0.0,
                    std::string("I2-A: the ") + g.name + " refusal must record itself, its workers and no combination");
            if (std::string(g.name) == "physical origins") {
                require(seen.counters.physical_origins > g.total - 1, "I2-A: the observed origins must show the crossing");
            }
        }
        // aggregate construction is charged in flight: a cap below what construction alone enumerates refuses before
        // the first aggregate is admitted, rather than after construction (or only at acceptance)
        {
            panvar::FragmentScoringOptions tight = base;
            tight.max_physical_origins = 10;
            panvar::IndexedUnionReport refused;
            require(adapter_refuses([&](panvar::FragmentScoringObservation*) {
                        (void)adapter_run(reference, f.fragments, segs, tight, &refused);
                    }) && b.origins_constructed > 10 && refused.batch.origin_entities_eligible > 0 &&
                        refused.batch.origin_entities_built == 0,
                    "I2-A: aggregate construction must be charged to the physical-origin cap as it runs");
        }
        // the store budget. I2-C: it bounds the allocated capacity held at once (the ledger's peak), not a cumulative charge
        require(!adapter_refuses([&](panvar::FragmentScoringObservation*) {
                    (void)adapter_run(reference, f.fragments, segs, base, nullptr, b.store_peak_bytes);
                }),
                "I2-A/C: a store budget equal to the ledger's peak must accept");
        require(adapter_refuses([&](panvar::FragmentScoringObservation*) {
                    (void)adapter_run(reference, f.fragments, segs, base, nullptr, b.store_peak_bytes - 1);
                }),
                "I2-A/C: a store budget one byte below the ledger's peak must refuse");
        // per-text proposed starts: no exact boundary is claimed (see the header), but the guard fires
        panvar::FragmentScoringOptions starts = base;
        starts.max_proposed_starts_per_haplotype = 1;
        require(adapter_refuses([&](panvar::FragmentScoringObservation*) { (void)adapter_run(reference, f.fragments, segs, starts); }),
                "I2-A: the per-text proposed-start cap must fire");
        std::cout << "  ok: I2-A guards accept at their exact executed totals and refuse one above (" << guards.size() + 1
                  << " caps), proposed-start cap fires\n";
    }
    // declared requirements are refused up front
    {
        const panvar::FragmentScoringOptions base = adapter_options(1);
        require(adapter_invalid([&] { (void)adapter_run(reference, f.fragments, segs, base, nullptr, 0); }),
                "I2-A: a zero store budget must be refused");
        panvar::FragmentScoringOptions o = base;
        o.max_union_cached_placements_consumed = 0;
        require(adapter_invalid([&] { (void)adapter_run(reference, f.fragments, segs, o); }),
                "I2-A: a zero cached-placement cap must be refused");
        o = base;
        o.max_union_cached_states_consumed = 0;
        require(adapter_invalid([&] { (void)adapter_run(reference, f.fragments, segs, o); }),
                "I2-A: origin reuse without both consumption caps must be refused");
        o = base;
        std::vector<panvar::HaplotypeWork> work;
        o.record_haplotype_work = &work;
        require(adapter_invalid([&] { (void)adapter_run(reference, f.fragments, segs, o); }),
                "I2-A: per-haplotype diagnostics must be refused");
        std::vector<std::vector<std::uint32_t>> short_segs(segs.begin(), segs.end() - 1);
        require(adapter_invalid([&] { (void)adapter_run(reference, f.fragments, short_segs, base); }),
                "I2-A: a segmentation per haplotype is required");
        // a one-shot batch refuses a second discovery with new units
        panvar::ContentPlacementIndex index(100, library().max_gap_bp, true);
        index.add_haplotype(reference.interned.sequences[0], segs[0]);
        panvar::ContentIndexBatch batch(index, f.fragments, library(), index_options(1, 0));
        batch.set_one_shot(true);
        batch.discover();
        index.add_haplotype(reference.interned.sequences[1], segs[1]);
        bool refused = false;
        try {
            batch.discover();
        } catch (const std::logic_error&) {
            refused = true;
        }
        require(refused, "I2-A: a one-shot batch must refuse to discover new units");
        std::cout << "  ok: I2-A declared requirements refused up front; one-shot lifecycle enforced\n";
    }
}


// ---- I2-C: the placement-store ledger charges allocated capacity -----------------------------------------------------
// The ledger must equal an independent recount of the placement containers' capacity (unused capacity included), the
// budget must accept at the exact peak and refuse one byte below at 1 and 4 threads, reconstruction must overlap less
// than the I2-A layout's, refusals return nothing and leave the ledger consistent, the incremental index keeps its fine
// rows, and a negative control that charges occupied rows only must fail the capacity checks.

// ---- I4-B: factorised placement views (coarse interior units served as views over the fine units' rows) ---------------
// Against the reconstructed layout and the reference: identical masses, exposures, components, ranking and winners;
// identical executed and served work (the same units are discovered, the same rows expanded, the same aggregates built);
// the store holds the distinct evidence only (no reconstructed rows); its ledger equals the recount; its budget accepts
// at the views' own peak and refuses one byte below; the materialised view rows equal the reconstructed rows.
struct ViewPair {
    Run plain, views;
    panvar::IndexedUnionReport plain_report, views_report;
};
Run view_adapter_run(const Run& ref, const std::vector<panvar::ReadPair>& fragments, const std::vector<std::vector<std::uint32_t>>& segs,
                     const panvar::FragmentScoringOptions& o, bool views, panvar::IndexedUnionReport* report,
                     std::size_t store = kAdapterStore) {
    Run run;
    run.interned = ref.interned;
    run.candidates = ref.candidates;
    run.result = panvar::score_diplotypes_indexed_union(run.interned, run.candidates, fragments, library(), o, segs, store, report,
                                                        nullptr, false, views);
    run.ranking = panvar::certify_scored_ranking(run.result.scores, fragments.size(), o.mismatch_band, run.result.counters,
                                                 panvar::CertificationOptions{});
    return run;
}
ViewPair view_pair(const Run& ref, const std::vector<panvar::ReadPair>& fragments, const std::vector<std::vector<std::uint32_t>>& segs,
                   const panvar::FragmentScoringOptions& o, const std::string& label, bool against_reference = true) {
    ViewPair out;
    out.plain = view_adapter_run(ref, fragments, segs, o, false, &out.plain_report);
    out.views = view_adapter_run(ref, fragments, segs, o, true, &out.views_report);
    if (against_reference) require_identical(ref, out.views, label + " (views vs reference)", false);
    require_identical(out.plain, out.views, label + " (views vs reconstructed)", false);
    require_same_executed(out.plain.result.counters, out.views.result.counters, label + " (views vs reconstructed)");
    const panvar::ContentBatchReport &a = out.plain_report.batch, &b = out.views_report.batch;
    require(a.placements_consumed_last == b.placements_consumed_last && a.origins_consumed_last == b.origins_consumed_last &&
                a.states_consumed_last == b.states_consumed_last && a.origins_constructed == b.origins_constructed &&
                a.origin_entities_admitted == b.origin_entities_admitted && a.origin_entities_declined == b.origin_entities_declined &&
                a.unit_placements_kept == b.unit_placements_kept && a.direct_physical_origins == b.direct_physical_origins &&
                a.served_physical_origins == b.served_physical_origins,
            label + ": served/consumed work or origin admission differs between layouts");
    require(b.evidence_layout == "fine_views" && b.reconstructed_placements == 0 && b.coarse_units_reconstructed == 0 &&
                b.view_rows_represented == a.reconstructed_placements && b.view_coarse_units == a.coarse_units_reconstructed,
            label + ": views must copy no coarse row and represent exactly the rows reconstruction copies");
    return out;
}

void content_index_view_gates(const Fixture& f, const Run& reference) {
    const std::vector<std::size_t> order = {0, 1, 2, 3, 4};
    const auto segs = segments_in_interned_order(f, reference, order);
    // the base fixture at 1 and 4 threads, against the reconstructed layout and the reference
    for (const std::size_t threads : {std::size_t{1}, std::size_t{4}}) {
        const ViewPair v = view_pair(reference, f.fragments, segs, adapter_options(threads), "I4-B base, " + std::to_string(threads) + " thr");
        const auto &a = v.plain_report.batch, &b = v.views_report.batch;
        std::cout << "    I4-B base " << threads << " thr: plain peak " << a.store_peak_bytes << " (kept rows " << a.unit_placements_kept
                  << ", reconstructed " << a.reconstructed_placements << ", fine units " << v.plain_report.index.fine_units
                  << ", fine occurrences " << v.plain_report.index.fine_occurrences << "); views peak " << b.store_peak_bytes
                  << " (table " << b.view_table_bytes << ", zone peak " << b.view_zone_temporary_peak_bytes << ", zone total "
                  << b.view_zone_temporary_bytes_total << ")\n";
        require(b.view_zone_temporary_bytes_total > 0 && a.origin_entities_admitted > 0,
                "I4-B base: origin zones must be gathered from views");
        if (threads == 1) {
            std::cout << "  ok: I4-B views = reconstructed = reference (1 and 4 threads); store peak " << b.store_peak_bytes
                      << " vs " << a.store_peak_bytes << " B; " << b.view_rows_represented << " coarse rows represented, none "
                      << "copied; zone gathers peak " << b.view_zone_temporary_peak_bytes << " B\n";
        }
    }
    // repeats and distinct physical copies within insert range (cross-copy mates; served and direct origins)
    {
        const Fixture c = cross_copy_fixture();
        Config three;
        three.order = {0, 1, 2};
        const Run ref = score(c, three);
        const ViewPair v = view_pair(ref, c.fragments, segments_in_interned_order(c, ref, three.order), adapter_options(2), "I4-B cross-copy");
        require(v.views_report.batch.served_physical_origins > 0 && v.views_report.batch.direct_physical_origins > 0,
                "I4-B cross-copy: both served and direct origins");
    }
    // both orientations and mixed read lengths (prefilter off, as the route requires), clipped/short mates
    {
        Fixture g = f;
        for (std::size_t i = 1; i < g.fragments.size(); i += 2) std::swap(g.fragments[i].first, g.fragments[i].second);
        for (std::size_t i = 0; i < g.fragments.size(); i += 7) {
            g.fragments[i].second.bases.resize(93);
            g.fragments[i].second.qualities.resize(93);
        }
        for (std::size_t i = 3; i < g.fragments.size(); i += 11) {   // clipped at the start
            g.fragments[i].first.bases.erase(0, 9);
            g.fragments[i].first.qualities.erase(0, 9);
        }
        const Run ref = score(g, Config{});
        panvar::FragmentScoringOptions o = adapter_options(4);
        o.union_prefilter_cache_bytes = 0;
        o.max_union_cached_prefilters_consumed = 0;
        (void)view_pair(ref, g.fragments, segments_in_interned_order(g, ref, order), o, "I4-B orientations, mixed and clipped lengths");
    }
    // identical reads in different fragments: every fragment stays its own observation
    {
        Fixture d = f;
        const std::size_t n = d.fragments.size();
        for (std::size_t i = 0; i < n; i += 5) {
            panvar::ReadPair copy = d.fragments[i];
            copy.first.name += "_dup";
            copy.second.name += "_dup";
            d.fragments.push_back(copy);
        }
        const Run ref = score(d, Config{});
        (void)view_pair(ref, d.fragments, segments_in_interned_order(d, ref, order), adapter_options(4), "I4-B identical reads");
    }
    // tiny/empty blocks and a novel junction mosaic
    {
        std::vector<std::string> mosaic;
        std::vector<std::uint32_t> mosaic_segments;
        Fixture t = tiny_fixture(mosaic, mosaic_segments);
        t.haplotypes.push_back(mosaic[0]);
        t.segments.push_back(mosaic_segments);
        Config five;
        five.order = {0, 1, 2, 3, 4};
        const Run ref = score(t, five);
        (void)view_pair(ref, t.fragments, segments_in_interned_order(t, ref, five.order), adapter_options(4), "I4-B tiny blocks, junction");
    }
    // forced hash collisions in the content index; cache exhaustion (a one-byte emission cache; an origin budget that
    // admits only some entities; none)
    {
        panvar::set_content_index_hash_collapse(true);
        (void)view_pair(reference, f.fragments, segs, adapter_options(4), "I4-B forced content-index collisions");
        panvar::set_content_index_hash_collapse(false);
        panvar::FragmentScoringOptions tiny_cache = adapter_options(4);
        tiny_cache.emission_cache_bytes = 1;
        (void)view_pair(reference, f.fragments, segs, tiny_cache, "I4-B one-byte emission cache");
        const ViewPair full = view_pair(reference, f.fragments, segs, adapter_options(1), "I4-B origin budget probe");
        panvar::FragmentScoringOptions partial = adapter_options(2);
        partial.union_origin_cache_bytes = std::max<std::uint64_t>(1, full.plain_report.batch.origin_held_bytes / 2);
        const ViewPair half = view_pair(reference, f.fragments, segs, partial, "I4-B partial origin admission");
        require(half.views_report.batch.origin_entities_declined > 0 && half.views_report.batch.origin_entities_admitted > 0,
                "fixture: the half budget must admit some entities and decline others");
        panvar::FragmentScoringOptions none = adapter_options(2);
        none.union_origin_cache_bytes = 0;
        none.max_union_cached_origins_consumed = 0;
        none.max_union_cached_states_consumed = 0;
        (void)view_pair(reference, f.fragments, segs, none, "I4-B origin reuse off");
    }
    // the store budget: accepted at the views' own peak, refused one byte below (single worker: the peak is exact);
    // the refusal returns no result
    {
        panvar::IndexedUnionReport probe;
        (void)view_adapter_run(reference, f.fragments, segs, adapter_options(1), true, &probe);
        const std::uint64_t peak = probe.batch.store_peak_bytes;
        panvar::IndexedUnionReport at;
        require_identical(reference, view_adapter_run(reference, f.fragments, segs, adapter_options(1), true, &at, peak),
                          "I4-B store budget at the views' peak", false);
        bool refused = false;
        panvar::IndexedUnionReport below;
        try {
            (void)view_adapter_run(reference, f.fragments, segs, adapter_options(1), true, &below, peak - 1);
        } catch (const panvar::ResourceExhausted& e) {
            refused = std::string(e.what()).find("placement store") != std::string::npos;
        }
        require(refused && below.refusal.find("placement store") != std::string::npos,
                "I4-B: one byte below the views' peak must refuse at the placement store");
        // at 4 threads the zone gathers are charged per origin batch before any worker allocates, so the peak (and the
        // verdict) is the same on every schedule
        std::uint64_t peak4 = 0;
        for (int rep = 0; rep < 3; ++rep) {
            panvar::IndexedUnionReport r4;
            (void)view_adapter_run(reference, f.fragments, segs, adapter_options(4), true, &r4);
            require(peak4 == 0 || r4.batch.store_peak_bytes == peak4, "I4-B: the 4-thread store peak must not depend on the schedule");
            peak4 = r4.batch.store_peak_bytes;
        }
        panvar::IndexedUnionReport at4;
        require_identical(reference, view_adapter_run(reference, f.fragments, segs, adapter_options(4), true, &at4, peak4),
                          "I4-B store budget at the 4-thread peak", false);
        bool refused4 = false;
        try {
            panvar::IndexedUnionReport below4;
            (void)view_adapter_run(reference, f.fragments, segs, adapter_options(4), true, &below4, peak4 - 1);
        } catch (const panvar::ResourceExhausted&) {
            refused4 = true;
        }
        require(refused4, "I4-B: one byte below the 4-thread peak must refuse");
        std::cout << "  ok: I4-B store budget accepts at the views' peak (" << peak << " B at 1 thread, " << peak4
                  << " B at 4, schedule-independent over 3 runs) and refuses one byte below\n";
    }
    // the batch ledger equals the recount; materialised view rows equal the reconstructed rows (1 and 4 threads)
    for (const std::size_t threads : {std::size_t{1}, std::size_t{4}}) {
        panvar::ContentPlacementIndex index(100, library().max_gap_bp, true);
        for (std::size_t h = 0; h < reference.interned.sequences.size(); ++h) index.add_haplotype(reference.interned.sequences[h], segs[h]);
        panvar::ContentIndexBatch plain(index, f.fragments, library(), index_options(threads, 0));
        panvar::ContentIndexBatch views(index, f.fragments, library(), index_options(threads, 0));
        plain.set_one_shot(false);   // keep the fine rows so both expose every coarse unit
        views.set_fine_views(true);
        views.set_one_shot(true);
        plain.discover();
        views.discover();
        require(views.report().store_live_bytes == views.store_capacity_recount(),
                "I4-B: the ledger must equal the recounted capacity (rows + container + view table)");
        require_same_rows(plain, views, "I4-B materialised views");
    }
    // a legacy-layout batch refuses views
    {
        panvar::ContentPlacementIndex index(100, library().max_gap_bp, true);
        for (std::size_t h = 0; h < reference.interned.sequences.size(); ++h) index.add_haplotype(reference.interned.sequences[h], segs[h]);
        panvar::ContentIndexBatch b(index, f.fragments, library(), index_options(1, 0));
        b.set_legacy_store_layout(true);
        b.set_fine_views(true);
        require(adapter_invalid([&] { b.discover(); }), "I4-B: views with the legacy store layout must be refused");
    }
    std::cout << "  ok: I4-B views on repeats/cross-copy mates, orientations + mixed/clipped lengths, identical reads, tiny "
                 "blocks + junction, forced collisions, a one-byte emission cache, partial/no origin admission; ledger = "
                 "recount; materialised views = reconstructed rows; legacy layout refused\n";
}

struct LedgerRun {
    panvar::ContentBatchReport report;
    std::uint64_t recount = 0;
    bool refused = false;
};
LedgerRun ledger_run(const Fixture& f, const Run& reference, const std::vector<std::vector<std::uint32_t>>& segs, std::size_t threads,
                     std::size_t budget, bool legacy, bool one_shot = true) {
    panvar::ContentPlacementIndex index(100, library().max_gap_bp, true);
    for (std::size_t h = 0; h < reference.interned.sequences.size(); ++h) index.add_haplotype(reference.interned.sequences[h], segs[h]);
    panvar::ContentIndexBatch batch(index, f.fragments, library(), index_options(threads, 0));
    batch.set_placement_store_budget(budget);
    batch.set_legacy_store_layout(legacy);
    batch.set_one_shot(one_shot);
    LedgerRun out;
    try {
        batch.discover();
    } catch (const panvar::ResourceExhausted&) {
        out.refused = true;
    }
    out.report = batch.report();
    out.recount = batch.store_capacity_recount();
    return out;
}

// The capacity gate: the ledger's live bytes are the recounted capacity, and the exact layout's budget boundary sits at
// its own peak, at 1 and 4 threads. Returns a description of the first failure ("" when it holds).
std::string capacity_gate(const Fixture& f, const Run& reference, const std::vector<std::vector<std::uint32_t>>& segs) {
    std::uint64_t peak = 0;
    for (const std::size_t threads : {std::size_t{1}, std::size_t{4}}) {
        const LedgerRun free_run = ledger_run(f, reference, segs, threads, 0, false);
        if (free_run.refused) return "an unbudgeted run refused";
        if (free_run.report.store_live_bytes != free_run.recount) {
            return "ledger live bytes " + std::to_string(free_run.report.store_live_bytes) + " differ from the recounted capacity " +
                   std::to_string(free_run.recount);
        }
        if (peak == 0) peak = free_run.report.store_peak_bytes;
        if (free_run.report.store_peak_bytes != peak) return "the peak depends on the thread count";
        const LedgerRun at = ledger_run(f, reference, segs, threads, peak, false);
        const LedgerRun below = ledger_run(f, reference, segs, threads, peak - 1, false);
        if (at.refused) return "a budget equal to the peak refused at " + std::to_string(threads) + " threads";
        if (!below.refused) return "a budget one byte below the peak accepted at " + std::to_string(threads) + " threads";
        if (below.report.store_live_bytes != below.recount) return "a refused run's ledger differs from what it holds";
        if (below.report.store_live_bytes > peak - 1) return "a refused run's ledger exceeds its budget";
    }
    return "";
}

void content_index_ledger_gates(const Fixture& f, const Run& reference) {
    const std::vector<std::size_t> order = {0, 1, 2, 3, 4};
    const auto segs = segments_in_interned_order(f, reference, order);
    // exact layout: ledger = capacity, boundary exact, threads invariant
    {
        const std::string failure = capacity_gate(f, reference, segs);
        require(failure.empty(), "I2-C capacity gate: " + failure);
    }
    // unused capacity: spare slots per stored vector are charged, and the boundary moves with them
    panvar::set_store_spare_capacity(7);
    {
        const std::string failure = capacity_gate(f, reference, segs);
        require(failure.empty(), "I2-C capacity gate with unused capacity: " + failure);
        const LedgerRun spare = ledger_run(f, reference, segs, 1, 0, false);
        const LedgerRun plain = [&] {
            panvar::set_store_spare_capacity(0);
            LedgerRun r = ledger_run(f, reference, segs, 1, 0, false);
            panvar::set_store_spare_capacity(7);
            return r;
        }();
        require(spare.report.store_live_bytes > plain.report.store_live_bytes &&
                    spare.report.store_peak_bytes > plain.report.store_peak_bytes,
                "unused capacity must be charged");
    }
    // negative control: charging occupied rows only must fail the capacity gate
    panvar::set_store_charge_rows_only(true);
    {
        const std::string failure = capacity_gate(f, reference, segs);
        require(!failure.empty(), "I2-C negative control: charging rows only must fail the capacity gate");
        std::cout << "  ok: I2-C negative control (rows-only charging) fails the capacity gate: " << failure << "\n";
    }
    panvar::set_store_charge_rows_only(false);
    panvar::set_store_spare_capacity(0);
    // reconstruction overlap: the exact layout releases fine rows as soon as their last coarse unit is built
    const LedgerRun exact = ledger_run(f, reference, segs, 4, 0, false);
    const LedgerRun legacy = ledger_run(f, reference, segs, 4, 0, true);
    require(legacy.report.store_live_bytes == legacy.recount, "the legacy layout's ledger must measure its capacity");
    require(exact.report.store_peak_bytes > exact.report.store_live_bytes && exact.report.store_released_bytes > 0,
            "fine and coarse rows overlap during reconstruction, then fine rows are released");
    require(exact.report.store_peak_bytes < legacy.report.store_peak_bytes &&
                exact.report.store_live_bytes <= legacy.report.store_live_bytes,
            "the exact layout must hold less at its peak and after discovery than the I2-A layout");
    require(exact.report.store_allocated_bytes - exact.report.store_released_bytes == exact.report.store_live_bytes,
            "allocated minus released must be live");
    require(exact.report.worker_scratch_bytes > 0, "the workers' temporary storage must be reported");
    std::cout << "  ok: I2-C ledger = recounted capacity (spare capacity charged); boundary exact at 1 and 4 threads; peak "
              << exact.report.store_peak_bytes << " vs I2-A layout " << legacy.report.store_peak_bytes << ", after discovery "
              << exact.report.store_live_bytes << " vs " << legacy.report.store_live_bytes << " bytes; worker scratch "
              << exact.report.worker_scratch_bytes << "\n";
    // incremental (not one-shot): fine rows are kept, the ledger stays exact across discover() calls
    {
        std::vector<std::string> mosaic;
        std::vector<std::uint32_t> mosaic_segments;
        Fixture t = tiny_fixture(mosaic, mosaic_segments);
        t.haplotypes.push_back(mosaic[0]);
        t.segments.push_back(mosaic_segments);
        Config five;
        five.order = {0, 1, 2, 3, 4};
        const Run ref = score(t, five);
        panvar::ContentPlacementIndex ti(100, library().max_gap_bp, true);
        panvar::ContentIndexBatch tb(ti, t.fragments, library(), index_options(4, 0));
        std::vector<std::size_t> tu;
        for (std::size_t i = 0; i < 4; ++i) { ti.add_haplotype(t.haplotypes[i], t.segments[i]); tu.push_back(ref.interned.index_of[i]); }
        tb.discover();
        require(tb.report().store_live_bytes == tb.store_capacity_recount() && tb.report().store_released_bytes == 0,
                "incremental: the ledger is exact and nothing is released");
        ti.add_haplotype(t.haplotypes[4], t.segments[4]);
        tu.push_back(ref.interned.index_of[4]);
        tb.discover();
        require(tb.report().store_live_bytes == tb.store_capacity_recount(), "incremental: the ledger stays exact after growth");
        std::vector<panvar::HaplotypeCandidate> every;
        for (std::size_t a = 0; a < 5; ++a) for (std::size_t b = a; b < 5; ++b) every.push_back({a, b});
        require_stage(tb.score(every), ref, tu, every, t.fragments.size(), "I2-C incremental after growth");
    }
    // refusal through the adapter: no result; then a fresh call in the same process scores exactly (nothing leaked)
    {
        const LedgerRun free_run = ledger_run(f, reference, segs, 1, 0, false);
        bool refused = false;
        try {
            (void)adapter_run(reference, f.fragments, segs, adapter_options(4), nullptr, free_run.report.store_peak_bytes - 1);
        } catch (const panvar::ResourceExhausted&) {
            refused = true;
        }
        require(refused, "the adapter must refuse below the exact layout's peak");
        require_identical(reference, adapter_run(reference, f.fragments, segs, adapter_options(4)), "I2-C after a refusal", false);
        panvar::IndexedUnionReport legacy_report;
        Run legacy_run;
        legacy_run.interned = reference.interned;
        legacy_run.candidates = reference.candidates;
        legacy_run.result = panvar::score_diplotypes_indexed_union(reference.interned, reference.candidates, f.fragments, library(),
                                                                   adapter_options(1), segs, kAdapterStore, &legacy_report,
                                                                   nullptr, true);
        legacy_run.ranking = panvar::certify_scored_ranking(legacy_run.result.scores, f.fragments.size(), 6,
                                                            legacy_run.result.counters, panvar::CertificationOptions{});
        require_identical(reference, legacy_run, "I2-C legacy layout through the adapter", false);
        require(legacy_report.batch.store_layout == "legacy_i2a", "the legacy layout must report itself");
    }
    std::cout << "  ok: I2-C incremental ledger exact across growth (nothing released); adapter refusal returns nothing and the "
                 "next call is exact; legacy layout identical\n";
}

// ---- I2-H: the fixed-geometry emitting recursion is invisible in placement ------------------------------------------
// Read lengths around the kernel's word boundaries (64, 65, 128, 129) and the production length; substitutions, indels and
// a tandem repeat; reads at and just above the mismatch band; clipped windows at both locus ends; both orientations; the
// emission cache off, on, tiny (most fits not retained) and with forced hash collisions; 1 and 4 threads; the reference
// pass and the indexed batch. The experimental backend must equal the reference in every accepted origin (start, insert,
// orientation, latent-mass bits, latent pair count and diagnostics, so origin multiplicities too), every mass cell,
// exposure, candidate component and the certified ranking, and in every counter, logical and (single-threaded) executed.
// Ambiguous bases keep their refusal. At 64-65 bp the default band 6 leaves pieces under the 10-bp minimum (the scorer's
// own refusal), so those lengths use 5.
std::size_t dp_band(std::size_t read_length) { return read_length < 70 ? 5 : 6; }
Fixture dp_fixture(std::uint64_t seed, std::size_t read_length) {
    std::mt19937_64 rng(seed);
    Fixture f;
    const std::string unit = random_bases(rng, 70);
    std::string base = random_bases(rng, 500);
    for (int i = 0; i < 4; ++i) base += unit;
    base += std::string(14, 'T') + random_bases(rng, 486);   // a homopolymer for indel ambiguity
    std::string subs = base, ins = base, del = base;
    for (const std::size_t at : {std::size_t{100}, std::size_t{400}, std::size_t{900}}) subs[at] = subs[at] == 'A' ? 'C' : 'A';
    ins.insert(650, random_bases(rng, 2));
    ins.insert(783, "T");
    del.erase(1100, 3);
    del.erase(790, 1);
    f.haplotypes = {base, subs, ins, del, base.substr(0, 1000) + random_bases(rng, 300)};
    for (const std::string& h : f.haplotypes) f.segments.push_back({static_cast<std::uint32_t>(h.size())});
    const std::size_t insert = 2 * read_length + 60;
    std::size_t n = 0;
    for (const std::string* src : {&base, &subs, &ins, &del}) {
        for (std::size_t start = 0; start + insert <= src->size(); start += 41, ++n) {
            panvar::ReadPair pair = fragment_at(*src, start, insert, read_length, "q" + std::to_string(n), rng);
            const std::size_t band = dp_band(read_length);
            const std::size_t edits = n % 9 == 0 ? band : n % 9 == 1 ? band + 1 : n % 9 == 2 ? 2 : 0;
            for (std::size_t e = 0; e < edits; ++e) {
                char& b = pair.first.bases[(e * 13 + 5) % read_length];
                b = b == 'A' ? 'C' : 'A';
            }
            if (n % 2 == 1) std::swap(pair.first, pair.second);   // both orientations
            f.fragments.push_back(pair);
        }
    }
    // templates flush against both locus ends: clipped windows
    f.fragments.push_back(fragment_at(base, 0, insert, read_length, "start", rng));
    f.fragments.push_back(fragment_at(base, base.size() - insert, insert, read_length, "end", rng));
    return f;
}

void require_same_origins(const Run& a, const Run& b, const std::string& label) {
    const auto& x = a.result.recorded_origins;
    const auto& y = b.result.recorded_origins;
    require(!x.empty(), label + ": no origins were recorded (vacuous)");
    require(x.size() == y.size(), label + ": recorded (fragment, haplotype) cells differ");
    for (const auto& cell : x) {
        const auto found = y.find(cell.first);
        require(found != y.end() && found->second.size() == cell.second.size(), label + ": origin multiplicity differs");
        for (std::size_t i = 0; i < cell.second.size(); ++i) {
            const panvar::FragmentPlacement& p = cell.second[i];
            const panvar::FragmentPlacement& q = found->second[i];
            require(p.forward_start == q.forward_start && p.insert_bp == q.insert_bp && p.orientation == q.orientation &&
                        same_double(p.log_latent_mass, q.log_latent_mass) && p.latent_alignments == q.latent_alignments &&
                        p.reverse_start == q.reverse_start && p.forward_consumed_bp == q.forward_consumed_bp &&
                        p.reverse_consumed_bp == q.reverse_consumed_bp && p.forward_mismatches == q.forward_mismatches &&
                        p.reverse_mismatches == q.reverse_mismatches,
                    label + ": an accepted origin differs");
        }
    }
}

void require_same_dp_counters(const Run& a, const Run& b, const std::string& label, bool executed) {
    const panvar::FragmentScoringCounters& x = a.result.counters;
    const panvar::FragmentScoringCounters& y = b.result.counters;
    require(x.piece_hits == y.piece_hits && x.proposed_starts == y.proposed_starts && x.unique_starts == y.unique_starts &&
                x.full_read_verifications == y.full_read_verifications && x.dp_evaluations == y.dp_evaluations &&
                x.prefilter_rejections == y.prefilter_rejections && x.accepted_mate_placements == y.accepted_mate_placements &&
                x.physical_origins == y.physical_origins && x.statistical_state_classes == y.statistical_state_classes,
            label + ": a logical counter differs");
    require(!executed || (a.result.emission_cache.executed_prefilters == b.result.emission_cache.executed_prefilters &&
                          a.result.emission_cache.executed_dps == b.result.emission_cache.executed_dps &&
                          a.result.emission_cache.dp_cache_hits == b.result.emission_cache.dp_cache_hits &&
                          a.result.emission_cache.entries == b.result.emission_cache.entries),
            label + ": an executed counter differs");
}

void gap_recursion_backend_gates() {
    using panvar::GapRecursionBackend;
    std::size_t compared = 0;
    const std::uint64_t fixed_before = panvar::gapped_dp_test::fixed_geometry_executions();
    for (const std::size_t length : {std::size_t{64}, std::size_t{65}, std::size_t{100}, std::size_t{128}, std::size_t{129}}) {
        const Fixture f = dp_fixture(900 + length, length);
        for (const std::size_t threads : {std::size_t{1}, std::size_t{4}}) {
            for (int mode = 0; mode < 4; ++mode) {   // emission cache: 0 off, 1 on, 2 tiny budget, 3 forced collisions
                Config reference;
                reference.threads = threads;
                reference.band = dp_band(length);
                reference.record_all_origins = true;
                reference.emission_cache_bytes = mode == 0 ? 0 : mode == 2 ? std::size_t{4096} : std::size_t{1} << 22;
                if (mode == 1) {   // the union prefilter-decision cache too
                    reference.prefilter_budget = std::size_t{1} << 26;
                    reference.prefilter_cap = std::uint64_t{1} << 40;
                }
                Config experimental = reference;
                experimental.backend = GapRecursionBackend::fixed_geometry_experimental;
                if (mode == 3) panvar::set_emission_cache_hash_collapse(true);
                const Run a = score(f, reference);
                const Run b = score(f, experimental);
                if (mode == 3) panvar::set_emission_cache_hash_collapse(false);
                const std::string label = "L=" + std::to_string(length) + ", emission cache mode " + std::to_string(mode) +
                                          ", " + std::to_string(threads) + " threads";
                require_identical(a, b, label);
                require_same_origins(a, b, label);
                require_same_dp_counters(a, b, label, threads == 1);
                require(a.result.emission_cache.executed_dps > 0, label + ": no recursion ran (vacuous)");
                if (mode == 2) {
                    require(a.result.emission_cache.entries_not_retained > 0,
                            label + ": the tiny budget must decline entries (limited-budget fallback)");
                }
                ++compared;
            }
        }
        // the indexed batch: experimental equals the reference scorer and the indexed reference backend
        Config banded;
        banded.band = dp_band(length);
        const Run ref = score(f, banded);
        for (const GapRecursionBackend backend : {GapRecursionBackend::scaled_reference, GapRecursionBackend::fixed_geometry_experimental}) {
            panvar::FragmentScoringOptions o = index_options(4, std::size_t{1} << 26);
            o.gap_recursion_backend = backend;
            o.mismatch_band = dp_band(length);
            panvar::ContentPlacementIndex index(length, library().max_gap_bp, true);
            const auto segs = segments_in_interned_order(f, ref, {0, 1, 2, 3, 4});
            for (std::size_t h = 0; h < ref.interned.sequences.size(); ++h) index.add_haplotype(ref.interned.sequences[h], segs[h]);
            panvar::ContentIndexBatch batch(index, f.fragments, library(), o);
            batch.discover();
            Run run;
            run.interned = ref.interned;
            run.candidates = ref.candidates;
            run.result = batch.score(run.candidates);
            run.ranking = panvar::certify_scored_ranking(run.result.scores, f.fragments.size(), dp_band(length),
                                                         ref.result.counters, panvar::CertificationOptions{});
            require_identical(ref, run, "indexed, L=" + std::to_string(length) +
                                            (backend == GapRecursionBackend::scaled_reference ? ", reference DP" : ", fixed DP"),
                              false);
            ++compared;
        }
    }
    // Worker threads' executions are not visible here; the single-threaded runs above are on this thread.
    require(panvar::gapped_dp_test::fixed_geometry_executions() > fixed_before,
            "the experimental runs must actually execute the fixed-geometry recursion");
    // ambiguous bases: a read carrying N is handled identically on both backends
    {
        Fixture f = dp_fixture(999, 100);
        f.fragments[3].first.bases[10] = 'N';
        Config reference;
        reference.record_all_origins = true;
        Config experimental = reference;
        experimental.backend = GapRecursionBackend::fixed_geometry_experimental;
        std::string a, b;
        bool a_ok = true, b_ok = true;
        Run ra, rb;
        try { ra = score(f, reference); } catch (const std::exception& e) { a = e.what(); a_ok = false; }
        try { rb = score(f, experimental); } catch (const std::exception& e) { b = e.what(); b_ok = false; }
        require(a_ok == b_ok && a == b, "an ambiguous base must be handled identically on both backends");
        if (a_ok) {
            require_identical(ra, rb, "ambiguous base");
            require_same_origins(ra, rb, "ambiguous base");
        }
    }
    std::cout << "  ok: I2-H fixed-geometry recursion invisible in " << compared << " comparisons (read lengths 64/65/100/128/129; "
                 "edits at/above the band; indels in a repeat and a homopolymer; clipped ends; both orientations; emission "
                 "cache off/on/tiny/collided; 1/4 threads; every accepted origin; indexed); ambiguous bases unchanged\n";
}

// ---- I2-K: the union pass's emission results handed to an adaptive scoring context ------------------------------------
// The handoff is detached from a real indexed union batch and adopted by a PreparedScoringContext. Against a cold
// context (same reads, library, options, haplotypes) every mass cell, exposure, candidate score and logical counter must
// be bit-identical; the handed-off context must actually serve hits from inherited entries and execute fewer
// recursions; content the union never saw (a novel junction) must still be computed. Refusals (changed library or band,
// disabled cache, a budget the inheritance does not fit) must leave the context exactly cold. Exceptions and
// sequential samples must not inherit anything.
struct ContextRun {
    std::vector<double> mass, exposure;
    panvar::FragmentScoringResult scores;
    panvar::FragmentScoringCounters counters;
    std::uint64_t executed_dps = 0, cache_hits = 0;
    panvar::EmissionHandoffReport report;
};

ContextRun context_run(const Fixture& f, const std::vector<std::string>& sequences, panvar::FragmentScoringOptions o,
                       const panvar::LibraryModel& lib, panvar::EmissionCacheHandoff* in) {
    panvar::reset_emission_handoff_report();
    o.emission_handoff_in = in;
    ContextRun run;
    panvar::PreparedScoringContext ctx(f.fragments, lib, o);
    const std::vector<std::size_t> slots = ctx.add_haplotypes(sequences);
    for (std::size_t frag = 0; frag < f.fragments.size(); ++frag) {
        for (const std::size_t s : slots) run.mass.push_back(ctx.mass_of(frag, s));
    }
    for (const std::size_t s : slots) run.exposure.push_back(ctx.exposure_of(s));
    std::vector<panvar::HaplotypeCandidate> candidates;
    for (std::size_t x = 0; x < slots.size(); ++x) {
        for (std::size_t y = x; y < slots.size(); ++y) candidates.push_back({slots[x], slots[y]});
    }
    run.scores = ctx.score_candidates(candidates);
    run.counters = ctx.cumulative_counters();
    const panvar::EmissionCacheCounters live = panvar::read_emission_cache_counters();
    (void)live;
    run.report = panvar::emission_handoff_report();
    return run;
}

void require_same_context_science(const ContextRun& a, const ContextRun& b, const std::string& label) {
    require(a.mass.size() == b.mass.size() && a.exposure.size() == b.exposure.size(), label + ": shapes differ");
    for (std::size_t i = 0; i < a.mass.size(); ++i) require(same_double(a.mass[i], b.mass[i]), label + ": a mass cell differs");
    for (std::size_t i = 0; i < a.exposure.size(); ++i) require(same_double(a.exposure[i], b.exposure[i]), label + ": an exposure differs");
    require(a.scores.scores.size() == b.scores.scores.size(), label + ": candidates differ");
    for (std::size_t c = 0; c < a.scores.scores.size(); ++c) {
        const panvar::DiplotypeScore& x = a.scores.scores[c];
        const panvar::DiplotypeScore& y = b.scores.scores[c];
        require(same_double(x.log_likelihood, y.log_likelihood) && same_double(x.fragment_term, y.fragment_term) &&
                    same_double(x.exposure_term, y.exposure_term) && x.fragments_with_mass == y.fragments_with_mass,
                label + ": a candidate score differs");
    }
    const panvar::FragmentScoringCounters& p = a.counters;
    const panvar::FragmentScoringCounters& q = b.counters;
    require(p.piece_hits == q.piece_hits && p.unique_starts == q.unique_starts && p.full_read_verifications == q.full_read_verifications &&
                p.dp_evaluations == q.dp_evaluations && p.prefilter_rejections == q.prefilter_rejections &&
                p.accepted_mate_placements == q.accepted_mate_placements && p.physical_origins == q.physical_origins &&
                p.statistical_state_classes == q.statistical_state_classes,
            label + ": a logical counter (placements, origins, states) differs");
}

// One handoff, detached from an indexed union batch over the fixture's haplotypes.
void detach_from_union(const Fixture& f, const Run& ref, std::size_t length, std::size_t band, std::size_t threads,
                       std::size_t union_cache_bytes, const panvar::LibraryModel& lib, panvar::EmissionCacheHandoff& out) {
    panvar::FragmentScoringOptions o = index_options(threads, std::size_t{1} << 26);
    o.mismatch_band = band;
    o.emission_cache_bytes = union_cache_bytes;
    o.emission_handoff_out = &out;
    panvar::ContentPlacementIndex index(length, lib.max_gap_bp, true);
    const auto segs = segments_in_interned_order(f, ref, {0, 1, 2, 3, 4});
    for (std::size_t h = 0; h < ref.interned.sequences.size(); ++h) index.add_haplotype(ref.interned.sequences[h], segs[h]);
    panvar::ContentIndexBatch batch(index, f.fragments, lib, o);
    batch.discover();
}

void emission_handoff_gates() {
    std::size_t compared = 0, adopted = 0, refused = 0;
    std::uint64_t inherited_hits = 0;
    for (const std::size_t length : {std::size_t{65}, std::size_t{100}, std::size_t{129}}) {
        const Fixture f = dp_fixture(900 + length, length);
        const std::size_t band = dp_band(length);
        Config banded;
        banded.band = band;
        const Run ref = score(f, banded);
        // the union's haplotypes, plus a novel junction (content the union never placed reads on)
        std::vector<std::string> sequences = ref.interned.sequences;
        sequences.push_back(f.haplotypes[0].substr(0, 700) + "ACGTTGCAGGCTTACG" + f.haplotypes[3].substr(820));
        for (const std::size_t threads : {std::size_t{1}, std::size_t{4}}) {
            for (int mode = 0; mode < 7; ++mode) {
                // 0 adopted; 1 adopted under forced hash collisions; 2 adopted from an exhausted union cache;
                // 3 refused: adaptive cache disabled; 4 refused: the inheritance exceeds the adaptive budget;
                // 5 refused: changed library; 6 refused: changed band
                panvar::LibraryModel union_lib = library();
                panvar::LibraryModel adaptive_lib = library();
                std::size_t union_band = band, adaptive_band = band;
                std::size_t union_cache = std::size_t{1} << 26, adaptive_cache = std::size_t{1} << 26;
                if (mode == 2) union_cache = std::size_t{1} << 16;
                if (mode == 3) adaptive_cache = 0;
                if (mode == 5) adaptive_lib.substitution_error = 0.03;
                if (mode == 6) adaptive_band = band - 1;
                if (mode == 1) panvar::set_emission_cache_hash_collapse(true);
                panvar::EmissionCacheHandoff handoff;
                detach_from_union(f, ref, length, union_band, threads, union_cache, union_lib, handoff);
                require(handoff.holds() && handoff.entries() > 0 && handoff.allocated_bytes() > 0,
                        "the union batch must detach its emission results into the handoff");
                if (mode == 4) adaptive_cache = static_cast<std::size_t>(handoff.budgeted_bytes() - 1);
                panvar::FragmentScoringOptions o = options(threads);
                o.mismatch_band = adaptive_band;
                o.emission_cache_bytes = adaptive_cache;
                const std::string label = "L=" + std::to_string(length) + ", mode " + std::to_string(mode) + ", " +
                                          std::to_string(threads) + " threads";
                const ContextRun cold = context_run(f, sequences, o, adaptive_lib, nullptr);
                const ContextRun warm = context_run(f, sequences, o, adaptive_lib, &handoff);
                if (mode == 1) panvar::set_emission_cache_hash_collapse(false);
                require(!handoff.holds(), label + ": adoption, or its refusal, consumes the handoff");
                require_same_context_science(cold, warm, label);
                if (mode <= 2) {
                    require(warm.report.adopted && warm.report.inherited_hits > 0,
                            label + ": the inherited results must actually serve hits");
                    require(warm.counters.dp_evaluations == cold.counters.dp_evaluations,
                            label + ": the logical recursion stage is the cold one's");
                    ++adopted;
                    inherited_hits += warm.report.inherited_hits;
                } else {
                    const char* reason = mode == 3 ? "cache_disabled" : mode == 4 ? "budget" : "context";
                    require(!warm.report.adopted && warm.report.refusal == reason && warm.report.inherited_hits == 0,
                            label + ": refused as " + reason);
                    ++refused;
                }
                ++compared;
            }
        }
    }
    // Executed work: inherited hits replace recursions; novel content is still computed; and the context's own entry
    // counter starts at zero (inherited entries are not stage work).
    {
        const Fixture f = dp_fixture(1000, 100);
        Config banded;
        banded.band = 6;
        const Run ref = score(f, banded);
        std::vector<std::string> sequences = ref.interned.sequences;
        sequences.push_back(f.haplotypes[0].substr(0, 700) + "ACGTTGCAGGCTTACG" + f.haplotypes[3].substr(820));
        panvar::FragmentScoringOptions o = options(1);
        o.emission_cache_bytes = std::size_t{1} << 26;
        const auto executed = [&](panvar::EmissionCacheHandoff* in, std::uint64_t* entries, std::uint64_t* hits) {
            panvar::reset_emission_handoff_report();
            panvar::FragmentScoringOptions oo = o;
            oo.emission_handoff_in = in;
            panvar::PreparedScoringContext ctx(f.fragments, library(), oo);
            (void)panvar::take_gapped_emission_counters();
            (void)ctx.add_haplotypes(sequences);
            *entries = panvar::read_emission_cache_counters().entries;
            *hits = panvar::emission_handoff_report().inherited_hits;
            return ctx.cumulative_emission_cache().executed_dps + panvar::take_gapped_emission_counters().executed_dps;
        };
        std::uint64_t cold_entries = 0, warm_entries = 0, cold_hits = 0, warm_hits = 0;
        const std::uint64_t cold = executed(nullptr, &cold_entries, &cold_hits);
        panvar::EmissionCacheHandoff handoff;
        detach_from_union(f, ref, 100, 6, 1, std::size_t{1} << 26, library(), handoff);
        const std::uint64_t inherited = handoff.entries();
        const std::uint64_t warm = executed(&handoff, &warm_entries, &warm_hits);
        require(warm < cold && warm > 0, "inherited results must replace recursions, and the novel junction must still compute");
        require(warm_entries < cold_entries && warm_entries > 0,
                "the adaptive cache's own entry count must exclude inherited entries (only its new stores)");
        require(cold_hits == 0 && warm_hits > 0, "only the handed-off context serves inherited hits");
        std::cout << "  ok: I2-K handoff executed recursions " << cold << " -> " << warm << " (" << inherited
                  << " inherited entries, " << warm_hits << " inherited hits; new stage entries " << cold_entries << " -> "
                  << warm_entries << ")\n";
    }
    // An exception before adoption leaves the handoff with its owner (freed with it); a later context is cold.
    {
        const Fixture f = dp_fixture(1100, 100);
        Config banded;
        banded.band = 6;
        const Run ref = score(f, banded);
        panvar::EmissionCacheHandoff handoff;
        detach_from_union(f, ref, 100, 6, 1, std::size_t{1} << 26, library(), handoff);
        std::vector<panvar::ReadPair> bad = f.fragments;
        bad[0].first.bases.clear();   // invalid: the context refuses at preparation, before the cache is armed
        bool threw = false;
        try {
            panvar::FragmentScoringOptions o = options(1);
            o.emission_handoff_in = &handoff;
            panvar::PreparedScoringContext ctx(bad, library(), o);
        } catch (const std::exception&) {
            threw = true;
        }
        require(threw && handoff.holds(), "a context that throws before arming leaves the handoff with its owner");
        handoff.clear();
        // the next sample: a cold context sees no inherited entry and serves no inherited hit
        const ContextRun next = context_run(f, ref.interned.sequences, options(1), library(), &handoff);
        require(!next.report.adopted && next.report.inherited_hits == 0 && next.report.refusal == "not_offered",
                "a sequential sample inherits nothing");
        const ContextRun cold = context_run(f, ref.interned.sequences, options(1), library(), nullptr);
        require_same_context_science(cold, next, "sequential sample");
    }
    std::cout << "  ok: I2-K emission handoff bit-identical to cold contexts in " << compared << " comparisons (" << adopted
              << " adopted, " << inherited_hits << " inherited hits; " << refused
              << " refused: cache disabled, budget, library, band; read lengths 65/100/129; edits at/above the band; "
                 "indels in a repeat and a homopolymer; clipped ends; both orientations; a novel junction; forced "
                 "collisions; an exhausted union cache; 1/4 threads); exception and sequential sample inherit nothing\n";
}


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
                // Above one thread the SERVED count is scheduling-dependent: a worker may miss a decision another
                // worker has not stored yet (measured under load: 11823..11905 served of 26638 requests, against 11905
                // at one thread), so no fixed cap near the served count has a schedule-independent outcome, and the
                // cap's exact boundary is asserted at one thread above. At four threads the gate asserts only what
                // holds in every schedule: with nothing declined, concurrency can only lower the served count, so the
                // one-thread count as the cap always accepts exactly; and one below it, a run either refuses with no
                // ranking or completes having served no more than its cap, with the reference's scores.
                require(p.decisions_not_retained == 0, "the four-thread bound needs every decision retained");
                Config four_at = with;
                four_at.threads = 4;
                four_at.prefilter_cap = p.hits;
                const Run four_ok = score(f, four_at);
                require(four_ok.result.union_prefilter_cache.hits <= p.hits,
                        "four threads never serve more decisions than one thread");
                require_prefilter_exact(plain, four_ok, "prefilter consumption cap at the one-thread count, 4 threads");
                below.threads = 4;
                below.prefilter_cap = p.hits - 1;
                for (int rep = 0; rep < 8; ++rep) {
                    bool refused = false;
                    try {
                        const Run r = score(f, below);
                        require(r.result.union_prefilter_cache.hits <= below.prefilter_cap,
                                "a completing four-thread run served no more than its prefilter cap");
                        require_prefilter_exact(plain, r, "prefilter cap one below, 4 threads, completing schedule");
                    } catch (const panvar::ResourceExhausted&) {
                        refused = true;
                    }
                    (void)refused;   // either outcome is sound; neither may return a partial ranking
                }
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

        // ---- out-of-place start sorter ---------------------------------------------------------------------
        {
            Config full_scratch = p2;
            full_scratch.prefilter_budget = std::size_t{1} << 24;
            full_scratch.prefilter_cap = many;
            for (const std::size_t budget : {std::size_t{64} << 20, std::size_t{1}}) {
                for (const std::size_t threads : {std::size_t{1}, std::size_t{2}, std::size_t{4}, std::size_t{8}}) {
                    Config c = full_scratch;
                    c.sort_scratch = budget;
                    c.threads = threads;
                    const Run sorted_run = score(f, c);
                    const panvar::StartSortReport& s = sorted_run.result.start_sort;
                    require(s.peak_worker_bytes <= s.share_bytes, "a worker's sort scratch exceeded its share");
                    if (budget == 1) {
                        require(s.out_of_place == 0 && s.fallbacks > 0 && s.fallbacks == s.in_place,
                                "a one-byte scratch budget must order every list in place, counted as fallbacks");
                    } else {
                        require(s.out_of_place > 0 && s.fallbacks == 0, "a full scratch budget must order out of place");
                    }
                    require_identical(reference, sorted_run, "start sorter, budget " + std::to_string(budget) + ", " +
                                                                 std::to_string(threads) + " threads");
                }
            }
            Config permuted_scratch = full_scratch;
            permuted_scratch.sort_scratch = std::size_t{64} << 20;
            permuted_scratch.order = {3, 0, 4, 2, 1};
            permuted_scratch.threads = 4;
            require_identical(reference, score(f, permuted_scratch), "start sorter, permuted haplotypes");
            permuted_scratch.reverse_candidates = true;
            require_identical(reference, score(f, permuted_scratch), "start sorter, reversed requests");
            std::cout << "  ok: out-of-place start sorter within its scratch share, and every list in place (counted as "
                         "fallbacks) under a one-byte budget, at 1/2/4/8 threads, permuted haplotypes and reversed "
                         "requests\n";
        }

        panvar::set_origin_cache_self_check(true);
        for (const std::size_t threads : {std::size_t{1}, std::size_t{4}}) {
            Config c = p2;
            c.threads = threads;
            require_identical(reference, score(f, c), "origin reuse self-check");
        }
        panvar::set_origin_cache_self_check(false);
        std::cout << "  ok: every served fragment census equals a fresh enumeration (self-check, 1 and 4 threads)\n";
        content_index_gates(f, reference);
        content_index_prefilter_gates(f, reference);
        content_index_origin_gates(f, reference);
        content_index_incremental_gates(f, reference);
        content_index_fine_gates(f, reference);
        content_index_adapter_gates(f, reference);
        content_index_view_gates(f, reference);
        content_index_ledger_gates(f, reference);
        gap_recursion_backend_gates();
        emission_handoff_gates();
        std::cout << "union placement cache gates passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "union placement cache gate FAILED: " << error.what() << "\n";
        return 1;
    }
}
