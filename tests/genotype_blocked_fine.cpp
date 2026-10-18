// S2-F (EXPERIMENTAL): fragment blocking combined with the I4-E fine-core origin aggregates. The combined route (every
// block's batch builds and owns its own whole and fine aggregates under the unchanged total budget, split half/half) is
// checked against ordinary blocking, the bulk indexed adapter (control and fine) and the exhaustive reference scorer.
//
// Fixtures: the I4-E fixtures (similar long alleles with and without shared flanks: single-occurrence entities, a
// duplicated stretch within insert range, a context variant; empty, consecutive empty, terminal empty and tiny segments;
// shifted identical windows; a novel junction), reads of 90/100/110 bp in both mate orders with substitutions and 2 bp
// read-side indels, two libraries; and the S2-E fixture (repeats, a distant copy, indels, a junction, an empty block
// joined across, a short haplotype, an identical labelling) with the prefilter cache.
// Required:
// - bit-identical masses (per block), score components, audit views and certified top sets;
// - exact-once origin accounting per block and cumulatively (executed - constructed + served = the reference's logical
//   origins; logical state classes equal), executed and served work kept distinct;
// - partial admission, complete fallback and the exact fine-partition boundary inside a block;
// - forced fine-key hash collisions; 1 vs 4 threads identical in science, work and admission;
// - a cumulative refusal in a later block (pairing, construction-time and discovery-time caps) returns no scores and
//   reports the completed blocks, the refused block's observed work and the whole attempt.
#include "panvar/fragment_likelihood.hpp"
#include "panvar/genotype_model.hpp"

#include <algorithm>
#include <cstring>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using panvar::BlockedUnionResult;
using panvar::DiplotypeScore;
using panvar::HaplotypeCandidate;
using panvar::ReadPair;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}
bool bits(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }
bool same(const DiplotypeScore& a, const DiplotypeScore& b) {
    return bits(a.log_likelihood, b.log_likelihood) && bits(a.fragment_term, b.fragment_term) &&
           bits(a.exposure_term, b.exposure_term) && bits(a.absolute_term_sum, b.absolute_term_sum) &&
           a.fragments_with_mass == b.fragments_with_mass;
}
std::string revcomp(const std::string& s) {
    std::string out(s.rbegin(), s.rend());
    for (char& c : out) c = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : 'A';
    return out;
}
std::string random_bases(std::mt19937_64& rng, std::size_t n) {
    std::string s(n, 'A');
    for (char& c : s) c = "ACGT"[rng() % 4];
    return s;
}

panvar::LibraryModel library(double mean, double sd) {
    panvar::LibraryModel m;
    m.insert_mean_bp = mean;
    m.insert_sd_bp = sd;
    m.support_sigmas = 3.0;
    m.substitution_error = 0.02;
    m.outlier_weight = 0.005;
    m.fragment_start_intensity = 0.01;
    m.provenance = "fixture: declared";
    return m;
}

// I4-E's options; `audits` adds two audit libraries and the prepared library terms (blocked vs bulk only)
panvar::FragmentScoringOptions options(std::size_t threads, std::size_t origin_bytes, const panvar::LibraryModel* audits = nullptr,
                                       bool prefilter = false) {
    panvar::FragmentScoringOptions o;
    o.mismatch_band = 6;
    o.minimum_piece_bp = 10;
    o.threads = threads;
    o.piece_hit_publish_batch = 1;
    o.emission_cache_bytes = 1u << 20;
    o.union_origin_cache_bytes = origin_bytes;
    o.max_union_cached_origins_consumed = std::uint64_t{1} << 40;
    o.max_union_cached_states_consumed = std::uint64_t{1} << 40;
    o.max_union_cached_placements_consumed = std::uint64_t{1} << 40;
    if (prefilter) {
        o.union_prefilter_cache_bytes = 1 << 22;
        o.max_union_cached_prefilters_consumed = std::uint64_t{1} << 40;
    }
    if (audits != nullptr) {
        for (const double factor : {0.9, 1.1}) {
            panvar::LibraryModel l = *audits;
            l.fragment_start_intensity *= factor;
            o.audit_libraries.push_back(l);
        }
        o.prepared_library_terms = true;
    }
    return o;
}

struct Fixture {
    std::string name;
    std::vector<std::string> haplotypes;
    std::vector<std::vector<std::uint32_t>> segments;
    std::vector<ReadPair> fragments;
    bool uniform_reads = false;
};

void add(Fixture& f, const std::vector<std::string>& pieces) {
    std::string h;
    std::vector<std::uint32_t> lengths;
    for (const std::string& p : pieces) {
        h += p;
        lengths.push_back(static_cast<std::uint32_t>(p.size()));
    }
    f.haplotypes.push_back(h);
    f.segments.push_back(lengths);
}

void simulate(Fixture& f, double mean, std::size_t step, std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::size_t n = 0;
    for (const std::string& h : f.haplotypes) {
        std::vector<std::size_t> starts;
        for (std::size_t s = 0; s + 500 <= h.size(); s += step) starts.push_back(s);
        starts.push_back(h.size() - static_cast<std::size_t>(mean) - 20);   // the haplotype's right end
        for (const std::size_t start : starts) {
            const std::size_t insert = static_cast<std::size_t>(mean) - 30 + rng() % 60;
            if (start + insert > h.size()) continue;
            const std::size_t length = f.uniform_reads ? 100 : 90 + 10 * (n % 3);
            std::string forward = h.substr(start, length);
            std::string reverse = revcomp(h.substr(start + insert - length, length));
            if (rng() % 3 == 0) forward[rng() % forward.size()] = "ACGT"[rng() % 4];
            if (rng() % 3 == 0) reverse[rng() % reverse.size()] = "ACGT"[rng() % 4];
            if (!f.uniform_reads && n % 7 == 3) forward.erase(40, 2);
            if (!f.uniform_reads && n % 11 == 5) reverse.insert(50, "TT");
            ReadPair p{{"r" + std::to_string(n) + "/1", forward, std::string(forward.size(), 'I')},
                       {"r" + std::to_string(n) + "/2", reverse, std::string(reverse.size(), 'I')}};
            if (n % 2 == 1) std::swap(p.first, p.second);
            f.fragments.push_back(std::move(p));
            ++n;
        }
    }
}

// I4-E: alleles of one long segment between shared contexts; alleles_only: no flanks (no whole aggregate, fine only)
Fixture similar_alleles(bool alleles_only) {
    std::mt19937_64 rng(20261017);
    const std::string A = random_bases(rng, 1800), B = random_bases(rng, 1800), X = random_bases(rng, 4200);
    std::vector<std::string> alleles;
    for (const std::size_t at : {300u, 1500u, 2900u, 3800u}) {
        std::string a = X;
        a[at] = a[at] == 'A' ? 'C' : 'A';
        alleles.push_back(a);
    }
    alleles.push_back(X.substr(0, 2000) + X.substr(2030));
    alleles.push_back(X.substr(0, 1000) + random_bases(rng, 25) + X.substr(1000));
    alleles.push_back(X.substr(0, 2400) + X.substr(2000, 400) + X.substr(2400));   // duplicated copy: cross-copy mates
    std::string A2 = A;
    A2[1300] = A2[1300] == 'G' ? 'T' : 'G';
    Fixture f;
    f.name = alleles_only ? "similar alleles, no shared flanks" : "similar alleles between shared contexts";
    for (std::size_t i = 0; i < alleles.size(); ++i) {
        if (alleles_only) add(f, {alleles[i]});
        else add(f, {i == 2 ? A2 : A, alleles[i], B});
    }
    return f;
}

// I4-E: empty, consecutive empty, terminal empty and tiny segments; shifted copies of one segment; a novel junction
Fixture segments_fixture() {
    std::mt19937_64 rng(77);
    const std::string A = random_bases(rng, 1500), B = random_bases(rng, 1500), X = random_bases(rng, 3600);
    const std::string P1 = random_bases(rng, 900), P2 = random_bases(rng, 1300), t4 = random_bases(rng, 4), t7 = random_bases(rng, 7);
    const auto variant = [&](std::size_t at) {
        std::string v = X;
        v[at] = v[at] == 'A' ? 'G' : 'A';
        return v;
    };
    Fixture f;
    f.name = "empty, tiny, shifted and split segments";
    add(f, {A, "", X, B});
    add(f, {A, "", "", variant(1800), t4, B, ""});
    add(f, {P1, variant(600), t7, B});
    add(f, {"", P2, X.substr(0, 1700), X.substr(1700), B});   // novel junction inside X
    add(f, {A, variant(3000), "", "", B});
    return f;
}

// S2-E's fixture shapes: a repeated copy, a distant copy, indels, a junction, an empty block joined across, a short
// haplotype, an identical labelling; segments with empty and tiny pieces; uniform reads (the prefilter cache applies)
Fixture s2e_fixture() {
    std::mt19937_64 rng(20261018);
    const std::string a = random_bases(rng, 1800);
    std::string b = a;
    for (const std::size_t at : {200u, 650u, 1100u, 1500u}) b[at] = b[at] == 'A' ? 'C' : 'A';
    std::string e = a;
    e.insert(700, "GATTA");
    e.erase(1300, 3);
    const std::vector<std::string> bytes = {a, b, a.substr(0, 900) + a.substr(500, 300) + a.substr(900),
                                            a + random_bases(rng, 400) + a.substr(100, 350), e,
                                            a.substr(0, 1000) + b.substr(1000), a.substr(0, 600) + a.substr(680), a.substr(400, 520), b};
    Fixture f;
    f.name = "S2-E shapes (repeats, indels, junction, short, identical labelling)";
    f.uniform_reads = true;
    for (const std::string& s : bytes) {
        std::vector<std::string> pieces;
        std::size_t at = 0;
        for (const std::size_t cut : {0u, 300u, 7u, 0u, 400u, 13u, 600u}) {
            const std::size_t take = std::min<std::size_t>(cut, s.size() - at);
            pieces.push_back(s.substr(at, take));
            at += take;
        }
        pieces.push_back(s.substr(at));
        add(f, pieces);
    }
    return f;
}

struct Prepared {
    panvar::InternedHaplotypes interned;
    std::vector<HaplotypeCandidate> candidates;
    std::vector<std::vector<std::uint32_t>> segments;
};

Prepared prepare(const Fixture& f) {
    Prepared o;
    std::vector<panvar::MosaicPath> paths;
    for (std::size_t i = 0; i < f.haplotypes.size(); ++i) paths.push_back(panvar::MosaicPath{{i}});
    o.interned = panvar::intern_haplotypes(paths, f.haplotypes);
    for (std::size_t x = 0; x < o.interned.sequences.size(); ++x) {
        for (std::size_t y = x; y < o.interned.sequences.size(); ++y) o.candidates.push_back({x, y});
    }
    o.segments.assign(o.interned.sequences.size(), {});
    for (std::size_t i = 0; i < f.haplotypes.size(); ++i) o.segments[o.interned.index_of[i]] = f.segments[i];
    return o;
}

const std::size_t kStore = std::size_t{1} << 30;

std::vector<std::size_t> winners(const std::vector<DiplotypeScore>& scores, std::size_t nf, const panvar::FragmentScoringCounters& c) {
    const panvar::CertifiedRanking r = panvar::certify_scored_ranking(scores, nf, 6, c, panvar::CertificationOptions{});
    std::vector<std::size_t> w = r.winner_set;
    std::sort(w.begin(), w.end());
    return w;
}

BlockedUnionResult blocked(const Fixture& f, const Prepared& p, const panvar::LibraryModel& lib, const panvar::FragmentScoringOptions& o,
                           std::size_t block, bool fine) {
    return panvar::score_diplotypes_indexed_blocked(p.interned, p.candidates, f.fragments, lib, o, p.segments, kStore, true, block, {},
                                                    fine);
}

void require_scores(const std::vector<DiplotypeScore>& a, const std::vector<DiplotypeScore>& b, const std::string& label) {
    require(a.size() == b.size(), label + ": candidate count differs");
    for (std::size_t c = 0; c < a.size(); ++c) require(same(a[c], b[c]), label + ": candidate " + std::to_string(c) + " components differ");
}
void require_audits(const BlockedUnionResult& r, const panvar::FragmentScoringResult& bulk, const std::string& label) {
    require(r.result.audit_scores.size() == bulk.audit_scores.size(), label + ": audit view count differs");
    for (std::size_t k = 0; k < bulk.audit_scores.size(); ++k) require_scores(r.result.audit_scores[k], bulk.audit_scores[k], label + " (audit view)");
}

// per block: executed - constructed (whole + fine, discarded builds included) + served = the logical origins; and the
// executed count is the cap's observed counter; cumulatively the same identity against the reference's physical origins
void require_exact_once(const BlockedUnionResult& r, std::uint64_t reference_origins, const std::string& label) {
    std::uint64_t executed = 0, constructed = 0, served = 0;
    for (const panvar::BlockedBlockReport& b : r.blocks) {
        require(b.stage == "completed", label + ": completed block stage");
        executed += b.observed_physical_origins;
        constructed += b.origins_constructed + b.fine_origins_constructed;
        served += b.served_physical_origins;
        require(b.fine_served_physical_origins <= b.served_physical_origins, label + ": fine served is part of served");
    }
    require(executed == r.result.counters.physical_origins, label + ": per-block observed origins must sum to the executed counter");
    require(served == r.result.counters.reused_physical_origins, label + ": served work is counted apart from executed work");
    require(executed - constructed + served == reference_origins,
            label + ": an origin was lost or represented twice (executed " + std::to_string(executed) + " - constructed " +
                std::to_string(constructed) + " + served " + std::to_string(served) + " != reference " + std::to_string(reference_origins) + ")");
}

// per-block logical identity: each block's batch (baseline seeded) against the same block scored without origin reuse
void require_block_identities(const Fixture& f, const panvar::LibraryModel& lib, std::size_t block, std::size_t threads,
                              const panvar::FragmentScoringResult& bulk, std::size_t unique, const std::string& label) {
    std::size_t max_read = 1;
    for (const ReadPair& p : f.fragments) max_read = std::max({max_read, p.first.bases.size(), p.second.bases.size()});
    const Prepared p = prepare(f);
    panvar::ContentPlacementIndex index(max_read, lib.max_gap_bp, true);
    for (std::size_t h = 0; h < unique; ++h) index.add_haplotype(p.interned.sequences[h], p.segments[h]);
    panvar::UnionCumulativeCharges cumulative;
    const std::size_t nf = f.fragments.size();
    for (std::size_t begin = 0; begin < nf; begin += block) {
        const std::size_t end = std::min(nf, begin + block);
        const std::vector<ReadPair> part(f.fragments.begin() + begin, f.fragments.begin() + end);
        panvar::ContentIndexBatch plain(index, part, lib, options(threads, 0));
        plain.set_fine_views(true);
        plain.set_one_shot(true);
        plain.set_placement_store_budget(kStore);
        plain.discover();
        plain.build_columns();
        panvar::ContentIndexBatch batch(index, part, lib, options(threads, std::size_t{64} << 20));
        batch.set_fine_views(true);
        batch.set_fine_origin_aggregates(true);
        batch.set_one_shot(true);
        batch.set_placement_store_budget(kStore);
        batch.set_cumulative_baseline(cumulative);
        batch.discover();
        batch.build_columns();
        const panvar::ContentBatchReport& b = batch.report();
        const panvar::UnionCumulativeCharges own = batch.charges(), seen = batch.observed_charges();
        require(own.physical_origins == seen.physical_origins && own.piece_hits == seen.piece_hits &&
                    own.origins_consumed == seen.origins_consumed && own.states_consumed == seen.states_consumed &&
                    own.placements_consumed == seen.placements_consumed,
                label + ": a completed block's observed charges equal its charges");
        require(own.physical_origins - b.origins_constructed - b.fine_origins_constructed + b.served_physical_origins ==
                    plain.charges().physical_origins,
                label + ": block " + std::to_string(begin / block) + " exact-once origin identity");
        require(b.state_classes == plain.report().state_classes, label + ": block logical state classes");
        for (std::size_t h = 0; h < unique; ++h) {
            for (std::size_t i = 0; i < end - begin; ++i) {
                require(bits(batch.mass_columns()[h][i], bulk.log_mass[(begin + i) * unique + h]) &&
                            bits(plain.mass_columns()[h][i], bulk.log_mass[(begin + i) * unique + h]),
                        label + ": per-fragment mass of a block differs from bulk");
            }
        }
        // the acceptance caps' executed counters (the fields they test), as the blocked driver accumulates them
        cumulative.executed.physical_origins += own.executed.physical_origins;
        cumulative.executed.statistical_state_classes += own.executed.statistical_state_classes;
        cumulative.executed.latent_alignment_pairs += own.executed.latent_alignment_pairs;
        cumulative.executed.piece_hits += own.executed.piece_hits;
        cumulative.piece_hits += own.piece_hits;
        cumulative.physical_origins += own.physical_origins;
        cumulative.placements_consumed += own.placements_consumed;
        cumulative.origins_consumed += own.origins_consumed;
        cumulative.states_consumed += own.states_consumed;
        cumulative.prefilters_consumed += own.prefilters_consumed;
    }
}

void require_same_work(const BlockedUnionResult& a, const BlockedUnionResult& b, const std::string& label) {
    require(a.blocks.size() == b.blocks.size() && a.result.counters.physical_origins == b.result.counters.physical_origins &&
                a.result.counters.reused_physical_origins == b.result.counters.reused_physical_origins &&
                a.result.counters.merged_state_classes == b.result.counters.merged_state_classes,
            label + ": executed work differs");
    for (std::size_t k = 0; k < a.blocks.size(); ++k) {
        const panvar::BlockedBlockReport &x = a.blocks[k], &y = b.blocks[k];
        require(x.fine_windows_admitted == y.fine_windows_admitted && x.fine_windows_built == y.fine_windows_built &&
                    x.fine_charged_bytes == y.fine_charged_bytes && x.fine_zones_applied == y.fine_zones_applied &&
                    x.origin_entities_admitted == y.origin_entities_admitted && x.origin_held_bytes == y.origin_held_bytes &&
                    x.observed_physical_origins == y.observed_physical_origins,
                label + ": block " + std::to_string(k) + " admission or work differs");
    }
}

std::uint64_t sum_fine(const BlockedUnionResult& r, std::uint64_t panvar::BlockedBlockReport::*field) {
    std::uint64_t n = 0;
    for (const panvar::BlockedBlockReport& b : r.blocks) n += b.*field;
    return n;
}

// a cumulative cap reached in a later block: no scores; completed totals, the refused block's observed work and the
// whole attempt reported
void require_refusal(const BlockedUnionResult& cut, std::size_t blocks_before, const char* stage, const std::string& label) {
    require(!cut.completed && !cut.refusal.empty() && cut.result.scores.empty() && cut.result.audit_scores.empty(),
            label + ": a refusal returns no scores");
    require(cut.blocks.size() == blocks_before, label + ": refused in block " + std::to_string(cut.blocks.size() + 1) +
                                                    ", expected " + std::to_string(blocks_before + 1));
    require(cut.has_refused_block && cut.refused_block.stage == stage &&
                cut.refused_block.first_fragment == (cut.blocks.empty() ? 0 : cut.blocks.back().first_fragment + cut.blocks.back().fragments),
            label + ": the refused block is recorded with its stage (" + cut.refused_block.stage + ")");
    std::uint64_t completed_origins = 0, completed_hits = 0, completed_served = 0;
    for (const panvar::BlockedBlockReport& b : cut.blocks) {
        completed_origins += b.observed_physical_origins;
        completed_hits += b.observed_piece_hits;
        completed_served += b.origins_consumed;
    }
    const panvar::BlockedBlockReport& x = cut.refused_block;
    require(cut.completed_totals.physical_origins == completed_origins && cut.completed_totals.piece_hits == completed_hits &&
                cut.completed_totals.origins_consumed == completed_served,
            label + ": completed-block totals");
    require(cut.attempt_totals.physical_origins == completed_origins + x.observed_physical_origins &&
                cut.attempt_totals.piece_hits == completed_hits + x.observed_piece_hits &&
                cut.attempt_totals.origins_consumed == completed_served + x.origins_consumed &&
                cut.attempt_totals.placements_consumed == cut.completed_totals.placements_consumed + x.placements_consumed,
            label + ": whole-attempt totals = completed + refused block");
    require(x.cumulative_physical_origins == cut.attempt_totals.physical_origins, label + ": refused block's cumulative view");
}

void gates(const Fixture& base, double mean, double sd, std::uint64_t seed, bool boundary) {
    Fixture f = base;
    simulate(f, mean, 37, seed);
    const panvar::LibraryModel lib = library(mean, sd);
    const std::string label = f.name + ", insert " + std::to_string(static_cast<int>(mean));
    const Prepared p = prepare(f);
    const std::size_t nf = f.fragments.size(), unique = p.interned.sequences.size();
    const std::size_t big = std::size_t{64} << 20;
    const bool pre = f.uniform_reads;
    // reference and bulk (I4-E's comparison): masses, components, ranking
    const panvar::FragmentScoringResult ref = panvar::score_diplotypes(p.interned, p.candidates, f.fragments, lib, options(1, 0));
    panvar::IndexedUnionReport control_report, fine_report;
    const panvar::FragmentScoringResult control = panvar::score_diplotypes_indexed_union(
        p.interned, p.candidates, f.fragments, lib, options(2, big, nullptr, pre), p.segments, kStore, &control_report, nullptr, false, true, false);
    const panvar::FragmentScoringResult bulk_fine = panvar::score_diplotypes_indexed_union(
        p.interned, p.candidates, f.fragments, lib, options(2, big, nullptr, pre), p.segments, kStore, &fine_report, nullptr, false, true, true);
    for (std::size_t h = 0; h < unique; ++h) {
        require(bits(ref.exposure[h], control.exposure[h]) && bits(ref.exposure[h], bulk_fine.exposure[h]), label + ": exposure differs");
        for (std::size_t i = 0; i < nf; ++i) {
            require(bits(ref.mass_of(i, h), control.mass_of(i, h)) && bits(ref.mass_of(i, h), bulk_fine.mass_of(i, h)),
                    label + ": bulk mass differs from the reference");
        }
    }
    require_scores(ref.scores, control.scores, label + " (bulk control vs reference)");
    require_scores(ref.scores, bulk_fine.scores, label + " (bulk fine vs reference)");
    const std::vector<std::size_t> top = winners(ref.scores, nf, ref.counters);
    const std::uint64_t reference_origins = ref.counters.physical_origins;

    // blocked and blocked+fine against bulk (audit views on) and the reference's top set
    const panvar::FragmentScoringOptions audited = options(2, big, &lib, pre);
    const panvar::FragmentScoringResult bulk_audited = panvar::score_diplotypes_indexed_union(
        p.interned, p.candidates, f.fragments, lib, audited, p.segments, kStore, nullptr, nullptr, false, true, false);
    std::uint64_t fine_admitted = 0, fine_served = 0;
    for (const std::size_t block : {std::size_t{37}, nf - 1, nf}) {   // a short final block of one pair; the whole run
        const std::string at = label + " (block " + std::to_string(block) + ")";
        const BlockedUnionResult ordinary = blocked(f, p, lib, audited, block, false);
        const BlockedUnionResult combined = blocked(f, p, lib, audited, block, true);
        require(ordinary.completed && combined.completed, at + ": both blocked routes complete");
        require(combined.blocks.size() == (nf + block - 1) / block, at + ": every block processed");
        require_scores(combined.result.scores, bulk_audited.scores, at + " (combined vs bulk)");
        require_scores(ordinary.result.scores, combined.result.scores, at + " (combined vs ordinary)");
        require_audits(combined, bulk_audited, at + " (combined)");
        require_audits(ordinary, bulk_audited, at + " (ordinary)");
        require(winners(combined.result.scores, nf, combined.result.counters) == top &&
                    winners(ordinary.result.scores, nf, ordinary.result.counters) == top,
                at + ": certified top set differs from the reference's");
        require_exact_once(combined, reference_origins, at + " (combined)");
        require_exact_once(ordinary, reference_origins, at + " (ordinary)");
        require(sum_fine(ordinary, &panvar::BlockedBlockReport::fine_windows_built) == 0, at + ": ordinary blocking builds no fine window");
        for (const panvar::BlockedBlockReport& b : combined.blocks) {
            require(b.origin_whole_budget_bytes == big / 2 && b.origin_fine_budget_bytes == big - big / 2,
                    at + ": every block splits the unchanged total half/half");
            require(b.fine_charged_bytes <= b.origin_fine_budget_bytes && b.origin_held_bytes <= b.origin_whole_budget_bytes &&
                        b.fine_transition_peak_bytes >= b.fine_charged_bytes && b.fine_hash_collisions == 0,
                    at + ": block fine accounting");
        }
        fine_admitted += sum_fine(combined, &panvar::BlockedBlockReport::fine_windows_admitted);
        fine_served += sum_fine(combined, &panvar::BlockedBlockReport::fine_served_physical_origins);
        if (block == 37) {
            // 1 vs 4 threads; forced fine-key hash collisions
            const panvar::FragmentScoringOptions one = options(1, big, &lib, pre), four = options(4, big, &lib, pre);
            const BlockedUnionResult c1 = blocked(f, p, lib, one, block, true), c4 = blocked(f, p, lib, four, block, true);
            require_scores(c1.result.scores, c4.result.scores, at + " (1 vs 4 threads)");
            require_same_work(c1, c4, at + " (1 vs 4 threads)");
            panvar::set_fine_window_hash_collapse(true);
            const BlockedUnionResult collided = blocked(f, p, lib, four, block, true);
            panvar::set_fine_window_hash_collapse(false);
            require_scores(c4.result.scores, collided.result.scores, at + " (forced fine-key collisions)");
            require_same_work(c4, collided, at + " (forced fine-key collisions)");
            require(sum_fine(collided, &panvar::BlockedBlockReport::fine_hash_collisions) > 0, at + ": the collision control must collide");
            require_block_identities(f, lib, block, 2, control, unique, at);
            std::cout << "    " << at << ": " << combined.blocks.size() << " blocks; executed origins "
                      << ordinary.result.counters.physical_origins << " (ordinary) -> " << combined.result.counters.physical_origins
                      << " (combined); fine windows admitted " << sum_fine(combined, &panvar::BlockedBlockReport::fine_windows_admitted)
                      << ", fine served " << sum_fine(combined, &panvar::BlockedBlockReport::fine_served_physical_origins)
                      << ", fine built " << sum_fine(combined, &panvar::BlockedBlockReport::fine_origins_constructed) << '\n';
        }
    }
    // the I4-E fixtures must exercise fine windows; the S2-E shapes' short segments admit windows with no pairs in them
    require(f.uniform_reads || (fine_admitted > 0 && fine_served > 0), label + ": the fixture must exercise fine windows in blocks");

    // cumulative refusals in a later block (block 37): pairing/construction, served-origin consumption, discovery
    {
        const std::size_t block = 37;
        const panvar::FragmentScoringOptions o = options(2, big, &lib, pre);
        const BlockedUnionResult full = blocked(f, p, lib, o, block, true);
        require(full.blocks.size() >= 3, label + ": the refusal fixture needs three blocks");
        panvar::FragmentScoringOptions capped = o;
        capped.max_physical_origins = full.blocks[1].cumulative_physical_origins + 1;   // blocks 1-2 fit, block 3 does not
        const BlockedUnionResult cut = blocked(f, p, lib, capped, block, true);
        require_refusal(cut, 2, "columns", label + " (physical-origin cap)");
        require(cut.attempt_totals.physical_origins > capped.max_physical_origins &&
                    cut.refused_block.observed_physical_origins > 1 &&
                    cut.refused_block.observed_physical_origins <= full.blocks[2].observed_physical_origins + o.threads,
                label + ": the refused block's observed origins cross the cap and stop within the block");
        const BlockedUnionResult cut_ordinary = blocked(f, p, lib, capped, block, false);
        require(!cut_ordinary.completed && cut_ordinary.has_refused_block && cut_ordinary.result.scores.empty(),
                label + ": ordinary blocking refuses too, with its refused block recorded");
        // construction-time refusal: a cap inside block 3's aggregate construction (whole + fine construction come before
        // its pairing); the observed counter still includes the aborted construction
        if (full.blocks[2].origins_constructed + full.blocks[2].fine_origins_constructed > 2) {
            panvar::FragmentScoringOptions in_build = o;
            in_build.max_physical_origins = full.blocks[1].cumulative_physical_origins + 2;
            const BlockedUnionResult b = blocked(f, p, lib, in_build, block, true);
            require_refusal(b, 2, "columns", label + " (cap during construction)");
            require(b.refused_block.observed_physical_origins >= 2, label + ": aborted construction is observed work");
        }
        // pairing-time refusal: block 3's construction (whole + fine) completes, then its pairing crosses the cap
        {
            panvar::FragmentScoringOptions in_pairing = o;
            const std::uint64_t built = full.blocks[2].origins_constructed + full.blocks[2].fine_origins_constructed;
            require(full.blocks[2].observed_physical_origins > built + 1, label + ": block 3 pairs directly after construction");
            in_pairing.max_physical_origins = full.blocks[1].cumulative_physical_origins + built + 1;
            const BlockedUnionResult b = blocked(f, p, lib, in_pairing, block, true);
            require_refusal(b, 2, "columns", label + " (cap during pairing)");
            require(b.refused_block.origins_constructed + b.refused_block.fine_origins_constructed == built &&
                        b.refused_block.observed_physical_origins > built + 1 &&
                        b.refused_block.observed_physical_origins <= built + 1 + o.threads &&
                        b.refused_block.fine_windows_admitted == full.blocks[2].fine_windows_admitted,
                    label + ": a pairing-time refusal reports the completed construction and the pairing work observed");
        }
        panvar::FragmentScoringOptions served_cap = o;
        served_cap.max_union_cached_origins_consumed = full.blocks[0].origins_consumed + full.blocks[1].origins_consumed + 1;
        if (full.blocks[2].origins_consumed > 1) {
            const BlockedUnionResult s = blocked(f, p, lib, served_cap, block, true);
            require_refusal(s, 2, "columns", label + " (served-origin consumption cap)");
            require(s.attempt_totals.origins_consumed > served_cap.max_union_cached_origins_consumed,
                    label + ": the consumption cap's observed counter crosses it");
        }
        panvar::FragmentScoringOptions hit_cap = o;
        hit_cap.max_piece_hits = full.blocks[1].cumulative_piece_hits + 1;
        const BlockedUnionResult h = blocked(f, p, lib, hit_cap, block, true);
        require_refusal(h, 2, "discovery", label + " (piece-hit cap)");
        require(h.attempt_totals.piece_hits > hit_cap.max_piece_hits && h.refused_block.observed_physical_origins == 0,
                label + ": a discovery refusal observes piece hits and no pairing");
        std::cout << "    refusal (block 3 of " << full.blocks.size() << "): completed " << cut.completed_totals.physical_origins
                  << " + refused block " << cut.refused_block.observed_physical_origins << " = attempt "
                  << cut.attempt_totals.physical_origins << " origins (cap " << capped.max_physical_origins << "); piece-hit cap: "
                  << h.completed_totals.piece_hits << " + " << h.refused_block.observed_piece_hits << " hits\n";
    }

    if (!boundary) return;
    // exact budget boundaries inside a block (no whole aggregate here): block k's fine partition = exactly its full
    // admission's charge admits it all; one byte less ends its prefix one window earlier; a 1 B partition falls back
    // completely; other blocks admit partially; every variant identical
    const std::size_t block = 37;
    const BlockedUnionResult full = blocked(f, p, lib, options(2, big, &lib, pre), block, true);
    require(sum_fine(full, &panvar::BlockedBlockReport::origin_entities_admitted) == 0, label + ": the boundary fixture has no whole aggregate");
    std::size_t k = 0;
    for (std::size_t i = 0; i < full.blocks.size(); ++i) {
        if (full.blocks[i].fine_windows_admitted > 0 && (full.blocks[k].fine_windows_admitted == 0 ||
                                                          full.blocks[i].fine_charged_bytes < full.blocks[k].fine_charged_bytes)) k = i;
    }
    const std::uint64_t C = full.blocks[k].fine_charged_bytes;
    require(C > 0, label + ": a block with fine windows");
    const BlockedUnionResult exact = blocked(f, p, lib, options(2, 2 * C, &lib, pre), block, true);
    const BlockedUnionResult below = blocked(f, p, lib, options(2, 2 * C - 2, &lib, pre), block, true);
    const BlockedUnionResult none = blocked(f, p, lib, options(2, 2, &lib, pre), block, true);
    for (const BlockedUnionResult* r : {&exact, &below, &none}) {
        require(r->completed, label + ": boundary variant completes");
        require_scores(r->result.scores, full.result.scores, label + " (budget boundary)");
        require_exact_once(*r, reference_origins, label + " (budget boundary)");
    }
    require(exact.blocks[k].origin_fine_budget_bytes == C && exact.blocks[k].fine_windows_admitted == full.blocks[k].fine_windows_admitted &&
                exact.blocks[k].fine_charged_bytes == C,
            label + ": the exact partition admits every window of its block");
    require(below.blocks[k].fine_windows_admitted + 1 == full.blocks[k].fine_windows_admitted && below.blocks[k].fine_charged_bytes <= C - 1 &&
                below.blocks[k].fine_windows_declined >= 1,
            label + ": one byte below ends the block's prefix one window earlier");
    bool partial = false;
    for (std::size_t i = 0; i < full.blocks.size(); ++i) {
        partial = partial || (exact.blocks[i].fine_windows_admitted < full.blocks[i].fine_windows_admitted &&
                              exact.blocks[i].fine_windows_admitted > 0);
        require(none.blocks[i].fine_windows_admitted == 0 && none.blocks[i].fine_zones_applied == 0 &&
                    none.blocks[i].fine_windows_built <= 16 && none.blocks[i].fine_windows_declined == none.blocks[i].fine_windows_eligible,
                label + ": complete fallback in every block");
    }
    std::cout << "    budget boundary (block " << k + 1 << " of " << full.blocks.size() << "): exact " << C << " B admits "
              << exact.blocks[k].fine_windows_admitted << "; one byte below admits " << below.blocks[k].fine_windows_admitted
              << "; partial admission elsewhere: " << (partial ? "yes" : "no") << "; 1 B partition admits 0\n";
}

bool invalid_blocked(const Fixture& base, bool views, std::size_t origin_bytes) {
    Fixture f = base;
    simulate(f, 300, 97, 3);
    const Prepared p = prepare(f);
    try {
        (void)panvar::score_diplotypes_indexed_blocked(p.interned, p.candidates, f.fragments, library(300, 30), options(1, origin_bytes),
                                                       p.segments, kStore, views, 16, {}, true);
    } catch (const std::invalid_argument&) {
        return true;
    }
    return false;
}

}  // namespace

int main(int argc, char** argv) {
    const bool quick = argc > 1 && std::string(argv[1]) == "quick";
    try {
        if (quick) {
            gates(similar_alleles(true), 300.0, 30.0, 7, true);
            std::cout << "genotype blocked fine (quick): PASS\n";
            return 0;
        }
        for (const auto& [mean, sd] : {std::pair<double, double>{300.0, 30.0}, std::pair<double, double>{450.0, 50.0}}) {
            std::cout << "  library " << mean << " +- " << sd << '\n';
            gates(similar_alleles(false), mean, sd, 5, false);
            gates(similar_alleles(true), mean, sd, 7, true);
            gates(segments_fixture(), mean, sd, 9, false);
        }
        std::cout << "  library 300 +- 30, prefilter cache\n";
        gates(s2e_fixture(), 300.0, 30.0, 11, false);
        require(invalid_blocked(segments_fixture(), false, 1 << 20), "blocked fine aggregates without fine views must be refused");
        require(invalid_blocked(segments_fixture(), true, 0), "blocked fine aggregates without origin reuse must be refused");
        std::cout << "genotype blocked fine: PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "genotype blocked fine: FAIL: " << e.what() << '\n';
        return 1;
    }
}
