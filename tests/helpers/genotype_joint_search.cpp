// M3-D gates: the bounded joint-block pass (genotype_joint_search) on synthetic fixtures, then the caller integration.
// Library gates: incremental vs cold scores (length changes; both homologues; the second-step deltas are themselves
// gated in both orders by the probe's two-block-lik-gate), a two-block valley single moves cannot cross, no improvement,
// sequence-equivalent alternatives, tied winners, bound skips, a cap refusal, thread and request-order invariance, and
// cleanup of temporary slots and censuses.
#include "panvar/fragment_likelihood.hpp"
#include "panvar/genotype_joint_search.hpp"
#include "panvar/genotype_model.hpp"
#include "panvar/ls_hmm.hpp"
#include "panvar/mosaic_spelling.hpp"

#include <cmath>
#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

std::string rnd(std::uint64_t& s, std::size_t n) {
    std::string out;
    for (std::size_t i = 0; i < n; ++i) {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        out += "ACGT"[(s >> 33) % 4];
    }
    return out;
}
std::string revcomp(const std::string& x) {
    std::string r(x.rbegin(), x.rend());
    for (char& c : r) c = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : 'A';
    return r;
}
char other_base(char c) { return c == 'A' ? 'C' : c == 'C' ? 'G' : c == 'G' ? 'T' : 'A'; }

panvar::LibraryModel library() {
    panvar::LibraryModel lib;
    lib.insert_mean_bp = 300;
    lib.insert_sd_bp = 30;
    lib.support_sigmas = 4;
    lib.substitution_error = 0.002;
    lib.outlier_weight = 0.001;
    lib.fragment_start_intensity = 0.1;
    lib.provenance = "declared: M3-D fixture";
    return lib;
}

panvar::FragmentScoringOptions scoring(std::size_t threads) {
    panvar::FragmentScoringOptions o;
    o.threads = threads;
    o.mismatch_band = 6;
    return o;
}

// The caller's refinement configuration: origin-cache path, planned retention.
void configure(panvar::PreparedScoringContext& ctx) {
    ctx.enable_origin_cache_path(true);
    ctx.set_origin_cache_limits(256u << 20, 0);
}

panvar::ReadPair pair_at(const std::string& hap, std::size_t start, std::size_t insert, const std::string& id) {
    panvar::ReadPair p;
    p.first = {id + "/1", hap.substr(start, 100), std::string(100, 'I')};
    p.second = {id + "/2", revcomp(hap.substr(start + insert - 100, 100)), std::string(100, 'I')};
    return p;
}

// Background fragments tiling a haplotype every `step` bp, skipping any fragment that overlaps an excluded window.
void tile(std::vector<panvar::ReadPair>& out, const std::string& hap, std::size_t step, const std::vector<std::pair<std::size_t, std::size_t>>& avoid,
          const std::string& tag) {
    for (std::size_t s = 0; s + 300 <= hap.size(); s += step) {
        bool hit = false;
        for (const auto& w : avoid) hit = hit || (s < w.second && w.first < s + 300);
        if (!hit) out.push_back(pair_at(hap, s, 300, tag + std::to_string(s)));
    }
}

// A four-SV-block model: each block's two alleles are 60 bp differing at one base (position 30). Templates: `common`
// copies of all-reference (A), one all-alternate (B), plus any extra rows.
struct Fixture {
    panvar::PreparedGenotypeModel model;
    std::vector<std::string> invariant;
    std::vector<std::vector<std::string>> alleles;
};
Fixture snp_model(std::size_t common, const std::vector<std::vector<std::size_t>>& extra_rows, std::size_t blocks = 4,
                  const std::vector<int>& length_change = {}) {
    std::uint64_t s = 99;
    panvar::StructuralLocusInput in;
    in.locus_id = "m3d";
    std::vector<std::vector<std::string>> alleles;
    for (std::size_t b = 0; b <= blocks; ++b) in.invariant_segments.push_back(rnd(s, b == 0 || b == blocks ? 700 : 450));
    for (std::size_t b = 0; b < blocks; ++b) {
        std::string a0 = rnd(s, 60), a1 = a0;
        a1[30] = other_base(a1[30]);
        if (b < length_change.size() && length_change[b] > 0) a1 += rnd(s, static_cast<std::size_t>(length_change[b]));
        if (b < length_change.size() && length_change[b] < 0) a1 = a1.substr(0, a1.size() - static_cast<std::size_t>(-length_change[b]));
        in.blocks.push_back({"sv:" + std::to_string(b), {{"r" + std::to_string(b), a0}, {"x" + std::to_string(b), a1}}});
        alleles.push_back({a0, a1});
    }
    const auto row = [&](const std::vector<std::size_t>& a) {
        std::vector<std::string> ids;
        for (std::size_t b = 0; b < blocks; ++b) ids.push_back((a[b] == 0 ? "r" : "x") + std::to_string(b));
        return ids;
    };
    for (std::size_t k = 0; k < common; ++k) in.templates.push_back({"A" + std::to_string(k), panvar::TemplateFrameStatus::Complete, row(std::vector<std::size_t>(blocks, 0))});
    in.templates.push_back({"B", panvar::TemplateFrameStatus::Complete, row(std::vector<std::size_t>(blocks, 1))});
    for (std::size_t k = 0; k < extra_rows.size(); ++k) in.templates.push_back({"E" + std::to_string(k), panvar::TemplateFrameStatus::Complete, row(extra_rows[k])});
    return Fixture{panvar::prepare_genotype_model(in), in.invariant_segments, alleles};
}

// Canonical allele index of the fixture's allele (models canonicalise alleles by bytes).
std::vector<std::size_t> path_of(const Fixture& f, const std::vector<std::size_t>& which) {
    std::vector<std::size_t> p;
    for (std::size_t b = 0; b < which.size(); ++b) {
        const auto& block = f.model.blocks()[b];
        std::size_t idx = block.alleles.size();
        for (std::size_t a = 0; a < block.alleles.size(); ++a) if (block.alleles[a].sequence == f.alleles[b][which[b]]) idx = a;
        require(idx != block.alleles.size(), "fixture allele not found");
        p.push_back(idx);
    }
    return p;
}

std::string spell(const Fixture& f, const std::vector<std::size_t>& path) { return panvar::spell_mosaic(f.model, panvar::MosaicPath{path}).bases; }

// SNP coordinate of block b on a spelled haplotype.
std::size_t snp_at(const Fixture& f, const std::vector<std::size_t>& path, std::size_t b) {
    return panvar::spell_mosaic(f.model, panvar::MosaicPath{path}).block_intervals[b].begin + 30;
}

struct Run {
    panvar::JointSearchResult result;
    std::vector<std::tuple<std::vector<std::size_t>, std::vector<std::size_t>, double, double, double>> scored;   // paths, fragment, exposure, prior
    std::uint64_t transient_before = 0, live_before = 0;
    double incumbent_total = 0.0;
};

Run run_pass(const Fixture& f, const std::vector<panvar::ReadPair>& reads, const std::vector<std::size_t>& p0,
             const std::vector<std::size_t>& p1, std::size_t threads, const panvar::JointSearchOptions& opt,
             std::uint64_t cap_state_classes = 0) {
    panvar::FragmentScoringOptions so = scoring(threads);
    if (cap_state_classes) so.max_state_classes = cap_state_classes;
    panvar::PreparedScoringContext ctx(reads, library(), so);
    configure(ctx);
    ctx.begin_transient_round();
    const std::vector<std::size_t> slots = ctx.add_haplotypes({spell(f, p0), spell(f, p1)}, nullptr, false, true);
    const panvar::FragmentScoringResult inc = ctx.score_candidates({{std::min(slots[0], slots[1]), std::max(slots[0], slots[1])}});
    const std::vector<double> sw(f.model.blocks().size() - 1, opt.switch_probability);
    Run run;
    run.incumbent_total = inc.scores[0].log_likelihood + panvar::log_unordered_sequence_pair_prior(f.model, spell(f, p0), spell(f, p1), sw);
    run.transient_before = ctx.transient_payload_bytes();
    run.live_before = ctx.live_transient_slots().size();
    run.result = panvar::run_joint_block_search(
        ctx, f.model, p0, p1, slots[0], slots[1], run.incumbent_total, reads.size(), 6, opt,
        [&](const panvar::JointCandidate& c, const std::string&, const std::string&) {
            run.scored.push_back({c.first, c.second, c.score.fragment_term, c.score.exposure_term, c.prior});
        },
        {});
    if (!run.result.moved) {
        require(ctx.transient_payload_bytes() == run.transient_before, "cleanup: transient payload must return to its level");
        require(ctx.live_transient_slots().size() == run.live_before, "cleanup: no pass slot may stay live");
    } else {
        require(ctx.live_transient_slots().size() <= run.live_before + 2, "cleanup: only the winner's slots may stay live");
        const auto& w = run.result.catalogue[run.result.winner];
        require(ctx.holds_census(w.slot_first) && ctx.holds_census(w.slot_second), "the winner must keep censuses for the continuation");
    }
    require(run.result.full_placements == 0, "the pass must not place any haplotype in full");
    return run;
}

// The valley: incumbent all-reference homozygote; truth hap1 = (alt, alt, ref, ref). One read per SNP. Each single
// change pays a donor switch (start-frequency and r/n) its one read cannot pay for; the joint move's two reads can.
std::vector<panvar::ReadPair> valley_reads(const Fixture& f, const std::vector<std::size_t>& truth, const std::vector<std::size_t>& other) {
    const std::string t = spell(f, truth), o = spell(f, other);
    std::vector<std::pair<std::size_t, std::size_t>> avoid;
    for (std::size_t b = 0; b < f.model.blocks().size(); ++b) {
        const std::size_t x = snp_at(f, truth, b);
        avoid.push_back({x - 5, x + 6});
    }
    std::vector<panvar::ReadPair> reads;
    tile(reads, t, 37, avoid, "t");
    tile(reads, o, 41, avoid, "o");
    for (std::size_t b = 0; b < f.model.blocks().size(); ++b) {
        const std::size_t x = snp_at(f, truth, b);
        reads.push_back(pair_at(t, x - 50, 300, "snp" + std::to_string(b)));   // first mate covers the SNP once
    }
    return reads;
}

void gate_valley_cold_threads_order() {
    const Fixture f = snp_model(18, {});
    const auto A = path_of(f, {0, 0, 0, 0});
    const auto truth = path_of(f, {1, 1, 0, 0});
    const auto reads = valley_reads(f, truth, A);
    panvar::JointSearchOptions opt;
    const Run one = run_pass(f, reads, A, A, 1, opt);
    require(one.result.status == "complete", "valley: the pass must complete");
    require(one.result.decision == "moved", "valley: the joint move must win, got " + one.result.decision);
    const auto& w = one.result.catalogue[one.result.winner];
    require(std::minmax(w.first, w.second) == std::minmax(truth, A), "valley: the winner must be the truth pair");
    // both single-block intermediates lose to the incumbent: single-block refinement cannot cross
    for (const auto& item : one.scored) {
        const auto& p0 = std::get<0>(item);
        const auto& p1 = std::get<1>(item);
        const double total = std::get<2>(item) + std::get<3>(item) + std::get<4>(item);
        const bool single = (p0 == path_of(f, {1, 0, 0, 0}) || p1 == path_of(f, {1, 0, 0, 0}) || p0 == path_of(f, {0, 1, 0, 0}) ||
                             p1 == path_of(f, {0, 1, 0, 0}));
        if (single) require(total < one.incumbent_total, "valley: each single-block intermediate must lose");
    }
    // incremental vs cold, every catalogue member, in a fresh context (the pass's context is gone)
    {
        panvar::PreparedScoringContext cold(reads, library(), scoring(4));
        const std::vector<double> sw(3, opt.switch_probability);
        for (const auto& item : one.scored) {
            const std::string a = spell(f, std::get<0>(item)), b = spell(f, std::get<1>(item));
            const auto slots = cold.add_haplotypes({a, b});
            const auto r = cold.score_candidates({{std::min(slots[0], slots[1]), std::max(slots[0], slots[1])}});
            require(std::fabs(r.scores[0].fragment_term - std::get<2>(item)) <= 1e-9 &&
                        std::fabs(r.scores[0].exposure_term - std::get<3>(item)) <= 1e-9,
                    "incremental and cold scores must agree");
            require(panvar::log_unordered_sequence_pair_prior(f.model, a, b, sw) == std::get<4>(item), "the pass prior must be the sequence-level prior");
        }
    }
    // four threads and a permuted request order: identical catalogue (by identity), scores and decision
    panvar::JointSearchOptions permuted = opt;
    permuted.request_order_seed = 12345;
    const Run four = run_pass(f, reads, A, A, 4, permuted);
    require(four.result.decision == one.result.decision && four.scored.size() == one.scored.size(), "thread/order: same decision and size");
    std::map<std::pair<std::vector<std::size_t>, std::vector<std::size_t>>, std::tuple<double, double, double>> a, b;
    for (const auto& x : one.scored) a[{std::get<0>(x), std::get<1>(x)}] = {std::get<2>(x), std::get<3>(x), std::get<4>(x)};
    for (const auto& x : four.scored) b[{std::get<0>(x), std::get<1>(x)}] = {std::get<2>(x), std::get<3>(x), std::get<4>(x)};
    require(a == b, "thread/order: identical scores for every catalogue member");
    std::cout << "valley: moved to truth; " << one.scored.size() << " catalogue members; cold, thread and order invariance PASS\n";
}

void gate_no_improvement_and_length_changes() {
    // length-changing alternatives (+7, -5, +11, 0 bp): incremental vs cold; truth = incumbent, so no move
    const Fixture f = snp_model(6, {{1, 0, 1, 0}}, 4, {7, -5, 11, 0});
    const auto A = path_of(f, {0, 0, 0, 0});
    const auto other = path_of(f, {1, 0, 1, 0});
    std::vector<panvar::ReadPair> reads;
    tile(reads, spell(f, A), 23, {}, "a");
    tile(reads, spell(f, other), 29, {}, "b");
    panvar::JointSearchOptions opt;
    const Run r = run_pass(f, reads, A, other, 4, opt);
    require(r.result.status == "complete" && r.result.decision == "no_improvement",
            "no improvement: the truth incumbent must stay, got " + r.result.decision);
    panvar::PreparedScoringContext cold(reads, library(), scoring(4));
    for (const auto& item : r.scored) {
        const auto slots = cold.add_haplotypes({spell(f, std::get<0>(item)), spell(f, std::get<1>(item))});
        const auto s = cold.score_candidates({{std::min(slots[0], slots[1]), std::max(slots[0], slots[1])}});
        require(std::fabs(s.scores[0].fragment_term - std::get<2>(item)) <= 1e-9 && std::fabs(s.scores[0].exposure_term - std::get<3>(item)) <= 1e-9,
                "length changes: incremental and cold scores must agree");
    }
    std::cout << "no improvement + length changes (both homologues edited): " << r.scored.size() << " members agree with cold; PASS\n";
}

void gate_tie_skip_refusal_equivalence() {
    // tie: two alternatives at an unobserved block carried symmetrically, each completing an intact template with the
    // observed alternative at block 1 -- the two joint moves tie exactly at the top
    {
        panvar::StructuralLocusInput in;
        std::uint64_t s = 7;
        in.locus_id = "tie";
        in.invariant_segments = {rnd(s, 700), rnd(s, 450), rnd(s, 700)};
        const std::string a0 = rnd(s, 60), b0 = rnd(s, 60);
        std::string a1 = a0, a2 = a0, b1 = b0;
        a1[30] = other_base(a1[30]);
        a2[31] = other_base(a2[31]);
        b1[30] = other_base(b1[30]);
        in.blocks = {{"sv:0", {{"a0", a0}, {"a1", a1}, {"a2", a2}}}, {"sv:1", {{"b0", b0}, {"b1", b1}}}};
        for (int k = 0; k < 10; ++k) in.templates.push_back({"A" + std::to_string(k), panvar::TemplateFrameStatus::Complete, {"a0", "b0"}});
        in.templates.push_back({"T1", panvar::TemplateFrameStatus::Complete, {"a1", "b1"}});
        in.templates.push_back({"T2", panvar::TemplateFrameStatus::Complete, {"a2", "b1"}});
        Fixture f{panvar::prepare_genotype_model(in), in.invariant_segments, {{a0, a1, a2}, {b0, b1}}};
        const auto A = path_of(f, {0, 0});
        const auto T1 = path_of(f, {1, 1});
        const std::string t = spell(f, T1), o = spell(f, A);
        const std::size_t xa = snp_at(f, T1, 0), xb = snp_at(f, T1, 1);
        std::vector<panvar::ReadPair> reads;
        tile(reads, t, 37, {{xa - 300, xa + 300}, {xb - 5, xb + 6}}, "t");   // block 0 unobserved on the truth haplotype
        tile(reads, o, 41, {{xa - 300, xa + 300}, {xb - 5, xb + 6}}, "o");
        for (int k = 0; k < 3; ++k) reads.push_back(pair_at(t, xb - 40 - 7 * k, 300, "b1_" + std::to_string(k)));
        panvar::JointSearchOptions opt;
        const Run r = run_pass(f, reads, A, A, 4, opt);
        // the repair round's own rule: an exact top-set tie is uncertified (rescore_incomplete) or, if certified, tie_at_best
        require(r.result.status == "complete" && (r.result.decision == "rescore_incomplete" || r.result.decision == "tie_at_best") &&
                    r.result.winner_set_size == 2 && !r.result.incumbent_in_top_set && !r.result.moved,
                "tie: two tied improving joint moves must stop without choosing, got " + r.result.decision + " (winner set " +
                    std::to_string(r.result.winner_set_size) + ")");
        std::cout << "tie: " << r.result.decision << " with 2 tied improving winners, no move; PASS\n";
    }
    // bound skip: nothing scored, the context untouched
    {
        const Fixture f = snp_model(4, {});
        const auto A = path_of(f, {0, 0, 0, 0});
        std::vector<panvar::ReadPair> reads;
        tile(reads, spell(f, A), 31, {}, "a");
        panvar::JointSearchOptions opt;
        opt.single_limit = 7;   // 8 single alternatives exist (4 blocks x 2 homologues x 1 alternative)
        panvar::PreparedScoringContext ctx(reads, library(), scoring(1));
        configure(ctx);
        ctx.begin_transient_round();
        const auto slots = ctx.add_haplotypes({spell(f, A)}, nullptr, false, true);
        const std::uint64_t placed = ctx.full_placement_count(), bytes = ctx.transient_payload_bytes();
        const auto res = panvar::run_joint_block_search(ctx, f.model, A, A, slots[0], slots[0], 0.0, reads.size(), 6, opt, nullptr, {});
        require(res.status == "skipped_bounds" && res.decision == "skipped" && res.bounds.single_alternatives == 8 && !res.bounds.within,
                "bound skip must be explicit");
        require(ctx.full_placement_count() == placed && ctx.transient_payload_bytes() == bytes, "bound skip must touch nothing");
        std::cout << "bound skip: skipped_bounds (8 > 7), context untouched; PASS\n";
    }
    // cap refusal: a tiny cumulative cap is crossed inside the pass -> ResourceExhausted (the caller's incomplete contract)
    {
        const Fixture f = snp_model(18, {});
        const auto A = path_of(f, {0, 0, 0, 0});
        const auto reads = valley_reads(f, path_of(f, {1, 1, 0, 0}), A);
        // the cap is set one class above what placing the incumbent consumes, so only the pass itself can cross it
        std::uint64_t incumbent_classes = 0;
        {
            panvar::PreparedScoringContext probe(reads, library(), scoring(1));
            (void)probe.add_haplotypes({spell(f, A)}, nullptr, false, true);
            incumbent_classes = probe.cumulative_counters().statistical_state_classes;
        }
        bool refused = false;
        try {
            (void)run_pass(f, reads, A, A, 1, panvar::JointSearchOptions{}, incumbent_classes + 1);
        } catch (const panvar::ResourceExhausted&) {
            refused = true;
        }
        require(refused, "a cap crossed inside the pass must refuse (ResourceExhausted), never return a partial ranking");
        std::cout << "cap refusal: ResourceExhausted; PASS\n";
    }
    // sequence-equivalent alternatives: two allele paths spell the incumbent's bytes; never a candidate, never duplicated
    {
        panvar::StructuralLocusInput in;
        in.locus_id = "equivalent";
        std::uint64_t s = 11;
        const std::string left = rnd(s, 700), right = rnd(s, 700), core = rnd(s, 70);
        in.blocks = {{"b0", {{"a", "A" + core}, {"at", "AT" + core}}}, {"b1", {{"tg", "TG" + core}, {"g", "G" + core}}}};
        in.invariant_segments = {left, "", right};
        in.templates = {{"h0", panvar::TemplateFrameStatus::Complete, {"a", "tg"}}, {"h1", panvar::TemplateFrameStatus::Complete, {"at", "g"}},
                        {"h2", panvar::TemplateFrameStatus::Complete, {"a", "g"}}, {"h3", panvar::TemplateFrameStatus::Complete, {"at", "tg"}}};
        Fixture f{panvar::prepare_genotype_model(in), in.invariant_segments, {}};
        std::vector<std::size_t> inc;
        for (const auto& block : f.model.blocks()) {
            for (std::size_t a = 0; a < block.alleles.size(); ++a) {
                if (block.alleles[a].sequence == "A" + core || block.alleles[a].sequence == "TG" + core) inc.push_back(a);
            }
        }
        require(inc.size() == 2, "equivalence fixture");
        std::vector<panvar::ReadPair> reads;
        tile(reads, spell(f, inc), 19, {}, "e");
        const Run r = run_pass(f, reads, inc, inc, 1, panvar::JointSearchOptions{});
        std::set<std::pair<std::string, std::string>> seen;
        const std::string ib = spell(f, inc);
        for (const auto& item : r.scored) {
            std::string a = spell(f, std::get<0>(item)), b = spell(f, std::get<1>(item));
            if (b < a) std::swap(a, b);
            require(!(a == ib && b == ib), "a sequence-equivalent spelling of the incumbent must never be a candidate");
            require(seen.insert({a, b}).second, "catalogue members must be unique sequence pairs");
        }
        std::cout << "sequence-equivalent alternatives: " << r.scored.size() << " unique members, incumbent spelling excluded; PASS\n";
    }
}

}  // namespace

int main() {
    try {
        gate_valley_cold_threads_order();
        gate_no_improvement_and_length_changes();
        gate_tie_skip_refusal_equivalence();
        std::cout << "genotype joint search: PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "genotype joint search: FAIL: " << e.what() << '\n';
        return 1;
    }
}
