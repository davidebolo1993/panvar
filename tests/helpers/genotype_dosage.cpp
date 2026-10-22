// Fixtures for the dosage audit (docs/algorithms/genotype.md, stages 11, 12 and 14): the frozen Python random
// stream and normal quantile, every undefined-estimator case, the robust-z rules at zero spread, the top-genotype and
// allocation rules, nearest-panel ties and panel-order invariance, display-id collisions, and thread identity.

#include "panvar/genotype_dosage.hpp"
#include "panvar/genotype_dosage_integrate.hpp"
#include "panvar/genotype_headline.hpp"
#include "panvar/genotype_dosage_repair.hpp"
#include "panvar/genotype_repair.hpp"
#include "panvar/fragment_likelihood.hpp"
#include "panvar/mosaic_spelling.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <map>
#include <unistd.h>

namespace {

int failures = 0;
void require(bool ok, const std::string& what) {
    if (!ok) { ++failures; std::cerr << "FAIL: " << what << '\n'; }
}

using namespace panvar;
using namespace panvar::dosage;

std::string random_dna(std::mt19937_64& rng, std::size_t n) {
    static const char b[] = "ACGT";
    std::string s(n, 'A');
    for (char& c : s) c = b[rng() % 4];
    return s;
}

// ---------------------------------------------------------------------------------------------- frozen numerics

void test_python_random_and_quantile() {
    // values printed by CPython 3.14: random.Random(20260925).randrange(...), statistics.NormalDist().inv_cdf(...)
    PythonRandom r(20260925);
    const std::vector<std::uint32_t> a = {5397, 2129, 8358, 1921, 7072};
    for (std::uint32_t x : a) require(r.randrange(12345) == x, "randrange(12345) follows the Python stream");
    for (int i = 0; i < 3; ++i) require(r.randrange(1) == 0, "randrange(1) is 0 and consumes the stream as Python does");
    const std::vector<std::uint32_t> b = {359797802u, 2222143676u, 1136046544u, 3843596277u, 36158164u};
    for (std::uint32_t x : b) require(r.randrange(4000000000u) == x, "randrange(4e9) follows the Python stream");
    require(normal_inverse_cdf(0.9) == 1.2815515655446006, "inv_cdf(0.9) bit-identical");
    require(normal_inverse_cdf(1 - 0.001 / 10) == 3.7190164854557084, "z_critical(5) bit-identical");
    require(normal_inverse_cdf(1 - 0.001 / 46) == 4.0881736736581198, "z_critical(23) bit-identical");
    require(normal_inverse_cdf(1e-12) == -7.0344838253011321, "inv_cdf tail bit-identical");
    require(molecule_span_bp(350, 40, DosageSettings{}) == std::optional<std::size_t>(510), "C4 molecule span is 510");
    require(!molecule_span_bp(0, 0, DosageSettings{}), "a library without an insert distribution is refused");
    require(!molecule_span_bp(std::numeric_limits<double>::quiet_NaN(), 40, DosageSettings{}), "a NaN insert mean is refused");
}

void test_kmers() {
    const std::vector<KmerAt> f = canonical_kmers("ACGTTGCAAC", 5), r = canonical_kmers("GTTGCAACGT", 5);
    std::vector<Kmer> fk, rk;
    for (const KmerAt& x : f) fk.push_back(x.kmer);
    for (const KmerAt& x : r) rk.push_back(x.kmer);
    std::sort(fk.begin(), fk.end()); std::sort(rk.begin(), rk.end());
    require(fk == rk, "a sequence and its reverse complement give the same canonical k-mers");
    const std::vector<KmerAt> n = canonical_kmers("ACGTANACGTA", 5);
    require(n.size() == 2 && n[0].position == 0 && n[1].position == 6, "a k-mer with a base outside ACGT is skipped");
}

// ---------------------------------------------------------------------------------------------- synthetic locus

struct Locus {
    std::string ctx0, ctx2, x, insert;
    std::optional<PreparedGenotypeModel> model;
};

// ctx:0 | sv:1 (short = x, long = x + insert) | ctx:2; two templates, short and long.
Locus make_locus(std::size_t insert_bp, bool periodic_context, std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    Locus L;
    if (periodic_context) {
        const std::string p = random_dna(rng, 50);
        for (int i = 0; i < 40; ++i) { L.ctx0 += p; }
        const std::string q = random_dna(rng, 50);
        for (int i = 0; i < 40; ++i) { L.ctx2 += q; }
    } else {
        L.ctx0 = random_dna(rng, 2000);
        L.ctx2 = random_dna(rng, 2000);
    }
    L.x = random_dna(rng, 600);
    L.insert = random_dna(rng, insert_bp);
    StructuralLocusInput in;
    in.locus_id = "dosage-fixture";
    in.invariant_segments = {"", "", "", ""};
    in.blocks = {{"ctx:0", {{"c0", L.ctx0}}, GenotypeBlockKind::Context},
                 {"sv:1", {{"short", L.x}, {"long", L.x + L.insert}}, GenotypeBlockKind::Variant},
                 {"ctx:2", {{"c2", L.ctx2}}, GenotypeBlockKind::Context}};
    in.templates = {{"t-short", TemplateFrameStatus::Complete, {"c0", "short", "c2"}, ChainOrientation::Forward, ""},
                    {"t-long", TemplateFrameStatus::Complete, {"c0", "long", "c2"}, ChainOrientation::Forward, ""}};
    L.model.emplace(prepare_genotype_model(in));
    return L;
}

std::size_t allele_id(const PreparedGenotypeModel& m, std::size_t block, const std::string& source) {
    return m.blocks()[block].allele_id_for_source(source);
}

std::vector<CensusBlock> census_of(const PreparedGenotypeModel& m) {
    std::vector<CensusBlock> c(3);
    for (std::size_t b = 0; b < 3; ++b) { c[b].block = b; c[b].block_id = m.blocks()[b].block_id; }
    c[1].route = CensusRoute::Pair;
    return c;
}

ReadPair pair_of(const std::string& a, const std::string& b) {
    ReadPair p;
    p.first.bases = a;
    p.second.bases = b;
    return p;
}

DosageAudit audit_of(const Locus& L, const std::string& a1, const std::string& a2, const std::vector<ReadPair>& reads,
                     std::size_t threads = 1) {
    const PreparedGenotypeModel& m = *L.model;
    const std::vector<CensusBlock> census = census_of(m);
    std::vector<BlockFeatureModel> features(3);
    features[1] = pooled_features(m.blocks()[1], DosageSettings{});
    const MosaicPath h1{{0, allele_id(m, 1, a1), 0}}, h2{{0, allele_id(m, 1, a2), 0}};
    return audit_dosage(m, census, features, reported_call(m, h1, h2), reads, 510, DosageSettings{}, threads);
}

void test_undefined_estimator_cases() {
    // no tested feature: the alleles differ by fewer than 500 k-mers; m = 0, so no z_critical
    {
        const Locus L = make_locus(200, false, 1);
        const DosageAudit a = audit_of(L, "short", "long", {pair_of(L.ctx0.substr(0, 100), L.ctx2.substr(0, 100))});
        require(a.blocks[1].status == DosageStatus::NotEvaluated && a.blocks[1].reason == "no_tested_feature", "no_tested_feature");
        require(a.m == 0 && !a.z_critical, "m = 0 computes no z_critical");
    }
    // no reference: periodic context blocks have no single-copy k-mer
    {
        const Locus L = make_locus(1500, true, 2);
        const DosageAudit a = audit_of(L, "short", "long", {pair_of(L.ctx0.substr(0, 100), L.insert.substr(0, 100))});
        require(a.blocks[1].status == DosageStatus::NotEvaluated && a.blocks[1].reason == "no_reference", "no_reference");
    }
    // no reference hits: reads hit the feature only
    {
        const Locus L = make_locus(1500, false, 3);
        std::vector<ReadPair> reads;
        for (int i = 0; i < 20; ++i) reads.push_back(pair_of(L.insert.substr(10 * i, 100), L.insert.substr(700 + 10 * i, 100)));
        const DosageAudit a = audit_of(L, "short", "long", reads);
        require(a.blocks[1].status == DosageStatus::NotEvaluated && a.blocks[1].reason == "no_reference_hits", "no_reference_hits");
    }
    // zero reference in a bootstrap replicate: one pair of 40 carries every reference hit
    {
        const Locus L = make_locus(1500, false, 4);
        std::vector<ReadPair> reads;
        reads.push_back(pair_of(L.ctx0.substr(100, 100), L.ctx2.substr(100, 100)));
        for (int i = 0; i < 39; ++i) reads.push_back(pair_of(L.insert.substr(10 * i, 100), L.insert.substr(700 + 10 * i, 100)));
        const DosageAudit a = audit_of(L, "short", "long", reads);
        require(a.blocks[1].status == DosageStatus::NotEvaluated && a.blocks[1].reason == "bootstrap_zero_reference",
                "bootstrap_zero_reference fails the block closed");
        require(a.blocks[1].features.empty(), "no replicate value reaches a quantile after a zero reference");
    }
}

void test_zero_spread() {
    // every pair hits the reference only: every replicate estimates the feature at exactly 0, so the IQR is 0
    const Locus L = make_locus(1500, false, 5);
    std::vector<ReadPair> reads;
    for (int i = 0; i < 30; ++i) reads.push_back(pair_of(L.ctx0.substr(20 * i, 100), L.ctx2.substr(20 * i, 100)));
    const DosageAudit hom = audit_of(L, "short", "short", reads);
    require(!hom.blocks[1].features.empty() && hom.blocks[1].features[0].robust_sd == 0.0, "zero spread");
    require(hom.blocks[1].features[0].robust_z == 0.0, "robust z is 0 when the call is at the median");
    const DosageAudit het = audit_of(L, "short", "long", reads);
    require(std::isinf(het.blocks[1].features[0].robust_z), "robust z is +inf when the call is off the median, by rule");
    require(het.blocks[1].status == DosageStatus::AutoRepair, "a call off a zero-spread median outside its interval triggers");
}

// I2-H's disclosed gap, closed: the fixed-geometry recursion must actually EXECUTE inside dosage, not merely leave the
// dosage tables unchanged. The zero-spread fixture's heterozygous call triggers an AUTO repair; repair_block scores
// its candidates exactly through the dosage stage's own scoring lifecycle. Single-threaded, so the existing per-thread
// execution hook, read around repair_block alone, sees every recursion of dosage's repair scoring. Both backends must
// choose the same repair from bit-identical candidate scores.
void test_fixed_recursion_inside_dosage() {
    const Locus L = make_locus(1500, false, 5);
    const PreparedGenotypeModel& m = *L.model;
    std::vector<ReadPair> reads;
    for (int i = 0; i < 30; ++i) {
        ReadPair p = pair_of(L.ctx0.substr(20 * i, 100), L.ctx2.substr(20 * i, 100));
        p.first.name = "r" + std::to_string(i) + "/1";
        p.second.name = "r" + std::to_string(i) + "/2";
        p.first.qualities = p.second.qualities = std::string(100, 'I');
        reads.push_back(p);
    }
    const DosageAudit audit = audit_of(L, "short", "long", reads);
    require(audit.blocks[1].status == DosageStatus::AutoRepair, "fixture: the call triggers an AUTO repair");
    const std::vector<CensusBlock> census = census_of(m);
    std::vector<BlockFeatureModel> features(3);
    features[1] = pooled_features(m.blocks()[1], DosageSettings{});
    const MosaicPath h1{{0, allele_id(m, 1, "short"), 0}}, h2{{0, allele_id(m, 1, "long"), 0}};
    const ReportedCall call = reported_call(m, h1, h2);
    LibraryModel library;
    library.insert_mean_bp = 350.0;
    library.insert_sd_bp = 40.0;
    library.support_sigmas = 4.0;
    library.substitution_error = 0.002;
    library.outlier_weight = 0.001;
    library.fragment_start_intensity = 0.05;
    library.provenance = "fixture: declared";
    RepairResult result[2];
    std::uint64_t executions[2] = {0, 0};
    for (int arm = 0; arm < 2; ++arm) {
        ScoringInputs in;
        in.model = &m;
        in.first = h1;
        in.second = h2;
        in.call = &call;
        in.reads = &reads;
        in.library = &library;
        in.options.threads = 1;
        in.options.mismatch_band = 6;
        in.options.minimum_piece_bp = 10;
        in.options.gap_recursion_backend = arm == 0 ? GapRecursionBackend::scaled_reference
                                                    : GapRecursionBackend::fixed_geometry_experimental;
        in.switch_probability = std::vector<double>(m.blocks().size() - 1, 0.01);
        double prior = 0.0;
        in.persisted = score_reported_diplotype(in, prior);   // what the run persisted for its reported call
        const std::uint64_t before = gapped_dp_test::fixed_geometry_executions();
        result[arm] = repair_block(audit.blocks[1], census[1], features[1], {}, in, DosageSettings{}, RepairSettings{});
        executions[arm] = gapped_dp_test::fixed_geometry_executions() - before;
    }
    require(result[1].outcome == "REPAIRED" && !result[1].candidates.empty(),
            "the repair must reproduce the reported call and score candidates (vacuous otherwise)");
    require(executions[0] == 0, "the reference backend runs no fixed-geometry recursion in dosage");
    require(executions[1] > 0, "the fixed-geometry recursion must execute inside dosage's repair scoring");
    const auto bits = [](double x) { std::uint64_t u; std::memcpy(&u, &x, sizeof u); return u; };
    require(result[0].outcome == result[1].outcome && result[0].chosen == result[1].chosen &&
                result[0].candidates.size() == result[1].candidates.size() && result[0].top_set_size == result[1].top_set_size,
            "both backends must reach the same repair");
    for (std::size_t c = 0; c < result[0].candidates.size(); ++c) {
        const BlockCandidate& a = result[0].candidates[c];
        const BlockCandidate& b = result[1].candidates[c];
        require(a.label == b.label && bits(a.fragment_term) == bits(b.fragment_term) &&
                    bits(a.exposure_term) == bits(b.exposure_term) && bits(a.likelihood) == bits(b.likelihood) &&
                    a.chosen == b.chosen && a.in_top_set == b.in_top_set,
                "every repair candidate's exact score must agree bit for bit");
    }
    std::cout << "  fixed-geometry recursion inside dosage: " << executions[1] << " executed in repair scoring of "
              << result[1].candidates.size() << " candidates; candidate scores and the repair bit-identical to the "
                 "reference\n";
}

// ---- I2-J: the interned census backend must equal the reference in every CensusBlock field ----------------------------
// A census outcome, or the exact refusal message.
struct CensusOutcome {
    std::vector<CensusBlock> rows;
    std::string error;
};

CensusOutcome census_with(const PreparedGenotypeModel& m, const NodeSequences& nodes, std::size_t span,
                          const DosageSettings& s, CensusBackend backend) {
    CensusOutcome out;
    try {
        out.rows = run_census(m, nodes, span, s, backend);
    } catch (const std::exception& e) {
        out.error = e.what();
    }
    return out;
}

bool same_census(const CensusOutcome& a, const CensusOutcome& b) {
    if (a.error != b.error || a.rows.size() != b.rows.size()) return false;
    for (std::size_t i = 0; i < a.rows.size(); ++i) {
        const CensusBlock& x = a.rows[i];
        const CensusBlock& y = b.rows[i];
        if (x.block != y.block || x.block_id != y.block_id || x.alleles != y.alleles ||
            x.max_node_multiplicity != y.max_node_multiplicity || x.dosage_only_pairs != y.dosage_only_pairs ||
            x.max_repeated_bp != y.max_repeated_bp || x.route != y.route) {
            return false;
        }
    }
    return true;
}

// One block of alleles, each a walk (empty: no walk source) with extra source walks aliasing the same sequence.
struct WalkAllele {
    std::string walk;                    // "a+,b-,..." or "" for an allele with no walk source
    std::vector<std::string> aliases;    // further source walks spelling the same sequence (the first walk is primary)
};

PreparedGenotypeModel walk_model(const std::vector<std::vector<WalkAllele>>& blocks, std::mt19937_64& rng) {
    StructuralLocusInput in;
    in.locus_id = "census-gate";
    std::size_t widest = 0;
    for (std::size_t b = 0; b < blocks.size(); ++b) {
        GenotypeBlockInput block;
        block.block_id = "sv:" + std::to_string(b);
        block.kind = GenotypeBlockKind::Variant;
        for (std::size_t a = 0; a < blocks[b].size(); ++a) {
            const std::string sequence = random_dna(rng, 8) + std::to_string(b) + "_" + std::to_string(a);   // distinct
            const std::string seq = sequence.substr(0, 8) + std::string(a + 1, "ACGT"[b % 4]);
            // A unique tag before "walk:" keeps duplicate walks apart as source ids; primary_walk reads what follows it.
            const std::string source = blocks[b][a].walk.empty() ? "nowalk" + std::to_string(a)
                                                                 : "a" + std::to_string(a) + "|walk:" + blocks[b][a].walk;
            block.alleles.push_back({source, seq});
            for (std::size_t k = 0; k < blocks[b][a].aliases.size(); ++k) {
                block.alleles.push_back({"a" + std::to_string(a) + "x" + std::to_string(k) + "|walk:" + blocks[b][a].aliases[k], seq});
            }
        }
        widest = std::max(widest, blocks[b].size());
        in.blocks.push_back(block);
    }
    in.invariant_segments.assign(blocks.size() + 1, "");
    for (std::size_t t = 0; t < widest; ++t) {
        PanelTemplateInput row;
        row.name = "t" + std::to_string(t);
        row.frame_status = TemplateFrameStatus::Complete;
        for (std::size_t b = 0; b < blocks.size(); ++b) {
            const std::size_t a = t % blocks[b].size();
            const WalkAllele& allele = blocks[b][a];
            row.allele_source_ids.push_back(allele.walk.empty() ? "nowalk" + std::to_string(a)
                                                                : "a" + std::to_string(a) + "|walk:" + allele.walk);
        }
        in.templates.push_back(row);
    }
    return prepare_genotype_model(in);
}

void test_census_backends_agree() {
    std::size_t compared = 0, errors = 0, routes[3] = {0, 0, 0}, boundary = 0;
    const auto check = [&](const PreparedGenotypeModel& m, const NodeSequences& nodes, std::size_t span,
                           const DosageSettings& s, const std::string& label) {
        const CensusOutcome ref = census_with(m, nodes, span, s, CensusBackend::reference);
        const CensusOutcome exp = census_with(m, nodes, span, s, CensusBackend::interned_experimental);
        require(same_census(ref, exp), "census backends differ: " + label + (ref.error.empty() ? "" : " (" + ref.error + ")"));
        ++compared;
        if (!ref.error.empty()) ++errors;
        for (const CensusBlock& c : ref.rows) ++routes[static_cast<int>(c.route)];
    };
    // existing fixtures: the review regression's repeat models and the dosage locus fixtures (no walks: no alleles)
    {
        const NodeSequences nodes = {{"x", std::string(600, 'A')}, {"y", std::string(600, 'A')}, {"z", "CCC"}};
        for (const bool unrolled : {false, true}) {
            StructuralLocusInput x;
            x.locus_id = "repeat";
            x.blocks = {{"sv:0", {{"walk:x+", std::string(600, 'A')}, {unrolled ? "walk:x+,y+" : "walk:x+,x+", std::string(1200, 'A')}},
                         GenotypeBlockKind::Variant},
                        {"ctx:1", {{"walk:z+", "CCC"}}, GenotypeBlockKind::Context}};
            x.invariant_segments = {"", "", ""};
            for (const auto& [name, first] : {std::pair<std::string, std::string>{"short", "walk:x+"},
                                              {"long", unrolled ? "walk:x+,y+" : "walk:x+,x+"}}) {
                PanelTemplateInput row;
                row.name = name;
                row.frame_status = TemplateFrameStatus::Complete;
                row.allele_source_ids = {first, "walk:z+"};
                x.templates.push_back(row);
            }
            const PreparedGenotypeModel m = prepare_genotype_model(x);
            for (const std::size_t span : {0, 499, 500, 501, 600, 601}) check(m, nodes, span, DosageSettings{}, "repeat model");
        }
        for (const std::uint64_t seed : {1, 2, 5, 12}) {
            const Locus L = make_locus(1500, seed == 2, seed);
            check(*L.model, {}, 510, DosageSettings{}, "dosage locus fixture");
        }
    }
    // randomized and adversarial walk models
    std::mt19937_64 rng(20261016);
    const std::vector<std::string> names = {"1", "2", "10", "s1", "node:7", "A+B", "x_y", "Δ", "seg-3", "11", "a.b", "Z"};
    for (int trial = 0; trial < 300; ++trial) {
        const std::size_t pool = 3 + rng() % (names.size() - 2);   // 3..12 names
        std::vector<std::vector<WalkAllele>> blocks(1 + rng() % 3);
        for (auto& block : blocks) {
            const std::size_t n = 1 + rng() % 7;
            for (std::size_t a = 0; a < n; ++a) {
                WalkAllele allele;
                const std::size_t steps = rng() % 9;   // 0: an allele with no walk
                for (std::size_t k = 0; k < steps; ++k) {
                    if (k > 0) allele.walk += ",";
                    allele.walk += names[rng() % pool] + (rng() % 2 ? "+" : "-");   // both orientations, repeats
                }
                if (a > 0 && rng() % 5 == 0) allele.walk = block[rng() % a].walk;   // a duplicate walk
                if (!allele.walk.empty() && rng() % 4 == 0) allele.aliases.push_back(names[rng() % pool] + "+");
                block.push_back(allele);
            }
        }
        PreparedGenotypeModel m = walk_model(blocks, rng);
        NodeSequences nodes;
        for (std::size_t i = 0; i < pool; ++i) nodes[names[i]] = std::string(1 + rng() % 120, "ACGT"[i % 4]);
        const bool drop = trial % 17 == 0;   // a missing node
        if (drop) nodes.erase(names[rng() % pool]);
        DosageSettings s;
        s.max_novel_bp = rng() % 3 == 0 ? 0 : rng() % 150;
        for (const std::size_t span : {std::size_t{0}, std::size_t{1 + rng() % 200}, std::size_t{510}}) {
            check(m, nodes, span, s, "random trial " + std::to_string(trial));
        }
        // thresholds immediately below, at and above a boundary this model actually has: the largest length gap and
        // the absent-node sum of one pair, computed here independently of both backends
        if (!drop) {
            for (const auto& block : m.blocks()) {
                std::vector<std::map<std::string, std::size_t>> counts;
                for (const auto& allele : block.alleles) {
                    const std::vector<NodeStep> walk = primary_walk(allele);
                    if (walk.empty()) continue;
                    std::map<std::string, std::size_t> c;
                    for (const NodeStep& st : walk) ++c[st.node];
                    counts.push_back(c);
                }
                if (counts.size() < 2) continue;
                std::size_t li = 0, lj = 0, novel = 0;
                for (const auto& [node, v] : counts[0]) li += v * nodes.at(node).size();
                for (const auto& [node, v] : counts[1]) lj += v * nodes.at(node).size();
                for (const auto& [node, v] : counts[0]) if (!counts[1].count(node)) novel += nodes.at(node).size();
                const std::size_t gap = lj > li ? lj - li : li - lj;
                for (const std::size_t span : {gap > 0 ? gap - 1 : 0, gap, gap + 1}) {
                    for (const std::size_t t : {novel > 0 ? novel - 1 : 0, novel, novel + 1}) {
                        DosageSettings at;
                        at.max_novel_bp = t;
                        check(m, nodes, span, at, "boundary, trial " + std::to_string(trial));
                        ++boundary;
                    }
                }
                break;
            }
        }
        // reordered inputs: the same alleles in reverse order
        std::vector<std::vector<WalkAllele>> reversed = blocks;
        for (auto& block : reversed) std::reverse(block.begin(), block.end());
        check(walk_model(reversed, rng), nodes, 510, s, "reordered trial " + std::to_string(trial));
    }
    // a malformed step in a later source walk refuses on both (primary_walk parses every source walk)
    {
        StructuralLocusInput in;
        in.locus_id = "malformed";
        GenotypeBlockInput block;
        block.block_id = "sv:0";
        block.alleles = {{"walk:s1+,s1+", "AAAA"}, {"walk:s1", "AAAA"}, {"walk:s1+", "CCCC"}};
        in.blocks = {block};
        in.invariant_segments = {"", ""};
        for (const std::string& src : {std::string("walk:s1+,s1+"), std::string("walk:s1+")}) {
            PanelTemplateInput row;
            row.name = src;
            row.frame_status = TemplateFrameStatus::Complete;
            row.allele_source_ids = {src};
            in.templates.push_back(row);
        }
        check(prepare_genotype_model(in), {{"s1", "ACGTACGT"}}, 0, DosageSettings{}, "malformed alias walk");
    }
    require(errors > 0 && routes[0] > 0 && routes[1] > 0 && routes[2] > 0 && boundary > 0,
            "the census gate must exercise refusals, every route and boundary thresholds");
    std::cout << "  census backends agree on every CensusBlock field in " << compared << " censuses (" << errors
              << " identical refusals; routes not_sensitive/pair/repeated_unit " << routes[0] << "/" << routes[1] << "/"
              << routes[2] << "; " << boundary << " boundary-threshold cases)\n";
}

void test_thread_identity() {
    const Locus L = make_locus(1500, false, 6);
    std::mt19937_64 rng(7);
    const std::string long_hap = L.ctx0 + L.x + L.insert + L.ctx2, short_hap = L.ctx0 + L.x + L.ctx2;
    std::vector<ReadPair> reads;
    for (int i = 0; i < 400; ++i) {
        const std::string& h = i % 2 ? long_hap : short_hap;
        const std::size_t at = rng() % (h.size() - 350);
        reads.push_back(pair_of(h.substr(at, 100), h.substr(at + 250, 100)));
    }
    const std::string one = dosage_features_table(audit_of(L, "short", "long", reads, 1)) + dosage_blocks_table(audit_of(L, "short", "long", reads, 1));
    for (std::size_t t : {2, 4, 8}) {
        const std::string other = dosage_features_table(audit_of(L, "short", "long", reads, t)) + dosage_blocks_table(audit_of(L, "short", "long", reads, t));
        require(one == other, "the audit tables are byte-identical at " + std::to_string(t) + " threads");
    }
}

void test_unsupported_library() {
    const Locus L = make_locus(1500, false, 8);
    LibraryModel lib;   // no insert distribution
    const std::string out = (std::filesystem::temp_directory_path() / "panvar_dosage_fixture").string();
    const DosageShadowSummary s = run_dosage_shadow(*L.model, {}, {}, lib, FragmentScoringOptions{}, {0.01, 0.01}, out, out, DosageShadowOptions{});
    require(s.audit_status == "REFUSED_LIBRARY", "unsupported_library refuses the audit");
    {
        std::ifstream in(out + ".dosage_blocks.tsv");
        std::string line;
        std::size_t rows = 0, refused = 0;
        std::getline(in, line);
        while (std::getline(in, line)) { ++rows; refused += line.find("\tNOT_EVALUATED\tunsupported_library") != std::string::npos; }
        require(rows == 3 && refused == 3, "every block is NOT_EVALUATED / unsupported_library");
        std::filesystem::remove(out + ".dosage_blocks.tsv");
    }
    LibraryModel ok;
    ok.insert_mean_bp = 350; ok.insert_sd_bp = 40;
    const DosageShadowSummary n = run_dosage_shadow(*L.model, {}, {}, ok, FragmentScoringOptions{}, {0.01, 0.01}, out + "_absent", out, DosageShadowOptions{});
    require(n.audit_status == "NOT_RUN", "a run without a persisted call is not audited");
    std::filesystem::remove(out + ".dosage_status.tsv");
}

// ---------------------------------------------------------------------------------------------- choice rules

BlockCandidate cand(const std::string& a, const std::string& b, double ll, const std::string& orient = "AB",
                    std::vector<std::size_t> ms = {}, double prior = 0.0, std::size_t rank = 0) {
    BlockCandidate c;
    c.ctx1_sequence = a; c.ctx2_sequence = b; c.likelihood = ll; c.orientation = orient; c.multiset = std::move(ms);
    c.panel_log_prior = prior; c.rank = rank; c.label = "c" + std::to_string(rank);
    return c;
}

void test_choice_rules() {
    {   // distinct genotypes tied at the top: keep the call, report every tied candidate
        RepairResult r;
        r.candidates = {cand("X", "Y", -10.0), cand("X", "Z", -10.0 + 5e-10), cand("X", "W", -12.0)};
        choose_repair(r, CensusRoute::Pair, 1e-9);
        require(r.outcome == "TRIGGERED_UNREPAIRED" && r.reason == "top_genotype_tie" && r.chosen == -1, "distinct tied genotypes refuse");
        require(r.candidates[0].in_top_set && r.candidates[1].in_top_set && !r.candidates[2].in_top_set, "every tied candidate is reported");
    }
    {   // duplicate spelled genotype: merged before tie handling
        RepairResult r;
        r.candidates = {cand("X", "Y", -10.0), cand("X", "Y", -10.0), cand("X", "Z", -11.0)};
        choose_repair(r, CensusRoute::Pair, 1e-9);
        require(r.outcome == "REPAIRED" && r.top_set_size == 1, "one spelled genotype under two labels is not a tie");
    }
    {   // orientation-only tie: AB is taken and the tie recorded
        RepairResult r;
        r.candidates = {cand("Y", "X", -10.0, "BA"), cand("X", "Y", -10.0, "AB")};
        choose_repair(r, CensusRoute::Pair, 1e-9);
        require(r.outcome == "REPAIRED" && r.orientation_tie && r.candidates[static_cast<std::size_t>(r.chosen)].orientation == "AB",
                "an orientation tie takes AB");
    }
    {   // same content multiset: the panel-structure prior chooses the allocation, even below the top likelihood
        RepairResult r;
        r.candidates = {cand("A", "B", -10.0, "AB", {1, 1, 2}, -3.0, 0), cand("C", "D", -10.5, "AB", {1, 1, 2}, -1.0, 1),
                        cand("E", "F", -20.0, "AB", {1, 2, 2}, 0.0, 2)};
        choose_repair(r, CensusRoute::RepeatedUnit, 1e-9);
        require(r.outcome == "REPAIRED" && r.chosen == 1 && r.allocation_tie == 1, "the prior chooses within the content multiset");
    }
    {   // exact (prior, likelihood) tie inside one multiset: recorded as an allocation tie
        RepairResult r;
        r.candidates = {cand("A", "B", -10.0, "AB", {1, 1, 2}, -1.0, 3), cand("C", "D", -10.0, "AB", {1, 1, 2}, -1.0, 2)};
        choose_repair(r, CensusRoute::RepeatedUnit, 1e-9);
        require(r.outcome == "REPAIRED" && r.allocation_tie == 2 && r.candidates[0].allocation_tied && r.candidates[1].allocation_tied,
                "an exact allocation tie is recorded with both genotypes");
    }
    {   // allocation-equivalence set within 1e-9 after the prior; the smallest canonical SHA-256 pair serialises it
        RepairResult r;
        r.candidates = {cand("A", "B", -10.0, "AB", {1, 1, 2}, -1.0, 0), cand("C", "D", -10.0 + 5e-10, "AB", {1, 1, 2}, -1.0, 1),
                        cand("E", "F", -10.0 - 5e-9, "AB", {1, 1, 2}, -1.0, 2)};
        r.candidates[0].ctx1_id = "ff"; r.candidates[0].ctx2_id = "ee";
        r.candidates[1].ctx1_id = "dd"; r.candidates[1].ctx2_id = "cc";
        r.candidates[2].ctx1_id = "00"; r.candidates[2].ctx2_id = "01";
        choose_repair(r, CensusRoute::RepeatedUnit, 1e-9);
        require(r.outcome == "REPAIRED" && r.allocation_tie == 2 && r.representative_only, "a within-tolerance allocation set is unresolved");
        require(r.chosen == 1, "the representative is the smallest canonical SHA-256 pair, not the rank or the likelihood");
        require(!r.candidates[2].allocation_tied, "an arrangement beyond the tolerance is not in the set");
    }
    {   // different content multisets tied: never resolved by the allocation rule
        RepairResult r;
        r.candidates = {cand("A", "B", -10.0, "AB", {1, 1, 2}, -1.0, 0), cand("C", "D", -10.0, "AB", {1, 2, 2}, 0.0, 1)};
        choose_repair(r, CensusRoute::RepeatedUnit, 1e-9);
        require(r.outcome == "TRIGGERED_UNREPAIRED" && r.reason == "top_genotype_tie", "a cross-multiset tie refuses");
    }
}

// ---------------------------------------------------------------------------------------------- provenance

std::pair<std::size_t, std::size_t> cigar_lengths(const std::string& cigar, std::size_t& edits, std::size_t& columns) {
    std::size_t q = 0, t = 0, n = 0;
    edits = columns = 0;
    for (char c : cigar) {
        if (std::isdigit(static_cast<unsigned char>(c))) { n = n * 10 + static_cast<std::size_t>(c - '0'); continue; }
        columns += n;
        if (c == '=' || c == 'X') { q += n; t += n; if (c == 'X') edits += n; }
        else if (c == 'I') { q += n; edits += n; }
        else if (c == 'D') { t += n; edits += n; }
        n = 0;
    }
    return {q, t};
}

void test_provenance() {
    const std::vector<PanelSequence> panel = {{"p3", "TTTTTTTT"}, {"p1", "ACGTACGA"}, {"p2", "ACGTACGC"}, {"p4", "ACGTACGA"}};
    const std::vector<NearestPanel> n = nearest_panel_sequences("ACGTACGT", panel);
    require(n.size() == 3, "every equally nearest panel sequence is reported (p1, p2 and the duplicate p4)");
    for (const NearestPanel& x : n) require(x.distance == 1, "ties share the minimal distance");
    for (std::size_t i = 1; i < n.size(); ++i) {
        require(n[i - 1].sha256 < n[i].sha256 || (n[i - 1].sha256 == n[i].sha256 && n[i - 1].id < n[i].id), "ties ordered by (SHA-256, id)");
    }
    std::vector<PanelSequence> shuffled = panel;
    std::reverse(shuffled.begin(), shuffled.end());
    std::swap(shuffled[0], shuffled[2]);
    const std::vector<NearestPanel> m = nearest_panel_sequences("ACGTACGT", shuffled);
    bool same = m.size() == n.size();
    for (std::size_t i = 0; same && i < n.size(); ++i) same = m[i].id == n[i].id && m[i].cigar == n[i].cigar && m[i].sha256 == n[i].sha256;
    require(same, "provenance is invariant to panel ordering");
    require(panel_digest(panel) == panel_digest(shuffled), "the searched-panel digest is invariant to panel ordering");
    // CIGAR consistency on a longer case with an indel
    std::mt19937_64 rng(11);
    const std::string t = random_dna(rng, 3000);
    std::string q = t;
    q.erase(1200, 7);
    q[2000] = q[2000] == 'A' ? 'C' : 'A';
    const std::vector<NearestPanel> c = nearest_panel_sequences(q, {{"t", t}});
    std::size_t edits = 0, cols = 0;
    const auto [ql, tl] = cigar_lengths(c.at(0).cigar, edits, cols);
    require(ql == q.size() && tl == t.size(), "the CIGAR spans both sequences");
    require(static_cast<long>(edits) == c[0].distance && cols == c[0].alignment_columns, "the CIGAR re-derives the distance and columns");
    require(std::fabs(c[0].normalized - static_cast<double>(c[0].distance) / static_cast<double>(cols)) < 1e-15, "normalised by alignment columns");
    // display-id collisions fail closed
    std::map<std::string, std::string> shown;
    require(register_display_id(shown, "abcdef0123456789aaaa"), "first display id registers");
    require(register_display_id(shown, "abcdef0123456789aaaa"), "the same full id may be registered again");
    require(!register_display_id(shown, "abcdef012345ffffbbbb"), "a forced display-id collision is detected");
}

} // namespace

// ---------------------------------------------------------------------------------------------- repair mode (section 22)

// ctx:0 (one allele) | sv:1 {a, b} | ctx:2 {c, d}; templates t1 = (c0, a, c) and t2 = (c0, b, d).
PreparedGenotypeModel two_template_model(std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    StructuralLocusInput in;
    in.locus_id = "integration-fixture";
    in.invariant_segments = {"", "", "", ""};
    in.blocks = {{"ctx:0", {{"c0", random_dna(rng, 300)}}, GenotypeBlockKind::Context},
                 {"sv:1", {{"a", random_dna(rng, 200)}, {"b", random_dna(rng, 260)}}, GenotypeBlockKind::Variant},
                 {"ctx:2", {{"c", random_dna(rng, 300)}, {"d", random_dna(rng, 310)}}, GenotypeBlockKind::Context}};
    in.templates = {{"t1", TemplateFrameStatus::Complete, {"c0", "a", "c"}, ChainOrientation::Forward, ""},
                    {"t2", TemplateFrameStatus::Complete, {"c0", "b", "d"}, ChainOrientation::Forward, ""}};
    return prepare_genotype_model(in);
}

void test_template_switches() {
    const PreparedGenotypeModel m = two_template_model(11);
    const std::size_t a = allele_id(m, 1, "a"), b = allele_id(m, 1, "b"), c = allele_id(m, 2, "c"), d = allele_id(m, 2, "d");
    const auto sw = [&](std::size_t x, std::size_t y) { return template_switches(m, MosaicPath{{0, x, y}}); };
    require(sw(a, c) && sw(a, c)->count == 0 && sw(a, c)->blocks.empty(), "a panel template's own path has no switch");
    require(sw(a, d) && sw(a, d)->count == 1 && sw(a, d)->blocks == std::vector<std::size_t>{2}, "one switch, at the block that starts the new segment");
    require(sw(b, c) && sw(b, c)->count == 1 && sw(b, c)->blocks == std::vector<std::size_t>{2}, "the switch is placed by the left-greedy rule");
    require(!template_switches(m, MosaicPath{{0, a, 99}}), "an allele no template carries has no decomposition");
}

std::string slurp_file(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    std::ostringstream s; s << in.rdbuf(); return s.str();
}

// A synthetic committed run: status, equivalence set, block tables, FASTA and a dosage block table, for a strict call
// (hap1, hap2) whose sv:1 the dosage stage classes as `sv1_status`.
void write_synthetic_run(const PreparedGenotypeModel& m, const std::string& prefix, const MosaicPath& h1, const MosaicPath& h2,
                         const std::string& run_status, const std::string& sv1_status, const std::string& sv1_class) {
    const std::string cls = run_status == "COMPLETE" ? "strict_call" : "genotype_equivalent_call";
    atomic_write(prefix + ".mosaic_status.tsv", "key\tvalue\nstatus\t" + run_status + "\ncall_class\t" + cls +
                                                    "\nno_call_reason\tnone\nphase_support\tREADS_WITHIN_CANDIDATES\n"
                                                    "library_provenance\tdeclared: synthetic\n");
    const auto path = [](const MosaicPath& p) {
        std::string o; for (std::size_t i = 0; i < p.alleles.size(); ++i) o += (i ? "," : "") + std::to_string(p.alleles[i]); return o;
    };
    atomic_write(prefix + ".equivalence_set.tsv", "rank\tcandidate_id\trepresentative\thap1_alleles\thap2_alleles\tscore\tfragment_term\t"
                                                  "exposure_term\tsequence_pair_prior\n1\t0\ttrue\t" + path(h1) + "\t" + path(h2) +
                                                  "\t-10.5\t-8\t-1.5\t-1\n");
    std::string bg = "block_id\tkind\tallele1\tallele2\tallele1_bp\tallele2_bp\tzygosity\tconditional_genotype_posterior\tconditional_GQ\tfilter\n";
    std::string gt = "block_id\tkind\tallele1\tallele2\tallele1_bp\tallele2_bp\tallele1_sources\tallele2_sources\tzygosity\t"
                     "conditional_genotype_posterior\tconditional_GQ\tfilter\n";
    std::string db = "block_id\tdosage_status\tdosage_reason\trepair_class\torientation_tie\n";
    for (std::size_t b = 0; b < m.blocks().size(); ++b) {
        const PreparedGenotypeBlock& k = m.blocks()[b];
        const std::size_t x = h1.alleles[b], y = h2.alleles[b];
        const std::string row = k.block_id + "\t" + (b == 1 ? "SV" : "context") + "\t" + std::to_string(x) + "\t" + std::to_string(y) + "\t" +
                                std::to_string(k.alleles[x].sequence.size()) + "\t" + std::to_string(k.alleles[y].sequence.size());
        bg += row + "\t" + (x == y ? "HOM" : "HET") + "\t0.9\t10\tPASS\n";
        if (b == 1) gt += row + "\tsrc\tsrc\t" + (x == y ? "HOM" : "HET") + "\t0.9\t10\tPASS\n";
        db += k.block_id + "\t" + (b == 1 ? sv1_status : "NOT_SENSITIVE") + "\tnone\t" + (b == 1 ? sv1_class : "NA") + "\tfalse\n";
    }
    atomic_write(prefix + ".block_genotypes.tsv", bg);
    atomic_write(prefix + ".genotypes.tsv", gt);
    atomic_write(prefix + ".dosage_blocks.tsv", db);
    std::string fa;
    for (const auto& [name, p] : {std::make_pair("hap1", &h1), std::make_pair("hap2", &h2)}) {
        const std::string s = spell_mosaic(m, *p).bases;
        fa += ">" + prefix + "_" + name + " fingerprint=f bp=" + std::to_string(s.size()) + " phase_support=READS_WITHIN_CANDIDATES\n" + s + "\n";
    }
    atomic_write(prefix + ".mosaic_haplotypes.fa", fa);
}

// One block's cells of a block table, by column name.
std::map<std::string, std::string> block_row(const std::string& path, const std::string& block_id) {
    std::istringstream in(slurp_file(path));
    std::vector<std::string> head;
    std::map<std::string, std::string> out;
    for (std::string line; std::getline(in, line);) {
        std::vector<std::string> f;
        std::stringstream ls(line);
        for (std::string x; std::getline(ls, x, '\t');) f.push_back(x);
        if (head.empty()) { head = f; continue; }
        if (f.empty() || f[0] != block_id) continue;
        for (std::size_t i = 0; i < head.size() && i < f.size(); ++i) out[head[i]] = f[i];
    }
    return out;
}

// The headline of one block after write_locus_genotype, against the expected filters and reasons.
void require_headline(const std::string& prefix, const std::string& block_id, const std::string& genotype, const std::string& genotype_reasons,
                      const std::string& dosage, const std::string& dosage_reasons, const std::string& what) {
    for (const char* table : {".block_genotypes.tsv", ".genotypes.tsv"}) {
        auto r = block_row(prefix + table, block_id);
        require(r["genotype_filter"] == genotype && r["genotype_filter_reasons"] == genotype_reasons && r["dosage_filter"] == dosage &&
                    r["dosage_filter_reasons"] == dosage_reasons,
                what + table + ": got " + r["genotype_filter"] + "/" + r["genotype_filter_reasons"] + " " + r["dosage_filter"] + "/" +
                    r["dosage_filter_reasons"]);
    }
}

void test_repair_mode_integration() {
    const PreparedGenotypeModel m = two_template_model(12);
    const std::size_t a = allele_id(m, 1, "a"), b = allele_id(m, 1, "b"), c = allele_id(m, 2, "c"), d = allele_id(m, 2, "d");
    const std::filesystem::path dir = std::filesystem::temp_directory_path() / ("panvar_dosage_integration_" + std::to_string(::getpid()));
    std::filesystem::create_directories(dir);
    const MosaicPath h1{{0, a, c}}, h2{{0, a, d}};
    const auto repair_to = [&](std::size_t x, std::size_t y) {
        CompleteRepair r;
        r.block = 1; r.block_id = "sv:1"; r.route = "pair";
        r.seq1 = m.blocks()[1].alleles[x].sequence; r.seq2 = m.blocks()[1].alleles[y].sequence;
        r.allele1 = x; r.allele2 = y;
        r.fragment_term = -7; r.exposure_term = -1; r.prior = -2; r.repair_separation = 0.5;
        return r;
    };
    const std::vector<std::string> legacy = {".genotypes.tsv", ".block_genotypes.tsv", ".mosaic_haplotypes.fa", ".mosaic_status.tsv",
                                             ".equivalence_set.tsv"};
    {   // an applied COMPLETE_DIPLOTYPE_REPAIR: sv:1 a/a -> a/b; the originals are the inputs byte for byte
        const std::string p = (dir / "applied").string();
        write_synthetic_run(m, p, h1, h2, "COMPLETE", "REPAIRED", "COMPLETE_DIPLOTYPE_REPAIR");
        std::map<std::string, std::string> before;
        for (const std::string& s : legacy) before[s] = slurp_file(p + s);
        DosageShadowSummary summary;
        summary.complete_repairs = {repair_to(a, b)};
        const DosageIntegration done = apply_dosage_mode(&m, "repair", p, p, &summary, "COMPLETE", "none");
        require(done.applied == std::vector<std::string>{"sv:1"} && done.not_applied.empty(), "an eligible complete repair is applied");
        for (const std::string& s : legacy) require(slurp_file(p + ".original" + s) == before[s], "original" + s + " is the caller's file");
        const std::string bg = slurp_file(p + ".block_genotypes.tsv");
        require(bg.find("sv:1\tSV\t" + std::to_string(a) + "\t" + std::to_string(b) + "\t") != std::string::npos, "the repaired pair is the block genotype");
        require(bg.find("\tNA\tNA\tNOT_EVALUATED_AFTER_DOSAGE_REPAIR\tREPAIRED\t") != std::string::npos,
                "the repaired block's posterior and GQ are withdrawn, and its filter is never the caller's without the reads");
        require(bg.find("\tPASS\tNA\tNA\tnot_evaluated_after_dosage_repair\tnot_evaluated_after_dosage_repair\n") != std::string::npos,
                "the caller's filter is kept as original_filter, with the basis");
        require(slurp_file(p + ".mosaic_status.tsv").find("dosage_filters_basis\tnot_evaluated_after_dosage_repair\n") != std::string::npos,
                "the status says the reported genotype's filters were not evaluated");
        require(slurp_file(p + ".equivalence_set.tsv") == before[".equivalence_set.tsv"], "the search record is not rewritten");
        std::string fa_bases;
        { std::istringstream fa(slurp_file(p + ".mosaic_haplotypes.fa")); for (std::string l; std::getline(fa, l);) if (l.empty() || l[0] != '>') fa_bases += l; }
        require(fa_bases == spell_mosaic(m, h1).bases + spell_mosaic(m, MosaicPath{{0, b, d}}).bases, "hap2 is respelled with the repaired allele");
        require(slurp_file(p + ".mosaic_status.tsv").rfind(before[".mosaic_status.tsv"], 0) == 0, "the status keeps every caller line, first");
        write_locus_genotype(&m, p, done);
        const std::string locus = slurp_file(p + ".locus_genotype.tsv");
        std::vector<std::string> head, row;
        { std::istringstream l(locus); std::string h, r; std::getline(l, h); std::getline(l, r);
          std::stringstream hs(h), rs(r); for (std::string x; std::getline(hs, x, '\t');) head.push_back(x); for (std::string x; std::getline(rs, x, '\t');) row.push_back(x); }
        require(head == locus_genotype_columns() && row.size() == head.size(), "the locus table has its declared columns");
        std::map<std::string, std::string> v;
        for (std::size_t i = 0; i < head.size() && i < row.size(); ++i) v[head[i]] = row[i];
        require(v["reported_genotype_source"] == "caller_with_dosage_repair" && v["hap2_block_path"] == "0," + std::to_string(b) + "," + std::to_string(d),
                "the locus paths carry the repair");
        require(v["hap2_sha256"] == sha256_bytes(spell_mosaic(m, MosaicPath{{0, b, d}}).bases), "the locus digest is the respelled haplotype's");
        require(v["hap2_switches"] == "0" && v["hap1_switches"] == "0", "both repaired homologues are panel templates here");
        require(v["total_score"] == "-10" && v["score_basis"] == "dosage_repaired_block_rescored", "the locus score is the repaired diplotype's");
        require(v["locus_repair_scope"] == "COMPLETE" && v["dosage_complete_applied_blocks"] == "1", "the whole-locus repair scope");
    }
    const auto refused = [&](const std::string& name, const std::string& run_status, std::vector<CompleteRepair> repairs,
                             const std::string& reason, bool corrupt_fasta) {
        const std::string p = (dir / name).string();
        write_synthetic_run(m, p, h1, h2, run_status, "REPAIRED", "COMPLETE_DIPLOTYPE_REPAIR");
        if (corrupt_fasta) {
            std::string fa = slurp_file(p + ".mosaic_haplotypes.fa");
            const std::size_t at = fa.find('\n') + 1;
            fa[at] = fa[at] == 'A' ? 'C' : 'A';
            atomic_write(p + ".mosaic_haplotypes.fa", fa);
        }
        std::map<std::string, std::string> before;
        for (const std::string& s : legacy) before[s] = slurp_file(p + s);
        DosageShadowSummary summary;
        summary.complete_repairs = repairs;
        const DosageIntegration done = apply_dosage_mode(&m, "repair", p, p, &summary, "COMPLETE", "none");
        require(done.applied.empty() && done.not_applied.count("sv:1") && done.not_applied.at("sv:1") == reason, name + ": refused as " + reason);
        require(slurp_file(p + ".mosaic_haplotypes.fa") == before[".mosaic_haplotypes.fa"], name + ": FASTA untouched");
        const std::string bg = slurp_file(p + ".block_genotypes.tsv");
        require(bg.find("sv:1\tcontext") == std::string::npos && bg.find("sv:1\tSV\t" + std::to_string(a) + "\t" + std::to_string(a) + "\t") != std::string::npos,
                name + ": the caller's genotype stays");
        require(bg.find("DOSAGE_REPAIR_NOT_APPLIED\t" + reason) != std::string::npos, name + ": the block says why");
        for (const std::string& s : legacy) require(slurp_file(p + ".original" + s) == before[s], name + ": original" + s);
        // headline: a guard refused the candidate; the unchanged genotype is dosage-inconsistent and the dosage content
        // unconfirmed
        write_locus_genotype(&m, p, done);
        require_headline(p, "sv:1", "LOW_QUALITY", "EXACT_DIPLOTYPE_UNRESOLVED", "LOW_QUALITY", "COMPLETE_REPAIR_NOT_APPLIED:" + reason,
                         name + ": repair-mode refusal headline ");
    };
    CompleteRepair second = repair_to(a, b);
    second.block = 2; second.block_id = "ctx:2";
    refused("two", "COMPLETE", {repair_to(a, b), second}, "several_complete_repairs_not_jointly_scored", false);
    refused("equivalent", "CONDITIONAL_GENOTYPE_EQUIVALENT", {repair_to(a, b)}, "not_a_strict_call", false);
    CompleteRepair novel = repair_to(a, b);
    novel.allele2.reset();
    refused("novel", "COMPLETE", {novel}, "repaired_allele_not_a_panel_allele", false);
    refused("fasta", "COMPLETE", {repair_to(a, b)}, "reported_haplotypes_not_reproduced", true);
    {   // audit never touches a file; off copies nothing onto itself
        const std::string p = (dir / "audit").string();
        write_synthetic_run(m, p, h1, h2, "COMPLETE", "REPAIRED", "COMPLETE_DIPLOTYPE_REPAIR");
        std::map<std::string, std::string> before;
        for (const std::string& s : legacy) before[s] = slurp_file(p + s);
        DosageShadowSummary summary;
        summary.complete_repairs = {repair_to(a, b)};
        const DosageIntegration done = apply_dosage_mode(&m, "audit", p, p, &summary, "COMPLETE", "none");
        require(done.applied.empty() && done.not_applied.empty(), "audit applies nothing");
        require(done.audit_not_applied == std::map<std::string, std::string>{{"sv:1", headline::kAuditModeNotApplied}},
                "audit records the eligible complete repair as available, not applied because of the mode");
        for (const std::string& s : legacy) require(slurp_file(p + s) == before[s], "audit leaves " + s);
        require(!std::filesystem::exists(p + ".original.genotypes.tsv"), "audit writes no originals");
        // headline: usable dosage content, an unchanged dosage-inconsistent genotype
        write_locus_genotype(&m, p, done);
        require_headline(p, "sv:1", "LOW_QUALITY", "EXACT_DIPLOTYPE_UNRESOLVED", "PASS", "COMPLETE_REPAIR_AVAILABLE_NOT_APPLIED",
                         "audit-mode available repair headline ");
        require(block_row(p + ".block_genotypes.tsv", "sv:1")["allele2"] == std::to_string(a), "audit keeps the caller's genotype");
    }
    {   // audit mode, where repair mode's guard would refuse the repair: the guard's reason, never PASS
        const std::string p = (dir / "audit_refused").string();
        write_synthetic_run(m, p, h1, h2, "COMPLETE", "REPAIRED", "COMPLETE_DIPLOTYPE_REPAIR");
        std::string fa = slurp_file(p + ".mosaic_haplotypes.fa");
        const std::size_t at = fa.find('\n') + 1;
        fa[at] = fa[at] == 'A' ? 'C' : 'A';
        atomic_write(p + ".mosaic_haplotypes.fa", fa);
        DosageShadowSummary summary;
        summary.complete_repairs = {repair_to(a, b)};
        const DosageIntegration done = apply_dosage_mode(&m, "audit", p, p, &summary, "COMPLETE", "none");
        require(done.applied.empty() && done.audit_not_applied == std::map<std::string, std::string>{{"sv:1", "reported_haplotypes_not_reproduced"}},
                "audit records the guard's reason");
        write_locus_genotype(&m, p, done);
        require_headline(p, "sv:1", "LOW_QUALITY", "EXACT_DIPLOTYPE_UNRESOLVED", "LOW_QUALITY",
                         "COMPLETE_REPAIR_NOT_APPLIED:reported_haplotypes_not_reproduced", "audit-mode refused repair headline ");
    }
    std::filesystem::remove_all(dir);
}

// The dosage block row the stage writes for sv:1 from a classified repair: the fields integration and the headline read.
void write_dosage_row(const PreparedGenotypeModel& m, const std::string& prefix, const RepairResult& r) {
    std::string db = "block_id\tdosage_status\tdosage_reason\trepair_class\torientation_tie\tdosage_content_scope\tdiplotype_scope\t"
                     "allocation_equivalence_scope\tserialization_representative\trepaired_alleles_role\tallocation_status\n";
    for (std::size_t b = 0; b < m.blocks().size(); ++b) {
        if (b != 1) { db += m.blocks()[b].block_id + "\tNOT_SENSITIVE\tnone\tNA\tfalse\tNA\tNA\tNA\tfalse\tNA\tNOT_APPLICABLE\n"; continue; }
        db += m.blocks()[b].block_id + "\t" + dosage_status_name(repair_dosage_status(r)) + "\t" + r.reason + "\t" + r.repair_class + "\tfalse\t" +
              r.dosage_content_scope + "\t" + r.diplotype_scope + "\t" + r.allocation_equivalence_scope + "\t" +
              (r.serialization_representative ? "true" : "false") + "\t" + repaired_alleles_role(r) + "\tUNRESOLVED\n";
    }
    atomic_write(prefix + ".dosage_blocks.tsv", db);
}

// An unresolved allocation-equivalence set (repeated-unit route, complete dosage content, complete candidate
// enumeration, two allocation-equivalent distinct unordered genotypes) is conditional content: never a resolved
// diplotype, never in the primary block tables or FASTA, labelled REPRESENTATIVE_ONLY_NOT_RESOLVED in the dosage
// diagnostics. Plus the integration guard against a misclassified COMPLETE_DIPLOTYPE_REPAIR representative.
void test_unresolved_allocation() {
    const PreparedGenotypeModel m = two_template_model(13);
    const std::size_t a = allele_id(m, 1, "a"), b = allele_id(m, 1, "b"), c = allele_id(m, 2, "c"), d = allele_id(m, 2, "d");
    const std::string sa = m.blocks()[1].alleles[a].sequence, sb = m.blocks()[1].alleles[b].sequence;
    // two arrangements of one content multiset, equal panel prior and likelihood: the set has two distinct genotypes
    const auto tied_result = [&](const std::string& scope) {
        RepairResult r;
        r.route = "repeated_unit";
        r.scope = scope;                 // COMPLETE: complete candidate enumeration
        r.content_complete = true;       // complete dosage content
        r.candidates = {cand(sa, sb, -10.0, "AB", {1, 1, 2}, -1.0, 0), cand(sb, sb, -10.0, "AB", {1, 1, 2}, -1.0, 1)};
        r.candidates[0].ctx1_id = "00"; r.candidates[0].ctx2_id = "01";   // the smallest canonical pair serialises the set
        r.candidates[1].ctx1_id = "01"; r.candidates[1].ctx2_id = "01";
        r.candidates[0].units_a = {0}; r.candidates[0].units_b = {1, 2};
        r.candidates[1].units_a = {1}; r.candidates[1].units_b = {0, 2};
        choose_repair(r, CensusRoute::RepeatedUnit, 1e-9);
        classify_repair_scope(r, CensusRoute::RepeatedUnit);
        return r;
    };
    const RepairResult r = tied_result("COMPLETE");
    require(r.outcome == "REPAIRED" && r.chosen == 0 && r.allocation_tie == 2 && r.representative_only,
            "fixture: a complete enumeration whose allocation-equivalence set has two distinct genotypes");
    require(r.repair_class == "CONDITIONAL_DOSAGE_CONTENT", "an unresolved allocation-equivalence set is CONDITIONAL_DOSAGE_CONTENT, got " + r.repair_class);
    require(r.serialization_representative && !resolved_diplotype_repair(r) &&
                std::string(repaired_alleles_role(r)) == "REPRESENTATIVE_ONLY_NOT_RESOLVED",
            "the representative is labelled REPRESENTATIVE_ONLY_NOT_RESOLVED");
    require(r.dosage_content_scope == "COMPLETE" && r.diplotype_scope == "COMPLETE" && r.allocation_equivalence_scope == "COMPLETE",
            "search and content coverage stay COMPLETE: coverage is not uniqueness");
    require(r.reason == "unresolved_allocation_equivalence_set", "the reason names the unresolved set");
    require(repair_dosage_status(r) == DosageStatus::ContentRepaired, "the block's dosage status is CONTENT_REPAIRED");
    {   // controls: a single arrangement, and the pair route, stay complete and resolved
        RepairResult one;
        one.route = "repeated_unit"; one.scope = "COMPLETE"; one.content_complete = true;
        one.candidates = {cand(sa, sb, -10.0, "AB", {1, 1, 2}, -1.0, 0), cand(sb, sb, -12.0, "AB", {1, 1, 2}, -1.0, 1)};
        choose_repair(one, CensusRoute::RepeatedUnit, 1e-9);
        classify_repair_scope(one, CensusRoute::RepeatedUnit);
        require(one.repair_class == "COMPLETE_DIPLOTYPE_REPAIR" && !one.serialization_representative && resolved_diplotype_repair(one) &&
                    std::string(repaired_alleles_role(one)) == "RESOLVED" && repair_dosage_status(one) == DosageStatus::Repaired,
                "a single allocation is a resolved complete repair");
        RepairResult pr;
        pr.route = "pair"; pr.scope = "COMPLETE";
        pr.candidates = {cand(sa, sb, -10.0, "AB"), cand(sb, sa, -10.0, "BA")};
        choose_repair(pr, CensusRoute::Pair, 1e-9);
        classify_repair_scope(pr, CensusRoute::Pair);
        require(pr.repair_class == "COMPLETE_DIPLOTYPE_REPAIR" && !pr.serialization_representative && pr.allocation_equivalence_scope == "NOT_APPLICABLE",
                "a pair-route repair with an orientation tie only is still a resolved complete repair");
        const RepairResult cond = tied_result("CONDITIONAL");
        require(cond.repair_class == "CONDITIONAL_DOSAGE_CONTENT" && cond.reason == "none" && cond.diplotype_scope == "CONDITIONAL",
                "a conditional enumeration keeps its class and reason");
    }
    // the dosage diagnostics carry the representative, and only with its label
    const std::vector<std::pair<std::string, RepairResult>> repairs = {{"sv:1", r}};
    const std::string bfa = repaired_block_alleles_fasta(repairs);
    require(bfa.find("role=REPRESENTATIVE_ONLY_NOT_RESOLVED") != std::string::npos && bfa.find("role=RESOLVED") == std::string::npos &&
                bfa.find("repair_class=CONDITIONAL_DOSAGE_CONTENT") != std::string::npos,
            "dosage_repaired_block_alleles.fa labels the representative REPRESENTATIVE_ONLY_NOT_RESOLVED");
    {
        std::string bases;
        std::istringstream in(bfa);
        for (std::string l; std::getline(in, l);) bases += l.empty() || l[0] == '>' ? "\n" : l;
        require(bases.find(sb) != std::string::npos, "the representative's sequence is in the dosage diagnostic FASTA");
    }
    require(repair_blocks_table(repairs).find("\tCONDITIONAL_DOSAGE_CONTENT\tCOMPLETE\tCOMPLETE\tCOMPLETE\ttrue\t") != std::string::npos,
            "dosage_repairs.tsv: conditional class, complete coverage scopes, serialization representative");
    require(repair_candidates_table(repairs).find("\ttrue\tREPRESENTATIVE_ONLY_NOT_RESOLVED\n") != std::string::npos,
            "dosage_repair_candidates.tsv: the chosen row is REPRESENTATIVE_ONLY_NOT_RESOLVED");

    const std::filesystem::path dir = std::filesystem::temp_directory_path() / ("panvar_dosage_unresolved_" + std::to_string(::getpid()));
    std::filesystem::create_directories(dir);
    const MosaicPath h1{{0, a, c}}, h2{{0, a, d}};   // the caller's sv:1 is a/a; the representative is a/b
    const std::vector<std::string> legacy = {".genotypes.tsv", ".block_genotypes.tsv", ".mosaic_haplotypes.fa", ".mosaic_status.tsv",
                                             ".equivalence_set.tsv"};
    const std::string represented = spell_mosaic(m, MosaicPath{{0, b, d}}).bases;
    for (const std::string mode : {"repair", "audit", "off"}) {
        const std::string p = (dir / mode).string();
        write_synthetic_run(m, p, h1, h2, "COMPLETE", "NA", "NA");
        write_dosage_row(m, p, r);
        std::map<std::string, std::string> before;
        for (const std::string& s : legacy) before[s] = slurp_file(p + s);
        // the stage records a CompleteRepair only for COMPLETE_DIPLOTYPE_REPAIR: this block contributes none
        DosageShadowSummary summary;
        require(r.repair_class != "COMPLETE_DIPLOTYPE_REPAIR", mode + ": no complete repair is recorded");
        const DosageIntegration done = apply_dosage_mode(&m, mode, p, p, &summary, "COMPLETE", "none");
        require(done.applied.empty() && done.not_applied.empty() && done.audit_not_applied.empty(), mode + ": nothing is applied or refused");
        const std::string fa = slurp_file(p + ".mosaic_haplotypes.fa");
        require(fa == before[".mosaic_haplotypes.fa"], mode + ": the primary haplotype FASTA is unchanged");
        require(fa.find("dosage_repair_role=RESOLVED") == std::string::npos && fa.find("RESOLVED") == std::string::npos,
                mode + ": no primary FASTA header claims a resolved repair");
        {
            std::string bases;
            std::istringstream in(fa);
            for (std::string l; std::getline(in, l);) if (l.empty() || l[0] != '>') bases += l;
            require(bases.find(represented) == std::string::npos, mode + ": the representative haplotype is not in the primary FASTA");
        }
        for (const char* table : {".block_genotypes.tsv", ".genotypes.tsv"}) {
            auto row = block_row(p + table, "sv:1");
            require(row["allele1"] == std::to_string(a) && row["allele2"] == std::to_string(a), mode + table + ": the caller's alleles stay");
            if (mode != "repair") continue;   // audit and off add no dosage columns to the caller's tables
            require(row["repair_class"] == "CONDITIONAL_DOSAGE_CONTENT" && row["dosage_repair_applied"] == "NA" &&
                        row["serialization_representative"] == "true" && row["dosage_content_scope"] == "COMPLETE" &&
                        row["dosage_stage_filter"] == "DOSAGE_CONTENT_CORRECTED" &&
                        row["dosage_stage_filter_reason"] == "unresolved_allocation_equivalence_set",
                    mode + table + ": the conditional-content dosage fields");
        }
        if (mode == "off") {
            for (const std::string& s : legacy) require(slurp_file(p + s) == before[s], "off leaves " + s + " byte-identical");
        } else if (mode == "audit") {
            for (const std::string& s : legacy) require(slurp_file(p + s) == before[s], "audit leaves " + s + " byte-identical");
        } else {
            // repair mode annotates the block tables and the status; the FASTA and the search record are the caller's
            for (const std::string& s : legacy) require(slurp_file(p + ".original" + s) == before[s], "repair: original" + s + " is the caller's file");
            require(slurp_file(p + ".equivalence_set.tsv") == before[".equivalence_set.tsv"], "repair: the search record is not rewritten");
            const std::string st = slurp_file(p + ".mosaic_status.tsv");
            require(st.find("dosage_reported_genotype\tcaller\n") != std::string::npos && st.find("dosage_content_corrected_blocks\tsv:1\n") != std::string::npos &&
                        st.find("dosage_repaired_blocks\tnone\n") != std::string::npos,
                    "repair: the status reports the caller's genotype and a content correction at sv:1");
        }
        write_locus_genotype(&m, p, done);
        if (mode == "off") {
            require_headline(p, "sv:1", "PASS", "none", "NOT_EVALUATED", "DOSAGE_NOT_EVALUATED:dosage_mode_off", "off: unresolved-allocation headline ");
            continue;
        }
        // conditional content: the content is authoritative within the model, the exact diplotype is not resolved
        require_headline(p, "sv:1", "LOW_QUALITY", "EXACT_DIPLOTYPE_UNRESOLVED", "PASS", "CONTENT_CORRECTED_EXACT_DIPLOTYPE_UNRESOLVED",
                         mode + ": unresolved-allocation headline ");
        std::map<std::string, std::string> v;
        {
            std::istringstream l(slurp_file(p + ".locus_genotype.tsv"));
            std::string h, row;
            std::getline(l, h); std::getline(l, row);
            std::stringstream hs(h), rs(row);
            std::vector<std::string> head, cells;
            for (std::string x; std::getline(hs, x, '\t');) head.push_back(x);
            for (std::string x; std::getline(rs, x, '\t');) cells.push_back(x);
            for (std::size_t i = 0; i < head.size() && i < cells.size(); ++i) v[head[i]] = cells[i];
        }
        require(v["reported_genotype_source"] == "caller" && v["dosage_complete_class_blocks"] == "0" &&
                    v["dosage_complete_applied_blocks"] == "0" && v["dosage_content_corrected_block_ids"] == "sv:1" &&
                    v["locus_repair_scope"] == "CONDITIONAL",
                mode + ": the locus reports the caller's genotype, no complete repair and a content correction at sv:1");
        require(v["hap2_sha256"] == sha256_bytes(spell_mosaic(m, h2).bases), mode + ": the locus haplotype is the caller's");
    }
    {   // defence in depth: a misclassified COMPLETE_DIPLOTYPE_REPAIR that is a serialization representative is refused
        for (const std::string mode : {"repair", "audit"}) {
            const std::string p = (dir / ("misclassified_" + mode)).string();
            write_synthetic_run(m, p, h1, h2, "COMPLETE", "REPAIRED", "COMPLETE_DIPLOTYPE_REPAIR");
            std::map<std::string, std::string> before;
            for (const std::string& s : legacy) before[s] = slurp_file(p + s);
            CompleteRepair cr;
            cr.block = 1; cr.block_id = "sv:1"; cr.route = "repeated_unit";
            cr.seq1 = sa; cr.seq2 = sb; cr.allele1 = a; cr.allele2 = b;   // panel alleles: every other guard would pass
            cr.fragment_term = -7; cr.exposure_term = -1; cr.prior = -2; cr.repair_separation = 0.5;
            cr.serialization_representative = true;
            DosageShadowSummary summary;
            summary.complete_repairs = {cr};
            const DosageIntegration done = apply_dosage_mode(&m, mode, p, p, &summary, "COMPLETE", "none");
            const auto& refused = mode == "repair" ? done.not_applied : done.audit_not_applied;
            require(done.applied.empty() && refused == std::map<std::string, std::string>{{"sv:1", "unresolved_allocation_equivalence_set"}},
                    mode + ": a misclassified representative is refused as unresolved_allocation_equivalence_set");
            const std::string fa = slurp_file(p + ".mosaic_haplotypes.fa");
            require(fa == before[".mosaic_haplotypes.fa"] && fa.find("dosage_repair_role=RESOLVED") == std::string::npos,
                    mode + ": the primary FASTA is unchanged and never says dosage_repair_role=RESOLVED");
            require(block_row(p + ".block_genotypes.tsv", "sv:1")["allele2"] == std::to_string(a), mode + ": the caller's genotype stays");
            write_locus_genotype(&m, p, done);
            require_headline(p, "sv:1", "LOW_QUALITY", "EXACT_DIPLOTYPE_UNRESOLVED", "LOW_QUALITY",
                             "COMPLETE_REPAIR_NOT_APPLIED:unresolved_allocation_equivalence_set", mode + ": misclassified representative headline ");
        }
        // the guard takes precedence over run-level eligibility too: the reason is the representative's
        const std::string p = (dir / "misclassified_equivalent").string();
        write_synthetic_run(m, p, h1, h2, "CONDITIONAL_GENOTYPE_EQUIVALENT", "REPAIRED", "COMPLETE_DIPLOTYPE_REPAIR");
        CompleteRepair cr;
        cr.block = 1; cr.block_id = "sv:1"; cr.seq1 = sa; cr.seq2 = sb; cr.allele1 = a; cr.allele2 = b; cr.serialization_representative = true;
        DosageShadowSummary summary;
        summary.complete_repairs = {cr};
        const DosageIntegration done = apply_dosage_mode(&m, "repair", p, p, &summary, "COMPLETE", "none");
        require(done.applied.empty() && done.not_applied.at("sv:1") == "unresolved_allocation_equivalence_set",
                "the representative guard is checked before run-level eligibility");
    }
    std::filesystem::remove_all(dir);
}

int main() {
    try {
        test_python_random_and_quantile();
        test_kmers();
        test_undefined_estimator_cases();
        test_zero_spread();
        test_thread_identity();
        test_census_backends_agree();
        test_fixed_recursion_inside_dosage();
        test_unsupported_library();
        test_choice_rules();
        test_provenance();
        test_template_switches();
        test_repair_mode_integration();
        test_unresolved_allocation();
    } catch (const std::exception& e) {
        std::cerr << "FAIL: exception " << e.what() << '\n';
        return 1;
    }
    if (failures) { std::cerr << failures << " failure(s)\n"; return 1; }
    std::cout << "genotype_dosage: all fixtures pass\n";
    return 0;
}
