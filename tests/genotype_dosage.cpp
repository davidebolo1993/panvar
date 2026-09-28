// Fixtures for the dosage audit (docs/algorithms/genotype_dosage.md sections 5, 8, 9, 14): the frozen Python random
// stream and normal quantile, every undefined-estimator case, the robust-z rules at zero spread, the top-genotype and
// allocation rules, nearest-panel ties and panel-order invariance, display-id collisions, and thread identity.

#include "panvar/genotype_dosage.hpp"
#include "panvar/genotype_dosage_repair.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

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

int main() {
    try {
        test_python_random_and_quantile();
        test_kmers();
        test_undefined_estimator_cases();
        test_zero_spread();
        test_thread_identity();
        test_unsupported_library();
        test_choice_rules();
        test_provenance();
    } catch (const std::exception& e) {
        std::cerr << "FAIL: exception " << e.what() << '\n';
        return 1;
    }
    if (failures) { std::cerr << failures << " failure(s)\n"; return 1; }
    std::cout << "genotype_dosage: all fixtures pass\n";
    return 0;
}
