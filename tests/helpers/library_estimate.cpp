// Fixtures for the automatic-library diagnostic (docs/algorithms/genotype.md, stage 3): recovery of a declared insert
// distribution, read lengths and error rate on a synthetic locus; the failure conditions (insufficient placements, a
// two-mode insert distribution, outward-facing pairs); and 1/2/4/8-thread identity.

#include "panvar/genotype_command.hpp"
#include "panvar/genotype_model.hpp"
#include "panvar/library_estimate.hpp"
#include "panvar/mosaic_spelling.hpp"

#include <cmath>
#include <functional>
#include <iostream>
#include <random>
#include <string>
#include <vector>

namespace {

int failures = 0;
void require(bool ok, const std::string& what) {
    if (!ok) { ++failures; std::cerr << "FAIL: " << what << '\n'; }
}

using namespace panvar;

std::string random_dna(std::mt19937_64& rng, std::size_t n) {
    static const char b[] = "ACGT";
    std::string s(n, 'A');
    for (char& c : s) c = b[rng() % 4];
    return s;
}

std::string revcomp(const std::string& s) {
    std::string o(s.rbegin(), s.rend());
    for (char& c : o) c = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : 'A';
    return o;
}

// ctx:0 (6.6 kb, with an 800 bp segment twice) | sv:1 {short 400 bp, long 400 bp + 1.2 kb} | ctx:2 (6 kb); two templates
PreparedGenotypeModel locus(std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    StructuralLocusInput in;
    in.locus_id = "library-fixture";
    in.invariant_segments = {"", "", "", ""};
    const std::string x = random_dna(rng, 400);
    // an invariant segment present twice per haplotype, so the invariant markers carry two copy numbers and the marker
    // depth fit can separate rate from background, as it does on a real locus
    const std::string repeat = random_dna(rng, 800);
    const std::string c0 = random_dna(rng, 3000) + repeat + random_dna(rng, 2000) + repeat;
    in.blocks = {{"ctx:0", {{"c0", c0}}, GenotypeBlockKind::Context},
                 {"sv:1", {{"short", x}, {"long", x + random_dna(rng, 1200)}}, GenotypeBlockKind::Variant},
                 {"ctx:2", {{"c2", random_dna(rng, 6000)}}, GenotypeBlockKind::Context}};
    in.templates = {{"t-short", TemplateFrameStatus::Complete, {"c0", "short", "c2"}, ChainOrientation::Forward, ""},
                    {"t-long", TemplateFrameStatus::Complete, {"c0", "long", "c2"}, ChainOrientation::Forward, ""}};
    return prepare_genotype_model(in);
}

// Paired reads of length L from both haplotypes; inserts from `insert()`; substitutions at rate e; outward pairs when
// asked (the mates swapped to face away from each other).
std::vector<ReadPair> simulate(const PreparedGenotypeModel& m, std::size_t pairs, const std::function<double(std::mt19937_64&)>& insert,
                               double e, std::uint64_t seed, bool outward = false, std::size_t L = 100) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    std::vector<std::string> hap = {spell_panel_template(m, 0).bases, spell_panel_template(m, 1).bases};
    std::vector<ReadPair> out;
    const auto mutate = [&](std::string s) {
        static const char b[] = "ACGT";
        for (char& c : s) if (u(rng) < e) { char n = c; while (n == c) n = b[rng() % 4]; c = n; }
        return s;
    };
    while (out.size() < pairs) {
        const std::string& h = hap[rng() % 2];
        const long len = std::lround(insert(rng));
        if (len < static_cast<long>(L) || len >= static_cast<long>(h.size())) continue;
        const std::size_t start = rng() % (h.size() - static_cast<std::size_t>(len));
        const std::string frag = h.substr(start, static_cast<std::size_t>(len));
        std::string a = mutate(frag.substr(0, L)), b = mutate(revcomp(frag.substr(frag.size() - L)));
        if (outward) { a = revcomp(a); b = revcomp(b); }
        const std::string name = "r" + std::to_string(out.size());
        out.push_back({{name + "/1", a, std::string(L, 'I')}, {name + "/2", b, std::string(L, 'I')}});
    }
    return out;
}

double value_of(const LibraryEstimate& e, const std::string& p) {
    for (const LibraryEstimateRow& r : e.rows) if (r.parameter == p) return std::atof(r.value.c_str());
    return std::nan("");
}
std::string status_of(const LibraryEstimate& e, const std::string& p) {
    for (const LibraryEstimateRow& r : e.rows) if (r.parameter == p) return r.status + ":" + r.reason;
    return "missing";
}

} // namespace

int main() {
    const PreparedGenotypeModel m = locus(7);
    LibraryModel declared;
    declared.insert_mean_bp = 350; declared.insert_sd_bp = 40; declared.support_sigmas = 4; declared.substitution_error = 0.002;
    declared.outlier_weight = 0.001; declared.fragment_start_intensity = 0.01;
    LibraryEstimateOptions opt;
    opt.threads = 4;
    std::normal_distribution<double> n350(350.0, 40.0);
    {   // recovery
        const std::vector<ReadPair> reads = simulate(m, 4000, [&](std::mt19937_64& r) { return n350(r); }, 0.002, 1);
        const LibraryEstimate e = estimate_library(m, reads, opt, &declared);
        require(std::fabs(e.insert_median - 350) < 6, "insert centre recovered (" + std::to_string(e.insert_median) + ")");
        require(std::fabs(e.insert_sd - 40) < 5, "insert spread recovered (" + std::to_string(e.insert_sd) + ")");
        require(status_of(e, "insert_mean_bp") == "ESTIMATED:none", "a clean library is ESTIMATED (" + status_of(e, "insert_mean_bp") + ")");
        require(e.read1_length_median == 100 && e.read2_length_median == 100, "read lengths measured");
        require(std::fabs(e.invariant_mismatch_rate - 0.002) < 0.001, "error rate at invariant sites (" + std::to_string(e.invariant_mismatch_rate) + ")");
        require(e.accepted > 1000 && e.orientation_discordant == 0, "most pairs place decisively, none discordant");
        require(!e.multimodal, "a normal insert is not multimodal");
        // thread identity: every row but the timing row
        for (std::size_t t : {1, 2, 8}) {
            LibraryEstimateOptions o = opt;
            o.threads = t;
            const LibraryEstimate f = estimate_library(m, reads, o, &declared);
            bool same = f.accepted_inserts == e.accepted_inserts && f.rows.size() == e.rows.size();
            for (std::size_t i = 0; same && i < f.rows.size(); ++i) {
                if (f.rows[i].parameter == "seconds_total") continue;
                same = f.rows[i].value == e.rows[i].value && f.rows[i].spread == e.rows[i].spread && f.rows[i].status == e.rows[i].status;
            }
            require(same, "identical at " + std::to_string(t) + " threads");
        }
    }
    {   // insufficient placements
        const std::vector<ReadPair> reads = simulate(m, 60, [&](std::mt19937_64& r) { return n350(r); }, 0.002, 2);
        const LibraryEstimate e = estimate_library(m, reads, opt, &declared);
        require(status_of(e, "insert_mean_bp").rfind("FALLBACK:insufficient_placements", 0) == 0, "too few placements fall back (" + status_of(e, "insert_mean_bp") + ")");
        require(value_of(e, "insert_mean_bp") == 350, "the fallback is the declared library's value");
        const LibraryEstimate none = estimate_library(m, reads, opt, nullptr);
        for (const LibraryEstimateRow& r : none.rows) {
            if (r.parameter == "insert_mean_bp") require(r.value == "NA" && r.source == "none_available", "no declared library: no fallback value is invented");
        }
    }
    {   // a two-mode insert distribution
        std::normal_distribution<double> a(300.0, 30.0), b(700.0, 30.0);
        const std::vector<ReadPair> reads = simulate(m, 4000, [&](std::mt19937_64& r) { return (r() % 2) ? a(r) : b(r); }, 0.002, 3);
        const LibraryEstimate e = estimate_library(m, reads, opt, &declared);
        require(e.multimodal && status_of(e, "insert_mean_bp").rfind("FALLBACK:multimodal", 0) == 0, "two insert modes fall back (" + status_of(e, "insert_mean_bp") + ")");
    }
    {   // outward-facing pairs are never accepted
        const std::vector<ReadPair> reads = simulate(m, 2000, [&](std::mt19937_64& r) { return n350(r); }, 0.002, 4, true);
        const LibraryEstimate e = estimate_library(m, reads, opt, &declared);
        require(e.accepted == 0 && e.orientation_discordant > 1000, "outward pairs are discordant, none accepted");
    }
    // --library auto: accepted estimates and each refusal (docs/modules/genotype.md, Required inputs)
    const LibraryModel preset = genotype_preset_library("short-read");
    {
        const std::vector<ReadPair> reads = simulate(m, 4000, [&](std::mt19937_64& r) { return n350(r); }, 0.002, 11);
        const LibraryEstimate e = estimate_library(m, reads, opt, nullptr);
        const AutoLibrary a = build_auto_library(e, preset, opt);
        require(a.accepted && a.refusal.empty(), "a clean estimate is accepted (" + a.refusal + ")");
        require(a.library.insert_mean_bp == e.insert_mean && a.library.insert_sd_bp == e.insert_sd, "auto takes the insert from the estimate");
        require(e.marker_fit_converged && a.library.fragment_start_intensity == e.marker_fit_intensity &&
                a.library.fragment_start_intensity != preset.fragment_start_intensity, "auto takes the marker-fit intensity, never the preset's");
        require(a.library.substitution_error == preset.substitution_error && a.library.outlier_weight == preset.outlier_weight &&
                a.library.support_sigmas == preset.support_sigmas, "every other value is the preset's");
        require(a.anchor_intensity == e.fragment_start_intensity, "the anchor-depth diagnostic is carried beside it");
        bool insert_placed = false, lambda_fit = false;
        for (const auto& [k, v] : a.sources) {
            if (k == "insert_mean_bp") insert_placed = v == "placement_estimate";
            if (k == "fragment_start_intensity") lambda_fit = v == "marker_fit_estimate";
        }
        require(insert_placed && lambda_fit, "the sources are recorded");
        require(library_used_rows(a).find("fragment_start_intensity_ratio_marker_to_anchor") != std::string::npos, "the ratio row is written");
        // a failed intensity estimate refuses, names --library FILE and carries no library out
        const auto intensity_refused = [&](const std::string& name, LibraryEstimate f) {
            const AutoLibrary r = build_auto_library(f, preset, opt);
            require(!r.accepted && r.refusal.find("fragment-start intensity") != std::string::npos && r.refusal.find("--library FILE") != std::string::npos,
                    name + ": refused, naming --library FILE (" + r.refusal + ")");
            require(r.library.fragment_start_intensity == 0.0 && r.sources.empty(), name + ": no library survives");
        };
        LibraryEstimate f = e; f.marker_fit_converged = false; intensity_refused("fit_not_converged", f);
        f = e; f.marker_fit_intensity = std::nan(""); intensity_refused("intensity_nan", f);
        f = e; f.marker_fit_intensity = 0.0; intensity_refused("intensity_zero", f);
        f = e; f.marker_fit_intensity = -1e-3; intensity_refused("intensity_negative", f);
    }
    const auto refused = [&](const std::string& name, const std::vector<ReadPair>& reads, const std::string& phrase) {
        const AutoLibrary a = build_auto_library(estimate_library(m, reads, opt, nullptr), preset, opt);
        require(!a.accepted, name + ": refused");
        require(a.refusal.find(phrase) != std::string::npos, name + ": says why (" + a.refusal + ")");
        require(a.refusal.find("--library FILE") != std::string::npos, name + ": names --library FILE");
        require(std::isnan(a.library.insert_mean_bp) || a.library.insert_mean_bp == 0.0, name + ": no insert is carried out of a refusal");
    };
    refused("insufficient", simulate(m, 300, [&](std::mt19937_64& r) { return n350(r); }, 0.002, 12), "decisively placed");
    {
        std::normal_distribution<double> a(300.0, 30.0), b(700.0, 30.0);
        refused("multimodal", simulate(m, 4000, [&](std::mt19937_64& r) { return (r() % 2) ? a(r) : b(r); }, 0.002, 13), "more than one mode");
    }
    {   // 10% of fragments far outside the library: contaminated inserts
        std::uniform_real_distribution<double> wide(1500.0, 4500.0);
        refused("contaminated_inserts", simulate(m, 4000, [&](std::mt19937_64& r) { return (r() % 10 == 0) ? wide(r) : n350(r); }, 0.002, 14),
                "contaminated");
    }
    {   // 20% outward-facing pairs mixed in: contaminated orientation
        std::vector<ReadPair> reads = simulate(m, 3200, [&](std::mt19937_64& r) { return n350(r); }, 0.002, 15);
        const std::vector<ReadPair> out = simulate(m, 800, [&](std::mt19937_64& r) { return n350(r); }, 0.002, 16, true);
        for (ReadPair p : out) { p.first.name = "o" + p.first.name; p.second.name = "o" + p.second.name; reads.push_back(p); }
        refused("contaminated_orientation", reads, "not inward-facing");
    }
    // the preset: technology only, long-read refused
    {
        bool threw = false;
        try { genotype_preset_settings("long-read"); } catch (const std::exception& error) { threw = std::string(error.what()).find("not implemented") != std::string::npos; }
        require(threw, "--preset long-read is refused as not implemented");
        threw = false;
        try { genotype_preset_settings("c4"); } catch (const std::exception&) { threw = true; }
        require(threw, "an unknown preset (no locus presets) is refused");
        // M15-D: the preset carries the release candidate's resource configuration (cache budgets, consumption
        // capacities, the footprint guard, the state-class cap and the model neighbourhood bound), but never a time
        // limit, a dosage stage, a candidate count, a reference or a locus-specific number.
        for (const auto& [k, v] : genotype_preset_settings("short-read")) {
            const bool forbidden = k == "time_limit_seconds" || k == "dosage_mode" || k == "dosage_audit" || k == "max_mosaic_candidates" ||
                                   k.find("reference") != std::string::npos || (k == "repair_max_neighbours" && v != "model");
            require(!forbidden, "the preset sets no time limit, dosage stage, candidate count, reference or locus-specific cap (" + k + ")");
        }
        require(std::isnan(genotype_preset_library("short-read").insert_mean_bp), "the preset declares no insert");
    }
    if (failures) { std::cerr << failures << " failure(s)\n"; return 1; }
    std::cout << "library_estimate: all fixtures pass\n";
    return 0;
}
