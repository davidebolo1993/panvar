// I3-A development driver (not a public CLI): successive samples in ONE process, either each built cold (the ordinary
// command, run_genotype_command) or all analysed against one prepared locus (prepare_genotype_locus once, then
// run_genotype_command_prepared). Each sample is a `panvar genotype` argument list, one argument per line, in a file.
//
//   genotype_prepared_locus_driver cold|prepared CEILING_BYTES CLEANUP(0|1) OUT_PREFIX ARGS_FILE [ARGS_FILE ...]
//   genotype_prepared_locus_driver write  CEILING_BYTES 0 OUT_PREFIX INDEX ARGS_FILE              (I3-C)
//   genotype_prepared_locus_driver write-panel CEILING_BYTES 0 OUT_PREFIX INDEX ARGS_FILE         (I3-D: with section 16)
//   genotype_prepared_locus_driver loaded CEILING_BYTES CLEANUP(0|1) OUT_PREFIX INDEX ARGS_FILE [ARGS_FILE ...]   (I3-C)
//
// I3-C: PANVAR_DRIVER_DIGEST=0 skips the I3-A convenience digest (an immutability check after loading/preparing and after
// every sample, ~0.4-1.9 s per call on these cells) so a process's wall time is the analysis alone; the note then says
// "digest_check_skipped". The e2e gate keeps the immutability check, and a load always verifies the index's complete digests.
//
// I3-C: `write` prepares from ARGS_FILE with its complete record, encodes and publishes INDEX (never overwriting), then
// verifies it: read + rewrite identity, a full load whose model, panel, fingerprint and report equal the prepared ones,
// and the live storage the structural input would hold if a PreparedLocus kept it. `loaded` loads INDEX (bounded,
// verified) and analyses every ARGS_FILE against it; the load's phase timings and accounting go in the "load" row.
//
// I3-A.1: every phase boundary (start, preparation, each sample, each between-sample cleanup, teardown) is recorded in a
// BatchResourceLedger against CEILING_BYTES. Once the process-lifetime peak exceeds the ceiling the batch starts no
// further sample; completed samples' outputs are left exactly as written. OUT_PREFIX.phases.tsv has one row per phase;
// OUT_PREFIX.batch_resource.tsv is the batch's resource verdict, separate from every sample's scientific result.
//
// Metrics, by platform. macOS: the physical footprint (current, and the kernel ledger's process-lifetime peak) and the
// malloc zones' live allocations. Elsewhere -- or with PANVAR_DRIVER_FORCE_PORTABLE -- resident size (current) and
// ru_maxrss (lifetime peak), labelled as RSS, and live allocations reported as unavailable ("NA"), never as zero.
#include "panvar/genotype_command.hpp"
#include "panvar/genotype_repair.hpp"
#include "panvar/marker_model.hpp"
#include "panvar/prepared_locus.hpp"
#include "panvar/prepared_locus_file.hpp"

#if defined(__APPLE__) && !defined(PANVAR_DRIVER_FORCE_PORTABLE)
#include <malloc/malloc.h>
#define PANVAR_DRIVER_MACOS_METRICS 1
#endif

#include <chrono>
#include <cstdlib>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace {
std::vector<std::string> read_args(const std::string& file) {
    std::ifstream in(file);
    if (!in) throw std::runtime_error("cannot read " + file);
    std::vector<std::string> args;
    for (std::string line; std::getline(in, line);) {
        if (!line.empty()) args.push_back(line);
    }
    return args;
}

double since(std::chrono::steady_clock::time_point t) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count();
}

// Every reading is optional: a platform query that fails (or reports 0, which no live process has) is unavailable,
// never a measured zero (I3-B).
struct Metrics {
    std::string label;
    std::optional<std::uint64_t> current, lifetime_peak, live_allocations;
};

std::optional<std::uint64_t> positive(std::uint64_t v) { return v > 0 ? std::optional<std::uint64_t>(v) : std::nullopt; }

Metrics measure() {
    Metrics m;
#if defined(PANVAR_DRIVER_MACOS_METRICS)
    m.label = "macos_phys_footprint";
    const panvar::FootprintSample f = panvar::physical_footprint();
    m.current = positive(f.current);
    m.lifetime_peak = positive(f.peak);
    malloc_statistics_t stats{};
    malloc_zone_statistics(nullptr, &stats);
    m.live_allocations = stats.size_in_use;
#else
    m.label = "rss";
    m.current = positive(panvar::current_rss_bytes());
    m.lifetime_peak = positive(static_cast<std::uint64_t>(panvar::peak_rss_gib() * 1073741824.0));
#endif
    return m;
}

std::string text(const std::optional<std::uint64_t>& v) { return v ? std::to_string(*v) : std::string("NA"); }
}  // namespace

int main(int argc, char** argv) {
    if (argc < 6) {
        std::cerr << "usage: genotype_prepared_locus_driver cold|prepared CEILING_BYTES CLEANUP(0|1) OUT_PREFIX ARGS_FILE...\n"
                     "       genotype_prepared_locus_driver write|write-panel|loaded CEILING_BYTES CLEANUP(0|1) OUT_PREFIX INDEX ARGS_FILE...\n";
        return 2;
    }
    const std::string mode = argv[1];
    const std::uint64_t ceiling = std::stoull(argv[2]);
    const bool cleanup = std::string(argv[3]) == "1";
    const std::string prefix = argv[4];
    const bool persist_panel = mode == "write-panel";   // I3-D
    const std::string effective_mode = persist_panel ? std::string("write") : mode;
    if (mode != "cold" && mode != "prepared" && mode != "write" && mode != "write-panel" && mode != "loaded") {
        std::cerr << "mode is cold, prepared, write or loaded\n";
        return 2;
    }
    // I3-C: the index path precedes the argument files in write and loaded modes
    const bool indexed = effective_mode == "write" || mode == "loaded";
    const std::string index_path = indexed ? std::string(argv[5]) : std::string();
    const int first_args = indexed ? 6 : 5;
    if (argc <= first_args) {
        std::cerr << "no argument file\n";
        return 2;
    }
    panvar::BatchResourceLedger ledger(ceiling, measure().label + "_lifetime_peak");
    std::ofstream out(prefix + ".phases.tsv");
    out << "phase\tsample\tseconds\texit\tmetric\tcurrent_bytes\tlifetime_peak_bytes\tpeak_increase_bytes\t"
           "live_allocation_bytes\tretained_process_storage_bytes\temission_cache_table_bytes\temission_cache_component_bytes\t"
           "prepared_builds\tmay_continue\tnote\n";
    const auto row = [&](const std::string& phase, const std::string& sample, double seconds, int code, const std::string& note) {
        const Metrics m = measure();
        const bool may = ledger.record(phase, m.current, m.lifetime_peak);
        const panvar::RetainedProcessStorage r = panvar::retained_process_storage();
        out << phase << '\t' << sample << '\t' << std::fixed << std::setprecision(3) << seconds << '\t' << code << '\t'
            << m.label << '\t' << text(m.current) << '\t' << text(m.lifetime_peak) << '\t' << text(ledger.phases().back().peak_increase) << '\t'
            << text(m.live_allocations) << '\t' << r.total() << '\t' << r.emission_cache_table_bytes << '\t'
            << r.emission_cache_component_bytes << '\t' << panvar::prepared_locus_builds() << '\t' << (may ? 1 : 0) << '\t'
            << note << '\n';
        out.flush();
        std::cout << phase << ' ' << sample << ": " << seconds << " s, exit " << code << ", " << m.label << ' ' << text(m.current)
                  << " (lifetime peak " << text(m.lifetime_peak) << "), live " << text(m.live_allocations) << ", retained "
                  << r.total() << (may ? "" : "  -> CEILING EXCEEDED, no further sample") << ' ' << note << '\n';
        return may;
    };
    const auto batch = std::chrono::steady_clock::now();
    bool may = row("start", "-", 0.0, 0, "");
    int worst = 0;
    std::size_t completed = 0, not_started = 0;
    double cleanup_seconds = 0.0;
    std::string error;
    try {
        std::shared_ptr<const panvar::PreparedLocus> locus;
        std::string digest;
        const char* digest_env = std::getenv("PANVAR_DRIVER_DIGEST");
        const bool check_digest = digest_env == nullptr || std::string(digest_env) != "0";
        if (effective_mode == "write" && may) {
            const auto live = [] { return measure().live_allocations; };
            const auto delta = [](std::optional<std::uint64_t> a, std::optional<std::uint64_t> b) {
                return a && b ? std::to_string(static_cast<long long>(*b) - static_cast<long long>(*a)) : std::string("NA");
            };
            const auto live_start = live();
            auto t = std::chrono::steady_clock::now();
            panvar::PreparedLocusPreparation prepared = panvar::prepare_genotype_locus_record(read_args(argv[first_args]), persist_panel);
            const double prepare_seconds = since(t);
            const auto live_prepared = live();
            t = std::chrono::steady_clock::now();
            const std::string bytes = panvar::encode_prepared_locus(prepared.record);
            const double encode_seconds = since(t);
            t = std::chrono::steady_clock::now();
            panvar::publish_prepared_locus(bytes, index_path);
            const double publish_seconds = since(t);
            may = row("write", "-", prepare_seconds + encode_seconds + publish_seconds, 0,
                      "prepare_with_record=" + std::to_string(prepare_seconds) + ";encode=" + std::to_string(encode_seconds) +
                          ";publish=" + std::to_string(publish_seconds) + ";file_bytes=" + std::to_string(bytes.size()) +
                          ";live_delta_prepared_with_record=" + delta(live_start, live_prepared));
            // verification (not part of the write cost): rewrite identity, full load equality, retention measurements
            t = std::chrono::steady_clock::now();
            panvar::PreparedLocusRecord read = panvar::read_prepared_locus_record(index_path, panvar::PreparedLocusLoadLimits{});
            const bool rewrite_identical = panvar::encode_prepared_locus(read) == bytes;
            const double read_seconds = since(t);
            const auto live_with_read_record = live();
            std::string input_live = "NA";
            {
                panvar::StructuralLocusInput input = std::move(read.input);
                read = panvar::PreparedLocusRecord{};
                const auto live_input_only = live();
                input = panvar::StructuralLocusInput{};
                input_live = delta(live(), live_input_only);
            }
            panvar::PreparedLocusLoadReport report;
            t = std::chrono::steady_clock::now();
            const auto loaded = panvar::load_prepared_locus(index_path, panvar::PreparedLocusLoadLimits{}, &report);
            const double load_seconds = since(t);
            panvar::PreparedLocusRecord with_loaded_report = prepared.record;
            with_loaded_report.report = loaded->report();
            const bool equal = loaded->fingerprint() == prepared.locus->fingerprint() &&
                               panvar::prepared_locus_model_content_sha256(loaded->model()) == prepared.record.model_content_sha256 &&
                               panvar::prepared_locus_panel_content_sha256(loaded->panel()) == prepared.record.panel_content_sha256 &&
                               loaded->node_sequences() == prepared.locus->node_sequences() &&
                               panvar::prepared_marker_panel_difference(loaded->panel(), prepared.locus->panel()).empty() &&
                               panvar::encode_prepared_locus(with_loaded_report) == bytes && loaded->digest() == prepared.locus->digest();
            prepared.record = panvar::PreparedLocusRecord{};
            const auto live_locus_only = live();
            std::string templates = std::to_string(loaded->model().panel().template_count());
            may = row("verify", "-", read_seconds + load_seconds, rewrite_identical && equal ? 0 : 5,
                      std::string("rewrite_identical=") + (rewrite_identical ? "1" : "0") + ";loaded_equals_prepared=" +
                          (equal ? "1" : "0") + ";read_seconds=" + std::to_string(read_seconds) + ";load_seconds=" +
                          std::to_string(load_seconds) + ";structural_input_live_bytes=" + input_live +
                          ";live_delta_two_loci_after_records_released=" + delta(live_start, live_locus_only) +
                          ";live_delta_with_read_record=" + delta(live_start, live_with_read_record) + ";templates=" + templates +
                          ";curated_exclusions=" + std::to_string(loaded->report().curated_exclusions.size()) +
                          ";file_checksum=" + report.file_checksum + ";panel_source=" + loaded->panel_source() +
                          ";panel_field_difference=" + panvar::prepared_marker_panel_difference(loaded->panel(), prepared.locus->panel()));
            if (!rewrite_identical || !equal) worst = 5;
        }
        if (mode == "loaded" && may) {
            panvar::PreparedLocusLoadReport report;
            const std::uint64_t builds0 = panvar::marker_panel_builds(), precounts0 = panvar::prepared_locus_panel_precounts();
            const auto t = std::chrono::steady_clock::now();
            locus = panvar::load_prepared_locus(index_path, panvar::PreparedLocusLoadLimits{}, &report);
            const double load_seconds = since(t);
            const std::uint64_t load_builds = panvar::marker_panel_builds() - builds0,
                                load_precounts = panvar::prepared_locus_panel_precounts() - precounts0;
            if (check_digest) digest = locus->digest();
            std::string note = "id=" + std::to_string(locus->id()) + ";file_bytes=" + std::to_string(report.file_bytes) +
                               ";accounted_peak=" + std::to_string(report.accounted_peak_bytes) + ";accounted_final=" +
                               std::to_string(report.accounted_final_bytes) + ";accounting_limit=" +
                               std::to_string(report.accounting_limit) + ";peak_phase=" + report.peak_phase +
                               ";model_projection=" + std::to_string(report.model_projection_bytes) + ";panel_projection=" +
                               std::to_string(report.panel_projection_bytes) + ";syncmer_occurrences=" +
                               std::to_string(report.panel_syncmer_occurrences) + ";structural_input_accounted=" +
                               std::to_string(report.structural_input_accounted_bytes) + ";panel_source=" + report.panel_source +
                               ";panel_rebuild_reason=" + report.panel_rebuild_reason + ";panel_decoded_accounted=" +
                               std::to_string(report.panel_decoded_accounted_bytes) + ";load_panel_builds=" +
                               std::to_string(load_builds) + ";load_precounts=" + std::to_string(load_precounts) + ";";
            for (const auto& [name, s] : report.seconds) note += name + "=" + std::to_string(s) + ";";
            may = row("load", "-", load_seconds, 0, note);
        }
        if (mode == "prepared" && may) {
            const auto t = std::chrono::steady_clock::now();
            locus = panvar::prepare_genotype_locus(read_args(argv[5]));
            const double prepare_seconds = since(t);
            digest = locus->digest();
            std::string seconds;
            for (const auto& [name, s] : locus->preparation_seconds()) seconds += name + "=" + std::to_string(s) + ";";
            may = row("prepare", "-", prepare_seconds, 0, "id=" + std::to_string(locus->id()) + ";digest=" + digest + ";" + seconds);
        }
        for (int i = effective_mode == "write" ? argc : first_args; i < argc; ++i) {
            if (!may) {
                ++not_started;
                continue;
            }
            const std::vector<std::string> args = read_args(argv[i]);
            const auto t = std::chrono::steady_clock::now();
            const int code = locus ? panvar::run_genotype_command_prepared(args, *locus) : panvar::run_genotype_command(args);
            const double seconds = since(t);
            worst = std::max(worst, code);
            ++completed;
            std::string note;
            if (locus) note = !check_digest ? "digest_check_skipped" : locus->digest() == digest ? "digest_unchanged" : "DIGEST_CHANGED";
            may = row("sample", argv[i], seconds, code, note);
            if (may && cleanup) {
                const panvar::BetweenSampleCleanup c = panvar::release_between_samples();
                cleanup_seconds += c.seconds;
                may = row("cleanup", argv[i], c.seconds, 0,
                          "released_process_storage=" + std::to_string(c.released.total()) +
                              ";free_pages_released=" + std::to_string(c.free_pages_released));
            }
        }
        if (locus) {
            const auto t = std::chrono::steady_clock::now();
            locus.reset();
            row("teardown", "-", since(t), 0, "");
        }
    } catch (const std::exception& e) {
        error = e.what();
        row("error", "-", 0.0, 3, error);
        worst = std::max(worst, 3);
    }
    row("batch", "-", since(batch), worst, "");
    std::ofstream verdict(prefix + ".batch_resource.tsv");
    verdict << "key\tvalue\nbatch_resource_verdict\t" << panvar::batch_resource_verdict_name(ledger.verdict())
            << "\nresource_metric\t" << ledger.metric() << "\nceiling_bytes\t" << ledger.ceiling() << "\nmax_lifetime_peak_bytes\t"
            << text(ledger.max_peak()) << "\nunverified_boundaries\t" << ledger.unverified_boundaries() << "\nexceeded_during_phase\t" << (ledger.exceeded_during().empty() ? "none" : ledger.exceeded_during())
            << "\nenforcement\tdetection_at_phase_boundaries_not_prevention\nsamples_requested\t" << (effective_mode == "write" ? 0 : argc - first_args)
            << "\nsamples_completed\t" << completed << "\nsamples_not_started\t" << not_started << "\nbetween_sample_cleanup\t"
            << (cleanup ? "on" : "off") << "\ncleanup_seconds_total\t" << std::fixed << std::setprecision(3) << cleanup_seconds
            << "\nworst_sample_exit\t" << worst << "\nerror\t" << (error.empty() ? "none" : error) << '\n';
    std::cout << "batch resource verdict: " << panvar::batch_resource_verdict_name(ledger.verdict()) << '\n';
    if (ledger.verdict() == panvar::BatchResourceVerdict::ceiling_exceeded) return 4;   // typed: distinct from sample exits
    return worst;
}
