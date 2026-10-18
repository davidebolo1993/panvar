// Phase 7. The command is orchestration only, so its gate is that the orchestration is right: the
// declared execution order runs, every output appears, an incomplete run says so instead of
// presenting a partial ranking, and the model fingerprint identifies the model rather than the run.

#include "panvar/genotype_command.hpp"
#include "panvar/prepared_locus.hpp"
#include "panvar/prepared_locus_file.hpp"
#include "panvar/integrated_snarls.hpp"
#include "panvar/graph_sort.hpp"
#include "panvar/gfa_io.hpp"
#include "panvar/gfa.hpp"
#include "panvar/genotype_model.hpp"
#include "panvar/genotype_locus.hpp"
#include "panvar/marker_model.hpp"
#include "panvar/bubbles.hpp"
// For the declared piece-hit publication batch: the cap test bounds the work a refused parallel run may
// perform, and that bound is cap + threads * batch. Reading the default here rather than repeating the
// number keeps the test honest if the batch is ever retuned.
#include "panvar/fragment_likelihood.hpp"
#include "panvar/genotype_repair.hpp"
#include "panvar/ls_hmm.hpp"

#include <cstdlib>
#include <sys/wait.h>
#include <unistd.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <cstdio>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <optional>
#include <random>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

std::string dir;

std::string path(const std::string& name) { return dir + "/" + name; }

// Every block table this test's runs wrote -- on every call path -- carries quality_filter and quality_reason equal to
// block_quality(filter, genotype_support) of the same row: the two columns depend on nothing else, so no truth
// input can move them.
void test_block_quality_columns() {
    std::size_t tables = 0, rows = 0, pass = 0, low = 0;
    for (const auto& entry : std::filesystem::directory_iterator(dir)) {
        const std::string name = entry.path().string();
        const std::string suffix = ".block_genotypes.tsv";
        if (name.size() < suffix.size() || name.compare(name.size() - suffix.size(), suffix.size(), suffix) != 0) continue;
        std::ifstream in(name);
        std::string line;
        std::getline(in, line);
        std::vector<std::string> header;
        {
            std::stringstream fields(line);
            for (std::string cell; std::getline(fields, cell, '\t');) header.push_back(cell);
        }
        const auto column = [&](const std::string& key) {
            const auto at = std::find(header.begin(), header.end(), key);
            require(at != header.end(), name + ": no column " + key);
            return static_cast<std::size_t>(at - header.begin());
        };
        const std::size_t filter = column("filter"), support = column("genotype_support");
        const std::size_t quality = column("quality_filter"), reason = column("quality_reason");
        require(quality + 2 == header.size() && reason + 1 == header.size(), name + ": the quality columns come last");
        ++tables;
        while (std::getline(in, line)) {
            std::vector<std::string> cells;
            std::stringstream fields(line);
            for (std::string cell; std::getline(fields, cell, '\t');) cells.push_back(cell);
            const panvar::BlockQuality expected = panvar::block_quality(cells.at(filter), cells.at(support));
            require(cells.at(quality) == expected.filter && cells.at(reason) == expected.reason,
                    name + ": a row's quality columns differ from block_quality(filter, genotype_support)");
            ++rows;
            (expected.filter == "PASS" ? pass : low) += 1;
        }
    }
    // PASS needs the guard on and read support; most runs here keep the guard off, so their blocks are
    // LOW_QUALITY / representation_not_evaluated. The rule's every outcome is gated in genotype_block_quality.
    require(tables > 10 && pass > 0 && low > 0,
            "the block-quality check must see many tables and both outcomes (" + std::to_string(tables) + " tables)");
    std::cout << "  block quality columns: " << tables << " tables, " << rows << " rows (" << pass << " PASS, " << low
              << " LOW_QUALITY) equal block_quality(filter, genotype_support)\n";
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

// A hand-built locus: a shared backbone with two SV bubbles, each carrying two alleles whose
// interior span clears the 50 bp SV threshold, plus a 10 bp variant between them that does not.
// The short one is the point of the fixture: it is exactly the sub-threshold variation that must end
// up inside a context block rather than being dropped, because `--min-variant-bp` filters on
// interior SPAN, not on divergence.
struct TinyLocus {
    std::string gfa_path;
    std::vector<std::string> haplotypes;
    std::string reference_name;
    std::map<std::string, std::string> node_sequences;
};

// `context_left` and `context_right` are the shared spans either side of the sub-threshold variant.
// At their defaults they are wider than any simulated insert, so no fragment can link neighbouring
// variants and the locus cannot be phased from reads -- which one gate below relies on.
TinyLocus write_tiny_locus(const std::string& file = "tiny.gfa", std::size_t context_left = 350,
                           std::size_t context_right = 340, bool duplicate_first_template = false) {
    std::mt19937_64 rng(4);
    TinyLocus locus;
    locus.gfa_path = path(file);

    struct Node { std::string id, seq; };
    const std::vector<Node> nodes = {
        {"1", random_bases(rng, 700)},     // left flank
        {"2", random_bases(rng, 300)},     // SV block one, reference allele
        {"3", random_bases(rng, 90)},      // SV block one, alternate allele
        {"4", random_bases(rng, context_left)},    // shared context
        {"5", random_bases(rng, 10)},      // sub-threshold variant, allele one
        {"6", random_bases(rng, 10)},      // sub-threshold variant, allele two
        {"7", random_bases(rng, context_right)},   // shared context
        {"8", random_bases(rng, 260)},     // SV block two, reference allele
        {"9", random_bases(rng, 80)},      // SV block two, alternate allele
        {"10", random_bases(rng, 700)},    // right flank
    };
    std::ofstream out(locus.gfa_path);
    out << "H\tVN:Z:1.0\n";
    for (const Node& node : nodes) out << "S\t" << node.id << '\t' << node.seq << '\n';
    for (const char* edge : {"1\t+\t2\t+", "1\t+\t3\t+", "2\t+\t4\t+", "3\t+\t4\t+",
                             "4\t+\t5\t+", "4\t+\t6\t+", "5\t+\t7\t+", "6\t+\t7\t+",
                             "7\t+\t8\t+", "7\t+\t9\t+", "8\t+\t10\t+", "9\t+\t10\t+"}) {
        out << "L\t" << edge << "\t0M\n";
    }
    // Four templates covering both states of all three sites, in combinations no single template
    // spells alone, so a mosaic genuinely has to recombine.
    const std::vector<std::pair<std::string, std::vector<std::string>>> templates = {
        {"sampleA#1#chr1:0-100", {"1", "2", "4", "5", "7", "8", "10"}},
        {"sampleB#1#chr1:0-100", {"1", "3", "4", "6", "7", "9", "10"}},
        {"sampleC#1#chr1:0-100", {"1", "2", "4", "6", "7", "9", "10"}},
        {"sampleD#1#chr1:0-100", {"1", "3", "4", "5", "7", "8", "10"}},
    };
    std::vector<std::pair<std::string, std::vector<std::string>>> all_templates = templates;
    // A second carrier of sampleA's haplotype raises its prior without adding any sequence.
    if (duplicate_first_template) all_templates.push_back({"sampleE#1#chr1:0-100", templates.front().second});
    std::map<std::string, std::string> sequence_of;
    for (const Node& node : nodes) sequence_of[node.id] = node.seq;
    locus.node_sequences = sequence_of;
    for (const auto& item : all_templates) {
        out << "P\t" << item.first << '\t';
        for (std::size_t i = 0; i < item.second.size(); ++i) {
            if (i != 0) out << ',';
            out << item.second[i] << '+';
        }
        out << "\t*\n";
        std::string spelled;
        for (const std::string& node : item.second) spelled += sequence_of[node];
        locus.haplotypes.push_back(spelled);
    }
    locus.reference_name = templates.front().first;
    return locus;
}

void write_library(const std::string& file) {
    std::ofstream out(file);
    out << "insert_mean_bp 300\ninsert_sd_bp 30\nsupport_sigmas 4\n"
           "substitution_error 0.002\noutlier_weight 0.001\nfragment_start_intensity 0.085\n"
           // Declared: every invariant anchor of this tiny locus carries one copy, so the marker
           // background is not identifiable from the data. test_refusals proves the run refuses
           // without this line.
           "marker_background 0.01\n";
}

// Deterministic paired reads from two chosen haplotypes.
void simulate(const std::string& first, const std::string& second,
              const std::string& r1_path, const std::string& r2_path, std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::ofstream r1(r1_path);
    std::ofstream r2(r2_path);
    std::size_t n = 0;
    for (const std::string* hap : {&first, &second}) {
        const std::size_t pairs = 30 * hap->size() / (2 * 2 * 100);
        for (std::size_t i = 0; i < pairs; ++i) {
            const std::size_t insert = 260 + (rng() % 80);
            if (hap->size() <= insert) continue;
            const std::size_t start = rng() % (hap->size() - insert);
            const std::string mate1 = hap->substr(start, 100);
            const std::string mate2 = revcomp(hap->substr(start + insert - 100, 100));
            ++n;
            r1 << "@f" << n << "/1\n" << mate1 << "\n+\n" << std::string(100, 'I') << '\n';
            r2 << "@f" << n << "/2\n" << mate2 << "\n+\n" << std::string(100, 'I') << '\n';
        }
    }
    require(n > 100, "the simulated library must contain enough fragments to be informative");
}

std::map<std::string, std::string> read_status(const std::string& file) {
    std::ifstream input(file);
    require(static_cast<bool>(input), "status file " + file + " must exist");
    std::map<std::string, std::string> values;
    std::string line;
    std::getline(input, line);   // header
    while (std::getline(input, line)) {
        const std::size_t tab = line.find('\t');
        if (tab == std::string::npos) continue;
        values[line.substr(0, tab)] = line.substr(tab + 1);
    }
    return values;
}

bool exists(const std::string& file) {
    std::ifstream input(file);
    return static_cast<bool>(input);
}

std::string read_fasta_record(const std::string& file, std::size_t which) {
    std::ifstream input(file);
    std::string line;
    std::string current;
    std::vector<std::string> records;
    while (std::getline(input, line)) {
        if (!line.empty() && line[0] == '>') {
            if (!current.empty()) records.push_back(current);
            current.clear();
        } else {
            current += line;
        }
    }
    if (!current.empty()) records.push_back(current);
    require(which < records.size(), "haplotype FASTA must carry both called haplotypes");
    return records[which];
}

void test_end_to_end() {
    const TinyLocus locus = write_tiny_locus();
    write_library(path("lib.txt"));
    // Genotype sampleB against sampleC: a heterozygote at both SV blocks and both context blocks.
    simulate(locus.haplotypes[1], locus.haplotypes[2], path("r1.fq"), path("r2.fq"), 99);

    panvar::GenotypeCommandOptions options;
    options.gfa_path = locus.gfa_path;
    options.reference_path = locus.reference_name;
    options.reads1_path = path("r1.fq");
    options.reads2_path = path("r2.fq");
    options.library_path = path("lib.txt");
    options.out_prefix = path("run");
    options.max_mosaic_candidates = 16;
    options.kmer_size = 21;
    options.syncmer_s = 11;
    options.quiet = true;

    const panvar::GenotypeRunStatus status = panvar::run_genotype(options);
    require(status == panvar::GenotypeRunStatus::Complete,
            std::string("the end-to-end run must complete, got ") +
                panvar::genotype_run_status_name(status));

    // Every declared output exists.
    for (const char* suffix : {".genotypes.tsv", ".block_genotypes.tsv", ".mosaic_candidates.tsv",
                               ".mosaic_scores.tsv", ".mosaic_status.tsv", ".mosaic_haplotypes.fa"}) {
        require(exists(path("run") + suffix),
                std::string("the run must emit ") + suffix);
    }

    const std::map<std::string, std::string> reported = read_status(path("run.mosaic_status.tsv"));
    require(reported.at("status") == "COMPLETE", "the status file must agree with the return value");
    require(reported.at("every_allele_scored") == "true",
            "every block allele must be scored with no cap");
    require(reported.at("templates_refused") == "0", "no panel template may be refused here");
    require(reported.count("model_fingerprint") == 1 && reported.at("model_fingerprint").size() == 16,
            "the run must record a model fingerprint");
    require(std::stod(reported.at("numeric_separation_gap")) > 0.0,
            "a COMPLETE run must carry a positive numerical separation gap");
    require(reported.at("threads_used_combination") == "1" &&
                reported.at("threads_used_placement") == "1",
            "a one-thread run must report one worker in every phase");
    {
        // Diploid combination is parallel over fixed contiguous candidate ranges, each worker owning
        // complete score entries and no value accumulated across threads. The thread count is therefore
        // a scheduling choice and must not be visible in any number the run produces. Compared against
        // the one-thread run above, byte for byte.
        const auto slurp_run = [](const std::string& file, bool status_file) {
            std::ifstream input(file);
            if (!input) return std::string("<absent>");
            std::string line;
            std::string text;
            while (std::getline(input, line)) {
                // Timings, the thread counters and process memory are exactly what a thread count may change.
                if (status_file && (line.rfind("seconds_", 0) == 0 || line.rfind("threads_", 0) == 0 ||
                                    line.rfind("resources_", 0) == 0)) {
                    continue;
                }
                if (!line.empty() && line[0] == '>') line = line.substr(line.find(' '));
                text += line + '\n';
            }
            return text;
        };
        for (const std::size_t workers : {std::size_t{2}, std::size_t{4}}) {
            panvar::GenotypeCommandOptions threaded = options;
            threaded.threads = workers;
            const std::string name = "threads" + std::to_string(workers);
            threaded.out_prefix = path(name);
            require(panvar::run_genotype(threaded) == panvar::GenotypeRunStatus::Complete,
                    name + " must complete");
            for (const char* suffix : {".genotypes.tsv", ".block_genotypes.tsv", ".mosaic_scores.tsv",
                                       ".mosaic_status.tsv", ".mosaic_haplotypes.fa"}) {
                const bool status_file = std::string(suffix) == ".mosaic_status.tsv";
                require(slurp_run(path("run") + suffix, status_file) ==
                            slurp_run(path(name) + suffix, status_file),
                        name + " must reproduce the one-thread run exactly: " + suffix);
            }
            const auto threaded_status = read_status(path(name) + ".mosaic_status.tsv");
            require(threaded_status.at("threads_requested") == std::to_string(workers),
                    name + " must record what was asked for");
            // Placement parallelises over haplotypes, so it may legitimately use fewer workers than
            // asked for -- but never more, and never one when it could use two. Full-origin recording
            // is the one case that forces serial, and this fixture does not request it.
            const std::size_t used_placement = std::stoul(threaded_status.at("threads_used_placement"));
            require(used_placement >= 1 && used_placement <= workers,
                    name + " must never claim more placement workers than were requested");
            if (workers == 2) {
                require(used_placement == 2,
                        "placement must actually engage both workers when two are requested");
            }
        }

        // Fewer candidates than workers: the partition must degrade to what exists, not spawn empty
        // workers or divide by zero, and must still report only the threads it actually used.
        panvar::GenotypeCommandOptions starved_threads = options;
        starved_threads.threads = 64;
        starved_threads.max_mosaic_candidates = 0;
        starved_threads.phase_audit_candidates = 0;
        starved_threads.out_prefix = path("threads_starved");
        require(panvar::run_genotype(starved_threads) == panvar::GenotypeRunStatus::Complete,
                "a run with more workers than candidates must still complete");
        const auto starved_status = read_status(path("threads_starved.mosaic_status.tsv"));
        require(std::stoul(starved_status.at("threads_used_combination")) <=
                    std::stoul(starved_status.at("candidate_fragment_evaluations")),
                "a run cannot use more combination workers than it has candidate work");
    }

    // The called haplotypes must be the simulated truth, in either chain orientation.
    const std::string hap1 = read_fasta_record(path("run.mosaic_haplotypes.fa"), 0);
    const std::string hap2 = read_fasta_record(path("run.mosaic_haplotypes.fa"), 1);
    const auto matches = [&](const std::string& called, const std::string& truth) {
        return called == truth || called == revcomp(truth);
    };
    const bool recovered =
        (matches(hap1, locus.haplotypes[1]) && matches(hap2, locus.haplotypes[2])) ||
        (matches(hap1, locus.haplotypes[2]) && matches(hap2, locus.haplotypes[1]));
    require(recovered, "the certified winner must spell the two simulated truth haplotypes exactly");

    // Only variant blocks are reporting units.
    std::ifstream genotypes(path("run.genotypes.tsv"));
    std::string line;
    std::getline(genotypes, line);
    std::size_t rows = 0;
    while (std::getline(genotypes, line)) {
        require(line.rfind("sv:", 0) == 0,
                "only SV blocks are reported as calls; context blocks are modelled, not called");
        ++rows;
    }
    require(rows == 2, "the fixture's two SV blocks must both be reported");

    // The block-evaluation table carries every statistical block, context blocks included, and marks
    // which are reportable SV calls.
    std::ifstream blocks(path("run.block_genotypes.tsv"));
    std::getline(blocks, line);
    std::size_t sv_rows = 0;
    std::size_t context_rows = 0;
    while (std::getline(blocks, line)) {
        if (line.rfind("sv:", 0) == 0) {
            require(line.find("\tSV\ttrue\t") != std::string::npos, "an SV block must be reportable");
            ++sv_rows;
        } else {
            require(line.find("\tcontext\tfalse\t") != std::string::npos,
                    "a context block must be present and marked non-reportable");
            ++context_rows;
        }
    }
    require(sv_rows == std::stoul(reported.at("variant_blocks")) &&
                context_rows == std::stoul(reported.at("context_blocks")) && context_rows > 0,
            "the block table must carry every statistical block, not only the SV bubbles");

    // The fingerprint identifies the model, not the run: a different output path and a different K
    // must not move it.
    panvar::GenotypeCommandOptions again = options;
    again.out_prefix = path("other");
    again.max_mosaic_candidates = 8;
    require(panvar::run_genotype(again) == panvar::GenotypeRunStatus::Complete,
            "the second run must also complete");
    require(read_status(path("other.mosaic_status.tsv")).at("model_fingerprint") ==
                reported.at("model_fingerprint"),
            "the model fingerprint must not depend on the output prefix or on K");
}

void test_refusals() {
    panvar::GenotypeCommandOptions missing;
    missing.out_prefix = path("refused");
    bool threw = false;
    try {
        panvar::run_genotype(missing);
    } catch (const std::exception&) {
        threw = true;
    }
    require(threw, "a run without its required inputs must refuse");
    require(!exists(path("refused.genotypes.tsv")),
            "a refused run must not leave a genotypes table behind");

    // Fail closed on an unidentifiable marker background. The tiny locus's invariant anchors all carry
    // one copy, so background and per-copy rate trade off exactly and must not be invented.
    const TinyLocus locus = write_tiny_locus();
    {
        std::ofstream out(path("lib_undeclared.txt"));
        out << "insert_mean_bp 300\ninsert_sd_bp 30\nsupport_sigmas 4\n"
               "substitution_error 0.002\noutlier_weight 0.001\nfragment_start_intensity 0.02\n";
    }
    panvar::GenotypeCommandOptions undeclared;
    undeclared.gfa_path = locus.gfa_path;
    undeclared.reference_path = locus.reference_name;
    undeclared.reads1_path = path("r1.fq");
    undeclared.reads2_path = path("r2.fq");
    undeclared.library_path = path("lib_undeclared.txt");
    undeclared.out_prefix = path("undeclared");
    undeclared.max_mosaic_candidates = 16;
    undeclared.kmer_size = 21;
    undeclared.syncmer_s = 11;
    undeclared.quiet = true;
    require(panvar::run_genotype(undeclared) == panvar::GenotypeRunStatus::CandidateIncomplete,
            "an unidentifiable marker background must end the run as CANDIDATE_INCOMPLETE");
    const std::map<std::string, std::string> status = read_status(path("undeclared.mosaic_status.tsv"));
    require(status.at("status") == "CANDIDATE_INCOMPLETE" &&
                status.at("note").find("not separately identifiable") != std::string::npos,
            "the status file must name why the marker depth could not be estimated");
    require(!exists(path("undeclared.genotypes.tsv")) && !exists(path("undeclared.mosaic_scores.tsv")),
            "a run that never ranked must not emit a ranking or a call");
}


panvar::GenotypeCommandOptions tiny_options(const std::string& prefix, const std::string& reads1,
                                            const std::string& reads2,
                                            const std::string& gfa = "tiny.gfa") {
    panvar::GenotypeCommandOptions options;
    options.gfa_path = path(gfa);
    options.reference_path = "sampleA#1#chr1:0-100";
    options.reads1_path = reads1;
    options.reads2_path = reads2;
    options.library_path = path("lib.txt");
    options.out_prefix = path(prefix);
    options.max_mosaic_candidates = 16;
    options.kmer_size = 21;
    options.syncmer_s = 11;
    options.quiet = true;
    return options;
}

bool no_ranking_or_call(const std::string& prefix) {
    return !exists(path(prefix + ".mosaic_scores.tsv")) && !exists(path(prefix + ".genotypes.tsv")) &&
           !exists(path(prefix + ".block_genotypes.tsv")) && !exists(path(prefix + ".mosaic_haplotypes.fa"));
}

// Every stage that can stop early ends with a named status, no ranking and no call.
void test_fail_closed_statuses() {
    const std::string r1 = path("r1.fq");
    const std::string r2 = path("r2.fq");

    panvar::GenotypeCommandOptions beam_budget = tiny_options("beam_budget", r1, r2);
    beam_budget.max_beam_extensions = 1;
    require(panvar::run_genotype(beam_budget) == panvar::GenotypeRunStatus::CandidateIncomplete,
            "an exhausted beam budget must end the run as CANDIDATE_INCOMPLETE");
    require(read_status(path("beam_budget.mosaic_status.tsv")).at("note").find("resource budget") !=
                std::string::npos,
            "the status must name the exhausted proposal budget");
    require(no_ranking_or_call("beam_budget"), "a run stopped in proposal must not rank or call");

    panvar::GenotypeCommandOptions scorer_cap = tiny_options("scorer_cap", r1, r2);
    scorer_cap.max_physical_origins = 1;
    require(panvar::run_genotype(scorer_cap) == panvar::GenotypeRunStatus::RescoreIncomplete,
            "an exhausted scoring cap must end the run as RESCORE_INCOMPLETE");
    require(read_status(path("scorer_cap.mosaic_status.tsv")).at("note").find("resource cap") !=
                std::string::npos,
            "the status must name the exhausted scoring cap");
    require(no_ranking_or_call("scorer_cap"), "a run stopped in rescoring must not rank or call");

    // Stale outputs: a complete run, then a refused run under the same prefix. The refused run must not
    // be left beside the earlier ranking.
    require(panvar::run_genotype(tiny_options("stale", r1, r2)) == panvar::GenotypeRunStatus::Complete,
            "fixture: the first run under the reused prefix must complete");
    require(exists(path("stale.genotypes.tsv")), "fixture: the first run must write its call");
    panvar::GenotypeCommandOptions reused = tiny_options("stale", r1, r2);
    reused.max_physical_origins = 1;
    require(panvar::run_genotype(reused) == panvar::GenotypeRunStatus::RescoreIncomplete,
            "the second run under the reused prefix must refuse");
    require(no_ranking_or_call("stale"),
            "a refused run must remove the earlier run's ranking and call from its prefix");

    // A structural failure writes MODEL_INCOMPLETE instead of escaping with no status. One panel path
    // gives no bubbles and therefore no block chain.
    {
        std::ofstream gfa(path("single.gfa"));
        gfa << "H\tVN:Z:1.0\nS\t1\tACGTACGTACGTACGTACGT\nS\t2\tTTGCATTGCATTGCA\nL\t1\t+\t2\t+\t0M\n"
               "P\tsolo#1#chr1:0-100\t1+,2+\t*\n";
    }
    panvar::GenotypeCommandOptions single = tiny_options("single", r1, r2);
    single.gfa_path = path("single.gfa");
    single.reference_path = "solo#1#chr1:0-100";
    require(panvar::run_genotype(single) == panvar::GenotypeRunStatus::ModelIncomplete,
            "a locus whose structural model cannot be built must end MODEL_INCOMPLETE");
    require(read_status(path("single.mosaic_status.tsv")).at("note").find("structural model") !=
                std::string::npos,
            "the status must say the structural model could not be built");
    require(no_ranking_or_call("single"), "a run stopped at the model must not rank or call");
}

std::vector<std::vector<std::string>> read_rows(const std::string& file) {
    std::ifstream input(file);
    require(static_cast<bool>(input), "cannot read " + file);
    std::vector<std::vector<std::string>> rows;
    std::string line;
    std::getline(input, line);
    while (std::getline(input, line)) {
        std::vector<std::string> cells;
        std::stringstream fields(line);
        for (std::string cell; std::getline(fields, cell, '\t');) cells.push_back(cell);
        rows.push_back(std::move(cells));
    }
    return rows;
}

std::vector<std::size_t> parse_path(const std::string& joined) {
    std::vector<std::size_t> out;
    std::stringstream stream(joined);
    for (std::string cell; std::getline(stream, cell, ',');) out.push_back(std::stoul(cell));
    return out;
}

struct TruthFrame {
    std::optional<panvar::PreparedGenotypeModel> model;
    panvar::PathProjection truth1;
    panvar::PathProjection truth2;
};

// Rebuilds the model exactly as the command did, refusing unless it reproduces the run's fingerprint,
// and projects both truth paths onto that model's block frame from a graph that also carries them.
// Sorting depends only on the reference path, so node identifiers name the same sequence in both.
TruthFrame rebuild_truth_frame(const panvar::GenotypeCommandOptions& options, const std::string& fingerprint,
                               const std::vector<std::string>& truth1_nodes,
                               const std::vector<std::string>& truth2_nodes) {
    TruthFrame frame;
    panvar::GfaModel gfa = panvar::read_gfa_model(options.gfa_path);
    panvar::GraphSortOptions sort_options;
    sort_options.reference_path = options.reference_path;
    sort_options.flip = true;
    const std::string resolved = panvar::sort_graph_reference(gfa, sort_options).resolved_reference;
    panvar::write_gfa_model(options.out_prefix + ".eval.sorted.gfa", gfa);
    panvar::ParseGfaOptions parse_options;
    const panvar::Graph graph = panvar::parse_gfa(options.out_prefix + ".eval.sorted.gfa", parse_options);
    panvar::BubbleCallOptions bubble_options;
    bubble_options.reference_path = resolved;
    bubble_options.quiet = true;
    bubble_options.min_variant_bp = options.min_variant_bp;
    bubble_options.snarl_pairs_override = panvar::find_top_level_snarls_cactus(panvar::snarl_input_from_model(gfa));
    bubble_options.snarl_source_supplied = true;
    const std::vector<panvar::Bubble> bubbles = panvar::call_bubbles(graph, bubble_options);
    panvar::LocusInputOptions locus_options;
    locus_options.locus_id = resolved;
    locus_options.reference_path_name = resolved;
    const panvar::LocusInputResult built = panvar::build_locus_input(graph, bubbles, locus_options);
    frame.model.emplace(panvar::prepare_genotype_model(built.input));
    require(panvar::model_fingerprint(*frame.model) == fingerprint,
            "the evaluation must rebuild the exact model the run scored");

    {
        std::ifstream in(options.gfa_path);
        std::ofstream out(options.out_prefix + ".truth.gfa");
        out << in.rdbuf();
        for (const auto& item : {std::make_pair("truth#1#x", &truth1_nodes),
                                 std::make_pair("truth#2#x", &truth2_nodes)}) {
            out << "P\t" << item.first << '\t';
            for (std::size_t i = 0; i < item.second->size(); ++i) {
                out << (i ? "," : "") << (*item.second)[i] << '+';
            }
            out << "\t*\n";
        }
    }
    panvar::GfaModel truth_gfa = panvar::read_gfa_model(options.out_prefix + ".truth.gfa");
    panvar::sort_graph_reference(truth_gfa, sort_options);
    panvar::write_gfa_model(options.out_prefix + ".truth.sorted.gfa", truth_gfa);
    const panvar::Graph truth_graph = panvar::parse_gfa(options.out_prefix + ".truth.sorted.gfa", parse_options);
    const auto t1 = panvar::project_path_onto_locus(truth_graph, bubbles, built.report, "truth#1#x");
    const auto t2 = panvar::project_path_onto_locus(truth_graph, bubbles, built.report, "truth#2#x");
    require(t1.has_value() && t2.has_value(), "both truth paths must project onto the block chain");
    frame.truth1 = *t1;
    frame.truth2 = *t2;
    const panvar::PreparedGenotypeModel& model = *frame.model;
    for (std::size_t b = 0; b < model.blocks().size(); ++b) {
        require(frame.truth1.blocks[b].has_value() && frame.truth2.blocks[b].has_value(),
                "every truth block must have a positional frame");
        bool available1 = false;
        bool available2 = false;
        for (const auto& allele : model.blocks()[b].alleles) {
            available1 = available1 || allele.sequence == *frame.truth1.blocks[b];
            available2 = available2 || allele.sequence == *frame.truth2.blocks[b];
        }
        require(available1 && available2, "fixture: every truth block allele must exist in the panel");
    }
    // Neither truth haplotype may be a panel row at every block: recombination is required.
    for (const panvar::PathProjection* truth : {&frame.truth1, &frame.truth2}) {
        for (std::size_t row = 0; row < model.panel().template_count(); ++row) {
            bool everywhere = true;
            for (std::size_t b = 0; b < model.blocks().size(); ++b) {
                const auto& allele = model.blocks()[b].alleles[model.panel().template_alleles[row][b]];
                everywhere = everywhere && allele.sequence == *truth->blocks[b];
            }
            require(!everywhere, "fixture: a truth haplotype must not be a panel row");
        }
    }
    return frame;
}

// True when some retained beam candidate spells both truth haplotypes, in either order, at every block.
bool truth_phasing_in_candidates(const TruthFrame& frame, const std::string& candidates_file) {
    const panvar::PreparedGenotypeModel& model = *frame.model;
    for (const auto& row : read_rows(candidates_file)) {
        const std::vector<std::size_t> hap1 = parse_path(row[2]);
        const std::vector<std::size_t> hap2 = parse_path(row[3]);
        const auto spells = [&](const std::vector<std::size_t>& path, const panvar::PathProjection& truth) {
            for (std::size_t b = 0; b < model.blocks().size(); ++b) {
                if (model.blocks()[b].alleles[path[b]].sequence != *truth.blocks[b]) return false;
            }
            return true;
        };
        if ((spells(hap1, frame.truth1) && spells(hap2, frame.truth2)) ||
            (spells(hap1, frame.truth2) && spells(hap2, frame.truth1))) {
            return true;
        }
    }
    return false;
}

const std::vector<std::string> kTruth1Nodes = {"1", "2", "4", "5", "7", "9", "10"};
const std::vector<std::string> kTruth2Nodes = {"1", "3", "4", "6", "7", "8", "10"};

std::pair<std::string, std::string> spell_truths(const TinyLocus& locus) {
    const auto spell = [&](const std::vector<std::string>& ids) {
        std::string out;
        for (const std::string& id : ids) out += locus.node_sequences.at(id);
        return out;
    };
    std::pair<std::string, std::string> truths = {spell(kTruth1Nodes), spell(kTruth2Nodes)};
    for (const std::string& row : locus.haplotypes) {
        require(row != truths.first && row != truths.second,
                "fixture: neither truth haplotype may be a panel row");
    }
    return truths;
}

// The first evidence the mosaic assembles a haplotype the panel does not contain. The contexts are
// short enough that fragments link neighbouring variants, so the reads can phase the truth. The call
// must be decided by the fragments, not by the prior, and match the truth at every block, phased.
void test_off_panel_truth_is_recovered() {
    const TinyLocus locus = write_tiny_locus("tiny_phaseable.gfa", 80, 70);
    const auto truths = spell_truths(locus);
    simulate(truths.first, truths.second, path("phaseable_r1.fq"), path("phaseable_r2.fq"), 123);
    // The crowding fixture. Before equivalence merging, a width of 64 filled with template histories of the
    // two zero-switch phasings and the truth was never offered. Under the same budget of 64 selection units
    // the truth must be proposed. Every panel diplotype is still scored, so the run calls and is marked
    // conditional (review round 6) with the truth among what it proposed.
    //
    // What this budget BUYS changed when production moved to the grouped proposer, and the change is a
    // unit change, not a tuning one: the state limit now counts allele-prefix PAIRS, where it used to
    // count (prefix, template) states. One prefix pair stands for every template attribution of itself,
    // so 64 is a far wider search than it was -- wide enough that this small locus is now searched
    // exhaustively, and its candidate set is conditional only because 26 mosaics were truncated to
    // K = 16. The starved case below keeps the pruned, missed-target route covered, since the caller
    // branches on it.
    {
        panvar::GenotypeCommandOptions budget = tiny_options(
            "phaseable_budget", path("phaseable_r1.fq"), path("phaseable_r2.fq"), "tiny_phaseable.gfa");
        budget.proposal_state_limit = 64;
        require(panvar::run_genotype(budget) == panvar::GenotypeRunStatus::Complete,
                "a pruned search short of K must still call from the complete panel baseline");
        const auto budget_status = read_status(path("phaseable_budget.mosaic_status.tsv"));
        require(budget_status.at("proposal_search_exhaustive") == "true" &&
                    budget_status.at("proposal_states_pruned") == "0",
                "fixture: 64 prefix pairs must search this locus exhaustively under the grouped proposer");
        require(budget_status.at("mosaic_target_reached") == "true" &&
                    std::stoul(budget_status.at("mosaic_candidates_found")) > 16,
                "fixture: 64 prefix pairs must find more mosaics than K=16, so the set is truncated");
        require(budget_status.at("panel_baseline_complete") == "true" &&
                    budget_status.at("candidate_scope") == "CONDITIONAL" && exists(path("phaseable_budget.genotypes.tsv")),
                "a mosaic-deficient call must be marked conditional and still emit its genotypes");
        const TruthFrame budget_frame = rebuild_truth_frame(
            budget, budget_status.at("model_fingerprint"), kTruth1Nodes, kTruth2Nodes);
        require(truth_phasing_in_candidates(budget_frame, path("phaseable_budget.mosaic_candidates.tsv")),
                "under the budget that used to crowd it out, the merged beam must propose the truth");
    }

    // The other route to a conditional set, which the 64-state budget used to take: a search too narrow
    // to reach K at all. run_genotype branches on mosaic_target_reached, so this path must stay covered.
    {
        panvar::GenotypeCommandOptions starved = tiny_options(
            "phaseable_starved", path("phaseable_r1.fq"), path("phaseable_r2.fq"), "tiny_phaseable.gfa");
        starved.proposal_state_limit = 4;
        require(panvar::run_genotype(starved) == panvar::GenotypeRunStatus::Complete,
                "a search short of K must still call from the complete panel baseline");
        const auto starved_status = read_status(path("phaseable_starved.mosaic_status.tsv"));
        require(starved_status.at("proposal_search_exhaustive") == "false" &&
                    starved_status.at("mosaic_target_reached") == "false",
                "fixture: 4 prefix pairs must prune and miss K=16");
        require(starved_status.at("panel_baseline_complete") == "true" &&
                    starved_status.at("candidate_scope") == "CONDITIONAL",
                "a search that misses K must be conditional and still carry the complete panel baseline");
    }

    // Recovery, with a budget above every reachable prefix pair.
    panvar::GenotypeCommandOptions options =
        tiny_options("phaseable", path("phaseable_r1.fq"), path("phaseable_r2.fq"), "tiny_phaseable.gfa");
    options.proposal_state_limit = 1048576;
    const panvar::GenotypeRunStatus status = panvar::run_genotype(options);
    const std::map<std::string, std::string> reported = read_status(path("phaseable.mosaic_status.tsv"));
    require(status == panvar::GenotypeRunStatus::Complete,
            std::string("the phaseable off-panel run must complete, got ") +
                panvar::genotype_run_status_name(status) + ": " +
                (reported.count("note") ? reported.at("note") : ""));
    require(reported.at("proposal_search_exhaustive") == "true" || reported.at("mosaic_target_reached") == "true",
            "an accuracy-eligible run must have searched exhaustively or reached its target");
    require(reported.at("phase_support") == "READS_WITHIN_CANDIDATES",
            "the reads, not the prior, must choose the phaseable off-panel truth within its candidates");

    const TruthFrame frame =
        rebuild_truth_frame(options, reported.at("model_fingerprint"), kTruth1Nodes, kTruth2Nodes);
    require(truth_phasing_in_candidates(frame, path("phaseable.mosaic_candidates.tsv")),
            "an exhaustive proposal must offer the off-panel truth");
    const panvar::PreparedGenotypeModel& model = *frame.model;
    const auto rows = read_rows(path("phaseable.block_genotypes.tsv"));
    require(rows.size() == model.blocks().size(), "the block table must carry every block");
    int pairing = 0;
    for (std::size_t b = 0; b < rows.size(); ++b) {
        const auto& alleles = model.blocks()[b].alleles;
        const std::string& c1 = alleles[std::stoul(rows[b][3])].sequence;
        const std::string& c2 = alleles[std::stoul(rows[b][4])].sequence;
        const std::string& a = *frame.truth1.blocks[b];
        const std::string& d = *frame.truth2.blocks[b];
        const bool direct = c1 == a && c2 == d;
        const bool swapped = c1 == d && c2 == a;
        require(direct || swapped, "block " + rows[b][0] + " must be called as the truth allele pair");
        if (a != d) {
            const int here = direct ? 1 : -1;
            require(pairing == 0 || pairing == here, "the call must be phased consistently with the truth");
            pairing = here;
        }
    }
    require(pairing != 0, "fixture: the truth must be heterozygous somewhere, or phase is untested");
    std::cout << "  off-panel truth (phaseable): all " << rows.size()
              << " blocks called and phased correctly; truth proposed within 64 prefix pairs; phase support READS_WITHIN_CANDIDATES, "
              << "fragment gap " << reported.at("fragment_term_gap_to_runner_up") << '\n';
}

// The same off-panel truth where the reads cannot phase it: every context span is wider than every
// insert, so no fragment links neighbouring variants. Under the declared prior the zero-switch template
// phasings that share the truth's genotype are then MORE probable than the truth, which needs a switch on
// both haplotypes -- so the gate is not "the truth wins" but that phase is never claimed from reads. At an
// exhaustive proposal width for this locus: every winning candidate carries the truth genotype at every
// block, the truth phasing is retained with a fragment term identical to the winner's, and the winner is
// never decided by fragments.
void test_off_panel_phase_is_not_claimed_from_reads() {
    const TinyLocus locus = write_tiny_locus();
    const auto truths = spell_truths(locus);
    simulate(truths.first, truths.second, path("unphaseable_r1.fq"), path("unphaseable_r2.fq"), 123);
    panvar::GenotypeCommandOptions options =
        tiny_options("unphaseable", path("unphaseable_r1.fq"), path("unphaseable_r2.fq"));
    options.proposal_state_limit = 1048576;   // far above every reachable prefix pair: exhaustive
    const panvar::GenotypeRunStatus status = panvar::run_genotype(options);
    const std::map<std::string, std::string> reported = read_status(path("unphaseable.mosaic_status.tsv"));
    require(reported.at("proposal_search_exhaustive") == "true", "fixture: the unphaseable search must be exhaustive");
    // Three admissible outcomes, and the invariant across all of them is that phase is never claimed
    // from reads that cannot carry it. BLOCK_CONSENSUS is the case this locus now usually reaches: the
    // two phasings of one genotype are inseparable, so no whole-locus winner is certified, every block
    // genotype is nonetheless agreed by both, and those blocks are called with phase UNRESOLVED.
    require((status == panvar::GenotypeRunStatus::RescoreIncomplete && reported.at("phase_support") == "UNRESOLVED") ||
                (status == panvar::GenotypeRunStatus::BlockConsensus && reported.at("phase_support") == "UNRESOLVED") ||
                (status == panvar::GenotypeRunStatus::Complete && reported.at("phase_support") == "PRIOR_ONLY"),
            std::string("an unphaseable truth must never be reported as read-phased, got ") +
                panvar::genotype_run_status_name(status) + " with phase support " +
                (reported.count("phase_support") ? reported.at("phase_support") : "NA"));
    if (status == panvar::GenotypeRunStatus::BlockConsensus) {
        require(!exists(path("unphaseable.mosaic_haplotypes.fa")),
                "a block-consensus run must emit no phased haplotype: the reads chose no phasing");
        require(std::stoul(reported.at("consensus_blocks_called")) > 0,
                "a block-consensus run must call the blocks its tied candidates agree on");
    }
    const TruthFrame frame =
        rebuild_truth_frame(options, reported.at("model_fingerprint"), kTruth1Nodes, kTruth2Nodes);
    const panvar::PreparedGenotypeModel& model = *frame.model;

    std::map<std::string, double> fragment_term;
    std::map<std::string, double> error_bound;
    std::set<std::string> winners;
    for (const auto& row : read_rows(path("unphaseable.mosaic_scores.tsv"))) {
        fragment_term[row[0]] = std::stod(row[4]);
        error_bound[row[0]] = std::stod(row[7]);
        if (row.back() == "true") winners.insert(row[0]);
    }
    require(!winners.empty(), "the ranking must report its winner set");

    std::string truth_candidate;
    std::size_t winning_members = 0;
    for (const auto& row : read_rows(path("unphaseable.mosaic_candidates.tsv"))) {
        const std::vector<std::size_t> hap1 = parse_path(row[2]);
        const std::vector<std::size_t> hap2 = parse_path(row[3]);
        bool genotype = true;
        bool direct_everywhere = true;
        bool swapped_everywhere = true;
        for (std::size_t b = 0; b < model.blocks().size(); ++b) {
            const auto& alleles = model.blocks()[b].alleles;
            const std::string& c1 = alleles[hap1[b]].sequence;
            const std::string& c2 = alleles[hap2[b]].sequence;
            const std::string& a = *frame.truth1.blocks[b];
            const std::string& d = *frame.truth2.blocks[b];
            const bool direct = c1 == a && c2 == d;
            const bool swapped = c1 == d && c2 == a;
            genotype = genotype && (direct || swapped);
            direct_everywhere = direct_everywhere && direct;
            swapped_everywhere = swapped_everywhere && swapped;
        }
        if (direct_everywhere || swapped_everywhere) truth_candidate = row[11];
        if (winners.count(row[11]) != 0) {
            ++winning_members;
            require(genotype, "every winning candidate must carry the truth genotype at every block");
        }
    }
    require(winning_members >= winners.size(), "every winner must appear in the candidate manifest");
    require(!truth_candidate.empty(), "an exhaustive proposal must retain the truth phasing");
    const std::string winner = *winners.begin();
    require(std::abs(fragment_term[truth_candidate] - fragment_term[winner]) <=
                error_bound[truth_candidate] + error_bound[winner],
            "with no phase-informative fragment the truth and the winner must have the same fragment term");
    std::cout << "  off-panel truth (unphaseable): " << panvar::genotype_run_status_name(status)
              << ", phase support " << reported.at("phase_support") << ", " << winners.size()
              << " winning phasing(s) with the truth genotype; truth phasing retained"
              << (winners.count(truth_candidate) ? " in the winner set" : " below the winner by prior only")
              << '\n';
}


// Candidate-set status rules (review round 6). Every panel diplotype is always scored, so a pruned search that
// finds fewer off-panel mosaics than K still calls, marked CONDITIONAL with the shortfall explicit; K = 0 scores
// the panel alone; a panel baseline over its budget is refused whole; an exhaustive search with fewer mosaics
// than requested is complete and unconditional.
// The candidate scope is ONE value, and the status and every block row must carry it. They disagreed:
// a panel-only run's status said PANEL_ONLY while every block row said CONDITIONAL, taken from an empty
// beam. That changes no genotype and corrupts anything that summarises the rows.
void require_scope_agrees(const std::string& prefix, const std::string& expected) {
    const auto status = read_status(path(prefix + ".mosaic_status.tsv"));
    require(status.count("candidate_scope") == 1 && status.at("candidate_scope") == expected,
            prefix + " must report candidate_scope " + expected + ", got " +
                (status.count("candidate_scope") ? status.at("candidate_scope") : std::string("(absent)")));
    std::size_t checked = 0;
    std::size_t scope_column = 0;
    {
        // By name: the filter columns follow candidate_scope, so it is no longer the last column.
        std::ifstream header_file(path(prefix + ".block_genotypes.tsv"));
        std::string header;
        std::getline(header_file, header);
        std::stringstream fields(header);
        std::size_t at = 0;
        bool found = false;
        for (std::string cell; std::getline(fields, cell, '\t'); ++at) {
            if (cell == "candidate_scope") {
                scope_column = at;
                found = true;
            }
        }
        require(found, prefix + ".block_genotypes.tsv has no candidate_scope column");
    }
    for (const auto& row : read_rows(path(prefix + ".block_genotypes.tsv"))) {
        if (row.empty() || row.front() == "block_id") continue;   // header, if read_rows returns one
        require(row.size() > scope_column && row[scope_column] == expected,
                prefix + " block row " + row.front() + " reports scope " + row[scope_column] +
                    " while its status reports " + expected);
        ++checked;
    }
    require(checked > 0, prefix + " must emit block rows for the scope check to mean anything");
}

void test_proposal_status_rules() {
    const std::string r1 = path("phaseable_r1.fq");
    const std::string r2 = path("phaseable_r2.fq");

    panvar::GenotypeCommandOptions starved = tiny_options("starved", r1, r2, "tiny_phaseable.gfa");
    starved.proposal_state_limit = 1;
    require(panvar::run_genotype(starved) == panvar::GenotypeRunStatus::Complete,
            "a pruned search short of its mosaic target must still call from the complete panel baseline");
    const auto starved_status = read_status(path("starved.mosaic_status.tsv"));
    require(starved_status.at("proposal_search_exhaustive") == "false" &&
                starved_status.at("mosaic_target_reached") == "false" &&
                std::stoull(starved_status.at("proposal_states_pruned")) > 0 &&
                starved_status.at("panel_baseline_complete") == "true" &&
                starved_status.at("candidate_scope") == "CONDITIONAL" && starved_status.at("global_certified") == "false",
            "a mosaic-deficient call must report the shortfall, a complete panel and a conditional scope");
    require(exists(path("starved.genotypes.tsv")) && exists(path("starved.block_genotypes.tsv")),
            "a mosaic-deficient call must emit its genotypes");

    // K = 0: the complete panel alone.
    panvar::GenotypeCommandOptions panel_only = tiny_options("panel_only", r1, r2, "tiny_phaseable.gfa");
    panel_only.max_mosaic_candidates = 0;
    panel_only.proposal_state_limit = 1048576;
    require(panvar::run_genotype(panel_only) == panvar::GenotypeRunStatus::Complete, "a K = 0 run must call");
    const auto panel_only_status = read_status(path("panel_only.mosaic_status.tsv"));
    std::size_t panel_only_mosaics = 0;
    for (const auto& row : read_rows(path("panel_only.mosaic_candidates.tsv"))) panel_only_mosaics += row[13] == "mosaic";
    // With K = 0, no audit and no positive prefix, the beam is not built at all: its diagnostics are
    // absent rather than zero, and the scope is neither COMPLETE nor CONDITIONAL, since both of those
    // describe what a search returned. The panel baseline is still scored whole and still calls.
    require(panel_only_mosaics == 0 && panel_only_status.at("panel_baseline_complete") == "true" &&
                panel_only_status.at("candidate_scope") == "PANEL_ONLY" &&
                std::stod(panel_only_status.at("seconds_proposal")) == 0.0 &&
                panel_only_status.count("mosaic_candidates") == 0 &&
                panel_only_status.count("mosaic_candidates_found") == 0 &&
                panel_only_status.count("beam_extensions") == 0 &&
                exists(path("panel_only.genotypes.tsv")),
            "K = 0 with no audit must rank the panel baseline alone, building no beam and reporting none");

    // A panel baseline over its declared budget is refused whole: no partial baseline, no ranking, no call.
    panvar::GenotypeCommandOptions over_budget = tiny_options("panel_budget", r1, r2, "tiny_phaseable.gfa");
    over_budget.max_panel_pair_fragment_evaluations = 1;
    require(panvar::run_genotype(over_budget) == panvar::GenotypeRunStatus::CandidateIncomplete,
            "a panel baseline over its budget must refuse");
    require(read_status(path("panel_budget.mosaic_status.tsv")).at("panel_baseline_complete") == "false" &&
                no_ranking_or_call("panel_budget"),
            "a refused panel baseline must not rank or call and must say the baseline is incomplete");

    panvar::GenotypeCommandOptions exhaustive = tiny_options("exhaustive", r1, r2, "tiny_phaseable.gfa");
    exhaustive.max_mosaic_candidates = 64;
    exhaustive.proposal_state_limit = 1048576;
    require(panvar::run_genotype(exhaustive) == panvar::GenotypeRunStatus::Complete,
            "an exhaustive search with fewer candidates than requested must complete");
    const auto exhaustive_status = read_status(path("exhaustive.mosaic_status.tsv"));
    require(exhaustive_status.at("proposal_search_exhaustive") == "true" &&
                exhaustive_status.at("mosaic_target_reached") == "false" &&
                exhaustive_status.at("candidate_set_conditional") == "false" &&
                std::stoul(exhaustive_status.at("mosaic_candidates_found")) < 64,
            "an exhaustive search below K must be reported complete and unconditional");
    require(exhaustive_status.at("candidate_scope") == "COMPLETE" && exhaustive_status.at("global_certified") == "true",
            "an exhaustive, untruncated, complete ranking is the only globally certified result");
    require(exhaustive_status.at("marker_depth_converged") == "true",
            "a completed run must rest on a converged marker-depth fit");
    // The protected panel baseline: every unordered panel pair is a candidate, reported complete, and the
    // manifest carries each once beside the off-panel mosaics that K counts.
    std::size_t panel_rows = 0;
    std::size_t mosaic_rows = 0;
    for (const auto& row : read_rows(path("exhaustive.mosaic_candidates.tsv"))) {
        require(row.size() == 14 && (row[13] == "panel" || row[13] == "mosaic"),
                "every manifest row must name its source");
        (row[13] == "panel" ? panel_rows : mosaic_rows) += 1;
    }
    require(exhaustive_status.at("panel_baseline_complete") == "true" &&
                std::stoul(exhaustive_status.at("panel_diplotypes")) == panel_rows && panel_rows > 0 &&
                std::stoul(exhaustive_status.at("mosaic_candidates")) == mosaic_rows &&
                std::stoul(exhaustive_status.at("candidates_total")) == panel_rows + mosaic_rows,
            "the status must report the complete panel baseline and the manifest must hold every panel pair");
    require(exhaustive_status.at("marker_background_effective") == exhaustive_status.at("marker_background") &&
                std::stod(exhaustive_status.at("marker_background_floor")) == panvar::kMarkerBackgroundFloor,
            "the status must report the model floor and the background proposal scoring actually used");

    // One scope, all three states, and every block row agreeing with the status beside it.
    require_scope_agrees("starved", "CONDITIONAL");
    require_scope_agrees("panel_only", "PANEL_ONLY");
    require_scope_agrees("exhaustive", "COMPLETE");
}

// Prior-only phase. The unphaseable locus with a second carrier of sampleA's haplotype: reads from
// sampleA and sampleB cannot link neighbouring variants, so the zero-switch phasings are fragment
// contenders, and the extra carrier lets the prior pick one. The call is made, but its phase must be
// reported PRIOR_ONLY, while every block genotype -- on which all contenders agree -- stays READS.
// The per-block conditional posterior, recomputed from the run's own tables: for each block, the summed
// candidate-set probability of every sequence candidate whose representative carries the winner's unordered
// genotype there; GQ is -10 log10 of the rest. Scope and separation must be stated on every row.
void check_block_posteriors(const std::string& prefix) {
    std::vector<std::vector<std::size_t>> hap1;
    std::vector<std::vector<std::size_t>> hap2;
    for (const auto& row : read_rows(path(prefix + ".mosaic_candidates.tsv"))) {
        hap1.push_back(parse_path(row[2]));
        hap2.push_back(parse_path(row[3]));
    }
    std::vector<std::pair<std::size_t, double>> representatives;   // (candidate index, set probability)
    std::size_t winner = hap1.size();
    for (const auto& row : read_rows(path(prefix + ".mosaic_scores.tsv"))) {
        representatives.emplace_back(std::stoul(row[1]), std::stod(row[10]));
        if (row[13] == "true") winner = std::stoul(row[1]);
    }
    require(winner < hap1.size(), "the scores must name a winner in the manifest");
    const auto rows = read_rows(path(prefix + ".block_genotypes.tsv"));
    for (std::size_t b = 0; b < rows.size(); ++b) {
        const auto called = std::minmax(hap1[winner][b], hap2[winner][b]);
        double same = 0.0;
        double other = 0.0;
        for (const auto& item : representatives) {
            (std::minmax(hap1[item.first][b], hap2[item.first][b]) == called ? same : other) += item.second;
        }
        const double posterior = std::stod(rows[b][11]);
        const double quality = std::stod(rows[b][12]);
        require(std::abs(posterior - same) <= 1e-9 && posterior > 0.0 && posterior <= 1.0,
                prefix + ": block " + rows[b][0] + " posterior " + rows[b][11] + " must be the winner genotype's set mass " +
                    std::to_string(same));
        if (other > 1e-9) {
            require(std::abs(quality + 10.0 * std::log10(other)) <= 1e-3 * std::max(1.0, quality),
                    prefix + ": block " + rows[b][0] + " GQ must be -10 log10 of the other genotypes' mass");
        } else {
            require(quality >= 80.0, prefix + ": block " + rows[b][0] + " with negligible rival mass must have a high GQ");
        }
        require((rows[b][13] == "true" || rows[b][13] == "false") &&
                    (rows[b][14] == "COMPLETE" || rows[b][14] == "CONDITIONAL"),
                "every block row must state its numerical separation and candidate scope");
    }
}

// The confidence calculation on sets where it matters; the C4-like fixtures above give posteriors of one.
void test_block_genotype_confidence() {
    using Row = std::vector<std::pair<std::size_t, std::size_t>>;
    // Winner 0.6 and a same-genotype rival 0.3 at block 0; a rival 0.1 with another genotype at block 0 only.
    const auto mixed = panvar::block_genotype_confidence(
        {std::log(0.6), std::log(0.3), std::log(0.1)}, {1e-9, 1e-9, 1e-9},
        {Row{{0, 1}, {2, 2}}, Row{{0, 1}, {2, 2}}, Row{{1, 1}, {2, 2}}});
    require(std::abs(mixed[0].posterior - 0.9) <= 1e-12 && std::abs(mixed[0].quality - 10.0) <= 1e-9 &&
                mixed[0].numerically_separated,
            "a 0.1 rival with another genotype must give posterior 0.9 and GQ 10");
    require(mixed[1].posterior == 1.0 && std::isinf(mixed[1].quality) && mixed[1].numerically_separated,
            "a block on which every candidate agrees must have posterior 1 and unbounded GQ");
    // A rival with another genotype within floating-point error of the winner is not numerically separated.
    const auto tied = panvar::block_genotype_confidence({-5.0, -5.0 - 1e-12}, {1e-9, 1e-9},
                                                        {Row{{0, 0}}, Row{{0, 1}}});
    require(!tied[0].numerically_separated && std::abs(tied[0].posterior - 0.5) <= 1e-9,
            "a genotype tied within numerical error must be reported as not separated, at posterior one half");
    // Rival mass of about 1e-300: in linear space the posterior is exactly one and the GQ infinite; in log space
    // the GQ is about 3000.
    const auto tiny = panvar::block_genotype_confidence({0.0, -300.0 * std::log(10.0)}, {1e-9, 1e-9},
                                                        {Row{{0, 0}}, Row{{1, 1}}});
    require(tiny[0].posterior == 1.0 && std::isfinite(tiny[0].quality) && std::abs(tiny[0].quality - 3000.0) <= 1e-6,
            "a negligible rival must still give a finite GQ, got " + std::to_string(tiny[0].quality));
    std::cout << "  block confidence: posterior " << mixed[0].posterior << " GQ " << mixed[0].quality
              << "; tie separated=" << (tied[0].numerically_separated ? "true" : "false") << "; 1e-300 rival GQ "
              << tiny[0].quality << "\n";
}

void test_prior_only_phase_is_reported() {
    const TinyLocus locus = write_tiny_locus("tiny_duplicate.gfa", 350, 340, true);
    simulate(locus.haplotypes[0], locus.haplotypes[1], path("duplicate_r1.fq"), path("duplicate_r2.fq"), 321);
    panvar::GenotypeCommandOptions options =
        tiny_options("prior_only", path("duplicate_r1.fq"), path("duplicate_r2.fq"), "tiny_duplicate.gfa");
    options.proposal_state_limit = 1048576;
    const panvar::GenotypeRunStatus status = panvar::run_genotype(options);
    const auto reported = read_status(path("prior_only.mosaic_status.tsv"));
    require(status == panvar::GenotypeRunStatus::Complete,
            std::string("the prior-only fixture must complete, got ") + panvar::genotype_run_status_name(status) +
                ": " + (reported.count("note") ? reported.at("note") : ""));
    require(reported.at("phase_support") == "PRIOR_ONLY" && std::stoul(reported.at("fragment_contenders")) > 0,
            "a winner chosen among fragment contenders by the prior must be reported PRIOR_ONLY, got " +
                reported.at("phase_support"));
    const auto rows = read_rows(path("prior_only.block_genotypes.tsv"));
    std::size_t supported = 0;
    for (const auto& row : rows) {
        require(row[10] == "PRIOR_ONLY", "every block row must carry the run's phase support");
        if (row[9] == "READS") ++supported;
    }
    require(supported == rows.size() && reported.at("blocks_genotype_read_supported") == std::to_string(rows.size()),
            "block genotypes on which every fragment contender agrees must remain read-supported");
    check_block_posteriors("prior_only");
    std::cout << "  prior-only phase: call made, phase PRIOR_ONLY, " << supported << "/" << rows.size()
              << " block genotypes read-supported, " << reported.at("fragment_contenders") << " contenders\n";
}

} // namespace

// K truncates one proposal ordering, so the development harness ranks every prefix of a single proposal
// instead of rerunning it. Each prefix must reproduce a standalone run with that K byte for byte -- its
// candidate manifest, scores, block calls, haplotypes and status apart from stage timings -- including
// prefixes the production rule refuses because a pruned search found fewer candidates than they ask for.
void test_candidate_prefixes_match_standalone() {
    const std::string r1 = path("phaseable_r1.fq");
    const std::string r2 = path("phaseable_r2.fq");
    const auto slurp = [](const std::string& file, bool status_file, bool drop_proposal = false,
                          bool drop_last_column = false) {
        std::ifstream input(file);
        if (!input) return std::string("<absent>");
        std::string line;
        std::string text;
        while (std::getline(input, line)) {
            if (status_file && line.rfind("seconds_", 0) == 0) continue;
            // Process memory and the worker gauge describe the process, not the ranking.
            if (status_file && (line.rfind("resources_", 0) == 0 || line.rfind("threads_peak_", 0) == 0)) continue;
            // Prefixes share one scoring pass over their union, so its work counters describe that pass.
            if (status_file) {
                const std::string key = line.substr(0, line.find('\t'));
                static const std::set<std::string> shared_work = {
                    "candidate_bases_scanned", "piece_hits", "proposed_starts", "unique_starts",
                    "full_read_verifications", "accepted_mate_placements", "valid_fr_states",
                    "physical_origins", "statistical_state_classes",
                    // The union origin reuse's served work and the merged census count describe the pass too.
                    "reused_physical_origins", "reused_state_aggregates", "reused_latent_alignment_pairs",
                    "reused_valid_fr_states", "merged_state_classes",
                    // So do the start sorter's list counts: the ladder's pass orders the union's lists.
                    "start_sort_out_of_place", "start_sort_in_place", "start_sort_fallbacks",
                    // The exact-path prior cache is shared by every ranked set, so a prefix view counts
                    // the paths interned while ranking the prefixes before it, where a standalone run
                    // starts empty. Deliberately cumulative: the invariants these prove are run-level --
                    // evaluations equal distinct paths, and an audit triggers no second full set of them
                    // -- and a per-set counter would hide exactly the reuse being tested.
                    "distinct_haplotype_paths", "haploid_prior_evaluations", "diploid_prior_lookups",
                    // candidates x fragments of the SHARED scoring pass, which covers the union of every
                    // prefix. A standalone run of one K scores only its own set, so this differs for
                    // every prefix, not just K = 0 -- it describes the pass, not the candidate set.
                    "candidate_fragment_evaluations",
                    // Only a refused run emits these, and they measure the work the shared budget
                    // allowed before it stopped -- a property of the pass, not of the candidate set.
                    "observed_piece_hits", "observed_physical_origins",
                    // How the shared pass's verifications split between the pair-HMM and the exact
                    // feasibility prefilter. Same reasoning as full_read_verifications above: the
                    // ladder's pass covers the union of every prefix, so it verifies roughly twice
                    // what a standalone K = 0 run does (26,412 against 13,186 on this fixture), and
                    // the split moves with it. The CALL must still be identical, which is what the
                    // rest of this comparison checks.
                    "dp_evaluations", "prefilter_rejections",
                    // Likewise the exact emission cache's work counters, whenever one is armed.
                    // emission_cache_bytes is deliberately NOT here: that is a declared option and
                    // both runs must agree on it.
                    "emission_cache_dp_requests", "emission_cache_dp_hits",
                    "emission_cache_dp_misses", "emission_cache_entries",
                    "emission_cache_payload_bytes", "emission_cache_executed_prefilters",
                    "emission_cache_executed_dps", "emission_cache_entries_not_retained",
                    // The explicit memory accounts describe the pass's own buffers, which the
                    // ladder sizes for the union of every prefix.
                    "scorer_mass_matrix_bytes", "scorer_haplotype_sequence_bytes",
                    "scorer_form_sequence_bytes", "scorer_piece_index_bytes",
                    "scorer_emission_cache_table_bytes", "scorer_emission_cache_arena_bytes",
                    "scorer_placement_rows_peak", "scorer_placement_capacity_bytes_peak",
                    "scorer_proposed_starts_peak", "scorer_proposed_starts_capacity_bytes_peak"};
                if (shared_work.count(key) != 0) continue;
                // K = 0 alone is asymmetric, and legitimately so: a standalone K = 0 run requests no
                // mosaics, no audit and no prefix, so it builds no beam at all, while the K = 0 VIEW of a
                // ladder shares the beam its larger prefixes needed. The call must still be identical --
                // that is what this gate is for -- but the proposal diagnostics describe a search that
                // happened in one run and not in the other. Dropped for K = 0 only.
                static const std::set<std::string> proposal_diagnostics = {
                    "candidate_scope", "candidate_set_conditional", "mosaic_candidates",
                    "mosaic_candidates_found", "mosaic_target_reached", "proposal_state_limit",
                    "proposal_final_states", "proposal_states_pruned", "proposal_search_exhaustive",
                    "beam_extensions", "candidates_total", "beam_diplotypes_found",
                    // The full-beam audit describes the same absent search: a standalone K = 0 run has
                    // no beam, so it has no prefix to truncate, nothing outside K to audit, and a scope
                    // of PANEL_ONLY. Dropped for K = 0 ONLY -- at every other K both runs build a beam
                    // and these must agree exactly, which is what keeps the audit inside this gate.
                    "primary_prefix_k", "full_beam_candidates", "primary_prefix_truncated",
                    "full_beam_audit_attempted", "candidates_exactly_audited", "audit_completed",
                    "final_candidate_scope",
                    // The proposal resource contract describes the same absent search. A standalone
                    // K = 0 run builds no proposer, so it has no predicted peak, no observed peak and
                    // nothing to gate against the ceiling. Dropped for K = 0 ONLY -- wherever both runs
                    // build a beam these must agree exactly, which is what keeps the contract inside
                    // this gate rather than outside it.
                    "proposal_payload_ceiling_bytes", "proposal_payload_predicted_bytes",
                    "proposal_payload_observed_bytes", "proposal_payload_predicted_rolling_bytes",
                    "proposal_payload_predicted_ancestry_bytes", "proposal_payload_predicted_pair_bytes",
                    "proposal_payload_predicted_index_bytes",
                    "proposal_payload_predicted_child_records",
                    "proposal_payload_observed_child_records"};
                if (drop_proposal && proposal_diagnostics.count(key) != 0) continue;
            }
            if (drop_last_column) {
                // The per-block candidate_scope column. A standalone K = 0 run is PANEL_ONLY while the
                // K = 0 view of a ladder shares the beam its larger prefixes needed, so this column
                // legitimately differs -- and only this column. Everything else on the row, the
                // genotype, posterior, GQ and phase, must still match byte for byte.
                // It is the 15th column; the representation-filter columns follow it and are compared.
                std::vector<std::string> cells;
                std::stringstream fields(line);
                for (std::string cell; std::getline(fields, cell, '\t');) cells.push_back(cell);
                if (cells.size() > 14) cells.erase(cells.begin() + 14);
                line.clear();
                for (std::size_t c = 0; c < cells.size(); ++c) line += (c == 0 ? "" : "\t") + cells[c];
            }
            if (!line.empty() && line[0] == '>') line = line.substr(line.find(' '));   // header names the prefix
            text += line + '\n';
        }
        return text;
    };
    const std::vector<std::size_t> prefixes = {0, 2, 4, 16, 64};
    for (const std::size_t state_limit : {std::size_t{64}, std::size_t{1048576}}) {
        const std::string base = "prefixes_" + std::to_string(state_limit);
        panvar::GenotypeCommandOptions ladder = tiny_options(base, r1, r2, "tiny_phaseable.gfa");
        ladder.proposal_state_limit = state_limit;
        ladder.max_mosaic_candidates = prefixes.back();
        ladder.development_candidate_prefixes = prefixes;
        panvar::run_genotype(ladder);
        std::size_t refused = 0;
        std::size_t complete = 0;
        for (const std::size_t k : prefixes) {
            const std::string name = base + "_standalone_K" + std::to_string(k);
            panvar::GenotypeCommandOptions standalone = tiny_options(name, r1, r2, "tiny_phaseable.gfa");
            standalone.proposal_state_limit = state_limit;
            standalone.max_mosaic_candidates = k;
            const panvar::GenotypeRunStatus status = panvar::run_genotype(standalone);
            if (status == panvar::GenotypeRunStatus::CandidateIncomplete) ++refused;
            if (status == panvar::GenotypeRunStatus::Complete) ++complete;
            const std::string extracted = path(base + ".K" + std::to_string(k));
            for (const std::string suffix : {".genotypes.tsv", ".block_genotypes.tsv", ".mosaic_candidates.tsv",
                                             ".mosaic_scores.tsv", ".mosaic_status.tsv", ".mosaic_haplotypes.fa"}) {
                const bool status_file = suffix == ".mosaic_status.tsv";
                const bool block_file = suffix == ".block_genotypes.tsv";
                require(slurp(extracted + suffix, status_file, k == 0, k == 0 && block_file) ==
                            slurp(path(name) + suffix, status_file, k == 0, k == 0 && block_file),
                        "prefix K=" + std::to_string(k) + " at state limit " + std::to_string(state_limit) +
                            " must reproduce the standalone " + suffix);
            }
            if (k == 0) {
                // Dropping the proposal diagnostics above is only sound if the two runs really did take
                // different paths. Assert that they did, so the relaxation cannot quietly widen into
                // hiding a difference that matters.
                const auto standalone_status = read_status(path(name) + ".mosaic_status.tsv");
                const auto view_status = read_status(extracted + ".mosaic_status.tsv");
                require(standalone_status.at("candidate_scope") == "PANEL_ONLY" &&
                            standalone_status.count("beam_extensions") == 0,
                        "a standalone K = 0 run requests no search and must build no beam");
                require(view_status.count("candidate_scope") == 1 &&
                            view_status.at("candidate_scope") != "PANEL_ONLY" &&
                            view_status.count("beam_extensions") == 1,
                        "the K = 0 view of a ladder shares the beam its larger prefixes required");
            }
        }
        require(refused == 0 && complete == prefixes.size(),
                "every prefix, K = 0 and mosaic-deficient ones included, must call from the complete panel");
        if (state_limit == 64) {
            require(read_status(path(base + ".K64.mosaic_status.tsv")).at("mosaic_target_reached") == "false",
                    "fixture: at 64 prefix pairs the largest prefix must be mosaic-deficient");
        }
        std::cout << "  candidate prefixes at state limit " << state_limit << ": " << complete << " complete, "
                  << refused << " refused, each identical to its standalone run\n";
    }
}

// Phase support from a declared audit set (review round 6). Without an audit, read-supported phase is only
// READS_WITHIN_CANDIDATES. With one, READS_AUDITED needs the same winner, the same block genotypes and
// read-supported phase in the larger set; a claim the larger set does not bear out is withdrawn.
void test_phase_audit() {
    panvar::GenotypeCommandOptions phaseable = tiny_options(
        "audit_phaseable", path("phaseable_r1.fq"), path("phaseable_r2.fq"), "tiny_phaseable.gfa");
    phaseable.proposal_state_limit = 1048576;
    phaseable.max_mosaic_candidates = 4;
    phaseable.phase_audit_candidates = 64;
    require(panvar::run_genotype(phaseable) == panvar::GenotypeRunStatus::Complete,
            "the audited phaseable run must complete");
    const auto confirmed = read_status(path("audit_phaseable.mosaic_status.tsv"));
    require(confirmed.at("phase_support_within_candidates") == "READS_WITHIN_CANDIDATES" &&
                confirmed.at("phase_audit_status") == "COMPLETE" && confirmed.at("phase_audit_same_winner") == "true" &&
                confirmed.at("phase_audit_same_block_genotypes") == "true" &&
                confirmed.at("phase_support") == "READS_AUDITED",
            "a phaseable truth's read-supported phase must be confirmed by the larger audit set, got " +
                confirmed.at("phase_support"));
    for (const auto& row : read_rows(path("audit_phaseable.block_genotypes.tsv"))) {
        require(row[10] == "READS_AUDITED", "block genotypes must carry the audited phase claim");
    }
    check_block_posteriors("audit_phaseable");

    // A claim the audit cannot check is withdrawn. Which width makes the audit unreachable changed with
    // the grouped proposer, and the new behaviour is better rather than broken: 64 prefix pairs now
    // searches this locus EXHAUSTIVELY, and an exhaustive search satisfies the audit even when it finds
    // fewer mosaics than were asked for, because it has proved there are no others to find. Only a
    // genuinely pruned width leaves the audit unreached, so phase is UNRESOLVED there -- while the
    // block genotypes are still emitted.
    panvar::GenotypeCommandOptions narrow =
        tiny_options("audit_unreachable", path("phaseable_r1.fq"), path("phaseable_r2.fq"), "tiny_phaseable.gfa");
    narrow.proposal_state_limit = 4;
    narrow.max_mosaic_candidates = 1;
    narrow.phase_audit_candidates = 64;
    require(panvar::run_genotype(narrow) == panvar::GenotypeRunStatus::Complete,
            "the primary set must still call when its audit set cannot be reached");
    const auto withdrawn = read_status(path("audit_unreachable.mosaic_status.tsv"));
    require(withdrawn.at("phase_support_within_candidates") == "READS_WITHIN_CANDIDATES" &&
                withdrawn.at("phase_audit_status") == "MOSAIC_TARGET_UNREACHED" &&
                withdrawn.at("phase_support") == "UNRESOLVED" && exists(path("audit_unreachable.genotypes.tsv")),
            "an unreachable audit must leave phase UNRESOLVED and still emit genotypes, got " +
                withdrawn.at("phase_support"));

    // Every branch of the decision, including the ones no small fixture reaches: an audit that completes with
    // the same call but a phase tie, and one that changes the winner or the block genotypes.
    const auto decide = panvar::audited_phase_support;
    require(decide("READS_WITHIN_CANDIDATES", true, true, true, "READS_WITHIN_CANDIDATES") == "READS_AUDITED" &&
                decide("READS_WITHIN_CANDIDATES", true, true, true, "PRIOR_ONLY") == "PRIOR_ONLY" &&
                decide("READS_WITHIN_CANDIDATES", true, false, true, "READS_WITHIN_CANDIDATES") == "UNRESOLVED" &&
                decide("READS_WITHIN_CANDIDATES", true, true, false, "READS_WITHIN_CANDIDATES") == "UNRESOLVED" &&
                decide("READS_WITHIN_CANDIDATES", false, false, false, "NA") == "UNRESOLVED" &&
                decide("PRIOR_ONLY", true, true, true, "READS_WITHIN_CANDIDATES") == "PRIOR_ONLY" &&
                decide("UNRESOLVED", true, true, true, "READS_WITHIN_CANDIDATES") == "UNRESOLVED",
            "the audited phase decision must confirm only an agreeing, read-supported audit");

    bool refused = false;
    try {
        panvar::GenotypeCommandOptions invalid = narrow;
        invalid.max_mosaic_candidates = 4;
        invalid.phase_audit_candidates = 2;
        panvar::run_genotype(invalid);
    } catch (const std::runtime_error&) {
        refused = true;
    }
    require(refused, "an audit set smaller than the primary set must be refused");
    std::cout << "  phase audit: exhaustive phaseable " << confirmed.at("phase_support") << "; K=1 at 64 states "
              << withdrawn.at("phase_support_within_candidates") << " -> " << withdrawn.at("phase_support")
              << " (audit " << withdrawn.at("phase_audit_status") << ")\n";
}

// A true panel-only run builds no beam at all, and must still make exactly the call the former K = 0
// path made with the beam present. Both paths exist in this binary: a run with development prefixes
// {0, 8} computes the beam, and its K0 outputs ARE the former path; a run with K = 0, no audit and no
// prefix takes the new beam-free path. Whole files cannot match -- timings, the proposal diagnostics and
// the scope column legitimately differ, and the FASTA headers carry the run's own prefix -- so this
// compares what must not change: the genotype call and the haplotypes it selected.
void test_panel_only_path() {
    const TinyLocus locus = write_tiny_locus();
    write_library(path("lib.txt"));
    simulate(locus.haplotypes[1], locus.haplotypes[2], path("po_1.fq"), path("po_2.fq"), 99);

    const auto slurp = [](const std::string& file) {
        std::ifstream in(file);
        std::string all, line;
        while (std::getline(in, line)) {
            all += line;
            all += '\n';
        }
        return all;
    };
    const auto sequences_only = [](const std::string& fasta) {
        std::string out;
        std::size_t at = 0;
        while (at < fasta.size()) {
            const std::size_t end = fasta.find('\n', at);
            const std::string line = fasta.substr(at, end == std::string::npos ? end : end - at);
            if (!line.empty() && line.front() != '>') {
                out += line;
                out += '\n';
            }
            if (end == std::string::npos) break;
            at = end + 1;
        }
        return out;
    };

    panvar::GenotypeCommandOptions with_beam = tiny_options("po_beam", path("po_1.fq"), path("po_2.fq"));
    with_beam.max_mosaic_candidates = 8;
    with_beam.development_candidate_prefixes = {0, 8};
    require(panvar::run_genotype(with_beam) == panvar::GenotypeRunStatus::Complete,
            "the reference run that builds a beam must complete");

    panvar::GenotypeCommandOptions panel_only = tiny_options("po_only", path("po_1.fq"), path("po_2.fq"));
    panel_only.max_mosaic_candidates = 0;
    panel_only.phase_audit_candidates = 0;
    require(panvar::run_genotype(panel_only) == panvar::GenotypeRunStatus::Complete,
            "a panel-only run must complete without building a beam");

    const std::map<std::string, std::string> status = read_status(path("po_only.mosaic_status.tsv"));
    require(status.at("candidate_scope") == "PANEL_ONLY",
            "a run with no search reports PANEL_ONLY, which is neither COMPLETE nor CONDITIONAL");
    require(status.at("phase_support") == "UNRESOLVED",
            "a panel-only run carries no audit and claims no read-supported phase");
    require(std::stod(status.at("seconds_proposal")) == 0.0,
            "no proposal ran, so the proposal took no time");
    require(status.count("beam_extensions") == 0 && status.count("mosaic_candidates") == 0,
            "a run with no beam must report no beam diagnostics at all");
    require(status.at("panel_baseline_complete") == "true",
            "every panel diplotype is still scored");

    require(slurp(path("po_beam.K0.genotypes.tsv")) == slurp(path("po_only.genotypes.tsv")),
            "panel-only must make exactly the call the former K=0 path made");
    require(sequences_only(slurp(path("po_beam.K0.mosaic_haplotypes.fa"))) ==
                sequences_only(slurp(path("po_only.mosaic_haplotypes.fa"))),
            "panel-only must select exactly the haplotypes the former K=0 path selected");
}

// The block-consensus rule, gated where it lives. When no whole-locus winner can be separated, a block
// is called only when every tied candidate agrees on its unordered allele pair.
//
// Tested as a pure function rather than through a contrived run: forcing two DIFFERENT genotypes to tie
// end to end would need their total scores to coincide within ~1e-9, which on simulated reads is a
// property of the tie threshold rather than of this rule. The cases that are reachable end to end --
// two phasings of one genotype (agreement everywhere, no phased haplotype), a resource refusal and an
// incomplete candidate search -- are asserted in the run-level tests instead:
// test_off_panel_phase_is_not_claimed_from_reads covers the first, and test_fail_closed_statuses covers
// the other two, where RescoreIncomplete and CandidateIncomplete exclude BLOCK_CONSENSUS by identity
// and no_ranking_or_call already requires that nothing was called.
void test_block_consensus_rule() {
    using Row = std::vector<std::pair<std::size_t, std::size_t>>;

    // Unanimous at every block: all called.
    const auto all_agree = panvar::block_consensus({Row{{0, 1}, {2, 2}, {3, 4}},
                                                    Row{{0, 1}, {2, 2}, {3, 4}}});
    require(all_agree.size() == 3, "consensus must report one entry per block");
    for (const auto& call : all_agree) require(call.agreed, "unanimous blocks must all be called");
    require(all_agree[0].allele1 == 0 && all_agree[0].allele2 == 1 && all_agree[2].allele1 == 3,
            "an agreed block must carry the pair every tied candidate holds");

    // Disagreement at exactly one block: only that block is no-called.
    const auto one_differs = panvar::block_consensus({Row{{0, 1}, {2, 2}, {3, 4}},
                                                      Row{{0, 1}, {2, 5}, {3, 4}}});
    require(one_differs[0].agreed && !one_differs[1].agreed && one_differs[2].agreed,
            "only the block on which the tied candidates disagree may be no-called");
    require(one_differs[1].allele1 == 0 && one_differs[1].allele2 == 0,
            "a no-called block must carry no genotype, not the first candidate's");

    // Order independence: reversing the candidate list cannot change the answer.
    const auto reversed = panvar::block_consensus({Row{{0, 1}, {2, 5}, {3, 4}},
                                                   Row{{0, 1}, {2, 2}, {3, 4}}});
    require(reversed.size() == one_differs.size(), "reversal must not change the block count");
    for (std::size_t b = 0; b < reversed.size(); ++b) {
        require(reversed[b].agreed == one_differs[b].agreed &&
                    reversed[b].allele1 == one_differs[b].allele1 &&
                    reversed[b].allele2 == one_differs[b].allele2,
                "reversing the tied candidate order must give identical consensus");
    }

    // A single winner agrees with itself: the unique-winner path is unchanged by this rule.
    const auto single = panvar::block_consensus({Row{{7, 9}, {1, 1}}});
    require(single.size() == 2 && single[0].agreed && single[1].agreed && single[0].allele1 == 7,
            "one candidate must agree at every block");

    // No candidates is no consensus, not vacuous agreement.
    require(panvar::block_consensus({}).empty(), "an empty tied set must report no consensus at all");

    // Three candidates, one dissenting at one block.
    const auto three = panvar::block_consensus({Row{{0, 0}, {1, 2}}, Row{{0, 0}, {1, 2}},
                                                Row{{0, 0}, {2, 2}}});
    require(three[0].agreed && !three[1].agreed,
            "a single dissenter among three must no-call its block and leave the rest called");

    std::cout << "  block consensus: unanimous called, one dissent no-called, order-independent\n";
}

// The tie-triggered full-beam exact audit. A K-truncated prefix is never grounds for consensus: the
// candidates below K were ordered by proposal mass but never exactly scored, so their order proves
// nothing about whether one of them wins or disagrees blockwise. When the primary ties, the search was
// exhaustive and the beam held candidates outside K, the complete retained beam is scored exactly and
// the OUTPUT COMES FROM THAT RANKING.
//
// This locus is the trigger: an exhaustive search finding 26 mosaics against the default K of 16, whose
// two phasings of one genotype cannot be separated.
//
// Directly gated here: the audit fires, the tables describe the audited set, a non-exhaustive proposal
// never audits, and a unique primary winner is untouched. The cases where an omitted candidate outscores
// the tie or disagrees at one block cannot be manufactured end to end without controlling exact scores;
// what makes them impossible is that the decision is derived from the audited set, which is asserted
// below by requiring the scored candidate count to equal the audited count.
void test_full_beam_audit() {
    panvar::GenotypeCommandOptions options =
        tiny_options("audit", path("unphaseable_r1.fq"), path("unphaseable_r2.fq"));
    options.proposal_state_limit = 1048576;   // exhaustive
    // max_mosaic_candidates stays at the fixture default of 16, below the 26 this locus yields.
    const panvar::GenotypeRunStatus status = panvar::run_genotype(options);
    const auto reported = read_status(path("audit.mosaic_status.tsv"));

    require(reported.at("proposal_search_exhaustive") == "true",
            "fixture: the audit only runs on an exhaustive search");
    require(reported.at("primary_prefix_truncated") == "true",
            "fixture: the beam must hold candidates outside K, or there is nothing to audit");
    require(reported.at("full_beam_audit_attempted") == "true",
            std::string("a tie on a truncated exhaustive beam must trigger the audit, got status ") +
                panvar::genotype_run_status_name(status));
    require(reported.at("audit_completed") == "true", "the audit must rank to completion here");

    // The eight concepts stay distinct: candidate_set_conditional still describes the PREFIX, while
    // final_candidate_scope describes the set that decided the call.
    require(reported.at("candidate_set_conditional") == "true",
            "the prefix was truncated, so its own scope is still conditional");
    require(reported.at("final_candidate_scope") == "COMPLETE",
            "the audited full beam of an exhaustive search omits nothing");

    // The tables must describe the audited set, not the prefix beside it.
    const std::size_t audited = std::stoul(reported.at("candidates_exactly_audited"));
    const std::size_t full_beam = std::stoul(reported.at("full_beam_candidates"));
    const std::size_t prefix_k = std::stoul(reported.at("primary_prefix_k"));
    const std::size_t panel = std::stoul(reported.at("panel_diplotypes"));
    require(full_beam > prefix_k, "fixture: the full beam must exceed K");
    // An identity, not an inequality: the audited set is exactly the panel baseline plus every mosaic
    // the beam retained. `>=` would also be satisfied by an audit that scored the panel twice, or that
    // covered only part of the beam and made the count up elsewhere.
    require(audited == panel + full_beam,
            "the audited set must be exactly the panel baseline plus the whole retained beam, got " +
                std::to_string(audited) + " against " + std::to_string(panel) + " + " +
                std::to_string(full_beam));
    require(read_rows(path("audit.mosaic_candidates.tsv")).size() == audited,
            "the candidate table must list the audited set, not the top-K prefix");

    // Whichever way the audited ranking falls, the contract holds.
    if (status == panvar::GenotypeRunStatus::BlockConsensus) {
        require(reported.at("phase_support") == "UNRESOLVED" &&
                    !exists(path("audit.mosaic_haplotypes.fa")),
                "an audited tie must report phase UNRESOLVED and emit no phased haplotype");
        const std::size_t called = std::stoul(reported.at("consensus_blocks_called"));
        const std::size_t no_called = std::stoul(reported.at("consensus_blocks_no_called"));
        require(called + no_called == read_rows(path("audit.block_genotypes.tsv")).size(),
                "every block must be called or explicitly no-called");
    } else {
        require(status == panvar::GenotypeRunStatus::Complete,
                std::string("an audited ranking must end COMPLETE or BLOCK_CONSENSUS, got ") +
                    panvar::genotype_run_status_name(status));
        require(exists(path("audit.mosaic_haplotypes.fa")),
                "a uniquely audited winner is an ordinary complete call, phased haplotype included");
    }

    // A unique primary winner is left alone. This is the case the audit must NOT touch: the search is
    // exhaustive and the prefix IS truncated, so both of the audit's other preconditions hold and only
    // separation keeps it from running. Without this fixture, "the audit changes nothing when the
    // primary winner is unique" would rest on reading the condition rather than on running it.
    panvar::GenotypeCommandOptions separated = tiny_options(
        "audit_separated", path("phaseable_r1.fq"), path("phaseable_r2.fq"), "tiny_phaseable.gfa");
    separated.proposal_state_limit = 1048576;   // exhaustive
    separated.max_mosaic_candidates = 1;        // below what this locus yields, so the prefix truncates
    const panvar::GenotypeRunStatus separated_status = panvar::run_genotype(separated);
    const auto separated_reported = read_status(path("audit_separated.mosaic_status.tsv"));
    require(separated_reported.at("proposal_search_exhaustive") == "true" &&
                separated_reported.at("primary_prefix_truncated") == "true",
            "fixture: this case is only informative while the search is exhaustive AND K truncates");
    require(separated_status == panvar::GenotypeRunStatus::Complete,
            "fixture: the separated locus must still produce a unique winner");
    require(separated_reported.at("full_beam_audit_attempted") == "false",
            "a uniquely separated primary winner must not be audited: there is no tie to resolve");
    require(separated_reported.at("candidates_exactly_audited") == "0" &&
                separated_reported.at("audit_completed") == "false",
            "an audit that never ran must report no audited candidates, not a silent full-beam pass");
    // ...and the original, pre-audit result stands: the PREFIX decided, so the tables describe the
    // prefix and the scope reported is the prefix's own.
    require(separated_reported.at("final_candidate_scope") == "CONDITIONAL",
            "with the prefix deciding, the final scope is the prefix's own conditional scope");
    require(read_rows(path("audit_separated.mosaic_candidates.tsv")).size() ==
                std::stoul(separated_reported.at("panel_diplotypes")) +
                    std::stoul(separated_reported.at("primary_prefix_k")),
            "an unaudited run must list the top-K prefix, exactly as it did before the audit existed");
    require(exists(path("audit_separated.mosaic_haplotypes.fa")),
            "a unique winner is an ordinary complete call, phased haplotype included");

    // A non-exhaustive proposal never audits, and never reaches consensus.
    panvar::GenotypeCommandOptions pruned =
        tiny_options("audit_pruned", path("unphaseable_r1.fq"), path("unphaseable_r2.fq"));
    pruned.proposal_state_limit = 2;   // forces pruning, so the search is not exhaustive
    const panvar::GenotypeRunStatus pruned_status = panvar::run_genotype(pruned);
    const auto pruned_reported = read_status(path("audit_pruned.mosaic_status.tsv"));
    if (pruned_reported.count("proposal_search_exhaustive") != 0) {
        require(pruned_reported.at("proposal_search_exhaustive") == "false",
                "fixture: a state limit of 2 must prune this locus");
        require(pruned_reported.at("full_beam_audit_attempted") == "false",
                "a non-exhaustive proposal must never be audited: it omits candidates it never proposed");
        require(pruned_status != panvar::GenotypeRunStatus::BlockConsensus,
                "a non-exhaustive proposal must never reach block consensus");
    }

    std::cout << "  full-beam audit: " << prefix_k << " of " << full_beam << " in the prefix, " << audited
              << " audited, " << panvar::genotype_run_status_name(status) << "\n";
}

// Block consensus, end to end. The pure-function gate above proves the rule; this proves a run can
// actually reach it, which nothing else does: every large-locus candidate set is CONDITIONAL, and the
// contract refuses consensus there, so without this fixture the whole path would ship undemonstrated.
//
// The precondition is a candidate set that is NOT conditional, which needs an exhaustive search whose
// every found mosaic is scored -- K above the number found, not the default 16 against 26. That is a
// real configuration, not a contrivance: it is what "score everything the search produced" means.
void test_block_consensus_end_to_end() {
    panvar::GenotypeCommandOptions options = tiny_options(
        "consensus", path("unphaseable_r1.fq"), path("unphaseable_r2.fq"));
    options.proposal_state_limit = 1048576;   // exhaustive
    options.max_mosaic_candidates = 64;       // above the 26 this locus yields, so nothing is truncated
    const panvar::GenotypeRunStatus status = panvar::run_genotype(options);
    const auto reported = read_status(path("consensus.mosaic_status.tsv"));

    require(reported.at("candidate_set_conditional") == "false",
            "fixture: the candidate set must not be conditional, or consensus is refused by contract");
    require(status == panvar::GenotypeRunStatus::BlockConsensus,
            std::string("two inseparable phasings of one genotype must reach block consensus, got ") +
                panvar::genotype_run_status_name(status) + ": " +
                (reported.count("note") ? reported.at("note") : ""));
    require(std::stoul(reported.at("winner_set_size")) > 1,
            "fixture: block consensus requires a genuine tie, not a single winner");

    const std::size_t called = std::stoul(reported.at("consensus_blocks_called"));
    const std::size_t no_called = std::stoul(reported.at("consensus_blocks_no_called"));
    require(called > 0, "the tied candidates agree somewhere, so some block must be called");
    require(reported.at("phase_support") == "UNRESOLVED",
            "block consensus must report phase UNRESOLVED: the reads chose no whole-locus history");

    // The output contract: genotypes yes, phased haplotype never.
    require(exists(path("consensus.genotypes.tsv")) && exists(path("consensus.block_genotypes.tsv")),
            "block consensus must emit its block genotypes");
    require(!exists(path("consensus.mosaic_haplotypes.fa")),
            "block consensus must emit no phased haplotype sequence");

    const auto rows = read_rows(path("consensus.block_genotypes.tsv"));
    require(rows.size() == called + no_called,
            "every block must appear exactly once, called or explicitly no-called");
    std::size_t nocall_rows = 0;
    for (const auto& row : rows) {
        require(row[10] == "UNRESOLVED", "every block row must carry phase UNRESOLVED under consensus");
        if (row[3] == ".") {
            ++nocall_rows;
            require(row[7] == "NOCALL", "a no-called block must be marked NOCALL, not left blank");
        }
    }
    require(nocall_rows == no_called, "the status and the rows must agree on how many blocks were no-called");

    std::cout << "  block consensus end to end: " << reported.at("winner_set_size") << " tied, " << called
              << " blocks called, " << no_called << " no-called, no phased haplotype\n";
}

// The production resource contract: one fixed ceiling over the proposer's whole predicted peak,
// refused BEFORE the proposal is constructed. This is the production path, not the library call --
// what matters is that a refusal leaves no ranking and no call behind it, because a partial output
// beside a status claiming refusal is indistinguishable from a successful run.
void test_proposal_payload_ceiling() {
    const std::string r1 = path("phaseable_r1.fq");
    const std::string r2 = path("phaseable_r2.fq");

    // What the contract predicts for this locus, read from a run under the frozen ceiling.
    panvar::GenotypeCommandOptions baseline =
        tiny_options("mem_baseline", r1, r2, "tiny_phaseable.gfa");
    require(panvar::run_genotype(baseline) == panvar::GenotypeRunStatus::Complete,
            "fixture: the baseline run must complete under the frozen ceiling");
    const auto reported = read_status(path("mem_baseline.mosaic_status.tsv"));
    const std::uint64_t predicted = std::stoull(reported.at("proposal_payload_predicted_bytes"));
    const std::uint64_t observed = std::stoull(reported.at("proposal_payload_observed_bytes"));
    require(predicted > 0 && observed > 0, "both predicted and observed bytes must be reported");
    require(observed <= predicted,
            "the observed peak must not exceed the prediction the run refused against");
    require(std::stoull(reported.at("proposal_payload_ceiling_bytes")) == 268435456ULL,
            "the frozen production ceiling must be reported as declared");
    require(std::stoull(reported.at("proposal_payload_predicted_index_bytes")) > 0,
            "the temporary and index storage must be part of the reported prediction");

    // Exactly at the predicted peak: admitted.
    panvar::GenotypeCommandOptions at_ceiling = tiny_options("mem_at", r1, r2, "tiny_phaseable.gfa");
    at_ceiling.max_proposal_payload_bytes = predicted;
    require(panvar::run_genotype(at_ceiling) == panvar::GenotypeRunStatus::Complete,
            "a ceiling exactly at the predicted peak must admit the run");

    // One byte below: refused, with no partial output.
    panvar::GenotypeCommandOptions below = tiny_options("mem_below", r1, r2, "tiny_phaseable.gfa");
    below.max_proposal_payload_bytes = predicted - 1;
    require(panvar::run_genotype(below) == panvar::GenotypeRunStatus::CandidateIncomplete,
            "a ceiling one byte below the predicted peak must refuse the run");
    const auto refused = read_status(path("mem_below.mosaic_status.tsv"));
    require(refused.at("note").find("payload ceiling") != std::string::npos,
            "the refusal must name the ceiling it exceeded, got: " + refused.at("note"));
    require(no_ranking_or_call("mem_below"),
            "a run refused at the memory ceiling must leave no ranking and no call");

    std::cout << "  proposal payload ceiling: predicted " << predicted << " B, observed " << observed
              << " B; admitted at the prediction, refused one byte below, no partial output\n";
}

// A panel template leaves the panel only by an explicit, recorded decision. The curated model is a
// DIFFERENT model and has to say so: a different fingerprint, the excluded template named with the
// reason it would otherwise have carried, and a digest of the decision. A manifest that no longer
// matches its graph is refused rather than quietly ignored -- that silence is exactly how a stale
// manifest would travel with a changed graph and curate nothing while appearing to curate.
void test_curated_panel() {
    const std::string r1 = path("phaseable_r1.fq");
    const std::string r2 = path("phaseable_r2.fq");

    panvar::GenotypeCommandOptions whole = tiny_options("curate_whole", r1, r2, "tiny_phaseable.gfa");
    require(panvar::run_genotype(whole) == panvar::GenotypeRunStatus::Complete,
            "the uncurated panel must complete");
    const auto whole_status = read_status(path("curate_whole.mosaic_status.tsv"));
    require(whole_status.at("templates_curated_out") == "0" &&
                whole_status.at("curation_digest") == "none",
            "a run with no manifest must report no curation, not an absent field");

    // sampleD, not sampleA: sampleA is this fixture's reference path, and curating it out would remove
    // the coordinate frame rather than a panel row.
    const std::string manifest = path("curation.tsv");
    {
        std::ofstream file(manifest);
        file << "# one template, excluded deliberately\nsampleD#1#chr1:0-100\n";
    }
    panvar::GenotypeCommandOptions curated = tiny_options("curate_one", r1, r2, "tiny_phaseable.gfa");
    curated.curation_manifest_path = manifest;
    require(panvar::run_genotype(curated) == panvar::GenotypeRunStatus::Complete,
            "the curated panel must complete");
    const auto curated_status = read_status(path("curate_one.mosaic_status.tsv"));
    require(curated_status.at("templates_curated_out") == "1",
            "the curated run must count the exclusion");
    require(curated_status.at("curated_exclusion").rfind("sampleD#1#chr1:0-100", 0) == 0,
            "the curated run must name the excluded template");
    require(curated_status.at("curation_digest") != "none" &&
                curated_status.at("curation_digest").size() == 16,
            "a curated run must carry a digest of the curation decision");
    require(curated_status.at("model_fingerprint") != whole_status.at("model_fingerprint"),
            "curating a template out must change the model fingerprint: it is a different panel, and a "
            "fingerprint that did not move would make the two runs indistinguishable");

    // Drift between a manifest and its graph is an error, not a no-op.
    const std::string stale = path("curation_stale.tsv");
    {
        std::ofstream file(stale);
        file << "sampleZ#1#chr1:0-100\n";
    }
    panvar::GenotypeCommandOptions drifted = tiny_options("curate_stale", r1, r2, "tiny_phaseable.gfa");
    drifted.curation_manifest_path = stale;
    require(panvar::run_genotype(drifted) == panvar::GenotypeRunStatus::ModelIncomplete,
            "a manifest naming a template the graph does not have must refuse");
    require(no_ranking_or_call("curate_stale"), "a refused curation must leave no ranking or call");

    // An empty manifest is indistinguishable from no curation at all, so passing one is an error.
    const std::string empty = path("curation_empty.tsv");
    {
        std::ofstream file(empty);
        file << "# nothing but a comment\n";
    }
    panvar::GenotypeCommandOptions blank = tiny_options("curate_empty", r1, r2, "tiny_phaseable.gfa");
    blank.curation_manifest_path = empty;
    require(panvar::run_genotype(blank) == panvar::GenotypeRunStatus::ModelIncomplete,
            "an empty curation manifest must refuse rather than score the uncurated panel");
}

// Resource caps must not depend on how the work was partitioned. Every cap here is a running total
// compared as `total > cap`, so the boundary sits exactly at the total a run reaches: a cap equal to it
// accepts, one below refuses.
//
// The three caps reach that verdict by genuinely different routes, which is what makes their agreement
// informative rather than a formality -- and the routes are distinguished by CAP, not by thread count.
// max_piece_hits and max_physical_origins are shared budgets checked in flight, so at every thread count
// the refusal is raised by whichever worker observes the global total cross; with a cap one below the
// reached total that happens at the final drain. max_state_classes is never checked in flight, so its
// refusal can only come from the deferred check after the ordered counter merge. Identical verdicts
// across both routes, at every thread count, are evidence the shared budget and the merge agree.
//
// The totals are measured here rather than written in: they belong to this fixture's graph and reads,
// and a hardcoded number would sit nowhere near the boundary if either ever changed.
void test_resource_caps_are_thread_invariant() {
    const std::string r1 = path("phaseable_r1.fq");
    const std::string r2 = path("phaseable_r2.fq");

    panvar::GenotypeCommandOptions calibrate =
        tiny_options("caps_calibrate", r1, r2, "tiny_phaseable.gfa");
    require(panvar::run_genotype(calibrate) == panvar::GenotypeRunStatus::Complete,
            "fixture: the calibration run must complete under the default caps");
    const auto base = read_status(path("caps_calibrate.mosaic_status.tsv"));
    const std::uint64_t totals[3] = {std::stoull(base.at("piece_hits")),
                                     std::stoull(base.at("physical_origins")),
                                     std::stoull(base.at("statistical_state_classes"))};
    const char* cap_names[3] = {"max_piece_hits", "max_physical_origins", "max_state_classes"};
    for (const std::uint64_t total : totals) {
        require(total > 1, "fixture: every calibrated counter must be non-trivial");
    }
    const auto apply_cap = [](panvar::GenotypeCommandOptions& options, int which, std::uint64_t value) {
        if (which == 0) options.max_piece_hits = value;
        else if (which == 1) options.max_physical_origins = value;
        else options.max_state_classes = value;
    };

    // 1 and 2 exercise the shared budget, 4 raises the contention, and 64 is deliberately excessive for a
    // fixture with a handful of haplotypes: it must be clamped to the work available, not honoured. The
    // verdict has to be the same at all four.
    const std::size_t thread_counts[4] = {1, 2, 4, 64};
    // The declared parallel overshoot. piece_hits is published in batches, so a worker can be up to one
    // batch past the cap when it publishes; physical_origins is an exact atomic, so a worker can be at
    // most one increment past. Either way the overshoot is per worker and bounded, which is the whole
    // difference between a shared budget and a per-worker one.
    const std::uint64_t publish_batch = panvar::FragmentScoringOptions{}.piece_hit_publish_batch;

    for (int which = 0; which < 3; ++which) {
        std::vector<std::string> accepted_counters;
        std::vector<std::uint64_t> work_before_refusal;
        for (const std::size_t workers : thread_counts) {
            // Just above: a cap equal to the total accepts, and must still call.
            const std::string at_name =
                "caps_at_" + std::to_string(which) + "_t" + std::to_string(workers);
            panvar::GenotypeCommandOptions at = tiny_options(at_name, r1, r2, "tiny_phaseable.gfa");
            at.threads = workers;
            apply_cap(at, which, totals[which]);
            require(panvar::run_genotype(at) == panvar::GenotypeRunStatus::Complete,
                    std::string(cap_names[which]) + " at exactly the reached total must accept at " +
                        std::to_string(workers) + " thread(s)");
            const auto at_status = read_status(path(at_name) + ".mosaic_status.tsv");
            accepted_counters.push_back(at_status.at("piece_hits") + "/" +
                                        at_status.at("physical_origins") + "/" +
                                        at_status.at("statistical_state_classes"));

            // Just below: one less refuses, ends RESCORE_INCOMPLETE, and leaves no ranking or call.
            const std::string below_name =
                "caps_below_" + std::to_string(which) + "_t" + std::to_string(workers);
            panvar::GenotypeCommandOptions below =
                tiny_options(below_name, r1, r2, "tiny_phaseable.gfa");
            below.threads = workers;
            apply_cap(below, which, totals[which] - 1);
            require(panvar::run_genotype(below) == panvar::GenotypeRunStatus::RescoreIncomplete,
                    std::string(cap_names[which]) + " one below the reached total must refuse at " +
                        std::to_string(workers) + " thread(s)");
            require(read_status(path(below_name) + ".mosaic_status.tsv").at("note").find("resource cap") !=
                        std::string::npos,
                    "a refusal must name the exhausted cap");
            require(no_ranking_or_call(below_name),
                    "a refused run must leave no partial ranking or call, at any thread count");

            // The verdict alone does not prove the budget is shared: a per-worker budget refuses on the
            // same merged totals, just after roughly threads x cap work. So assert the work performed.
            const auto refused = read_status(path(below_name) + ".mosaic_status.tsv");
            if (which == 2) {
                // max_state_classes is a post-merge acceptance cap, not a breaker: classes are counted
                // after the per-fragment merge, so the pass runs to completion and is judged afterwards.
                // A refusal here therefore reports a whole run's work. That is the documented difference
                // from the two caps above, and this is where it is pinned.
                require(std::stoull(refused.at("observed_piece_hits")) == totals[0],
                        "a statistical-class refusal must report a complete pass, not a truncated one");
            } else {
                const std::string key =
                    which == 0 ? "observed_piece_hits" : "observed_physical_origins";
                const std::uint64_t overshoot =
                    which == 0 ? workers * publish_batch : static_cast<std::uint64_t>(workers);
                const std::uint64_t observed = std::stoull(refused.at(key));
                // An upper bound alone is satisfied by a silent zero, which is what a refusal reports if
                // the plumbing is wrong. The work performed must be real.
                require(observed > 0, "a refusal must report the work it performed, not a zero");
                require(observed <= (totals[which] - 1) + overshoot,
                        std::string("work before refusal must not exceed the cap plus the declared "
                                    "parallel overshoot for ") +
                            cap_names[which] + " at " + std::to_string(workers) +
                            " thread(s): " + std::to_string(observed) + " > " +
                            std::to_string((totals[which] - 1) + overshoot));
                work_before_refusal.push_back(observed);
            }
        }
        for (std::size_t i = 1; i < accepted_counters.size(); ++i) {
            require(accepted_counters[i] == accepted_counters.front(),
                    std::string("accepted runs must report identical counters at every thread count "
                                "for ") +
                        cap_names[which] + ": " + accepted_counters.front() + " vs " +
                        accepted_counters[i]);
        }
        // The shared-budget property itself. Work before refusal must not grow with the thread count
        // beyond the declared per-worker overshoot; under a per-worker budget it grows toward
        // threads x cap, and the 4- and 64-thread rows are what catch that.
        for (std::size_t i = 1; i < work_before_refusal.size(); ++i) {
            const std::uint64_t overshoot =
                thread_counts[i] * (which == 0 ? publish_batch : std::uint64_t{1});
            require(work_before_refusal[i] <= work_before_refusal.front() + overshoot,
                    std::string("work before refusal must not grow with the thread count for ") +
                        cap_names[which] + ": " + std::to_string(work_before_refusal.front()) +
                        " serial vs " + std::to_string(work_before_refusal[i]) + " at " +
                        std::to_string(thread_counts[i]) + " thread(s)");
        }

        // The runs above cannot actually tell a shared budget from a per-worker one. With the cap one
        // below the total no single worker ever reaches it, so every arrangement refuses at the
        // post-merge check having already done the whole pass -- 28,264 piece hits at 1, 2 and 64
        // threads alike. Both semantics produce that number.
        //
        // Setting the cap far below the total separates them. A shared budget stops near the cap
        // wherever the work sits; T per-worker budgets each spend a full cap, so the work grows toward
        // T x cap, which for this fixture is the entire pass.
        if (which != 2) {
            const std::uint64_t deep = totals[which] / 4;
            require(deep > 1, "fixture: the calibrated total must leave room for a quarter cap");
            for (const std::size_t workers : {std::size_t{1}, std::size_t{4}}) {
                const std::string name =
                    "caps_deep_" + std::to_string(which) + "_t" + std::to_string(workers);
                panvar::GenotypeCommandOptions run = tiny_options(name, r1, r2, "tiny_phaseable.gfa");
                run.threads = workers;
                apply_cap(run, which, deep);
                require(panvar::run_genotype(run) == panvar::GenotypeRunStatus::RescoreIncomplete,
                        std::string("a cap far below the reached total must refuse for ") +
                            cap_names[which] + " at " + std::to_string(workers) + " thread(s)");
                const auto status = read_status(path(name) + ".mosaic_status.tsv");
                const std::string key =
                    which == 0 ? "observed_piece_hits" : "observed_physical_origins";
                const std::uint64_t unit = which == 0 ? publish_batch : std::uint64_t{1};
                const std::uint64_t observed = std::stoull(status.at(key));
                require(observed > 0, "a refusal must report the work it performed, not a zero");
                require(observed <= deep + workers * unit,
                        std::string("work before refusal must stay within the cap plus the declared "
                                    "parallel overshoot for ") +
                            cap_names[which] + " at " + std::to_string(workers) +
                            " thread(s): " + std::to_string(observed) + " > " +
                            std::to_string(deep + workers * unit));

                // A refusal must also report WHO ran, not only how much work was done. Both fields were
                // previously left at their success-path defaults, so a refused run claimed a single
                // placement thread however many it actually used -- which is how a 30x ACOT refusal came
                // to report one worker while its own CPU time proved it had used two.
                const std::size_t reported_threads =
                    std::stoul(status.at("observed_threads_used_placement"));
                require(reported_threads >= 1 && reported_threads <= workers,
                        "a refusal must report a placement thread count it could actually have used");
                if (workers > 1) {
                    require(reported_threads > 1,
                            "a parallel refusal must not report the serial default: reporting one worker "
                            "for a run that used several is the defect this field replaces");
                }
                require(status.at("observed_threads_used_combination") == "0",
                        "a placement refusal must report zero combination workers rather than a default "
                        "of one: that phase never ran");
            }
        }
    }
}

// ---- the frozen final decision contract, as the pure function production calls ----------------------
// Synthetic pools: every candidate is its own sequence pair, scored as the caller would hand it over.
panvar::PoolCandidate pooled(std::vector<std::size_t> a, std::vector<std::size_t> b, double fragment, double prior) {
    const auto spell = [](const std::vector<std::size_t>& v) {
        std::string text = "h";
        for (const std::size_t x : v) text += std::to_string(x) + ".";
        return text;
    };
    panvar::PoolCandidate c;
    c.first_sequence = panvar::sequence_digest(spell(a));
    c.second_sequence = panvar::sequence_digest(spell(b));
    c.first = std::move(a);
    c.second = std::move(b);
    c.score.fragment_term = fragment;
    c.score.exposure_term = -100.0;
    c.score.log_likelihood = fragment - 100.0;
    c.score.absolute_term_sum = std::fabs(fragment) + 100.0;
    c.score.fragments_with_mass = 10;
    c.sequence_prior = prior;
    c.path_prior = prior;
    return c;
}

// Everything the decision says, independent of the order the pool was handed over in.
std::string render_decision(const panvar::PooledDecision& d) {
    std::ostringstream out;
    out << std::setprecision(17) << panvar::top_set_decision_name(d.decision) << '|' << d.phase_support << '|';
    for (const std::size_t g : d.top_groups) {
        for (const std::size_t m : d.groups[g]) {
            for (const std::size_t x : d.candidates[m].first) out << x << ',';
            out << '/';
            for (const std::size_t x : d.candidates[m].second) out << x << ',';
            out << ';';
        }
    }
    out << '|';
    for (const auto& pair : d.genotype) out << pair.first << '/' << pair.second << ',';
    out << '|';
    for (const std::size_t b : d.disagreeing_blocks) out << b << ',';
    const panvar::PoolCandidate& rep = d.candidates[d.representative[d.ranking.candidates.front().candidate_id]];
    out << "|rep:";
    for (const std::size_t x : rep.first) out << x << ',';
    out << '/';
    for (const std::size_t x : rep.second) out << x << ',';
    out << '|' << d.genotype_separation_gap << '|' << d.has_same_genotype_alternative << d.same_genotype_fragment_identical
        << '|' << d.same_genotype_prior_gap;
    return out.str();
}

panvar::PooledDecision decide_both_ways(std::vector<panvar::PoolCandidate> pool, std::size_t blocks,
                                        const std::string& label) {
    const panvar::PooledDecision forward = panvar::decide_pooled_genotype(pool, blocks, 10, 6);
    std::reverse(pool.begin(), pool.end());
    for (panvar::PoolCandidate& c : pool) {   // and each identity presented with its homologues exchanged
        std::swap(c.first, c.second);
        std::swap(c.first_sequence, c.second_sequence);
    }
    const panvar::PooledDecision reversed = panvar::decide_pooled_genotype(pool, blocks, 10, 6);
    require(render_decision(forward) == render_decision(reversed),
            label + ": reversing the candidate order must give the identical decision:\n  " +
                render_decision(forward) + "\n  " + render_decision(reversed));
    return forward;
}

void test_pooled_decision_rule() {
    using panvar::TopSetDecision;
    // 1. Two tied candidates differing only in the phase between their two heterozygous blocks.
    {
        const auto d = decide_both_ways({pooled({0, 0, 2}, {1, 1, 2}, -1000.0, -5.0),
                                         pooled({0, 1, 2}, {1, 0, 2}, -1000.0, -5.0),
                                         pooled({0, 0, 2}, {0, 0, 2}, -1100.0, -5.0)},
                                        3, "global phase");
        require(d.decision == TopSetDecision::GenotypeEquivalent && d.top_groups.size() == 2,
                "two tied phase alternatives of one genotype must be genotype-equivalent");
        require(d.phase_support == "UNRESOLVED", "a tied top set must leave phase UNRESOLVED");
        require(d.genotype == std::vector<std::pair<std::size_t, std::size_t>>{{0, 1}, {0, 1}, {2, 2}},
                "the equivalent set must emit its shared unordered genotype");
        require(d.has_other_genotype && d.genotype_separation_certified,
                "the genotype must be separated from the best different genotype");
    }
    // 2. Two tied candidates with blockwise phase switches but identical unordered genotypes.
    {
        const auto d = decide_both_ways({pooled({0, 0, 0, 0}, {1, 1, 1, 1}, -1000.0, -5.0),
                                         pooled({0, 1, 0, 1}, {1, 0, 1, 0}, -1000.0, -5.0),
                                         pooled({0, 0, 0, 0}, {0, 1, 1, 1}, -1100.0, -5.0)},
                                        4, "blockwise switches");
        require(d.decision == TopSetDecision::GenotypeEquivalent && d.phase_support == "UNRESOLVED",
                "tied blockwise phase switches of one genotype must be genotype-equivalent and unphased");
    }
    // 3. A tie differing at one block: refused, and the disagreeing block is named.
    {
        const auto d = decide_both_ways({pooled({0, 0, 0}, {1, 1, 1}, -1000.0, -5.0),
                                         pooled({0, 0, 0}, {1, 1, 0}, -1000.0, -5.0),
                                         pooled({0, 0, 0}, {0, 0, 0}, -1100.0, -5.0)},
                                        3, "one-block disagreement");
        require(d.decision == TopSetDecision::DistinctTie, "a tie that disagrees at any block must refuse");
        require(d.disagreeing_blocks == std::vector<std::size_t>{2}, "the refusal must name the disagreeing block");
    }
    // 4. A unique winner. A same-genotype phase alternative with an IDENTICAL fragment term, lower only
    // through the prior, does not make it a tie: it makes its phase PRIOR_ONLY.
    {
        const auto d = decide_both_ways({pooled({0, 0, 0}, {1, 1, 1}, -1000.0, -5.0),
                                         pooled({0, 1, 1}, {1, 0, 0}, -1000.0, -12.0),
                                         pooled({0, 0, 0}, {1, 1, 0}, -1004.0, -5.0)},
                                        3, "unique winner");
        require(d.decision == TopSetDecision::Unique && d.top_groups.size() == 1,
                "a unique numerical winner must follow the unique path whatever lower same-genotype candidates exist");
        require(d.phase_support == "PRIOR_ONLY", "a phase chosen only by the prior must be PRIOR_ONLY");
        require(d.has_same_genotype_alternative && d.same_genotype_fragment_identical &&
                    std::fabs(d.same_genotype_prior_gap - 7.0) < 1e-12,
                "the best same-genotype alternative must be reported with its identical fragment term and prior gap");
        require(std::fabs(d.genotype_separation_gap - 4.0) < 1e-9 && d.genotype_separation_certified,
                "separation must be measured against the best DIFFERENT genotype, not the phase alternative");
        const auto reads = decide_both_ways({pooled({0, 0, 0}, {1, 1, 1}, -1000.0, -5.0),
                                             pooled({0, 1, 1}, {1, 0, 0}, -1010.0, -5.0),
                                             pooled({0, 0, 0}, {1, 1, 0}, -1004.0, -5.0)},
                                            3, "read-supported unique winner");
        require(reads.decision == TopSetDecision::Unique && reads.phase_support == "READS_WITHIN_CANDIDATES",
                "a phase the fragments choose must be read-supported within the candidates");
    }
    std::cout << "  pooled decision rule: phase-only and blockwise-switch ties unphased, one-block tie refused, "
                 "unique winner unique (PRIOR_ONLY / READS), all order-invariant\n";
}

// ---- phase-switch symmetry: exact equivalence of two tied branches, and the adversarial case -------------
panvar::SymmetryTraceItem traced(std::size_t round, std::vector<std::size_t> a, std::vector<std::size_t> b,
                                 double fragment, double prior) {
    panvar::SymmetryTraceItem item;
    item.round = round;
    item.first = std::move(a);
    item.second = std::move(b);
    item.score.fragment_term = fragment;
    item.score.exposure_term = -100.0;
    item.score.log_likelihood = fragment - 100.0;
    item.score.absolute_term_sum = std::fabs(fragment) + 100.0;
    item.score.fragments_with_mass = 10;
    item.prior = prior;
    return item;
}

// The local phase switch at `block`: the two haplotypes' alleles there exchanged.
panvar::SymmetryTraceItem switched(panvar::SymmetryTraceItem item, std::size_t block) {
    std::swap(item.first[block], item.second[block]);
    return item;
}

void test_phase_switch_symmetry() {
    // Branch A from a0 = (0,0,2 | 1,1,3); branch B from a1 = the local switch at block 1 = (0,1,2 | 1,0,3).
    const std::size_t block = 1;
    std::vector<panvar::SymmetryTraceItem> a = {
        traced(1, {0, 0, 2}, {1, 1, 3}, -1000.0, -5.0),    // incumbent a0
        traced(1, {4, 0, 2}, {1, 1, 3}, -990.0, -6.0),     // improving move at block 0
        traced(1, {0, 0, 5}, {1, 1, 3}, -1004.0, -6.0),
        traced(2, {4, 0, 2}, {1, 1, 3}, -990.0, -6.0),     // incumbent after the move
        traced(2, {4, 0, 2}, {1, 1, 6}, -995.0, -7.0),
    };
    std::vector<panvar::SymmetryTraceItem> b;
    for (const auto& item : a) b.push_back(switched(item, block));
    const panvar::SymmetryReport same = panvar::check_phase_switch_symmetry(a, b, block, 10);
    require(same.equivalent && same.first_break_round == 0 && same.depths.size() == 2,
            "branches that are local switches of each other, bit for bit, must be equivalent");
    require(same.numerically_equivalent, "bit-identical switched branches are numerically equivalent too");

    // Summation order: the same switched branches with one counterpart's exposure term an ulp away (its
    // total within the certified bound) are not bit-exact, but are equivalent under the certify rule's equality.
    std::vector<panvar::SymmetryTraceItem> b_ulp = b;
    b_ulp[3].score.exposure_term = std::nextafter(b_ulp[3].score.exposure_term, 0.0);
    b_ulp[3].score.log_likelihood = std::nextafter(b_ulp[3].score.log_likelihood, 0.0);
    const panvar::SymmetryReport ulp = panvar::check_phase_switch_symmetry(a, b_ulp, block, 10);
    require(!ulp.equivalent && ulp.first_break_round == 2 && ulp.numerically_equivalent,
            "an ulp of summation order must break bit identity but not the certified numerical equality");
    // ... while a difference beyond the bound (1e-6 nats here, far above ~1e-9) breaks both.
    std::vector<panvar::SymmetryTraceItem> b_far = b;
    b_far[3].score.log_likelihood += 1e-6;
    const panvar::SymmetryReport far = panvar::check_phase_switch_symmetry(a, b_far, block, 10);
    require(!far.numerically_equivalent && far.numeric_first_break_round == 2,
            "a difference beyond the certified bound must break the numerical equality");

    // ADVERSARIAL: locally genotype-equivalent at depth 1, but at depth 2 the switched candidate carries a
    // different prior (a switch the Li-Stephens path can take on one haplotype and not the other), so branch B
    // moves elsewhere and ends on a different genotype. Equivalence must break at depth 2, and the pooled rule
    // over the two terminals must refuse rather than pick one.
    std::vector<panvar::SymmetryTraceItem> b2 = b;
    b2[4].prior = std::nextafter(b2[4].prior, 0.0);       // one bit
    b2.push_back(traced(2, {0, 1, 2}, {4, 0, 7}, -990.0, -6.0));   // a move only branch B makes
    const panvar::SymmetryReport broken = panvar::check_phase_switch_symmetry(a, b2, block, 10);
    require(!broken.equivalent && broken.first_break_round == 2,
            "a single-bit score difference at depth 2 must break equivalence at depth 2, not be inferred away");
    require(!broken.numerically_equivalent && broken.numeric_first_break_round == 2,
            "a move only one branch makes breaks the numerical equality too");
    std::vector<panvar::PoolCandidate> terminals = {pooled({4, 0, 2}, {1, 1, 3}, -990.0, -6.0),
                                                    pooled({0, 1, 2}, {4, 0, 7}, -990.0, -6.0)};
    require(panvar::decide_pooled_genotype(terminals, 3, 10, 6).decision == panvar::TopSetDecision::DistinctTie,
            "branches that end on different genotypes, tied, must be refused");

    // A candidate one branch scores and the other does not (their visited sets differ) also breaks it.
    std::vector<panvar::SymmetryTraceItem> b3 = b;
    b3.erase(b3.begin() + 2);
    const panvar::SymmetryReport visited = panvar::check_phase_switch_symmetry(a, b3, block, 10);
    require(!visited.equivalent && visited.first_break_round == 1,
            "a switched neighbourhood missing a candidate must break equivalence at that depth");
    require(!visited.numerically_equivalent && visited.numeric_first_break_round == 1,
            "a visited-set asymmetry breaks the numerical equality too");
    std::cout << "  phase-switch symmetry: exact switched branches equivalent; an ulp breaks bit identity only; a gap "
                 "beyond the bound, one-bit prior break, divergent terminals and visited-set asymmetry all detected\n";
}

// ---- improving-move tie branching: the engine and the calling contract, on a scripted landscape --------
panvar::TieItem tie_item(std::size_t round, std::vector<std::size_t> a, std::vector<std::size_t> b, double fragment,
                         double prior) {
    const panvar::PoolCandidate c = pooled(a, b, fragment, prior);
    panvar::TieItem item;
    item.round = round;
    item.first = c.first;
    item.second = c.second;
    item.first_sequence = c.first_sequence;
    item.second_sequence = c.second_sequence;
    item.score = c.score;
    item.prior = prior;
    return item;
}

struct ScriptedTies {
    std::map<panvar::DiplotypeIdentity, panvar::BranchResult> script;
    std::vector<std::string> ran;              // start identities, in the order branches ran
    std::vector<std::string> ties_seen;
    std::vector<panvar::PoolCandidate> pool;
    std::set<panvar::DiplotypeIdentity> pooled_ids;
    std::string throw_for;                     // a start that throws
    void add(const panvar::TieItem& item) {
        const panvar::DiplotypeIdentity key = panvar::diplotype_identity(item.first, item.second);
        if (!pooled_ids.insert(key).second) return;
        panvar::PoolCandidate c;
        c.first = item.first;
        c.second = item.second;
        c.first_sequence = item.first_sequence;
        c.second_sequence = item.second_sequence;
        c.score = item.score;
        c.sequence_prior = item.prior;
        c.path_prior = item.prior;
        pool.push_back(c);
    }
};

panvar::BranchResult converged_at(std::size_t round, const panvar::TieItem& start, const panvar::TieItem& terminal) {
    panvar::BranchResult r;
    r.outcome = panvar::BranchOutcome::Converged;
    r.items = {start, terminal};
    r.items[0].round = round;
    r.items[1].round = round;
    r.rounds = 2;
    r.moves = 1;
    r.last_round = 2;
    r.stop = "no_improvement";
    r.moved_to = {panvar::diplotype_identity(terminal.first, terminal.second)};
    r.terminal_first = terminal.first;
    r.terminal_second = terminal.second;
    r.terminal_total = terminal.score.log_likelihood + terminal.prior;
    return r;
}

std::string identity_text(const panvar::TieItem& item) {
    std::string text;
    const panvar::DiplotypeIdentity key = panvar::diplotype_identity(item.first, item.second);
    for (const auto* side : {&key.first, &key.second}) {
        for (const std::size_t x : *side) text += std::to_string(x);
        text += "|";
    }
    return text;
}

struct TieRun {
    panvar::TieBranchingReport report;
    panvar::PooledDecision decided;
    panvar::TieBranchVerdict verdict;
    std::vector<std::string> ran;
};

TieRun run_scripted(ScriptedTies scripted, const std::vector<panvar::TieItem>& root_items, std::size_t max_branches,
                    bool reverse) {
    for (const panvar::TieItem& item : root_items) scripted.add(item);
    const panvar::TieSet root = panvar::terminal_tie(root_items, root_items.front().round, 10, 6);
    const auto run_branch = [&](const panvar::BranchStart& start) {
        if (identity_text(start.start) == scripted.throw_for) throw std::runtime_error("scripted branch failure");
        scripted.ran.push_back(identity_text(start.start));
        const panvar::BranchResult result =
            scripted.script.at(panvar::diplotype_identity(start.start.first, start.start.second));
        for (const panvar::TieItem& item : result.items) scripted.add(item);
        return result;
    };
    const auto on_tie = [&](const panvar::TieSet& tie) { scripted.ties_seen.push_back(std::to_string(tie.tied.size())); };
    TieRun out;
    out.report = panvar::run_tie_branching(root, {}, max_branches, 10, 6, run_branch, on_tie, reverse);
    out.ran = scripted.ran;
    out.decided = panvar::decide_pooled_genotype(scripted.pool, 3, 10, 6);
    out.verdict = panvar::tie_branch_verdict(!out.report.incomplete.empty(), out.decided.decision);
    return out;
}

// The branch lifecycle points run_tie_branching gives its caller: after the root tie is recorded and seeded (-1), and
// after each branch is fully processed -- a split only once its own tie is recorded and its children seeded. At each,
// the seeds still waiting are exactly the branches not yet run, so a caller releasing everything else keeps the
// frontier and nothing a finished branch owned.
void test_tie_branch_settle_points() {
    const panvar::TieItem I = tie_item(2, {0, 0, 0}, {1, 1, 1}, -1000.0, -5.0);
    const panvar::TieItem a0 = tie_item(2, {0, 2, 0}, {1, 3, 1}, -990.0, -5.0);
    const panvar::TieItem a1 = tie_item(2, {0, 3, 0}, {1, 2, 1}, -990.0, -5.0);
    const std::vector<panvar::TieItem> root_items = {I, a0, a1};
    const auto key = [](const panvar::TieItem& x) { return panvar::diplotype_identity(x.first, x.second); };
    std::map<panvar::DiplotypeIdentity, panvar::BranchResult> script;
    panvar::BranchResult split;
    split.outcome = panvar::BranchOutcome::Split;
    split.items = {tie_item(1, {0, 2, 0}, {1, 3, 1}, -990.0, -5.0), tie_item(1, {0, 2, 8}, {1, 3, 1}, -985.0, -5.0),
                   tie_item(1, {0, 2, 1}, {1, 3, 8}, -985.0, -5.0)};
    split.rounds = 1;
    split.last_round = 1;
    split.stop = "rescore_incomplete";
    split.terminal_first = {0, 2, 0};
    split.terminal_second = {1, 3, 1};
    script[key(a0)] = split;
    script[key(a1)] = converged_at(1, a1, tie_item(0, {0, 3, 0}, {1, 2, 9}, -995.0, -5.0));
    script[key(split.items[1])] = converged_at(1, split.items[1], tie_item(0, {0, 2, 8}, {1, 3, 6}, -975.0, -5.0));
    script[key(split.items[2])] = converged_at(1, split.items[2], tie_item(0, {0, 2, 6}, {1, 3, 8}, -975.0, -5.0));
    std::set<std::set<std::string>> terminal_sets;
    for (const bool reverse : {false, true}) {
        std::vector<std::string> events;
        std::set<std::string> waiting;   // seeds built and not yet run
        const panvar::TieSet root = panvar::terminal_tie(root_items, 2, 10, 6);
        const auto on_tie = [&](const panvar::TieSet& tie) {
            for (const panvar::TieItem& n : tie.tied) waiting.insert(identity_text(n));
            events.push_back("tie" + std::to_string(tie.id));
        };
        const auto run_branch = [&](const panvar::BranchStart& start) {
            require(waiting.count(identity_text(start.start)) == 1, "a branch must start from a built, waiting seed");
            waiting.erase(identity_text(start.start));
            events.push_back("run" + std::to_string(start.id));
            return script.at(key(start.start));
        };
        std::vector<std::size_t> frontier;
        const auto settled = [&](long branch) {
            events.push_back("settle" + std::to_string(branch));
            frontier.push_back(waiting.size());
        };
        const panvar::TieBranchingReport report =
            panvar::run_tie_branching(root, {}, 16, 10, 6, run_branch, on_tie, reverse, settled);
        require(report.incomplete.empty() && report.branches.size() == 4, "four branches, complete");
        require(events.front() == "tie0" && events[1] == "settle-1", "the root settles after its seeds are built");
        for (std::size_t i = 0; i < events.size(); ++i) {
            if (events[i] != "tie1") continue;
            require(i + 1 < events.size() && events[i + 1].rfind("settle", 0) == 0,
                    "a split branch settles only after its own tie is recorded and its children seeded");
        }
        // Canonical: the split a0 runs first and adds two seeds; reversed: a1 converges first, then a0 splits.
        const std::vector<std::size_t> want = reverse ? std::vector<std::size_t>{2, 1, 2, 1, 0}
                                                      : std::vector<std::size_t>{2, 3, 2, 1, 0};
        require(frontier == want, "the waiting seeds at each settle point are the branches not yet run");
        require(frontier.back() == 0, "nothing waits once every branch has run");
        std::set<std::string> terminals(report.terminal_genotypes.begin(), report.terminal_genotypes.end());
        terminal_sets.insert(terminals);
    }
    require(terminal_sets.size() == 1, "canonical and reversed branch orders reach the same terminals");
    std::cout << "  tie-branch settle points: root, then each branch after processing (a split after its children are "
                 "seeded); the waiting seeds are the frontier; both orders agree\n";
}

void test_tie_branching_engine() {
    // The production tie: incumbent I, two improving phase assignments a0 / a1 of one pair, tied.
    const panvar::TieItem I = tie_item(2, {0, 0, 0}, {1, 1, 1}, -1000.0, -5.0);
    const panvar::TieItem a0 = tie_item(2, {0, 2, 0}, {1, 3, 1}, -990.0, -5.0);
    const panvar::TieItem a1 = tie_item(2, {0, 3, 0}, {1, 2, 1}, -990.0, -5.0);
    const panvar::TieItem other = tie_item(2, {4, 0, 0}, {1, 1, 1}, -1003.0, -5.0);
    const std::vector<panvar::TieItem> root = {I, a0, a1, other};
    const panvar::TieSet tie = panvar::terminal_tie(root, 2, 10, 6);
    require(tie.tied.size() == 2 && identity_text(tie.tied[0]) == identity_text(a0),
            "the certify rule finds both improving assignments, in canonical identity order");
    const auto key = [](const panvar::TieItem& x) { return panvar::diplotype_identity(x.first, x.second); };

    // 1. unique: the canonical-FIRST branch ends lower; the pooled rule takes the other branch's terminal
    {
        ScriptedTies s;
        s.script[key(a0)] = converged_at(1, a0, tie_item(0, {0, 2, 5}, {1, 3, 1}, -980.0, -5.0));
        const panvar::TieItem y = tie_item(0, {0, 3, 0}, {1, 2, 6}, -970.0, -5.0);
        s.script[key(a1)] = converged_at(1, a1, y);
        const TieRun r = run_scripted(s, root, 16, false);
        require(r.report.incomplete.empty() && r.report.branches.size() == 2, "two branches, complete");
        require(r.decided.decision == panvar::TopSetDecision::Unique, "one pooled winner");
        const panvar::PoolCandidate& best = r.decided.candidates[r.decided.representative[r.decided.ranking.candidates.front().candidate_id]];
        require(panvar::diplotype_identity(best.first, best.second) == key(y),
                "the winner is the second branch's terminal: never the first branch because it is first");
        require(r.verdict.status == "COMPLETE" && r.verdict.search_conditional == "tie_branched" &&
                    r.verdict.no_call_reason == "none",
                "unique -> COMPLETE, search_conditional=tie_branched");
    }
    // 2. genotype-equivalent: both branches end on one unordered genotype, phased differently, tied
    {
        ScriptedTies s;
        s.script[key(a0)] = converged_at(1, a0, tie_item(0, {0, 2, 5}, {1, 3, 1}, -980.0, -5.0));
        s.script[key(a1)] = converged_at(1, a1, tie_item(0, {0, 2, 1}, {1, 3, 5}, -980.0, -5.0));
        const TieRun r = run_scripted(s, root, 16, false);
        require(r.decided.decision == panvar::TopSetDecision::GenotypeEquivalent && r.decided.phase_support == "UNRESOLVED",
                "one unordered genotype, phase unresolved");
        require(r.verdict.status == "CONDITIONAL_GENOTYPE_EQUIVALENT", "genotype-equivalent -> CONDITIONAL_GENOTYPE_EQUIVALENT");
    }
    // 3. divergent: tied terminals that disagree at a block
    {
        ScriptedTies s;
        s.script[key(a0)] = converged_at(1, a0, tie_item(0, {0, 2, 5}, {1, 3, 1}, -980.0, -5.0));
        s.script[key(a1)] = converged_at(1, a1, tie_item(0, {0, 3, 7}, {1, 2, 1}, -980.0, -5.0));
        const TieRun r = run_scripted(s, root, 16, false);
        require(r.report.terminal_genotypes.size() == 2 && r.decided.decision == panvar::TopSetDecision::DistinctTie,
                "two terminal genotypes, tied");
        require(r.verdict.no_call_reason == "divergent_tie" && r.verdict.call_class == "no_call", "divergent -> no-call divergent_tie");
    }
    // 4. nested: branch a0 stops at a tie of its own; its two assignments become branches 2 and 3
    ScriptedTies nested;
    {
        panvar::BranchResult split;
        split.outcome = panvar::BranchOutcome::Split;
        split.items = {tie_item(1, {0, 2, 0}, {1, 3, 1}, -990.0, -5.0), tie_item(1, {0, 2, 8}, {1, 3, 1}, -985.0, -5.0),
                       tie_item(1, {0, 2, 1}, {1, 3, 8}, -985.0, -5.0)};
        split.rounds = 1;
        split.last_round = 1;
        split.stop = "rescore_incomplete";
        split.terminal_first = {0, 2, 0};
        split.terminal_second = {1, 3, 1};
        nested.script[key(a0)] = split;
        nested.script[key(a1)] = converged_at(1, a1, tie_item(0, {0, 3, 0}, {1, 2, 9}, -995.0, -5.0));
        const panvar::TieItem b0 = split.items[1];
        const panvar::TieItem b1 = split.items[2];
        nested.script[key(b0)] = converged_at(1, b0, tie_item(0, {0, 2, 8}, {1, 3, 6}, -975.0, -5.0));
        nested.script[key(b1)] = converged_at(1, b1, tie_item(0, {0, 2, 6}, {1, 3, 8}, -975.0, -5.0));
    }
    const TieRun deep = run_scripted(nested, root, 4, false);
    require(deep.report.incomplete.empty() && deep.report.ties.size() == 2 && deep.report.branches.size() == 4,
            "a nested tie adds its branches: two ties, four branches");
    require(deep.report.branches[2].start.history.size() == 1 && deep.report.branches[2].start.parent_tie == 1,
            "a nested branch inherits its parent's start in its history, and names its tie");
    require(deep.ran.size() == 4 && deep.ran[0] == identity_text(a0) && deep.ran[1] == identity_text(a1),
            "breadth-first, canonical order: the production tie's branches run first");
    require(deep.verdict.status == "CONDITIONAL_GENOTYPE_EQUIVALENT", "the nested terminals pool to one genotype");
    // 5. the branch cap: exactly the four it needs succeeds; three refuses, with no decision to make
    const TieRun capped = run_scripted(nested, root, 3, false);
    require(!capped.report.incomplete.empty() && capped.report.incomplete.find("branch cap 3") != std::string::npos,
            "one branch below what the ties need refuses");
    require(capped.verdict.status == "REPAIR_INCOMPLETE" && capped.verdict.no_call_reason == "tie_branching_incomplete",
            "incompleteness -> REPAIR_INCOMPLETE tie_branching_incomplete");
    // 5b. the cap is proven insufficient, not exhausted: once a0 splits the worklist holds four unique starts against
    // a cap of three, so the search stops right there -- no further branch runs and no further tie is seen.
    require(capped.ran.size() == 1 && capped.ran[0] == identity_text(a0) && capped.report.branches.size() == 1 &&
                capped.report.ties.size() == 2,
            "the cap-bound search stops as soon as its unique starts exceed the cap, running nothing after that");
    require(capped.report.incomplete.find("at least 4 branches required; proven after 1 branches") != std::string::npos,
            "the refusal names the proven lower bound and where it was proven");
    {
        // A root tie alone wider than the cap: refused before any branch runs.
        const TieRun root_only = run_scripted(nested, root, 1, false);
        require(root_only.ran.empty() && root_only.report.branches.empty() &&
                    root_only.report.incomplete.find("branch cap 1 binds (at least 2") != std::string::npos,
                "a root tie wider than the cap is refused before any branch runs");
        // Reversed order reaches the same verdict (it may prove it at a different branch).
        const TieRun capped_reversed = run_scripted(nested, root, 3, true);
        require(capped_reversed.verdict.status == capped.verdict.status &&
                    capped_reversed.verdict.no_call_reason == capped.verdict.no_call_reason &&
                    capped_reversed.report.incomplete.find("branch cap 3 binds") != std::string::npos,
                "forward and reversed enumeration reach the same cap-bound verdict");
        // The semantics are those of any other incomplete branching (a throwing branch below).
        ScriptedTies s = nested;
        s.throw_for = identity_text(a1);
        const TieRun thrown = run_scripted(s, root, 16, false);
        require(thrown.verdict.status == capped.verdict.status && thrown.verdict.call_class == capped.verdict.call_class &&
                    thrown.verdict.no_call_reason == capped.verdict.no_call_reason &&
                    thrown.verdict.search_conditional == capped.verdict.search_conditional,
                "a cap-bound stop has exactly the status, call class, reason and conditional of any incomplete branching");
    }
    // 6. an exception inside a branch fails closed
    {
        ScriptedTies s = nested;
        s.throw_for = identity_text(a1);
        const TieRun r = run_scripted(s, root, 16, false);
        require(r.report.incomplete.find("scripted branch failure") != std::string::npos &&
                    r.verdict.no_call_reason == "tie_branching_incomplete",
                "a throwing branch makes the run incomplete");
    }
    // 7. reversed order within every tie: the same branches, the same pooled decision and representative
    {
        const TieRun reversed = run_scripted(nested, root, 4, true);
        require(reversed.ran[0] == identity_text(a1), "the gate hook really reverses the order");
        require(render_decision(reversed.decided) == render_decision(deep.decided) &&
                    reversed.report.branches.size() == deep.report.branches.size(),
                "the pooled result does not depend on the branch order");
    }
    // 8. a split whose last round holds no improving tie is an error, never a silent stop
    {
        ScriptedTies s = nested;
        panvar::BranchResult bad = s.script[key(a0)];
        bad.items = {tie_item(1, {0, 2, 0}, {1, 3, 1}, -990.0, -5.0), tie_item(1, {0, 2, 8}, {1, 3, 1}, -995.0, -5.0)};
        s.script[key(a0)] = bad;
        const TieRun r = run_scripted(s, root, 16, false);
        require(r.report.incomplete.find("without an improving tie") != std::string::npos, "a split with no tie refuses");
    }
    std::cout << "  tie branching: unique (not the first branch), genotype-equivalent, divergent, nested BFS, cap exact "
                 "and minus one, exceptions, reversed order and tieless splits all behave as contracted\n";
}

// ---- representation-conflict guard: the rule, at its bound and under every symmetry --------------------
panvar::DiplotypeScore component_score(double fragment, double exposure) {
    panvar::DiplotypeScore s;
    s.fragment_term = fragment;
    s.exposure_term = exposure;
    s.log_likelihood = fragment + exposure;
    s.absolute_term_sum = std::fabs(fragment) + std::fabs(exposure);
    return s;
}

void test_representation_guard() {
    const double B = panvar::kRepresentationBoundNats;
    require(std::fabs(B - std::log(1000.0)) < 1e-15, "the bound is ln(1000)");
    // The sv:7 shape: fragments +253 for the call, exposure -164 against it, marker 186 against it.
    const panvar::DiplotypeScore called = component_score(-30172.49, -2323.65);
    const panvar::DiplotypeScore comp = component_score(-30425.36, -2159.96);
    const panvar::RepresentationEvidence sv7 =
        panvar::representation_conflict(186.3, called, -9.2, {{comp, -45.2}, {comp, -45.6}});
    require(sv7.conflict && std::fabs(sv7.statistic - 163.69) < 1e-6, "two depth readings against the call: flagged");
    // The exact-mosaic shape: the marker disagrees by 63 nats, the exposure term does not -- never flagged.
    require(!panvar::representation_conflict(63.3, component_score(-1000.0, -200.0), 0.0,
                                             {{component_score(-1032.4, -200.6), -18.2}}).conflict,
            "a marker disagreement without an exposure-scale difference is not a conflict");
    // A depth channel that agrees with the call can never flag it, whatever the components.
    require(!panvar::representation_conflict(0.0, called, -9.2, {{comp, -45.2}}).conflict, "marker margin 0: PASS");
    // At the bound exactly: flagged; one ulp below on any channel: PASS.
    const panvar::DiplotypeScore at = component_score(-100.0 + B, -50.0 - B);
    const panvar::DiplotypeScore rival = component_score(-100.0, -50.0);
    require(panvar::representation_conflict(B, at, 0.0, {{rival, 0.0}}).conflict, "exactly at the bound: flagged");
    require(!panvar::representation_conflict(std::nextafter(B, 0.0), at, 0.0, {{rival, 0.0}}).conflict,
            "one ulp below the bound on the marker: PASS");
    // Assignment order and homologue swap: the competitor is the better assignment by a total order, so the
    // order they are offered in, or which homologue carries which allele, cannot change the evidence.
    const panvar::DiplotypeScore a0 = component_score(-30425.36, -2159.96);
    const panvar::DiplotypeScore a1 = component_score(-30425.36, -2159.90);   // same total after its prior
    const panvar::RepresentationEvidence forward = panvar::representation_conflict(186.3, called, -9.2, {{a0, -45.20}, {a1, -45.26}});
    const panvar::RepresentationEvidence backward = panvar::representation_conflict(186.3, called, -9.2, {{a1, -45.26}, {a0, -45.20}});
    require(forward.conflict == backward.conflict && forward.statistic == backward.statistic &&
                forward.exposure_diff == backward.exposure_diff && forward.likelihood_margin == backward.likelihood_margin,
            "competitor assignment order (and so homologue swap) cannot change the evidence");
    // The marker-best pair: emission descending, the first in index order on a tie.
    panvar::BlockProposalEmissions emission;
    emission.allele_count = 3;
    emission.log_score.assign(9, -10.0);
    emission.log_score[0 * 3 + 2] = -1.0;
    emission.log_score[2 * 3 + 0] = -1.0;
    emission.log_score[1 * 3 + 2] = -1.0;
    emission.log_score[2 * 3 + 1] = -1.0;
    require(panvar::marker_best_pair(emission) == std::make_pair(std::size_t{0}, std::size_t{2}),
            "a marker tie goes to the first pair in index order");
    std::cout << "  representation guard: flags the two-depth-readings conflict, not a marker-only disagreement; "
                 "bound exact and one ulp below; order and homologue swap invariant; marker ties canonical\n";
}

// ---- panel-anchored adaptive completion, end to end ------------------------------------------------
panvar::GenotypeCommandOptions repair_options(const std::string& prefix, const std::string& r1,
                                              const std::string& r2, const std::string& gfa) {
    panvar::GenotypeCommandOptions options = tiny_options(prefix, r1, r2, gfa);
    options.repair_rounds = 11;
    options.repair_local_pairs = 32;
    options.repair_origin_cache_bytes = 1 << 24;
    options.threads = 2;
    return options;
}

std::string slurp(const std::string& file) {
    std::ifstream input(file);
    std::stringstream buffer;
    buffer << input.rdbuf();
    return buffer.str();
}

// The scientific columns of a rounds table: everything up to the invariance flags, no timing.
std::string rounds_science(const std::string& file) {
    std::string out;
    for (const auto& row : read_rows(file)) {
        for (std::size_t c = 0; c < 28 && c < row.size(); ++c) out += row[c] + '\t';
        out += '\n';
    }
    return out;
}

bool no_final_outputs(const std::string& prefix) {
    for (const char* suffix : {".genotypes.tsv", ".block_genotypes.tsv", ".mosaic_candidates.tsv", ".mosaic_scores.tsv",
                               ".mosaic_haplotypes.fa", ".equivalence_set.tsv", ".repair_rounds.tsv"}) {
        if (exists(path(prefix + suffix))) return false;
    }
    return true;
}

void require_truth_genotype(const panvar::GenotypeCommandOptions& options, const std::string& prefix) {
    const auto reported = read_status(path(prefix + ".mosaic_status.tsv"));
    const TruthFrame frame = rebuild_truth_frame(options, reported.at("model_fingerprint"), kTruth1Nodes, kTruth2Nodes);
    const auto rows = read_rows(path(prefix + ".block_genotypes.tsv"));
    require(rows.size() == frame.model->blocks().size(), "the block table must carry every block");
    for (std::size_t b = 0; b < rows.size(); ++b) {
        const auto& alleles = frame.model->blocks()[b].alleles;
        const std::string& c1 = alleles[std::stoul(rows[b][3])].sequence;
        const std::string& c2 = alleles[std::stoul(rows[b][4])].sequence;
        const std::string& a = *frame.truth1.blocks[b];
        const std::string& d = *frame.truth2.blocks[b];
        require((c1 == a && c2 == d) || (c1 == d && c2 == a),
                prefix + ": block " + rows[b][0] + " must be called as the truth's unordered allele pair");
    }
}

// ---- the representation guard's outside-pool path, end to end on the tiny locus ------------------------
// Gate hooks score every marker competitor in the live context (never from the pool) and accept any margin,
// so the path Sol asked to see is exercised on a cheap locus. The hooks also flag every scored block, which is
// what lets the cap-minus-one case keep a LowQual and an unevaluated block side by side.
std::vector<std::string> guard_prefilter(const std::string& prefix) {
    std::vector<std::string> out;
    for (const char* table : {".genotypes.tsv", ".block_genotypes.tsv"}) {
        for (const auto& row : read_rows(path(prefix + table))) {
            std::string line;
            // The guard's eight filter columns and the two combined-quality columns built on them come last.
            for (std::size_t c = 0; c + 10 < row.size(); ++c) line += row[c] + '\t';
            out.push_back(line);
        }
    }
    for (const char* table : {".mosaic_candidates.tsv", ".mosaic_scores.tsv", ".equivalence_set.tsv", ".repair_rounds.tsv"}) {
        for (const auto& row : read_rows(path(prefix + table))) {
            std::string line;
            for (std::size_t c = 0; c < row.size() && (std::string(table) != ".repair_rounds.tsv" || c < 28); ++c) line += row[c] + '\t';
            out.push_back(line);
        }
    }
    std::ifstream fasta(path(prefix + ".mosaic_haplotypes.fa"));
    for (std::string line; std::getline(fasta, line);) {
        if (!line.empty() && line[0] != '>') out.push_back(line);
    }
    // The status apart from timing, process gauges, the guard's own keys and the cumulative repair work,
    // which legitimately includes the guard's scoring when it is on.
    for (const auto& [key, value] : read_status(path(prefix + ".mosaic_status.tsv"))) {
        static const std::vector<std::string> skip = {"seconds_", "resources_", "threads_", "representation_",
                                                      "locus_filter", "lowqual_blocks", "not_evaluated_blocks",
                                                      "quality_",
                                                      "repair_cumulative_", "repair_margin_", "repair_origin_cache_bytes_end",
                                                      "settings_overridden", "phys_", "rss_", "emission_cache_",
                                                      "scorer_",
                                                      // the declared cap itself, which the cap cases set
                                                      "repair_max_physical_origins"};
        bool skipped = false;
        for (const std::string& prefix_key : skip) skipped = skipped || key.rfind(prefix_key, 0) == 0;
        if (!skipped) out.push_back(key + "=" + value);
    }
    return out;
}

std::vector<std::string> guard_filters(const std::string& prefix) {
    std::vector<std::string> out;
    for (const auto& row : read_rows(path(prefix + ".block_genotypes.tsv"))) {
        std::string line = row[0];
        for (std::size_t c = row.size() - 10; c < row.size(); ++c) line += '\t' + row[c];
        out.push_back(line);
    }
    return out;
}

void test_representation_guard_outside_pool() {
    const std::string r1 = path("phaseable_r1.fq");
    const std::string r2 = path("phaseable_r2.fq");
    const auto options_for = [&](const std::string& name, bool guard, bool force, std::size_t threads) {
        panvar::GenotypeCommandOptions o = repair_options(name, r1, r2, "tiny_phaseable.gfa");
        o.max_mosaic_candidates = 0;
        o.threads = threads;
        o.representation_guard = guard;
        o.test_representation_guard_force_context = force;
        if (guard) {
            o.test_representation_guard_bound_nats = -std::numeric_limits<double>::max();
            o.test_representation_guard_runner_up = true;   // the marker agrees here: take its runner-up
        }
        return o;
    };
    const auto status_of = [&](const std::string& name) { return read_status(path(name + ".mosaic_status.tsv")); };
    require(panvar::run_genotype(options_for("guard_off", false, false, 2)) == panvar::GenotypeRunStatus::Complete,
            "fixture: the guard-off run calls");
    require(panvar::run_genotype(options_for("guard_pool", true, false, 2)) == panvar::GenotypeRunStatus::Complete &&
                panvar::run_genotype(options_for("guard_ctx", true, true, 2)) == panvar::GenotypeRunStatus::Complete,
            "the guard never changes the call status");
    const auto off = status_of("guard_off");
    const auto pooled_run = status_of("guard_pool");
    const auto ctx = status_of("guard_ctx");
    const std::size_t scored = std::stoul(ctx.at("representation_competitors_scored"));
    require(scored > 0, "the fixture must reach the outside-pool path (competitors scored in the live context)");
    require(guard_prefilter("guard_ctx") == guard_prefilter("guard_off") &&
                guard_prefilter("guard_pool") == guard_prefilter("guard_off"),
            "call, ranking, tables and FASTA are identical with the guard on (pool or live context) and off");
    require(guard_filters("guard_ctx") == guard_filters("guard_pool"),
            "a competitor scored outside the pool carries exactly the evidence the pool holds for it");
    // The guard's own work reconciles exactly with the refreshed cumulative counters.
    const std::uint64_t guard_work = std::stoull(ctx.at("representation_guard_physical_origins"));
    const std::uint64_t total_on = std::stoull(ctx.at("repair_cumulative_physical_origins"));
    const std::uint64_t total_off = std::stoull(off.at("repair_cumulative_physical_origins"));
    require(guard_work > 0 && total_on == total_off + guard_work,
            "cumulative physical origins = production work + the guard's reported work, exactly");
    require(std::stoull(pooled_run.at("representation_guard_physical_origins")) == 0 &&
                pooled_run.at("repair_cumulative_physical_origins") == off.at("repair_cumulative_physical_origins"),
            "competitors held in the pool cost no scoring");
    require(ctx.at("representation_guard_complete") == "true" && ctx.at("not_evaluated_blocks") == "none",
            "every block completed the rule");
    // Exact-cap admission, then one below: the guard's last scoring is refused -> NOT_EVALUATED, never PASS,
    // and the genotype untouched; the flags the guard did reach are kept beside it.
    panvar::GenotypeCommandOptions exact = options_for("guard_cap_exact", true, true, 2);
    exact.repair_max_physical_origins = total_on;
    require(panvar::run_genotype(exact) == panvar::GenotypeRunStatus::Complete &&
                guard_filters("guard_cap_exact") == guard_filters("guard_ctx") &&
                status_of("guard_cap_exact").at("representation_guard_complete") == "true",
            "at exactly the cumulative cap the guard completes with the same evidence");
    panvar::GenotypeCommandOptions short_by_one = options_for("guard_cap_minus1", true, true, 2);
    short_by_one.repair_max_physical_origins = total_on - 1;
    require(panvar::run_genotype(short_by_one) == panvar::GenotypeRunStatus::Complete,
            "a guard refusal never changes the call status");
    const auto minus = status_of("guard_cap_minus1");
    require(minus.at("representation_guard_complete") == "false" && minus.at("not_evaluated_blocks") != "none",
            "one origin short: a block is NOT_EVALUATED and the guard is marked incomplete");
    require(minus.at("locus_filter") != "PASS", "an incomplete guard is never PASS");
    bool unevaluated_row = false;
    for (const std::string& line : guard_filters("guard_cap_minus1")) {
        require(line.find("\tPASS\t") == std::string::npos || line.find("NOT_EVALUATED") == std::string::npos,
                "a row is never both");
        unevaluated_row = unevaluated_row || line.find("\tNOT_EVALUATED\tcompetitor_not_scored") != std::string::npos;
    }
    require(unevaluated_row, "the refused block says why");
    require(guard_prefilter("guard_cap_minus1") == guard_prefilter("guard_off"),
            "a guard scoring refusal leaves the genotype, ranking, tables and FASTA exactly as they were");
    if (minus.at("lowqual_blocks") != "none") {
        require(minus.at("locus_filter") == "LOWQUAL_REPRESENTATION_CONFLICT",
                "a flag and an unevaluated block coexist: the locus is LowQual and the guard is incomplete");
    }
    // Thread count cannot change the filter.
    for (const std::size_t threads : {std::size_t{1}, std::size_t{4}, std::size_t{8}}) {
        const std::string name = "guard_ctx_t" + std::to_string(threads);
        require(panvar::run_genotype(options_for(name, true, true, threads)) == panvar::GenotypeRunStatus::Complete &&
                    guard_filters(name) == guard_filters("guard_ctx"),
                "the filter is identical at " + std::to_string(threads) + " threads");
    }
    std::cout << "  representation guard outside the pool: " << scored << " competitors scored in the live context, "
              << "call and ranking unchanged, evidence equal to the pool's, " << guard_work
              << " guard origins reconciled, cap exact / minus one (NOT_EVALUATED, never PASS), 1/4/8 threads\n";
}

void test_adaptive_completion() {
    const std::string r1 = path("phaseable_r1.fq");
    const std::string r2 = path("phaseable_r2.fq");

    // The panel alone as the union (K = 0): the truth is off-panel, so only the repair can reach it.
    panvar::GenotypeCommandOptions anchored = repair_options("adaptive_k0", r1, r2, "tiny_phaseable.gfa");
    anchored.max_mosaic_candidates = 0;
    const panvar::GenotypeRunStatus status = panvar::run_genotype(anchored);
    const auto reported = read_status(path("adaptive_k0.mosaic_status.tsv"));
    require(status == panvar::GenotypeRunStatus::Complete,
            std::string("the panel-anchored completion must call the phaseable locus, got ") +
                panvar::genotype_run_status_name(status) + ": " + (reported.count("note") ? reported.at("note") : ""));
    require(reported.at("algorithm") == "panel_anchored_adaptive_completion" &&
                reported.at("repair_converged") == "true" && std::stoul(reported.at("repair_moves")) > 0,
            "the repair must run, move off the panel anchor and converge");
    require(reported.at("final_decision") == "UNIQUE" && reported.at("phase_support") == "READS_WITHIN_CANDIDATES",
            "the phaseable truth must be a unique, read-phased winner of the pool");
    require(reported.at("final_candidate_scope") == "CONDITIONAL" && reported.at("global_certified") == "false" &&
                reported.at("search_globally_exhaustive") == "false",
            "a local search from one anchor must never claim a complete candidate set");
    require(reported.at("pool_score_mismatches") == "0", "the union and the repair must score shared candidates identically");
    require(reported.at("haplotype_representation") == "read_supported_within_candidates",
            "a read-phased unique call must say what its FASTA is");
    require(reported.count("sha256_gfa") && reported.at("sha256_gfa").size() == 64 &&
                reported.at("sha256_reads1").size() == 64,
            "the inputs must be recorded by content");
    require(reported.at("conditional_gq_calibration") == "uncalibrated", "GQ must be labelled uncalibrated");
    require(reported.at("call_class") == "strict_call" && reported.at("no_call_reason") == "none",
            "a unique call must be classed a strict call");
    require(reported.at("union_top_hap1") != "NA" && reported.at("union_top_hap2") != "NA",
            "the union's top candidate must be recorded as the diagnostic fallback");
    require_truth_genotype(anchored, "adaptive_k0");
    require(exists(path("adaptive_k0.repair_rounds.tsv")) && exists(path("adaptive_k0.equivalence_set.tsv")) &&
                exists(path("adaptive_k0.mosaic_haplotypes.fa")),
            "a completed repair must leave its trajectory, its equivalence set and its haplotypes");
    require(!exists(path("adaptive_k0.repair_rounds.partial.tsv")) && !exists(path("adaptive_k0.sorted.gfa")),
            "a completed run must leave no partial rounds and no staged graph");
    require(slurp(path("adaptive_k0.mosaic_haplotypes.fa")).find("phase_support=READS_WITHIN_CANDIDATES") !=
                std::string::npos,
            "the FASTA header must carry the phase support");

    // The frozen union (K = 16) and the thread count: identical trajectory and identical outputs.
    std::string reference_rounds;
    std::string reference_tables;
    for (const std::size_t threads : {std::size_t{1}, std::size_t{4}}) {
        const std::string name = "adaptive_t" + std::to_string(threads);
        panvar::GenotypeCommandOptions run = repair_options(name, r1, r2, "tiny_phaseable.gfa");
        run.threads = threads;
        require(panvar::run_genotype(run) == panvar::GenotypeRunStatus::Complete, "the K = 16 completion must call");
        require_truth_genotype(run, name);
        const std::string rounds = rounds_science(path(name + ".repair_rounds.tsv"));
        const std::string tables = slurp(path(name + ".genotypes.tsv")) + slurp(path(name + ".block_genotypes.tsv")) +
                                   slurp(path(name + ".mosaic_scores.tsv")) + slurp(path(name + ".equivalence_set.tsv")) +
                                   slurp(path(name + ".mosaic_candidates.tsv"));
        if (reference_rounds.empty()) {
            reference_rounds = rounds;
            reference_tables = tables;
        } else {
            require(rounds == reference_rounds, "the repair trajectory must not depend on the thread count");
            require(tables == reference_tables, "the pooled decision and its tables must not depend on the thread count");
        }
    }

    // The unphaseable locus: no fragment links neighbouring variants.
    {
        panvar::GenotypeCommandOptions run =
            repair_options("adaptive_unphaseable", path("unphaseable_r1.fq"), path("unphaseable_r2.fq"), "tiny.gfa");
        const panvar::GenotypeRunStatus result = panvar::run_genotype(run);
        const auto unphased = read_status(path("adaptive_unphaseable.mosaic_status.tsv"));
        require(result == panvar::GenotypeRunStatus::ConditionalGenotypeEquivalent,
                std::string("tied phase histories of one genotype must be CONDITIONAL_GENOTYPE_EQUIVALENT, got ") +
                    panvar::genotype_run_status_name(result) + ": " + (unphased.count("note") ? unphased.at("note") : ""));
        require(std::stoul(unphased.at("top_set_size")) > 1 && unphased.at("phase_support") == "UNRESOLVED",
                "the equivalent set must hold several candidates with phase UNRESOLVED");
        require(!exists(path("adaptive_unphaseable.mosaic_haplotypes.fa")),
                "a genotype-equivalent call must emit no haplotype FASTA");
        require(read_rows(path("adaptive_unphaseable.equivalence_set.tsv")).size() >= 2,
                "the complete equivalence set must be reported");
        for (const auto& row : read_rows(path("adaptive_unphaseable.block_genotypes.tsv"))) {
            require(row[3] != "." && row[10] == "UNRESOLVED", "every block is genotyped, every block unphased");
        }
        require(unphased.at("haplotype_representation") == "none", "no haplotype representation is emitted");
        require(unphased.at("call_class") == "genotype_equivalent_call", "the call must be classed genotype-equivalent");
    }

    // A round cap: the search is stopped before it converges, so no genotype.
    {
        panvar::GenotypeCommandOptions capped = repair_options("adaptive_cap", r1, r2, "tiny_phaseable.gfa");
        capped.max_mosaic_candidates = 0;
        capped.repair_rounds = 1;
        require(panvar::run_genotype(capped) == panvar::GenotypeRunStatus::RepairIncomplete,
                "a search stopped at its round cap must not call");
        const auto s = read_status(path("adaptive_cap.mosaic_status.tsv"));
        require(s.at("repair_stop_reason") == "round_cap" && s.at("repair_converged") == "false",
                "the round cap must be named as the stop");
        require(s.at("call_class") == "no_call" && s.at("no_call_reason") == "round_cap",
                "a round-capped search is a no-call for the round cap");
        require(!exists(path("adaptive_cap.genotypes.tsv")) && !exists(path("adaptive_cap.mosaic_haplotypes.fa")) &&
                    !exists(path("adaptive_cap.mosaic_scores.tsv")),
                "a round-capped search must emit no genotype and no ranking");
        require(exists(path("adaptive_cap.repair_rounds.tsv")), "the trajectory it did run is kept");
    }
    // A cumulative cap: refused in flight, no genotype, no ranking.
    {
        panvar::GenotypeCommandOptions refused = repair_options("adaptive_refused", r1, r2, "tiny_phaseable.gfa");
        refused.repair_max_physical_origins = 1;
        require(panvar::run_genotype(refused) == panvar::GenotypeRunStatus::RescoreIncomplete,
                "a repair past its cumulative cap must refuse");
        require(no_final_outputs("adaptive_refused") && !exists(path("adaptive_refused.sorted.gfa")),
                "a refused repair must leave no final output and no staged graph");
    }
    std::cout << "  adaptive completion: off-panel truth reached from the panel anchor ("
              << reported.at("repair_moves") << " moves), thread-invariant; unphaseable locus genotype-equivalent "
                 "with no FASTA; round cap and cumulative cap emit no genotype\n";
}

// The watchdog stops the whole process, so it is exercised in a child.
int run_in_child(const panvar::GenotypeCommandOptions& options) {
    std::cout.flush();
    const pid_t child = fork();
    require(child >= 0, "fork failed");
    if (child == 0) {
        try {
            panvar::run_genotype(options);
        } catch (...) {
            std::_Exit(3);
        }
        std::_Exit(0);
    }
    int wait_status = 0;
    waitpid(child, &wait_status, 0);
    return WIFEXITED(wait_status) ? WEXITSTATUS(wait_status) : -1;
}

void test_watchdog_and_interruption() {
    const std::string r1 = path("phaseable_r1.fq");
    const std::string r2 = path("phaseable_r2.fq");
    {
        panvar::GenotypeCommandOptions interrupted = repair_options("adaptive_interrupted", r1, r2, "tiny_phaseable.gfa");
        interrupted.max_mosaic_candidates = 0;
        interrupted.test_interrupt_after_repair_rounds = 1;
        require(run_in_child(interrupted) == 124, "an interrupted run must stop through the time-limit path");
        const auto s = read_status(path("adaptive_interrupted.mosaic_status.tsv"));
        require(s.at("status") == "TIME_LIMIT" && s.at("repair_rounds_committed") == "1",
                "the status must say TIME_LIMIT and how many rounds were committed");
        require(read_rows(path("adaptive_interrupted.repair_rounds.partial.tsv")).size() == 1,
                "the committed round must be preserved as a partial row");
        require(no_final_outputs("adaptive_interrupted") && !exists(path("adaptive_interrupted.sorted.gfa")),
                "an interrupted run must leave no final output and no staged graph");
        // The preserved row is the complete run's first round.
        const auto partial = read_rows(path("adaptive_interrupted.repair_rounds.partial.tsv"));
        const auto complete = read_rows(path("adaptive_k0.repair_rounds.tsv"));
        for (std::size_t c = 0; c < 28; ++c) {
            require(partial[0][c] == complete[0][c], "a preserved partial round must equal the uninterrupted run's");
        }
    }
    {
        panvar::GenotypeCommandOptions breached = repair_options("adaptive_footprint", r1, r2, "tiny_phaseable.gfa");
        breached.footprint_limit_bytes = 1;
        require(run_in_child(breached) == 125, "a footprint breach must stop the process as a resource refusal");
        const auto s = read_status(path("adaptive_footprint.mosaic_status.tsv"));
        require(s.at("status") == "RESOURCE_REFUSAL" && s.at("no_call_reason") == "footprint_limit",
                "the status must say RESOURCE_REFUSAL, for the footprint limit");
        require(no_final_outputs("adaptive_footprint") && !exists(path("adaptive_footprint.sorted.gfa")),
                "a refused run must leave no final output and no staged graph");
    }
    std::cout << "  watchdog: interruption keeps its committed round and nothing final; a footprint breach refuses\n";
}

void test_profile(const std::string& source_dir) {
    {
        std::ofstream abc(path("abc.txt"));
        abc << "abc";
    }
    require(panvar::sha256_file(path("abc.txt")) ==
                "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad",
            "SHA-256 must match the FIPS 180-4 test vector");
    panvar::GenotypeCommandOptions options;
    panvar::apply_genotype_profile(source_dir + "/profiles/genotype_c4.tsv", options);
    require(options.profile_name == "c4-panel-anchored-v1" && options.proposal_state_limit == 256 &&
                options.max_mosaic_candidates == 128 && options.phase_audit_candidates == 0 && options.threads == 8 &&
                options.repair_rounds == 11 && options.repair_local_pairs == 32 &&
                options.repair_origin_cache_bytes == 402653184 && options.footprint_limit_bytes == 2684354560ULL &&
                options.emission_cache_bytes == 536870912 && options.max_physical_origins == 700000000 &&
                options.max_state_classes == 600000000 && options.repair_max_neighbours == 670,
            "the committed C4 profile must set the frozen production settings");
    require(options.profile_sha256.size() == 64, "the profile must be recorded by content");
    const auto refuses = [&](const std::string& body, const std::string& label) {
        {
            std::ofstream bad(path("bad_profile.tsv"));
            bad << body;
        }
        panvar::GenotypeCommandOptions ignored;
        bool threw = false;
        try {
            panvar::apply_genotype_profile(path("bad_profile.tsv"), ignored);
        } catch (const std::exception&) {
            threw = true;
        }
        require(threw, "a profile with " + label + " must be refused");
    };
    refuses("profile\tx\nno_such_setting\t1\n", "an unknown setting");
    refuses("profile\tx\nrepair_rounds\t11\nrepair_rounds\t12\n", "a repeated setting");
    refuses("repair_rounds\t11\n", "no name");
    refuses("profile\tx\nrepair_rounds\t11x\n", "an unparsable value");
    std::cout << "  profile: committed C4 profile applied and hashed; malformed profiles refused\n";
}


// ---- I2-A (EXPERIMENTAL): the indexed union backend at the caller -------------------------------------------------
// Default off and invisible; when selected, every scientific output equals the reference backend's in the same binary
// (byte for byte, files and status) at 1 and 4 threads, through the adaptive repair; only backend provenance, executed
// work, cache, memory and time keys may differ. A refusal on the indexed route ends RESCORE_INCOMPLETE with no ranking.
void with_union_caches(panvar::GenotypeCommandOptions& o) {
    o.union_placement_cache_bytes = std::size_t{64} << 20;
    o.max_union_cached_placements_consumed = std::uint64_t{1} << 40;
    o.union_origin_cache_bytes = std::size_t{64} << 20;
    o.max_union_cached_origins_consumed = std::uint64_t{1} << 40;
    o.max_union_cached_states_consumed = std::uint64_t{1} << 40;
    o.union_prefilter_cache_bytes = std::size_t{16} << 20;
    o.max_union_cached_prefilters_consumed = std::uint64_t{1} << 40;
}

bool indexed_allowed_key(const std::string& key) {
    // perf_compare.py's DIAGNOSTIC class (time, memory, cache, executed work, profile identity), plus the indexed
    // backend's own provenance block and the executed accepted-placement count
    static const std::regex diagnostic(
        "(^seconds|_seconds$|seconds_|wall|rss|footprint|phys_|_gib|timing|cache|cached|scratch|_hits$|_misses$|"
        "^profile$|^profile_name$|^profile_sha256$|^physical_origins$|^statistical_state_classes$|^merged_state_classes$|"
        "^piece_hits$|^reused_|^scorer_|^threads_peak|^resources_|^repair_cumulative_|^repair_margin_|_peak_bytes$|"
        "^max_union_|^union_|^emission_|^start_sort|^repair_origin_cache|^tie_seed_|^tie_branch_physical|"
        "^representation_guard_physical_origins$|^marker_lookup|^seeding_|"
        "^candidate_bases_scanned$|^dp_evaluations$|^full_read_verifications$|^prefilter_rejections$|^proposed_starts$|"
        "^unique_starts$|^valid_fr_states$|^threads_start_sort_|"
        "^indexed_union_|^threads_indexed_union_|^union_scoring_backend|^accepted_mate_placements$|"
        // I2-B: the adaptive window-discovery block (provenance, work, memory, time)
        "^adaptive_window_|^adaptive_fine_|^worker_seconds_adaptive_|"
        // I2-D: the fused-seed provenance key and its full-placement work counter
        "^adaptive_fused_seed$|^adaptive_full_placements$|"
        // I2-E: the marker-proposal reuse block (work and storage)
        "^marker_proposal_|"
        // I2-F: the prepared-combination block (provenance and work)
        "^combination_|"
        // I2-H: the recursion-backend provenance keys
        "^gap_recursion_backend|"
        // I2-K: the emission-handoff provenance and work block
        "^adaptive_emission_handoff|"
        // I3-A: the prepared-locus provenance
        "^prepared_locus)");
    return std::regex_search(key, diagnostic);
}

void require_same_science(const std::string& a, const std::string& b, const std::string& label) {
    for (const char* suffix : {".genotypes.tsv", ".block_genotypes.tsv", ".mosaic_candidates.tsv", ".mosaic_scores.tsv",
                               ".mosaic_haplotypes.fa", ".equivalence_set.tsv", ".locus_genotype.tsv",
                               ".library_sensitivity.tsv"}) {
        require(exists(path(a + suffix)) == exists(path(b + suffix)), label + ": " + suffix + " present in one run only");
        if (exists(path(a + suffix))) {
            // the run's own prefix (FASTA record names, haplotype ids) differs by name only
            const auto placeheld = [](std::string text, const std::string& prefix) {
                for (std::size_t at = text.find(prefix); at != std::string::npos; at = text.find(prefix, at + 5)) {
                    text.replace(at, prefix.size(), "<RUN>");
                }
                return text;
            };
            require(placeheld(slurp(path(a + suffix)), path(a)) == placeheld(slurp(path(b + suffix)), path(b)),
                    label + ": " + suffix + " differs");
        }
    }
    const auto sa = read_status(path(a + ".mosaic_status.tsv")), sb = read_status(path(b + ".mosaic_status.tsv"));
    std::set<std::string> keys;
    for (const auto& kv : sa) keys.insert(kv.first);
    for (const auto& kv : sb) keys.insert(kv.first);
    for (const std::string& key : keys) {
        if (indexed_allowed_key(key)) continue;
        const auto x = sa.find(key), y = sb.find(key);
        require(x != sa.end() && y != sb.end(), label + ": status key '" + key + "' present in one run only");
        // the run's own output prefix differs by name only
        std::string vx = x->second, vy = y->second;
        for (std::string* v : {&vx, &vy}) {
            for (const std::string& name : {a, b}) {
                const std::size_t at = v->find(path(name));
                if (at != std::string::npos) v->replace(at, path(name).size(), "<RUN>");
            }
        }
        require(vx == vy, label + ": status key '" + key + "' differs: " + x->second + " vs " + y->second);
    }
}

void test_indexed_union_backend() {
    const std::string r1 = path("r1.fq"), r2 = path("r2.fq");
    // default off: no backend key in the status or the settings table
    panvar::GenotypeCommandOptions plain = tiny_options("ix_default", r1, r2);
    with_union_caches(plain);
    require(panvar::run_genotype(plain) == panvar::GenotypeRunStatus::Complete, "fixture: the default run must call");
    const auto plain_status = read_status(path("ix_default.mosaic_status.tsv"));
    for (const auto& kv : plain_status) {
        require(kv.first.rfind("indexed_union_", 0) != 0 && kv.first.rfind("union_scoring_backend", 0) != 0 &&
                    kv.first.rfind("threads_indexed_union_", 0) != 0,
                "a default run must not report the experimental backend (" + kv.first + ")");
    }
    require(slurp(path("ix_default.settings.tsv")).find("union_scoring_backend") == std::string::npos &&
                slurp(path("ix_default.settings.tsv")).find("indexed_union_store_bytes") == std::string::npos,
            "a default run's settings table must not list the experimental keys");

    // selected: same science at 1 and 4 threads, without and with the adaptive repair
    const auto indexed = [&](panvar::GenotypeCommandOptions o) {
        o.union_scoring_backend = "indexed_fine_experimental";
        o.indexed_union_store_bytes = std::size_t{256} << 20;
        o.setting_sources["union_scoring_backend"] = "profile test";
        o.setting_sources["indexed_union_store_bytes"] = "profile test";
        return o;
    };
    for (const std::size_t threads : {std::size_t{1}, std::size_t{4}}) {
        const std::string name = "ix_indexed_t" + std::to_string(threads);
        panvar::GenotypeCommandOptions o = tiny_options(name, r1, r2);
        with_union_caches(o);
        o.threads = threads;
        panvar::GenotypeCommandOptions ref = tiny_options("ix_ref_t" + std::to_string(threads), r1, r2);
        with_union_caches(ref);
        ref.threads = threads;
        require(panvar::run_genotype(ref) == panvar::GenotypeRunStatus::Complete, "fixture: the reference must call");
        require(panvar::run_genotype(indexed(o)) == panvar::GenotypeRunStatus::Complete, "the indexed backend must call");
        require_same_science("ix_ref_t" + std::to_string(threads), name,
                             "indexed backend, " + std::to_string(threads) + " thread(s)");
        const auto st = read_status(path(name + ".mosaic_status.tsv"));
        require(st.at("union_scoring_backend") == "indexed_fine_experimental" &&
                    st.at("union_scoring_backend_stages") == "initial_union_pass" && st.at("indexed_union_refusal") == "none" &&
                    std::stoull(st.at("indexed_union_store_charged_bytes")) > 0,
                "the indexed run must record its backend, the stage it served and its store");
        const std::string settings = slurp(path(name + ".settings.tsv"));
        require(settings.find("union_scoring_backend\tindexed_fine_experimental\tprofile test") != std::string::npos &&
                    settings.find("indexed_union_store_bytes\t268435456\tprofile test") != std::string::npos,
                "the experimental keys must be recorded with their source");
    }
    {
        const std::string pr1 = path("phaseable_r1.fq"), pr2 = path("phaseable_r2.fq");
        panvar::GenotypeCommandOptions ref = repair_options("ix_repair_ref", pr1, pr2, "tiny_phaseable.gfa");
        ref.max_mosaic_candidates = 0;
        with_union_caches(ref);
        panvar::GenotypeCommandOptions ix = ref;
        ix.out_prefix = path("ix_repair_idx");
        require(panvar::run_genotype(ref) == panvar::GenotypeRunStatus::Complete &&
                    panvar::run_genotype(indexed(ix)) == panvar::GenotypeRunStatus::Complete,
                "fixture: both adaptive runs must call");
        require_same_science("ix_repair_ref", "ix_repair_idx", "indexed backend through the adaptive repair");
        // the repair trajectory's scientific columns (perf_compare.py's SCIENTIFIC_ROUND_COLUMNS: the first 28)
        const auto science = [](const std::string& file) {
            std::ifstream in(file);
            std::string line, out;
            while (std::getline(in, line)) {
                std::size_t at = 0;
                for (int k = 0; k < 28 && at != std::string::npos; ++k) at = line.find('\t', at + 1);
                out += line.substr(0, at) + "\n";
            }
            return out;
        };
        require(science(path("ix_repair_ref.repair_rounds.tsv")) == science(path("ix_repair_idx.repair_rounds.tsv")),
                "the repair trajectory must be the reference's");
        require(read_status(path("ix_repair_idx.mosaic_status.tsv")).at("pool_score_mismatches") == "0",
                "the repair must score the indexed union's candidates identically");
    }
    // refusals on the indexed route end RESCORE_INCOMPLETE, with no ranking and the refusal recorded
    {
        panvar::GenotypeCommandOptions capped = tiny_options("ix_refused_origins", r1, r2);
        with_union_caches(capped);
        capped.max_physical_origins = 1;
        require(panvar::run_genotype(indexed(capped)) == panvar::GenotypeRunStatus::RescoreIncomplete,
                "an exhausted physical-origin cap on the indexed route must end RESCORE_INCOMPLETE");
        const auto st = read_status(path("ix_refused_origins.mosaic_status.tsv"));
        require(st.at("note").find("resource cap") != std::string::npos &&
                    st.at("indexed_union_refusal").find("physical-origin cap") != std::string::npos &&
                    std::stoull(st.at("observed_physical_origins")) > 1 && st.count("observed_threads_used_placement"),
                "the refusal must be named and its observed work reported");
        require(no_ranking_or_call("ix_refused_origins") && no_final_outputs("ix_refused_origins"),
                "a refused indexed pass must leave no ranking, no call and no final table");

        panvar::GenotypeCommandOptions store = tiny_options("ix_refused_store", r1, r2);
        with_union_caches(store);
        store = indexed(store);
        store.indexed_union_store_bytes = 1;
        require(panvar::run_genotype(store) == panvar::GenotypeRunStatus::RescoreIncomplete,
                "an exhausted placement-store budget must end RESCORE_INCOMPLETE");
        require(read_status(path("ix_refused_store.mosaic_status.tsv")).at("indexed_union_refusal").find("placement store") !=
                        std::string::npos &&
                    no_ranking_or_call("ix_refused_store"),
                "the store refusal must be named, with no ranking or call");
    }
    // declared requirements
    {
        panvar::GenotypeCommandOptions bad;
        bool threw = false;
        try {
            panvar::apply_genotype_setting("profile test", "union_scoring_backend", "fast", bad);
        } catch (const std::exception&) {
            threw = true;
        }
        require(threw, "an unknown union backend must be refused");
        for (int which = 0; which < 3; ++which) {
            panvar::GenotypeCommandOptions o = tiny_options("ix_invalid", r1, r2);
            with_union_caches(o);
            if (which == 0) o = indexed(o), o.indexed_union_store_bytes = 0;
            if (which == 1) o.indexed_union_store_bytes = 1;   // without the indexed backend
            if (which == 2) o = indexed(o), o.max_union_cached_placements_consumed = 0;
            threw = false;
            try {
                (void)panvar::run_genotype(o);
            } catch (const std::exception&) {
                threw = true;
            }
            require(threw, "an inconsistent experimental backend setting must be refused before any work");
        }
    }
    std::cout << "  indexed union backend: default off and unrecorded; same science at 1 and 4 threads and through the "
                 "repair; refusals RESCORE_INCOMPLETE with no ranking; inconsistent settings refused\n";
}


// ---- I2-B (EXPERIMENTAL): fine window discovery in the adaptive search, at the caller -----------------------------
// Control: the I2-A indexed union with the reference adaptive search (key set explicitly to "reference"); fine: the same
// with adaptive_window_discovery fine_experimental. Same science, the complete repair trajectory's scientific columns,
// at 2 and 4 threads; a served-cap refusal stops the adaptive search with no genotype; inconsistent settings refused.
void test_adaptive_fine_windows() {
    const std::string pr1 = path("phaseable_r1.fq"), pr2 = path("phaseable_r2.fq");
    // default off: no adaptive_window block
    {
        const auto st = read_status(path("ix_repair_ref.mosaic_status.tsv"));
        for (const auto& kv : st) {
            require(kv.first.rfind("adaptive_window_", 0) != 0 && kv.first.rfind("adaptive_fine_", 0) != 0,
                    "a default run must not report the adaptive window block (" + kv.first + ")");
        }
    }
    const auto arm = [&](const std::string& name, bool fine, std::size_t threads) {
        panvar::GenotypeCommandOptions o = repair_options(name, pr1, pr2, "tiny_phaseable.gfa");
        o.max_mosaic_candidates = 0;
        o.threads = threads;
        with_union_caches(o);
        o.union_scoring_backend = "indexed_fine_experimental";
        o.indexed_union_store_bytes = std::size_t{256} << 20;
        o.adaptive_window_discovery = fine ? "fine_experimental" : "reference";
        o.setting_sources["adaptive_window_discovery"] = "profile test";
        if (fine) {
            o.adaptive_fine_store_bytes = std::size_t{64} << 20;
            o.adaptive_fine_max_served_placements = std::uint64_t{1} << 40;
        }
        return o;
    };
    const auto science = [](const std::string& file) {
        std::ifstream in(file);
        std::string line, out;
        while (std::getline(in, line)) {
            std::size_t at = 0;
            for (int k = 0; k < 28 && at != std::string::npos; ++k) at = line.find('\t', at + 1);
            out += line.substr(0, at) + "\n";
        }
        return out;
    };
    for (const std::size_t threads : {std::size_t{2}, std::size_t{4}}) {
        const std::string c = "fw_control_t" + std::to_string(threads), f = "fw_fine_t" + std::to_string(threads);
        require(panvar::run_genotype(arm(c, false, threads)) == panvar::GenotypeRunStatus::Complete &&
                    panvar::run_genotype(arm(f, true, threads)) == panvar::GenotypeRunStatus::Complete,
                "fixture: both adaptive arms must call");
        require_same_science(c, f, "fine adaptive windows, " + std::to_string(threads) + " threads");
        require(science(path(c + ".repair_rounds.tsv")) == science(path(f + ".repair_rounds.tsv")),
                "the complete adaptive trajectory must be the control's");
        const auto sc = read_status(path(c + ".mosaic_status.tsv")), sf = read_status(path(f + ".mosaic_status.tsv"));
        require(sf.at("pool_score_mismatches") == "0" && sc.at("adaptive_window_discovery") == "reference" &&
                    sf.at("adaptive_window_discovery") == "fine_experimental" && sc.at("adaptive_fine_units_served") == "0" &&
                    std::stoull(sf.at("adaptive_fine_units_discovered")) > 0 &&
                    sc.at("adaptive_window_enumerations") == sf.at("adaptive_window_enumerations"),
                "both arms must report their path; the fine arm must discover; both enumerate the same windows");
    }
    // a served-placement cap refusal: typed, the adaptive search stops, no genotype is written
    {
        panvar::GenotypeCommandOptions o = arm("fw_refused", true, 2);
        o.adaptive_fine_max_served_placements = 1;
        const panvar::GenotypeRunStatus status = panvar::run_genotype(o);
        const auto st = read_status(path("fw_refused.mosaic_status.tsv"));
        require(status != panvar::GenotypeRunStatus::Complete && no_final_outputs("fw_refused") &&
                    st.at("no_call_reason") == "cumulative_cap",
                "a served-placement refusal must stop the adaptive search with no genotype");
    }
    // declared requirements
    {
        bool threw = false;
        panvar::GenotypeCommandOptions bad;
        try {
            panvar::apply_genotype_setting("profile test", "adaptive_window_discovery", "fast", bad);
        } catch (const std::exception&) {
            threw = true;
        }
        require(threw, "an unknown adaptive window path must be refused");
        for (int which = 0; which < 2; ++which) {
            panvar::GenotypeCommandOptions o = arm("fw_invalid", which == 0, 2);
            if (which == 0) o.adaptive_fine_store_bytes = 0;            // fine without a store budget
            if (which == 1) o.adaptive_fine_store_bytes = 1 << 20;      // a budget without the fine path
            threw = false;
            try {
                (void)panvar::run_genotype(o);
            } catch (const std::exception&) {
                threw = true;
            }
            require(threw, "an inconsistent adaptive window setting must be refused before any work");
        }
    }
    std::cout << "  adaptive fine windows: default off; same science and trajectory at 2 and 4 threads; served-cap refusal "
                 "stops with no genotype; inconsistent settings refused\n";
}


// ---- I2-C: the union store layouts and the stage attribution block, at the caller ------------------------------------
void test_store_layout_and_timing_detail() {
    const std::string pr1 = path("phaseable_r1.fq"), pr2 = path("phaseable_r2.fq");
    const auto arm = [&](const std::string& name, const std::string& layout, bool timing) {
        panvar::GenotypeCommandOptions o = repair_options(name, pr1, pr2, "tiny_phaseable.gfa");
        o.max_mosaic_candidates = 0;
        with_union_caches(o);
        o.union_scoring_backend = "indexed_fine_experimental";
        o.indexed_union_store_bytes = std::size_t{256} << 20;
        o.indexed_union_store_layout = layout;
        o.adaptive_window_discovery = "fine_experimental";
        o.setting_sources["adaptive_window_discovery"] = "profile test";
        o.adaptive_fine_store_bytes = std::size_t{64} << 20;
        o.adaptive_fine_max_served_placements = std::uint64_t{1} << 40;
        o.timing_detail = timing;
        if (timing) o.setting_sources["timing_detail"] = "profile test";
        return o;
    };
    require(panvar::run_genotype(arm("sl_legacy", "legacy_i2a", true)) == panvar::GenotypeRunStatus::Complete &&
                panvar::run_genotype(arm("sl_exact", "exact", true)) == panvar::GenotypeRunStatus::Complete &&
                panvar::run_genotype(arm("sl_plain", "exact", false)) == panvar::GenotypeRunStatus::Complete,
            "fixture: every layout arm must call");
    require_same_science("sl_legacy", "sl_exact", "exact against legacy store layout");
    const auto legacy = read_status(path("sl_legacy.mosaic_status.tsv")), exact = read_status(path("sl_exact.mosaic_status.tsv"));
    require(legacy.at("indexed_union_store_layout") == "legacy_i2a" && exact.at("indexed_union_store_layout") == "exact" &&
                exact.at("indexed_union_store_budget_basis") == "allocated_capacity" &&
                std::stoull(exact.at("indexed_union_store_peak_bytes")) < std::stoull(legacy.at("indexed_union_store_peak_bytes")),
            "both layouts report themselves; the exact layout peaks lower");
    // the attribution block: present only when asked for, disjoint parts, remainders reported (never negative)
    for (const auto& kv : read_status(path("sl_plain.mosaic_status.tsv"))) {
        require(kv.first.rfind("seconds_detail_", 0) != 0, "the attribution block must be absent without timing_detail");
    }
    for (const char* key : {"seconds_detail_model_graph_load_sort_parse", "seconds_detail_markers_read_fastq",
                            "seconds_detail_model_unattributed", "seconds_detail_markers_unattributed",
                            "seconds_detail_adaptive_seed_placement_census", "seconds_detail_adaptive_neighbour_priors",
                            "seconds_detail_adaptive_pool_hook", "seconds_detail_adaptive_round_unattributed",
                            "seconds_detail_adaptive_post_search"}) {
        require(exact.count(key) != 0, std::string("the attribution block must report ") + key);
    }
    for (const char* key : {"seconds_detail_model_unattributed", "seconds_detail_markers_unattributed",
                            "seconds_detail_adaptive_round_unattributed", "seconds_detail_adaptive_run_outside_rounds",
                            "seconds_detail_adaptive_post_search"}) {
        require(std::stod(exact.at(key)) > -0.005, std::string("a remainder cannot be materially negative: ") + key);
    }
    bool threw = false;
    try {
        panvar::GenotypeCommandOptions o = tiny_options("sl_invalid", path("r1.fq"), path("r2.fq"));
        o.indexed_union_store_layout = "legacy_i2a";   // without the indexed backend
        (void)panvar::run_genotype(o);
    } catch (const std::exception&) {
        threw = true;
    }
    require(threw, "a store layout without the indexed backend must be refused");
    std::cout << "  store layouts and attribution: exact = legacy science, exact peaks lower; detail block only when asked, "
                 "remainders reported; a stray layout key refused\n";
}


// ---- I2-D: fused anchor seeding at the caller --------------------------------------------------------------------------
void test_fused_anchor_seeding() {
    const std::string pr1 = path("phaseable_r1.fq"), pr2 = path("phaseable_r2.fq");
    const auto arm = [&](const std::string& name, bool fused, std::size_t threads) {
        panvar::GenotypeCommandOptions o = repair_options(name, pr1, pr2, "tiny_phaseable.gfa");
        o.max_mosaic_candidates = 0;
        o.threads = threads;
        with_union_caches(o);
        o.union_scoring_backend = "indexed_fine_experimental";
        o.indexed_union_store_bytes = std::size_t{256} << 20;
        o.adaptive_window_discovery = "fine_experimental";
        o.setting_sources["adaptive_window_discovery"] = "profile test";
        o.adaptive_fine_store_bytes = std::size_t{64} << 20;
        o.adaptive_fine_max_served_placements = std::uint64_t{1} << 40;
        o.adaptive_fused_seed = fused;
        o.setting_sources["adaptive_fused_seed"] = "profile test";
        return o;
    };
    const auto science = [](const std::string& file) {
        std::ifstream in(file);
        std::string line, out;
        while (std::getline(in, line)) {
            std::size_t at = 0;
            for (int k = 0; k < 28 && at != std::string::npos; ++k) at = line.find('\t', at + 1);
            out += line.substr(0, at) + "\n";
        }
        return out;
    };
    for (const std::size_t threads : {std::size_t{2}, std::size_t{4}}) {
        const std::string c = "fs_control_t" + std::to_string(threads), f = "fs_fused_t" + std::to_string(threads);
        require(panvar::run_genotype(arm(c, false, threads)) == panvar::GenotypeRunStatus::Complete &&
                    panvar::run_genotype(arm(f, true, threads)) == panvar::GenotypeRunStatus::Complete,
                "fixture: both seeding arms must call");
        require_same_science(c, f, "fused seeding, " + std::to_string(threads) + " threads");
        require(science(path(c + ".repair_rounds.tsv")) == science(path(f + ".repair_rounds.tsv")),
                "the complete adaptive trajectory must be the control's");
        const auto sc = read_status(path(c + ".mosaic_status.tsv")), sf = read_status(path(f + ".mosaic_status.tsv"));
        const std::uint64_t pc = std::stoull(sc.at("adaptive_full_placements")), pf = std::stoull(sf.at("adaptive_full_placements"));
        const std::uint64_t seeded = sc.at("repair_anchor_hap1") == sc.at("repair_anchor_hap2") ? 1 : 2;
        require(sf.at("adaptive_fused_seed") == "1" && sc.at("adaptive_fused_seed") == "0" && pc == pf + seeded &&
                    sf.at("pool_score_mismatches") == "0",
                "the fused arm must place each distinct anchor once fewer, with nothing else changed");
    }
    std::cout << "  fused anchor seeding: same science and trajectory at 2 and 4 threads; one fewer full placement per "
                 "distinct anchor\n";
}


// ---- I2-E: marker-proposal contribution reuse at the caller --------------------------------------------------------------
void test_marker_proposal_reuse() {
    const std::string pr1 = path("phaseable_r1.fq"), pr2 = path("phaseable_r2.fq");
    const auto arm = [&](const std::string& name, std::size_t budget) {
        panvar::GenotypeCommandOptions o = repair_options(name, pr1, pr2, "tiny_phaseable.gfa");
        o.max_mosaic_candidates = 8;   // the beam runs, so its candidates come from these emissions
        o.marker_proposal_reuse_table_bytes = budget;
        o.setting_sources["marker_proposal_reuse_table_bytes"] = "profile test";
        return o;
    };
    require(panvar::run_genotype(arm("mp_off", 0)) == panvar::GenotypeRunStatus::Complete &&
                panvar::run_genotype(arm("mp_on", std::size_t{4} << 20)) == panvar::GenotypeRunStatus::Complete &&
                panvar::run_genotype(arm("mp_tiny", 16)) == panvar::GenotypeRunStatus::Complete,
            "fixture: every reuse arm must call");
    require_same_science("mp_off", "mp_on", "marker-proposal reuse");
    require_same_science("mp_off", "mp_tiny", "marker-proposal reuse, tiny budget");
    const auto off = read_status(path("mp_off.mosaic_status.tsv")), on = read_status(path("mp_on.mosaic_status.tsv")),
               tiny = read_status(path("mp_tiny.mosaic_status.tsv"));
    // (this fixture's terms are each needed once, so reuse serves none here; the unit gates prove reuse itself)
    require(off.at("marker_proposal_reuse") == "off" && on.at("marker_proposal_reuse") == "on" &&
                std::stoull(on.at("marker_proposal_contributions_executed")) +
                        std::stoull(on.at("marker_proposal_contributions_reused")) ==
                    std::stoull(on.at("marker_proposal_contribution_visits")) &&
                on.at("marker_proposal_contribution_visits") == off.at("marker_proposal_contribution_visits") &&
                tiny.at("marker_proposal_contributions_fallback") == tiny.at("marker_proposal_contribution_visits") &&
                std::stoull(tiny.at("marker_proposal_reuse_peak_table_bytes")) <= 16,
            "the reuse block must account for every visit; a tiny budget falls back for all and allocates within it");
    for (const auto& kv : read_status(path("ix_default.mosaic_status.tsv"))) {
        require(kv.first.rfind("marker_proposal_", 0) != 0, "a default run must not report the reuse block");
    }
    std::cout << "  marker-proposal reuse: same science (candidate tables included) on, off and at a tiny budget; work "
                 "reported; default runs unchanged\n";
}


// ---- I2-F: prepared main-library combination terms at the caller --------------------------------------------------------
void test_prepared_combination_terms() {
    const std::string pr1 = path("phaseable_r1.fq"), pr2 = path("phaseable_r2.fq");
    const auto science = [](const std::string& file) {
        std::ifstream in(file);
        std::string line, out;
        while (std::getline(in, line)) {
            std::size_t at = 0;
            for (int k = 0; k < 28 && at != std::string::npos; ++k) at = line.find('\t', at + 1);
            out += line.substr(0, at) + "\n";
        }
        return out;
    };
    for (const std::size_t threads : {std::size_t{1}, std::size_t{4}}) {
        const auto arm = [&](const std::string& name, bool prepared) {
            panvar::GenotypeCommandOptions o = repair_options(name, pr1, pr2, "tiny_phaseable.gfa");
            o.max_mosaic_candidates = 8;
            o.threads = threads;
            o.combination_prepared_library_terms = prepared;
            o.setting_sources["combination_prepared_library_terms"] = "profile test";
            return o;
        };
        const std::string c = "pc_scalar_t" + std::to_string(threads), p = "pc_prepared_t" + std::to_string(threads);
        require(panvar::run_genotype(arm(c, false)) == panvar::GenotypeRunStatus::Complete &&
                    panvar::run_genotype(arm(p, true)) == panvar::GenotypeRunStatus::Complete,
                "fixture: both combination arms must call");
        require_same_science(c, p, "prepared combination terms, " + std::to_string(threads) + " threads");
        require(science(path(c + ".repair_rounds.tsv")) == science(path(p + ".repair_rounds.tsv")),
                "the adaptive trajectory must be the scalar path's");
        const auto sc = read_status(path(c + ".mosaic_status.tsv")), sp = read_status(path(p + ".mosaic_status.tsv"));
        require(sc.at("combination_prepared_calls") == "0" && std::stoull(sp.at("combination_prepared_calls")) >= 2 &&
                    std::stoll(sp.at("combination_net_validations_avoided")) > 0 &&
                    std::stoll(sp.at("combination_net_logs_avoided")) > 0 &&
                    std::stoull(sp.at("combination_hypothetical_scalar_logs")) >= std::stoull(sp.at("combination_prepared_logs")),
                "the prepared arm must report the union and adaptive combinations it prepared");
    }
    for (const auto& kv : read_status(path("ix_default.mosaic_status.tsv"))) {
        require(kv.first.rfind("combination_", 0) != 0, "a default run must not report the combination block");
    }
    std::cout << "  prepared combination terms: same science and trajectory at 1 and 4 threads; union and adaptive "
                 "combinations prepared; default runs unchanged\n";
}


// ---- I2-H: the fixed-geometry emitting recursion at the caller ----------------------------------------------------------
void test_gap_recursion_backend() {
    const std::string pr1 = path("phaseable_r1.fq"), pr2 = path("phaseable_r2.fq");
    const auto science = [](const std::string& file) {
        std::ifstream in(file);
        std::string line, out;
        while (std::getline(in, line)) {
            std::size_t at = 0;
            for (int k = 0; k < 28 && at != std::string::npos; ++k) at = line.find('\t', at + 1);
            out += line.substr(0, at) + "\n";
        }
        return out;
    };
    for (const std::size_t threads : {std::size_t{1}, std::size_t{4}}) {
        const auto arm = [&](const std::string& name, bool fixed) {
            panvar::GenotypeCommandOptions o = repair_options(name, pr1, pr2, "tiny_phaseable.gfa");
            with_union_caches(o);
            o.max_mosaic_candidates = 8;
            o.threads = threads;
            o.gap_recursion_fixed_geometry = fixed;
            o.setting_sources["gap_recursion_backend"] = "profile test";
            return o;
        };
        const std::string c = "gr_scaled_t" + std::to_string(threads), e = "gr_fixed_t" + std::to_string(threads);
        require(panvar::run_genotype(arm(c, false)) == panvar::GenotypeRunStatus::Complete &&
                    panvar::run_genotype(arm(e, true)) == panvar::GenotypeRunStatus::Complete,
                "fixture: both recursion arms must call");
        require_same_science(c, e, "gap recursion backend, " + std::to_string(threads) + " threads");
        require(science(path(c + ".repair_rounds.tsv")) == science(path(e + ".repair_rounds.tsv")),
                "the adaptive trajectory must be the reference recursion's");
        const auto sc = read_status(path(c + ".mosaic_status.tsv")), se = read_status(path(e + ".mosaic_status.tsv"));
        require(sc.at("gap_recursion_backend") == "scaled_reference" &&
                    se.at("gap_recursion_backend") == "fixed_geometry_experimental" &&
                    se.at("gap_recursion_backend_scope") == "union_adaptive_dosage",
                "each arm must record its recursion backend");
        require(std::stoull(se.at("dp_evaluations")) > 0, "the fixture must reach the recursion (vacuous otherwise)");
    }
    for (const auto& kv : read_status(path("ix_default.mosaic_status.tsv"))) {
        require(kv.first.rfind("gap_recursion_", 0) != 0, "a default run must not report the recursion backend");
    }
    std::cout << "  gap recursion backend: same science and trajectory at 1 and 4 threads; recorded; default runs unchanged\n";
}


// ---- I2-K: the union pass's emission results handed to the adaptive search ---------------------------------------------
void test_emission_handoff() {
    const std::string pr1 = path("phaseable_r1.fq"), pr2 = path("phaseable_r2.fq");
    const auto science = [](const std::string& file) {
        std::ifstream in(file);
        std::string line, out;
        while (std::getline(in, line)) {
            std::size_t at = 0;
            for (int k = 0; k < 28 && at != std::string::npos; ++k) at = line.find('\t', at + 1);
            out += line.substr(0, at) + "\n";
        }
        return out;
    };
    for (const std::size_t threads : {std::size_t{1}, std::size_t{4}}) {
        const auto arm = [&](const std::string& name, bool handoff) {
            panvar::GenotypeCommandOptions o = repair_options(name, pr1, pr2, "tiny_phaseable.gfa");
            with_union_caches(o);
            o.max_mosaic_candidates = 8;
            o.threads = threads;
            o.union_scoring_backend = "indexed_fine_experimental";   // the batch that detaches its emission results
            o.indexed_union_store_bytes = std::size_t{256} << 20;
            o.setting_sources["union_scoring_backend"] = "profile test";
            o.setting_sources["indexed_union_store_bytes"] = "profile test";
            o.adaptive_window_discovery = "fine_experimental";
            o.adaptive_fine_store_bytes = std::size_t{64} << 20;
            o.adaptive_fine_max_served_placements = 2000000000;
            o.setting_sources["adaptive_window_discovery"] = "profile test";
            o.setting_sources["adaptive_fine_store_bytes"] = "profile test";
            o.setting_sources["adaptive_fine_max_served_placements"] = "profile test";
            o.adaptive_fused_seed = true;
            o.setting_sources["adaptive_fused_seed"] = "profile test";
            o.repair_tie_branching = true;   // any branch context starts empty: the handoff is adopted once
            o.emission_cache_bytes = std::size_t{64} << 20;   // both arms: the cache whose results are handed over
            o.adaptive_emission_handoff = handoff;
            o.setting_sources["adaptive_emission_handoff"] = "profile test";
            return o;
        };
        const std::string c = "eh_off_t" + std::to_string(threads), e = "eh_on_t" + std::to_string(threads);
        require(panvar::run_genotype(arm(c, false)) == panvar::GenotypeRunStatus::Complete &&
                    panvar::run_genotype(arm(e, true)) == panvar::GenotypeRunStatus::Complete,
                "fixture: both handoff arms must call");
        require_same_science(c, e, "emission handoff, " + std::to_string(threads) + " threads");
        require(science(path(c + ".repair_rounds.tsv")) == science(path(e + ".repair_rounds.tsv")),
                "the adaptive trajectory must be the cold search's");
        const auto sc = read_status(path(c + ".mosaic_status.tsv")), se = read_status(path(e + ".mosaic_status.tsv"));
        require(sc.at("adaptive_emission_handoff") == "0" && sc.at("adaptive_emission_handoff_detached") == "0" &&
                    sc.at("adaptive_emission_handoff_inherited_hits") == "0",
                "the off arm detaches and inherits nothing");
        require(se.at("adaptive_emission_handoff_detached") == "1" && se.at("adaptive_emission_handoff_adopted") == "1" &&
                    se.at("adaptive_emission_handoff_refusal") == "none" &&
                    std::stoull(se.at("adaptive_emission_handoff_entries")) > 0 &&
                    std::stoull(se.at("adaptive_emission_handoff_inherited_hits")) > 0,
                "the on arm must adopt the union's results and serve hits from them");
        // a sequential sample in the same process: a fresh run with the handoff off reports and inherits nothing
        require(panvar::run_genotype(arm(c + "_again", false)) == panvar::GenotypeRunStatus::Complete,
                "fixture: a sequential run must call");
        require_same_science(c, c + "_again", "sequential sample after a handoff run");
        require(read_status(path(c + "_again.mosaic_status.tsv")).at("adaptive_emission_handoff_inherited_hits") == "0",
                "a sequential sample inherits nothing");
    }
    std::cout << "  emission handoff: same science and trajectory at 1 and 4 threads (rounds, fused seed, fine windows, tie "
                 "branching on); adopted with inherited hits; a sequential sample inherits nothing\n";
}


// ---- I3-A: one prepared locus, successive samples -----------------------------------------------------------------------
// Two genuinely different samples (an off-panel phaseable truth, and a pair of panel haplotypes from another seed) on one
// graph, analysed against ONE prepared locus in A-B-A, reversed and repeated orders, at 1 and 4 threads, with two
// declared libraries and with --library auto. Every prepared run must equal the ordinary (cold) run of the same sample
// and settings; the object must be prepared once and stay byte-identical (digest); incompatible graph, reference,
// bubbles, curation and marker or bubble settings must be refused before any output; a failing sample must not affect
// the next; and the command path (dosage, locus table) must agree too.
void test_prepared_locus() {
    const TinyLocus locus = write_tiny_locus("prep_locus.gfa", 80, 70);
    const auto truths = spell_truths(locus);
    simulate(truths.first, truths.second, path("prep_A_r1.fq"), path("prep_A_r2.fq"), 123);
    simulate(locus.haplotypes.front(), locus.haplotypes.back(), path("prep_B_r1.fq"), path("prep_B_r2.fq"), 456);
    {
        std::ofstream lib2(path("prep_lib2.txt"));
        lib2 << "insert_mean_bp 290\ninsert_sd_bp 35\nsupport_sigmas 4\nsubstitution_error 0.004\noutlier_weight 0.002\n"
                "fragment_start_intensity 0.08\nmarker_background 0.01\n";
    }
    const auto sample = [&](const std::string& name, char which, std::size_t threads, const std::string& library) {
        const std::string r1 = path(std::string("prep_") + which + "_r1.fq"), r2 = path(std::string("prep_") + which + "_r2.fq");
        panvar::GenotypeCommandOptions o = repair_options(name, r1, r2, "prep_locus.gfa");
        o.max_mosaic_candidates = 8;
        o.threads = threads;
        o.library_path = library;
        return o;
    };
    const std::string lib1 = path("lib.txt"), lib2 = path("prep_lib2.txt");
    const auto trajectory = [](const std::string& file) {   // the rounds' scientific columns (timing and memory follow)
        std::ifstream in(file);
        std::string line, out;
        while (std::getline(in, line)) {
            std::size_t at = 0;
            for (int k = 0; k < 28 && at != std::string::npos; ++k) at = line.find('\t', at + 1);
            out += line.substr(0, at) + "\n";
        }
        return out;
    };
    // the cold (ordinary) runs, one per distinct sample configuration
    struct Config { char which; std::size_t threads; std::string library, tag; };
    const std::vector<Config> configs = {{'A', 1, lib1, "A_t1"}, {'B', 1, lib1, "B_t1"}, {'A', 4, lib1, "A_t4"},
                                         {'B', 4, lib1, "B_t4"}, {'B', 1, lib2, "B_lib2"}, {'A', 1, "auto", "A_auto"}};
    std::map<std::string, panvar::GenotypeRunStatus> cold_status;
    for (const Config& c : configs) {
        cold_status[c.tag] = panvar::run_genotype(sample("prep_cold_" + c.tag, c.which, c.threads, c.library));
        require(read_status(path("prep_cold_" + c.tag + ".mosaic_status.tsv")).count("prepared_locus") == 0,
                "an ordinary run reports no prepared locus");
    }
    require(cold_status["A_t1"] == panvar::GenotypeRunStatus::Complete && cold_status["B_t1"] == panvar::GenotypeRunStatus::Complete,
            "fixture: both samples must call");
    require(slurp(path("prep_cold_A_t1.genotypes.tsv")) != slurp(path("prep_cold_B_t1.genotypes.tsv")),
            "fixture: the two samples must be genuinely different (different calls)");

    // prepare once
    const std::uint64_t builds_before = panvar::prepared_locus_builds();
    const std::shared_ptr<const panvar::PreparedLocus> prepared = panvar::prepare_genotype_locus(sample("prep_prepare", 'A', 1, lib1));
    const std::string digest = prepared->digest();
    require(panvar::prepared_locus_builds() == builds_before + 1, "one preparation");
    require(prepared->fingerprint() == read_status(path("prep_cold_A_t1.mosaic_status.tsv")).at("model_fingerprint"),
            "the prepared model is the ordinary run's model");
    // A-B-A, then reversed, then repeated, then libraries; every one against the same object
    const std::vector<std::string> order = {"A_t1", "B_t1", "A_t1", "B_t4", "A_t4", "A_t4", "B_t1", "B_lib2", "A_auto", "B_lib2"};
    std::size_t runs = 0;
    for (const std::string& tag : order) {
        const Config& c = *std::find_if(configs.begin(), configs.end(), [&](const Config& x) { return x.tag == tag; });
        const std::string name = "prep_run" + std::to_string(runs++) + "_" + tag;
        const panvar::GenotypeRunStatus status = panvar::run_genotype_prepared(sample(name, c.which, c.threads, c.library), *prepared);
        require(status == cold_status[tag], name + ": the run status must be the cold run's");
        require_same_science("prep_cold_" + tag, name, "prepared locus, " + name);
        require(exists(path("prep_cold_" + tag + ".repair_rounds.tsv")) == exists(path(name + ".repair_rounds.tsv")),
                name + ": the adaptive trajectory must exist in both or neither");
        if (exists(path("prep_cold_" + tag + ".repair_rounds.tsv"))) {
            require(trajectory(path("prep_cold_" + tag + ".repair_rounds.tsv")) == trajectory(path(name + ".repair_rounds.tsv")),
                    name + ": the adaptive trajectory must be the cold run's (scientific columns)");
        }
        const auto st = read_status(path(name + ".mosaic_status.tsv"));
        require(st.at("prepared_locus") == "1" && st.at("prepared_locus_id") == std::to_string(prepared->id()),
                name + ": the status must record the prepared locus it used");
        require(prepared->digest() == digest, name + ": the prepared locus must be unchanged by an analysis");
    }
    // a failing sample (unreadable reads), then a valid one
    {
        panvar::GenotypeCommandOptions broken = sample("prep_broken", 'A', 1, lib1);
        broken.reads1_path = path("prep_missing_r1.fq");
        bool failed = false;
        try {
            failed = panvar::run_genotype_prepared(broken, *prepared) != panvar::GenotypeRunStatus::Complete;
        } catch (const std::exception&) {
            failed = true;
        }
        require(failed, "fixture: a sample with unreadable reads must fail");
        require(panvar::run_genotype_prepared(sample("prep_after_broken", 'B', 1, lib1), *prepared) == cold_status["B_t1"],
                "the sample after a failure must run");
        require_same_science("prep_cold_B_t1", "prep_after_broken", "after a failing sample");
    }
    // incompatible inputs and settings: refused before any output, never analysed against the object
    {
        write_tiny_locus("prep_other.gfa", 80, 71);
        {
            std::ofstream curation(path("prep_curation.txt"));
            curation << locus.reference_name << '\n';
        }
        {
            std::ofstream bubbles(path("prep_bubbles.csv"));
            bubbles << "not,used\n";
        }
        const std::vector<std::pair<std::string, std::function<void(panvar::GenotypeCommandOptions&)>>> changes = {
            {"graph", [&](panvar::GenotypeCommandOptions& o) { o.gfa_path = path("prep_other.gfa"); }},
            {"reference", [&](panvar::GenotypeCommandOptions& o) { o.reference_path = "sampleA#1#chr1:0-101"; }},
            {"bubbles", [&](panvar::GenotypeCommandOptions& o) { o.bubbles_path = path("prep_bubbles.csv"); }},
            {"min_variant_bp", [&](panvar::GenotypeCommandOptions& o) { o.min_variant_bp = 40; }},
            {"curation", [&](panvar::GenotypeCommandOptions& o) { o.curation_manifest_path = path("prep_curation.txt"); }},
            {"kmer", [&](panvar::GenotypeCommandOptions& o) { o.kmer_size = 23; }},
            {"syncmer", [&](panvar::GenotypeCommandOptions& o) { o.syncmer_s = 9; }},
        };
        for (const auto& [what, change] : changes) {
            panvar::GenotypeCommandOptions o = sample("prep_changed_" + what, 'A', 1, lib1);
            change(o);
            std::filesystem::remove(path("prep_changed_" + what + ".mosaic_status.tsv"));   // a scratch dir may hold an old run
            bool refused = false;
            try {
                (void)panvar::run_genotype_prepared(o, *prepared);
            } catch (const panvar::PreparedLocusMismatch&) {
                refused = true;
            }
            require(refused, "a changed " + what + " must refuse the prepared locus");
            require(!exists(path("prep_changed_" + what + ".mosaic_status.tsv")), "a refused sample writes no output (" + what + ")");
        }
        // after the refusals, the object still serves the valid samples
        require(panvar::run_genotype_prepared(sample("prep_after_refusals", 'A', 1, lib1), *prepared) == cold_status["A_t1"],
                "the object still serves after refusals");
        require_same_science("prep_cold_A_t1", "prep_after_refusals", "after refusals");
        // a changed graph needs its own preparation (one more build), and then runs
        panvar::GenotypeCommandOptions other = sample("prep_other_cold", 'A', 1, lib1);
        other.gfa_path = path("prep_other.gfa");
        const panvar::GenotypeRunStatus other_cold = panvar::run_genotype(other);
        const std::uint64_t before_other = panvar::prepared_locus_builds();
        const auto other_prepared = panvar::prepare_genotype_locus(other);
        require(panvar::prepared_locus_builds() == before_other + 1, "a changed graph is prepared again");
        other.out_prefix = path("prep_other_prepared");
        require(panvar::run_genotype_prepared(other, *other_prepared) == other_cold, "the re-prepared graph runs");
        require_same_science("prep_other_cold", "prep_other_prepared", "re-prepared graph");
    }
    // the command path (genotype, dosage, locus table) against the prepared locus, at the CLI's own defaults. Its outputs
    // go in a subdirectory: the command appends the headline filter columns, which test_block_quality_columns (a scan of
    // the scratch directory's run_genotype tables) does not expect.
    {
        std::filesystem::create_directories(path("prep_cmd"));
        const auto args = [&](const std::string& prefix, char which) {
            return std::vector<std::string>{"--graph", path("prep_locus.gfa"), "--reference-path", "sampleA#1#chr1:0-100",
                                            "--reads1", path(std::string("prep_") + which + "_r1.fq"), "--reads2",
                                            path(std::string("prep_") + which + "_r2.fq"), "--library", lib1, "--kmer-size",
                                            "21", "--syncmer-s", "11", "--max-mosaic-candidates", "8", "--threads", "1",
                                            "--quiet", "--out", path(prefix)};
        };
        const std::uint64_t before_cmd = panvar::prepared_locus_builds();
        const auto command_locus = panvar::prepare_genotype_locus(args("prep_cmd/prepare", 'A'));
        const std::string command_digest = command_locus->digest();
        for (const char which : {'A', 'B', 'A'}) {
            const std::string cold = std::string("prep_cmd/cold_") + which, warm = std::string("prep_cmd/warm_") + which;
            const int a = panvar::run_genotype_command(args(cold, which));
            const int b = panvar::run_genotype_command_prepared(args(warm, which), *command_locus);
            require(a == b, std::string("command exit codes must agree for sample ") + which);
            const auto sa = read_status(path(cold + ".mosaic_status.tsv")), sb = read_status(path(warm + ".mosaic_status.tsv"));
            std::set<std::string> keys;
            for (const auto& kv : sa) keys.insert(kv.first);
            for (const auto& kv : sb) keys.insert(kv.first);
            for (const std::string& k : keys) {
                const auto x = sa.find(k), y = sb.find(k);
                const bool same = x != sa.end() && y != sb.end() && x->second == y->second;
                require(same || indexed_allowed_key(k) || k == "prefix" || k == "out_prefix",
                        std::string("command path, sample ") + which + ": status key " + k + " differs");
            }
            for (const char* suffix : {".genotypes.tsv", ".block_genotypes.tsv", ".mosaic_candidates.tsv", ".mosaic_scores.tsv",
                                       ".mosaic_haplotypes.fa", ".equivalence_set.tsv", ".library_sensitivity.tsv",
                                       ".dosage_status.tsv", ".dosage_blocks.tsv", ".dosage_census.tsv", ".dosage_features.tsv",
                                       ".dosage_repairs.tsv", ".locus_genotype.tsv"}) {
                require(exists(path(cold + suffix)) == exists(path(warm + suffix)), std::string("command path: ") + suffix + " presence");
                if (!exists(path(cold + suffix))) continue;
                std::string x = slurp(path(cold + suffix)), y = slurp(path(warm + suffix));
                if (std::string(suffix) == ".dosage_status.tsv") {   // its timing rows differ by run
                    const auto strip = [](const std::string& t) {
                        std::string out;
                        std::istringstream in(t);
                        for (std::string line; std::getline(in, line);) if (line.rfind("seconds_", 0) != 0) out += line + "\n";
                        return out;
                    };
                    x = strip(x);
                    y = strip(y);
                }
                const auto placeheld = [](std::string text, const std::string& prefix) {   // the run's own name only
                    for (std::size_t at = text.find(prefix); at != std::string::npos; at = text.find(prefix, at + 5)) {
                        text.replace(at, prefix.size(), "<RUN>");
                    }
                    return text;
                };
                const auto base = [](const std::string& prefix) { return prefix.substr(prefix.rfind('/') + 1); };
                x = placeheld(placeheld(placeheld(x, path(cold)), cold), base(cold));
                y = placeheld(placeheld(placeheld(y, path(warm)), warm), base(warm));
                require(x == y, std::string("command path: ") + suffix + " differs for sample " + which);
            }
        }
        require(panvar::prepared_locus_builds() == before_cmd + 1 && command_locus->digest() == command_digest,
                "the command route prepares once and leaves the object unchanged");
    }
    require(prepared->digest() == digest, "the prepared locus is unchanged after every analysis");
    std::cout << "  prepared locus: " << runs << " prepared runs (A-B-A, reversed, repeated; 1 and 4 threads; two declared "
              << "libraries and auto) equal their cold runs; prepared once; digest unchanged; 7 incompatible changes refused "
                 "before output; a failing sample then a valid one; a changed graph re-prepared; command path (dosage, "
                 "locus table) identical\n";
}


// ---- I3-A.1: between-sample cleanup is scientifically invisible ---------------------------------------------------------
// The same distinct samples as test_prepared_locus (run after it, reusing its cold references made WITHOUT any cleanup):
// A-B-A against one prepared locus with release_between_samples() after every sample, a failing sample followed by a
// valid one, 1 and 4 threads, and cold runs with the release between them. Every result must equal its cold reference
// (science, status, trajectory); the release must empty the process-global residue; the prepared locus must be
// untouched (digest) and not rebuilt.
void test_prepared_locus_cleanup() {
    const std::string lib1 = path("lib.txt");
    const auto sample = [&](const std::string& name, char which, std::size_t threads) {
        panvar::GenotypeCommandOptions o = repair_options(name, path(std::string("prep_") + which + "_r1.fq"),
                                                          path(std::string("prep_") + which + "_r2.fq"), "prep_locus.gfa");
        o.max_mosaic_candidates = 8;
        o.threads = threads;
        o.library_path = lib1;
        return o;
    };
    const auto trajectory = [](const std::string& file) {
        std::ifstream in(file);
        std::string line, out;
        while (std::getline(in, line)) {
            std::size_t at = 0;
            for (int k = 0; k < 28 && at != std::string::npos; ++k) at = line.find('\t', at + 1);
            out += line.substr(0, at) + "\n";
        }
        return out;
    };
    const auto same_as_cold = [&](const std::string& cold, const std::string& run) {
        require_same_science(cold, run, "cleanup, " + run);
        require(exists(path(cold + ".repair_rounds.tsv")) == exists(path(run + ".repair_rounds.tsv")) &&
                    (!exists(path(cold + ".repair_rounds.tsv")) ||
                     trajectory(path(cold + ".repair_rounds.tsv")) == trajectory(path(run + ".repair_rounds.tsv"))),
                run + ": the adaptive trajectory must be the cold run's");
    };
    const auto release = [&](const std::string& after) {
        const panvar::BetweenSampleCleanup c = panvar::release_between_samples();
        require(panvar::retained_process_storage().total() == 0, after + ": the release must empty the process-global residue");
        return c.released.total();
    };
    const std::shared_ptr<const panvar::PreparedLocus> prepared = panvar::prepare_genotype_locus(sample("clean_prepare", 'A', 1));
    const std::string digest = prepared->digest();
    const std::uint64_t builds = panvar::prepared_locus_builds();
    std::uint64_t released = 0;
    const std::vector<std::pair<std::string, std::string>> order = {{"A_t1", "clean_run0"}, {"B_t1", "clean_run1"},
                                                                    {"A_t1", "clean_run2"}, {"B_t4", "clean_run3"},
                                                                    {"A_t4", "clean_run4"}};
    for (const auto& [tag, name] : order) {
        const char which = tag[0];
        const std::size_t threads = tag.back() == '4' ? 4 : 1;
        (void)panvar::run_genotype_prepared(sample(name, which, threads), *prepared);
        same_as_cold("prep_cold_" + tag, name);
        released += release(name);
        require(prepared->digest() == digest, name + ": the prepared locus is untouched by the release");
    }
    // a failing sample, the release, then a valid one
    {
        panvar::GenotypeCommandOptions broken = sample("clean_broken", 'A', 1);
        broken.reads1_path = path("prep_missing_r1.fq");
        try {
            (void)panvar::run_genotype_prepared(broken, *prepared);
        } catch (const std::exception&) {
        }
        released += release("clean_broken");
        (void)panvar::run_genotype_prepared(sample("clean_after_broken", 'B', 1), *prepared);
        same_as_cold("prep_cold_B_t1", "clean_after_broken");
    }
    // cold runs with the release between them
    for (const std::string tag : {"A_t1", "B_t1", "A_t4"}) {
        const std::string name = "clean_cold_" + tag;
        (void)panvar::run_genotype(sample(name, tag[0], tag.back() == '4' ? 4 : 1));
        same_as_cold("prep_cold_" + tag, name);
        released += release(name);
    }
    require(panvar::prepared_locus_builds() == builds && prepared->digest() == digest, "prepared once, untouched");
    require(released > 0, "the samples must leave process-global residue for the release to remove (vacuous otherwise)");
    std::cout << "  prepared locus cleanup: A-B-A, threads 1/4, a failure then success and cold runs, each with the "
              << "between-sample release (" << released << " bytes of residue released in total), equal their cold "
                 "references; the prepared locus untouched and not rebuilt\n";
}

// ---- I3-C: a persistent prepared-locus index, written once and loaded ------------------------------------------------
// The same distinct samples (cold references from test_prepared_locus): two independent preparations write identical
// bytes; the written index reads, rewrites byte-identically and loads to an object equal to the prepared one; samples
// A-B-A, 1 and 4 threads, a second declared library and --library auto analysed against the LOADED object equal their
// cold runs (science, status, trajectory); the status records the origin and the file checksum; a changed graph is
// refused before output.
void test_prepared_locus_index() {
    const std::string lib1 = path("lib.txt"), lib2 = path("prep_lib2.txt");
    const auto sample = [&](const std::string& name, char which, std::size_t threads, const std::string& library) {
        panvar::GenotypeCommandOptions o = repair_options(name, path(std::string("prep_") + which + "_r1.fq"),
                                                          path(std::string("prep_") + which + "_r2.fq"), "prep_locus.gfa");
        o.max_mosaic_candidates = 8;
        o.threads = threads;
        o.library_path = library;
        return o;
    };
    const auto trajectory = [](const std::string& file) {
        std::ifstream in(file);
        std::string line, out;
        while (std::getline(in, line)) {
            std::size_t at = 0;
            for (int k = 0; k < 28 && at != std::string::npos; ++k) at = line.find('\t', at + 1);
            out += line.substr(0, at) + "\n";
        }
        return out;
    };
    const panvar::PreparedLocusPreparation one = panvar::prepare_genotype_locus_record(sample("index_prepare1", 'A', 1, lib1));
    const panvar::PreparedLocusPreparation two = panvar::prepare_genotype_locus_record(sample("index_prepare2", 'B', 4, lib2));
    const std::string bytes = panvar::encode_prepared_locus(one.record);
    require(panvar::encode_prepared_locus(two.record) == bytes, "independent preparations (other sample options) write identical bytes");
    const std::string file = path("prep_index.pvpl"), rewrite = path("prep_index_rewrite.pvpl");
    std::filesystem::remove(file);      // one named file each: an earlier run of this test may have left them
    std::filesystem::remove(rewrite);
    panvar::write_prepared_locus(one.record, file);
    panvar::write_prepared_locus(panvar::read_prepared_locus_record(file, panvar::PreparedLocusLoadLimits{}), rewrite);
    require(slurp(rewrite) == bytes && slurp(file) == bytes, "write, read, rewrite: byte-identical");
    panvar::PreparedLocusLoadReport report;
    const std::uint64_t builds = panvar::prepared_locus_builds();
    const std::shared_ptr<const panvar::PreparedLocus> loaded = panvar::load_prepared_locus(file, panvar::PreparedLocusLoadLimits{}, &report);
    require(panvar::prepared_locus_builds() == builds, "loading is not preparing");
    require(loaded->digest() == one.locus->digest() && loaded->fingerprint() == one.locus->fingerprint() &&
                panvar::prepared_locus_model_content_sha256(loaded->model()) == one.record.model_content_sha256 &&
                panvar::prepared_locus_panel_content_sha256(loaded->panel()) == one.record.panel_content_sha256,
            "the loaded object equals the prepared one");
    const std::string digest = loaded->digest();
    struct Run { std::string tag; char which; std::size_t threads; std::string library; };
    const std::vector<Run> order = {{"A_t1", 'A', 1, lib1}, {"B_t1", 'B', 1, lib1}, {"A_t1", 'A', 1, lib1}, {"B_t4", 'B', 4, lib1},
                                    {"A_t4", 'A', 4, lib1}, {"B_lib2", 'B', 1, lib2}, {"A_auto", 'A', 1, "auto"}};
    std::size_t runs = 0;
    for (const Run& r : order) {
        const std::string name = "index_run" + std::to_string(runs++) + "_" + r.tag;
        (void)panvar::run_genotype_prepared(sample(name, r.which, r.threads, r.library), *loaded);
        require_same_science("prep_cold_" + r.tag, name, "loaded index, " + name);
        require(exists(path("prep_cold_" + r.tag + ".repair_rounds.tsv")) == exists(path(name + ".repair_rounds.tsv")) &&
                    (!exists(path(name + ".repair_rounds.tsv")) ||
                     trajectory(path("prep_cold_" + r.tag + ".repair_rounds.tsv")) == trajectory(path(name + ".repair_rounds.tsv"))),
                name + ": the adaptive trajectory must be the cold run's");
        const auto st = read_status(path(name + ".mosaic_status.tsv"));
        require(st.at("prepared_locus") == "1" && st.at("prepared_locus_origin") == "loaded" &&
                    st.at("prepared_locus_file_sha256") == report.file_checksum &&
                    st.at("prepared_locus_id") == std::to_string(loaded->id()),
                name + ": the status records a loaded index, its checksum and the runtime id");
        require(loaded->digest() == digest, name + ": the loaded locus is unchanged by an analysis");
    }
    {   // a changed graph: refused before any output
        panvar::GenotypeCommandOptions o = sample("index_changed_graph", 'A', 1, lib1);
        o.gfa_path = path("prep_other.gfa");
        std::filesystem::remove(path("index_changed_graph.mosaic_status.tsv"));
        bool refused = false;
        try {
            (void)panvar::run_genotype_prepared(o, *loaded);
        } catch (const panvar::PreparedLocusMismatch&) {
            refused = true;
        }
        require(refused && !exists(path("index_changed_graph.mosaic_status.tsv")), "a changed graph refuses the loaded index before output");
    }
    // I3-D: the same samples against an index carrying the persisted panel: loaded without a pre-count or a build, every
    // panel field equal to the prepared (freshly built) panel, and every analysis equal to its cold run
    const panvar::PreparedLocusPreparation with_panel = panvar::prepare_genotype_locus_record(sample("index_prepare_panel", 'A', 1, lib1), true);
    const std::string panel_file = path("prep_index_panel.pvpl");
    std::filesystem::remove(panel_file);
    panvar::write_prepared_locus(with_panel.record, panel_file);
    const std::uint64_t panel_builds = panvar::marker_panel_builds(), precounts = panvar::prepared_locus_panel_precounts();
    panvar::PreparedLocusLoadReport panel_report;
    const auto panel_loaded = panvar::load_prepared_locus(panel_file, panvar::PreparedLocusLoadLimits{}, &panel_report);
    require(panvar::marker_panel_builds() == panel_builds && panvar::prepared_locus_panel_precounts() == precounts &&
                panel_loaded->panel_source() == "loaded",
            "the panel route loads without a pre-count or a panel build");
    require(panvar::prepared_marker_panel_difference(panel_loaded->panel(), with_panel.locus->panel()).empty() &&
                panel_loaded->digest() == one.locus->digest(),
            "the persisted panel equals the freshly built one, field by field");
    require(loaded->panel_source() == "rebuilt", "the v1 index's panel is rebuilt");
    std::size_t panel_runs = 0;
    for (const Run& r : order) {
        const std::string name = "index_panel_run" + std::to_string(panel_runs++) + "_" + r.tag;
        (void)panvar::run_genotype_prepared(sample(name, r.which, r.threads, r.library), *panel_loaded);
        require_same_science("prep_cold_" + r.tag, name, "panel-loaded index, " + name);
        require(exists(path("prep_cold_" + r.tag + ".repair_rounds.tsv")) == exists(path(name + ".repair_rounds.tsv")) &&
                    (!exists(path(name + ".repair_rounds.tsv")) ||
                     trajectory(path("prep_cold_" + r.tag + ".repair_rounds.tsv")) == trajectory(path(name + ".repair_rounds.tsv"))),
                name + ": the adaptive trajectory must be the cold run's");
        const auto st = read_status(path(name + ".mosaic_status.tsv"));
        require(st.at("prepared_locus_origin") == "loaded" && st.at("prepared_locus_panel_source") == "loaded" &&
                    st.at("prepared_locus_file_sha256") == panel_report.file_checksum,
                name + ": the status records the panel-loaded index");
    }
    require(read_status(path("index_run0_A_t1.mosaic_status.tsv")).at("prepared_locus_panel_source") == "rebuilt",
            "the v1-loaded runs record a rebuilt panel");
    std::cout << "  prepared locus index: independent preparations byte-identical (" << bytes.size() << " B); write/read/"
              << "rewrite identical; " << runs << " analyses against the LOADED index (A-B-A, threads 1/4, two declared "
              << "libraries and auto) equal their cold runs; origin and file checksum in the status; a changed graph refused; "
              << panel_runs << " analyses against the PANEL-loaded index (no pre-count, no panel build; panel field-equal to a "
              << "fresh build) equal their cold runs\n";
}

// ---- S1-B: output_schema (opt-in sequence-first reporting) ---------------------------------------------------------
// Three scenarios (strict call, genotype-equivalent call, round-cap stop), each run with the schema unset, explicitly
// legacy and sequence_first_v1 in its own directory under the same prefix name. Required: (1) every legacy output is
// identical across the three arms (directory normalised; status-like files key by key, timing/memory values aside;
// the trajectory's scientific columns), so scores, rankings, group weights and the trajectory are unchanged; (2) the
// declared provenance only: explicit legacy adds exactly the settings row and the status row output_schema;
// sequence_first_v1 adds those plus four sequence status rows and the <prefix>.sequence_* files; (3) the exit codes:
// legacy 0/0/1, sequence_first_v1 0/0/3 (PROVISIONAL is not success). Then a watchdog stop under sequence_first_v1
// writes no sequence output, and the prepared-locus route writes the same sequence outputs as the cold route.
namespace sf {
std::string normalised(const std::string& file, const std::string& dir) {
    std::string text = slurp(file);
    for (std::size_t at = text.find(dir); at != std::string::npos; at = text.find(dir, at)) text.replace(at, dir.size(), "<D>");
    return text;
}
std::map<std::string, std::string> kv(const std::string& file) {
    std::map<std::string, std::string> out;
    std::istringstream in(slurp(file));
    for (std::string line; std::getline(in, line);) {
        const std::size_t tab = line.find('\t');
        if (tab != std::string::npos) out[line.substr(0, tab)] = line.substr(tab + 1);
    }
    return out;
}
bool volatile_key(const std::string& k) {
    static const std::regex pattern("seconds|footprint|rss|peak|_bytes|timing|worker|threads_|profile_sha256");
    return std::regex_search(k, pattern);
}
std::set<std::string> files(const std::string& dir) {
    std::set<std::string> out;
    for (const auto& e : std::filesystem::directory_iterator(dir)) out.insert(e.path().filename().string());
    return out;
}
}  // namespace sf

void test_sequence_first_schema() {
    namespace fs = std::filesystem;
    struct Scenario {
        std::string name, r1, r2, gfa;
        std::function<void(panvar::GenotypeCommandOptions&)> configure;
        int legacy_exit, v1_exit;
        std::string prediction, status;
    };
    const std::vector<Scenario> scenarios = {
        {"strict", path("phaseable_r1.fq"), path("phaseable_r2.fq"), "tiny_phaseable.gfa", [](auto&) {}, 0, 0,
         "PREDICTION", "RESOLVED"},
        {"ge", path("unphaseable_r1.fq"), path("unphaseable_r2.fq"), "tiny.gfa", [](auto&) {}, 0, 0, "PREDICTION",
         "GENOTYPE_RESOLVED_SEQUENCE_SET"},
        {"cap", path("phaseable_r1.fq"), path("phaseable_r2.fq"), "tiny_phaseable.gfa",
         [](auto& o) { o.max_mosaic_candidates = 0; o.repair_rounds = 1; }, 1, 3, "PROVISIONAL", "NOT_DETERMINED"},
    };
    const std::vector<std::string> sequence_files = {"call.sequence_locus.tsv", "call.sequence_blocks.tsv",
                                                     "call.representation_intervals.tsv", "call.sequence_alternatives.tsv",
                                                     "call.sequence_haplotypes.fa"};
    for (const Scenario& sc : scenarios) {
        std::map<std::string, std::string> dirs;
        std::map<std::string, int> exits;
        for (const std::string arm : {"unset", "legacy", "v1"}) {
            const std::string rel = "sf_" + sc.name + "_" + arm;
            fs::remove_all(path(rel));
            fs::create_directories(path(rel));
            panvar::GenotypeCommandOptions o = repair_options(rel + "/call", sc.r1, sc.r2, sc.gfa);
            sc.configure(o);
            if (arm != "unset") {
                o.output_schema = arm == "v1" ? "sequence_first_v1" : "legacy";
                o.setting_sources["output_schema"] = "test";
            }
            exits[arm] = panvar::execute_genotype_command(o, nullptr);
            dirs[arm] = path(rel);
        }
        require(exits["unset"] == sc.legacy_exit && exits["legacy"] == sc.legacy_exit,
                sc.name + ": the legacy exit code is unchanged (unset and explicit legacy)");
        require(exits["v1"] == sc.v1_exit, sc.name + ": the sequence_first_v1 exit code (" + std::to_string(exits["v1"]) + ")");
        const std::set<std::string> base = sf::files(dirs["unset"]);
        require(sf::files(dirs["legacy"]) == base, sc.name + ": explicit legacy writes exactly the same files");
        std::set<std::string> v1 = sf::files(dirs["v1"]);
        for (const std::string& f : sequence_files) {
            const bool expected = f != "call.sequence_haplotypes.fa" || sc.prediction != "NO_PREDICTION";
            require(v1.count(f) == (expected ? 1u : 0u) && base.count(f) == 0, sc.name + ": sequence file " + f);
            v1.erase(f);
        }
        require(v1 == base, sc.name + ": sequence_first_v1 adds only the sequence files");
        for (const std::string& f : base) {
            const std::string a = dirs["unset"] + "/" + f;
            for (const std::string arm : {"legacy", "v1"}) {
                const std::string b = dirs[arm] + "/" + f;
                const std::string what = sc.name + " " + arm + " " + f;
                if (f.size() > 11 && f.compare(f.size() - 11, 11, "_status.tsv") == 0) {
                    auto ka = sf::kv(a), kb = sf::kv(b);
                    if (f == "call.mosaic_status.tsv") {
                        require(kb.at("output_schema") == (arm == "v1" ? "sequence_first_v1" : "legacy"), what + ": schema row");
                        kb.erase("output_schema");
                        if (arm == "v1") {
                            require(kb.at("sequence_prediction") == sc.prediction && kb.at("sequence_status") == sc.status &&
                                        kb.at("sequence_exit_code") == std::to_string(sc.v1_exit),
                                    what + ": the sequence status rows");
                            require(kb.at("sequence_publication") == "complete", what + ": publication complete (S1-C row)");
                            for (const char* k : {"sequence_prediction", "sequence_status", "sequence_completion",
                                                  "sequence_exit_code", "sequence_publication"}) kb.erase(k);
                        }
                    }
                    std::set<std::string> keys_a, keys_b;
                    for (const auto& item : ka) keys_a.insert(item.first);
                    for (const auto& item : kb) keys_b.insert(item.first);
                    require(keys_a == keys_b, what + ": no other status row is added or removed");
                    const auto strip = [](std::string text, const std::string& dir) {
                        for (std::size_t at = text.find(dir); at != std::string::npos; at = text.find(dir, at)) {
                            text.replace(at, dir.size(), "<D>");
                        }
                        return text;
                    };
                    for (const auto& [k, v] : ka) {
                        if (sf::volatile_key(k)) continue;
                        require(strip(v, dirs["unset"]) == strip(kb.at(k), dirs[arm]),
                                what + ": status value " + k + " (" + v + " vs " + kb.at(k) + ")");
                    }
                } else if (f == "call.settings.tsv") {
                    const std::string extra = std::string("output_schema\t") + (arm == "v1" ? "sequence_first_v1" : "legacy") +
                                              "\ttest\n";
                    std::string text = slurp(b);
                    const std::size_t at = text.find(extra);
                    require(at != std::string::npos, what + ": exactly one output_schema settings row");
                    text.erase(at, extra.size());
                    require(text == slurp(a), what + ": every other settings row unchanged");
                } else if (f == "call.repair_rounds.tsv") {
                    require(rounds_science(a) == rounds_science(b), what + ": the adaptive trajectory is unchanged");
                } else if (f == "call.repair_memory.tsv" || f == "call.tie_branch_memory.tsv") {
                    continue;   // memory ledgers: measurements, not results
                } else {
                    require(sf::normalised(a, dirs["unset"]) == sf::normalised(b, dirs[arm]), what + ": byte-identical");
                }
            }
        }
        const auto locus = sf::kv(dirs["v1"] + "/call.sequence_locus.tsv");
        require(locus.at("prediction") == sc.prediction && locus.at("sequence_status") == sc.status &&
                    locus.at("structure_confidence") == "NA" && locus.at("structure_confidence_reason") == "UNCALIBRATED" &&
                    locus.at("projection_coverage") == "EXPLORED",
                sc.name + ": the sequence locus record");
        if (sc.name == "ge") {
            require(!base.count("call.mosaic_haplotypes.fa") &&
                        slurp(dirs["v1"] + "/call.sequence_haplotypes.fa").find("representative_of_set=yes") != std::string::npos,
                    "ge: the legacy call writes no FASTA; the sequence report retains a marked representative");
        }
        if (sc.name == "cap") {
            require(locus.at("completion") == "round_cap" && locus.at("decision_scope") == "NOT_DETERMINED" &&
                        locus.at("legacy_no_call_reason") == "round_cap",
                    "cap: PROVISIONAL from the round cap, no decision set claimed");
        }
    }
    {   // a watchdog stop under sequence_first_v1: no prediction is written
        fs::remove_all(path("sf_watchdog"));
        fs::create_directories(path("sf_watchdog"));
        panvar::GenotypeCommandOptions o =
            repair_options("sf_watchdog/call", path("phaseable_r1.fq"), path("phaseable_r2.fq"), "tiny_phaseable.gfa");
        o.max_mosaic_candidates = 0;
        o.test_interrupt_after_repair_rounds = 1;
        o.output_schema = "sequence_first_v1";
        o.setting_sources["output_schema"] = "test";
        std::cout.flush();
        const pid_t child = fork();
        require(child >= 0, "fork failed");
        if (child == 0) {
            try {
                std::_Exit(panvar::execute_genotype_command(o, nullptr));
            } catch (...) {
                std::_Exit(3);
            }
        }
        int wait_status = 0;
        waitpid(child, &wait_status, 0);
        require(WIFEXITED(wait_status) && WEXITSTATUS(wait_status) == 124, "watchdog: the time-limit exit");
        require(read_status(path("sf_watchdog/call.mosaic_status.tsv")).at("status") == "TIME_LIMIT" &&
                    !exists(path("sf_watchdog/call.sequence_locus.tsv")) &&
                    !exists(path("sf_watchdog/call.sequence_haplotypes.fa")),
                "watchdog: TIME_LIMIT and no sequence output (no prediction from a stopped computation)");
    }
    {   // the prepared-locus route: the same sequence outputs as the cold route
        std::map<std::string, std::string> d;
        for (const std::string route : {"cold", "prepared"}) {
            const std::string rel = "sf_route_" + route;
            fs::remove_all(path(rel));
            fs::create_directories(path(rel));
            panvar::GenotypeCommandOptions o =
                repair_options(rel + "/call", path("prep_A_r1.fq"), path("prep_A_r2.fq"), "prep_locus.gfa");
            o.max_mosaic_candidates = 8;
            o.output_schema = "sequence_first_v1";
            o.setting_sources["output_schema"] = "test";
            if (route == "cold") {
                require(panvar::execute_genotype_command(o, nullptr) == 0, "route: the cold run predicts");
            } else {
                const std::shared_ptr<const panvar::PreparedLocus> locus = panvar::prepare_genotype_locus(o);
                require(panvar::execute_genotype_command(o, locus.get()) == 0, "route: the prepared run predicts");
            }
            d[route] = path(rel);
        }
        for (const std::string& f : sequence_files) {
            require(sf::normalised(d["cold"] + "/" + f, d["cold"]) == sf::normalised(d["prepared"] + "/" + f, d["prepared"]),
                    "route: " + f + " identical on the prepared-locus route");
        }
    }
    std::cout << "  output_schema: legacy outputs identical across unset/legacy/sequence_first_v1 (strict, genotype-"
                 "equivalent, round cap); declared provenance only; exits 0/0/1 vs 0/0/3; watchdog writes no prediction; "
                 "prepared-locus route identical\n";
}

// ---- S1-C: the sequence-first reports are run-owned (same-prefix lifecycle) ------------------------------------------
// Each sequence first publishes a PREDICTION at one prefix, then runs the next case at the SAME prefix. Required: no
// earlier report survives beside a later outcome (legacy, provisional, watchdog stop, normal no-prediction); a prepared-
// locus compatibility refusal, raised before the run owns its prefix, leaves every existing output byte-identical; a
// report that cannot be fully written leaves no record claiming a prediction (exit 2, status says FAILED).
void test_sequence_first_lifecycle() {
    namespace fs = std::filesystem;
    const std::string rel = "sf_life", prefix = rel + "/call";
    const std::vector<std::string> owned = {".sequence_locus.tsv", ".sequence_blocks.tsv", ".sequence_haplotypes.fa",
                                            ".representation_intervals.tsv", ".sequence_alternatives.tsv"};
    const auto options = [&](bool v1) {
        panvar::GenotypeCommandOptions o =
            repair_options(prefix, path("phaseable_r1.fq"), path("phaseable_r2.fq"), "tiny_phaseable.gfa");
        if (v1) {
            o.output_schema = "sequence_first_v1";
            o.setting_sources["output_schema"] = "test";
        }
        return o;
    };
    const auto predict = [&]() {
        fs::remove_all(path(rel));
        fs::create_directories(path(rel));
        require(panvar::execute_genotype_command(options(true), nullptr) == 0, "lifecycle: the first run predicts");
        require(sf::kv(path(prefix + ".sequence_locus.tsv")).at("prediction") == "PREDICTION" &&
                    exists(path(prefix + ".sequence_haplotypes.fa")),
                "lifecycle: a published PREDICTION with its FASTA");
    };
    // prediction -> legacy
    predict();
    require(panvar::execute_genotype_command(options(false), nullptr) == 0, "lifecycle: the legacy run calls");
    for (const std::string& f : owned) {
        require(!exists(path(prefix + f)), "prediction -> legacy: no sequence report survives (" + f + ")");
    }
    // prediction -> provisional
    predict();
    {
        panvar::GenotypeCommandOptions o = options(true);
        o.max_mosaic_candidates = 0;
        o.repair_rounds = 1;
        require(panvar::execute_genotype_command(o, nullptr) == 3, "prediction -> provisional: exit 3");
        const auto locus = sf::kv(path(prefix + ".sequence_locus.tsv"));
        require(locus.at("prediction") == "PROVISIONAL" &&
                    slurp(path(prefix + ".sequence_haplotypes.fa")).find("prediction=PROVISIONAL") != std::string::npos &&
                    slurp(path(prefix + ".sequence_blocks.tsv")).find("IDENTIFIED") == std::string::npos &&
                    read_rows(path(prefix + ".representation_intervals.tsv")).empty(),
                "prediction -> provisional: only the provisional report, nothing from the earlier prediction");
    }
    // prediction -> normal no-prediction outcome (a declared neighbourhood bound refuses before any search)
    predict();
    {
        panvar::GenotypeCommandOptions o = options(true);
        o.repair_max_neighbours = 1;
        require(panvar::execute_genotype_command(o, nullptr) == 1, "prediction -> no prediction: exit 1");
        require(sf::kv(path(prefix + ".sequence_locus.tsv")).at("prediction") == "NO_PREDICTION" &&
                    !exists(path(prefix + ".sequence_haplotypes.fa")) &&
                    read_rows(path(prefix + ".sequence_blocks.tsv")).empty(),
                "prediction -> no prediction: a NO_PREDICTION record, no FASTA, no earlier block rows");
    }
    // prediction -> watchdog stop
    predict();
    {
        panvar::GenotypeCommandOptions o = options(true);
        o.max_mosaic_candidates = 0;
        o.test_interrupt_after_repair_rounds = 1;
        std::cout.flush();
        const pid_t child = fork();
        require(child >= 0, "fork failed");
        if (child == 0) {
            try {
                std::_Exit(panvar::execute_genotype_command(o, nullptr));
            } catch (...) {
                std::_Exit(3);
            }
        }
        int wait_status = 0;
        waitpid(child, &wait_status, 0);
        require(WIFEXITED(wait_status) && WEXITSTATUS(wait_status) == 124, "prediction -> watchdog: exit 124");
        require(read_status(path(prefix + ".mosaic_status.tsv")).at("status") == "TIME_LIMIT",
                "prediction -> watchdog: TIME_LIMIT status");
        for (const std::string& f : owned) {
            require(!exists(path(prefix + f)), "prediction -> watchdog: no earlier report beside TIME_LIMIT (" + f + ")");
        }
    }
    // prediction -> prepared-locus compatibility refusal: nothing at the prefix changes
    predict();
    {
        std::map<std::string, std::string> before;
        for (const auto& e : fs::directory_iterator(path(rel))) before[e.path().filename().string()] = slurp(e.path().string());
        panvar::GenotypeCommandOptions other = options(true);
        other.gfa_path = path("tiny.gfa");   // same reference path, other graph content
        const std::shared_ptr<const panvar::PreparedLocus> locus = panvar::prepare_genotype_locus(other);
        bool refused = false;
        try {
            (void)panvar::execute_genotype_command(options(true), locus.get());
        } catch (const panvar::PreparedLocusMismatch&) {
            refused = true;
        }
        require(refused, "prepared refusal: the mismatched prepared locus is refused");
        std::map<std::string, std::string> after;
        for (const auto& e : fs::directory_iterator(path(rel))) after[e.path().filename().string()] = slurp(e.path().string());
        require(after == before, "prepared refusal: every existing output, the sequence reports included, is unchanged");
    }
    // a report that cannot be fully written
    predict();
    {
        fs::create_directories(path(prefix + ".sequence_alternatives.tsv.writing/blocker"));
        const int code = panvar::execute_genotype_command(options(true), nullptr);
        require(code == 2, "publication failure: exit 2 (got " + std::to_string(code) + ")");
        const auto status = read_status(path(prefix + ".mosaic_status.tsv"));
        require(status.at("sequence_prediction") == "NO_PREDICTION" && status.at("sequence_publication") == "FAILED" &&
                    status.at("sequence_exit_code") == "2",
                "publication failure: the status claims no prediction and names the failure");
        for (const char* f : {".sequence_blocks.tsv", ".sequence_haplotypes.fa", ".representation_intervals.tsv",
                              ".sequence_alternatives.tsv"}) {
            require(!exists(path(prefix + f)), std::string("publication failure: no partial ") + f);
        }
        if (exists(path(prefix + ".sequence_locus.tsv"))) {
            require(sf::kv(path(prefix + ".sequence_locus.tsv")).at("prediction") == "NO_PREDICTION",
                    "publication failure: no locus record claims a prediction");
        }
        fs::remove_all(path(prefix + ".sequence_alternatives.tsv.writing"));   // the test's own blocker
    }
    std::cout << "  sequence-first lifecycle: prediction -> legacy / provisional / no prediction / watchdog leave no "
                 "earlier report; a prepared refusal changes nothing; a failed publication leaves no prediction claim\n";
}

// ---- S1-D: output_schema sequence_first_v2 (panel references, membership FASTA rule, diagnostic dosage) -------------
// (1) Membership on real caller output: every haplotype row's origin is PANEL exactly when an active panel template of
//     the same model spells its sequence (independently spelled here from a prepared locus), and its listed paths are
//     exactly those templates. (2) The dosage audit is diagnostic: audit and off give identical predictions, hashes,
//     tables and caller outputs (only the recorded dosage mode differs). (3) Explicit repair is refused under v2 (exit 2,
//     no outputs) and still runs under v1. (4) The prepared-locus route writes the same v2 outputs as the cold route.
void test_sequence_first_v2() {
    namespace fs = std::filesystem;
    const auto options = [&](const std::string& rel, const std::string& dosage) {
        fs::remove_all(path(rel));
        fs::create_directories(path(rel));
        panvar::GenotypeCommandOptions o =
            repair_options(rel + "/call", path("phaseable_r1.fq"), path("phaseable_r2.fq"), "tiny_phaseable.gfa");
        o.output_schema = "sequence_first_v2";
        o.setting_sources["output_schema"] = "test";
        o.dosage_mode = dosage;
        o.dosage_mode_explicit = true;
        return o;
    };
    const auto rows = [&](const std::string& file) {
        std::vector<std::map<std::string, std::string>> out;
        const auto table = read_rows(file);
        std::istringstream in(slurp(file));
        std::string head;
        std::getline(in, head);
        std::vector<std::string> cols;
        std::istringstream h(head);
        for (std::string c; std::getline(h, c, '\t');) cols.push_back(c);
        for (std::string line; std::getline(in, line);) {
            std::map<std::string, std::string> row;
            std::istringstream l(line);
            std::size_t k = 0;
            for (std::string c; std::getline(l, c, '\t'); ++k) row[cols[k]] = c;
            out.push_back(row);
        }
        return out;
    };
    // (1) membership, against the model's active panel spelled independently
    const panvar::GenotypeCommandOptions audit = options("sf2_audit", "audit");
    require(panvar::execute_genotype_command(audit, nullptr) == 0, "v2: the strict fixture predicts");
    const std::shared_ptr<const panvar::PreparedLocus> prepared = panvar::prepare_genotype_locus(audit);
    const panvar::PreparedGenotypeModel& model = prepared->model();
    std::map<std::string, std::vector<std::string>> panel_by_sha;
    for (std::size_t t = 0; t < model.panel().template_count(); ++t) {
        panel_by_sha[panvar::sha256_bytes(panvar::spell_mosaic(model, panvar::MosaicPath{model.panel().template_alleles[t]}).bases)]
            .push_back(model.panel().template_names[t]);
    }
    const auto haps = rows(path("sf2_audit/call.sequence_haplotypes.tsv"));
    require(haps.size() == 2, "v2: two haplotype rows");
    bool any_off = false;
    for (const auto& h : haps) {
        const auto found = panel_by_sha.find(h.at("sha256"));
        std::vector<std::string> want = found == panel_by_sha.end() ? std::vector<std::string>{} : found->second;
        std::sort(want.begin(), want.end());
        std::string joined;
        for (const auto& n : want) joined += (joined.empty() ? "" : ",") + n;
        require(h.at("origin") == (want.empty() ? "OFF_PANEL" : "PANEL") && h.at("panel_paths") == (want.empty() ? "NA" : joined),
                "v2: membership equals exact identity with the active panel (" + h.at("haplotype") + ")");
        any_off = any_off || want.empty();
    }
    const auto locus = sf::kv(path("sf2_audit/call.sequence_locus.tsv"));
    require(locus.at("fasta_emitted") == (any_off ? "yes" : "no") &&
                exists(path("sf2_audit/call.sequence_haplotypes.fa")) == any_off &&
                locus.at("fasta_reason") == (any_off ? "OFF_PANEL_PAIR" : "PANEL_REPRESENTED_PAIR") &&
                locus.at("dosage_mode") == "audit",
            "v2: the FASTA decision follows membership");
    // (2) audit vs off: nothing but the recorded dosage mode changes
    const panvar::GenotypeCommandOptions off = options("sf2_off", "off");
    require(panvar::execute_genotype_command(off, nullptr) == 0, "v2 off: predicts");
    for (const char* f : {"call.sequence_haplotypes.tsv", "call.sequence_blocks.tsv", "call.representation_intervals.tsv",
                          "call.sequence_alternatives.tsv", "call.mosaic_scores.tsv", "call.mosaic_candidates.tsv",
                          "call.equivalence_set.tsv"}) {
        require(sf::normalised(path("sf2_audit/") + f, path("sf2_audit")) == sf::normalised(path("sf2_off/") + f, path("sf2_off")),
                std::string("audit vs off: ") + f + " identical");
    }
    {   // the legacy block table: every column identical except the dosage filter pair, which records the mode (off:
        // DOSAGE_NOT_EVALUATED:dosage_mode_off; audit: its finding) -- a diagnostic annotation, never a genotype change
        const auto table = [&](const std::string& file) {
            std::vector<std::vector<std::string>> out;
            std::istringstream in(slurp(file));
            for (std::string line; std::getline(in, line);) {
                std::vector<std::string> cells;
                std::istringstream l(line);
                for (std::string c; std::getline(l, c, '\t');) cells.push_back(c);
                out.push_back(cells);
            }
            return out;
        };
        const auto a = table(path("sf2_audit/call.block_genotypes.tsv")), o = table(path("sf2_off/call.block_genotypes.tsv"));
        require(a.size() == o.size() && !a.empty() && a.front() == o.front(), "audit vs off: same block rows and columns");
        for (std::size_t r = 1; r < a.size(); ++r) {
            for (std::size_t c = 0; c < a.front().size(); ++c) {
                const std::string& col = a.front()[c];
                if (col == "dosage_filter" || col == "dosage_filter_reasons") continue;
                require(a[r][c] == o[r][c], "audit vs off: block column " + col + " unchanged by the audit");
            }
            require(o[r][std::find(o.front().begin(), o.front().end(), "dosage_filter_reasons") - o.front().begin()] ==
                        "DOSAGE_NOT_EVALUATED:dosage_mode_off",
                    "explicit off is recorded as not evaluated");
        }
    }
    auto la = sf::kv(path("sf2_audit/call.sequence_locus.tsv")), lo = sf::kv(path("sf2_off/call.sequence_locus.tsv"));
    require(la.at("dosage_mode") == "audit" && lo.at("dosage_mode") == "off", "audit vs off: recorded modes");
    la.erase("dosage_mode");
    lo.erase("dosage_mode");
    require(la == lo, "audit vs off: the locus record is otherwise identical (prediction, hashes, origins)");
    require(exists(path("sf2_audit/call.dosage_status.tsv")) &&
                read_status(path("sf2_audit/call.dosage_status.tsv")).at("dosage_audit") == "shadow",
            "audit: the dosage stage ran as an audit and wrote its own status");
    // (3) repair refused under v2 only
    {
        panvar::GenotypeCommandOptions rep = options("sf2_repair", "repair");
        require(panvar::execute_genotype_command(rep, nullptr) == 2 && !exists(path("sf2_repair/call.mosaic_status.tsv")),
                "v2: explicit repair is refused before any output");
        rep.output_schema = "sequence_first_v1";
        require(panvar::execute_genotype_command(rep, nullptr) == 0, "v1: repair still runs");
    }
    // (4) the prepared-locus route
    {
        const panvar::GenotypeCommandOptions route = options("sf2_prepared", "audit");
        require(panvar::execute_genotype_command(route, prepared.get()) == 0, "v2 prepared: predicts");
        for (const char* f : {"call.sequence_haplotypes.tsv", "call.sequence_locus.tsv", "call.sequence_blocks.tsv",
                              "call.representation_intervals.tsv", "call.sequence_alternatives.tsv"}) {
            require(sf::normalised(path("sf2_audit/") + f, path("sf2_audit")) ==
                        sf::normalised(path("sf2_prepared/") + f, path("sf2_prepared")),
                    std::string("v2 prepared route: ") + f + " identical to the cold route");
        }
    }
    std::cout << "  sequence_first_v2: membership = exact identity with the active panel (" << (any_off ? "off-panel" : "panel")
              << " pair here); audit leaves prediction and hashes unchanged; repair refused only under v2; prepared route identical\n";
}

// ---- S1-D.1: the v2 output contract across the whole prefix, and one dosage contract at every command entry --------
namespace s1d1 {
// Every FASTA a run left under its prefix (any file of the prefix's directory named <base>.*.fa or <base>.*.fasta).
std::set<std::string> fasta_inventory(const std::string& dir) {
    std::set<std::string> out;
    for (const std::string& f : sf::files(dir)) {
        const bool fa = f.size() > 3 && (f.compare(f.size() - 3, 3, ".fa") == 0 ||
                                         (f.size() > 6 && f.compare(f.size() - 6, 6, ".fasta") == 0));
        if (fa && f.rfind("call.", 0) == 0) out.insert(f);
    }
    return out;
}
// Every reference an output makes to a prefix file (call.<suffix>.tsv|.fa, also inside absolute paths and before a
// '#fragment') names a file that exists: nothing points to a suppressed or removed output.
std::vector<std::string> dangling_references(const std::string& dir) {
    // (a staging name in an error message, <file>.writing, is not a reference to an output)
    static const std::regex ref(R"(call\.[A-Za-z0-9_.]+?\.(tsv|fa|fasta)\b(?!\.writing))");
    std::vector<std::string> out;
    for (const std::string& f : sf::files(dir)) {
        if (f.rfind("call.", 0) != 0) continue;
        const std::string text = slurp(dir + "/" + f);
        for (auto it = std::sregex_iterator(text.begin(), text.end(), ref); it != std::sregex_iterator(); ++it) {
            const std::string name = it->str();
            if (!std::filesystem::exists(dir + "/" + name)) out.push_back(f + " -> " + name);
        }
    }
    return out;
}
std::map<std::string, std::string> locus_row(const std::string& file) {
    std::istringstream in(slurp(file));
    std::string head, row;
    std::getline(in, head);
    std::getline(in, row);
    std::vector<std::string> cols, cells;
    std::istringstream h(head), r(row);
    for (std::string c; std::getline(h, c, '\t');) cols.push_back(c);
    for (std::string c; std::getline(r, c, '\t');) cells.push_back(c);
    std::map<std::string, std::string> out;
    for (std::size_t k = 0; k < cols.size() && k < cells.size(); ++k) out[cols[k]] = cells[k];
    return out;
}
}  // namespace s1d1

// (1) Direct API, default-constructed dosage fields (Sol's S1-D review reproducer): audit runs, and its artefacts
//     exist. (2) Declared off: not run, no audit artefacts. (3) Repair -- declared, or set by value -- is refused with
//     exit 2 before the prefix is touched (a sentinel file survives). (4) The prepared-locus command route agrees with
//     (1)-(3). (5) A caller-only API run (no dosage stage) says so. (6) The FASTA inventory is the policy's: v2 = the
//     sequence FASTA only for an off-panel pair, nothing else; legacy and v1 keep their FASTA unchanged. (7) No output
//     references a missing prefix file.
void test_sequence_first_v2_contract() {
    namespace fs = std::filesystem;
    const auto fresh = [&](const std::string& rel, const std::string& schema) {
        fs::remove_all(path(rel));
        fs::create_directories(path(rel));
        panvar::GenotypeCommandOptions o =
            repair_options(rel + "/call", path("phaseable_r1.fq"), path("phaseable_r2.fq"), "tiny_phaseable.gfa");
        o.output_schema = schema;   // dosage_mode / dosage_mode_explicit / dosage_mode_source left as constructed
        return o;
    };
    const auto audit_ran = [&](const std::string& rel, const std::string& source, const std::string& what) {
        const auto locus = sf::kv(path(rel + "/call.sequence_locus.tsv"));
        require(locus.at("dosage_mode") == "audit" && locus.at("dosage_mode_source") == source,
                what + ": recorded dosage_mode audit, source " + source + " (got " + locus.at("dosage_mode") + ", " +
                    locus.at("dosage_mode_source") + ")");
        require(exists(path(rel + "/call.dosage_status.tsv")) &&
                    read_status(path(rel + "/call.dosage_status.tsv")).at("dosage_audit_status") == "COMPLETE" &&
                    exists(path(rel + "/call.dosage_blocks.tsv")),
                what + ": the audit ran (dosage_status COMPLETE, dosage_blocks written)");
        const auto row = s1d1::locus_row(path(rel + "/call.locus_genotype.tsv"));
        require(row.at("dosage_mode") == "audit" && row.at("dosage_audit_status") == "COMPLETE",
                what + ": the locus table records the audit actually run");
    };
    const auto v2_inventory = [&](const std::string& rel, const std::string& what) {
        const auto locus = sf::kv(path(rel + "/call.sequence_locus.tsv"));
        const bool off_panel = locus.at("fasta_emitted") == "yes";
        const std::set<std::string> want =
            off_panel ? std::set<std::string>{"call.sequence_haplotypes.fa"} : std::set<std::string>{};
        std::set<std::string> got = s1d1::fasta_inventory(path(rel));
        for (auto it = got.begin(); it != got.end();) {   // a non-empty dosage FASTA is an audit diagnostic (kept)
            const bool diagnostic = it->rfind("call.dosage_", 0) == 0 && fs::file_size(path(rel + "/" + *it)) > 0;
            it = diagnostic ? got.erase(it) : std::next(it);
        }
        require(got == want, what + ": v2 FASTA inventory is the policy's (" + std::to_string(got.size()) + " files)");
        const auto row = s1d1::locus_row(path(rel + "/call.locus_genotype.tsv"));
        require(row.at("main_fasta") == (off_panel ? "call.sequence_haplotypes.fa" : "NA") && row.at("original_fasta") == "NA" &&
                    row.at("hap1_id") == "NA" && row.at("hap2_id") == "NA",
                what + ": the locus table links only the authoritative FASTA");
        const auto dangling = s1d1::dangling_references(path(rel));
        require(dangling.empty(), what + ": no dangling output reference" + (dangling.empty() ? "" : " (" + dangling.front() + ")"));
    };
    // (1) direct API, dosage untouched
    const panvar::GenotypeCommandOptions api = fresh("sf21_api", "sequence_first_v2");
    require(api.dosage_mode == "off" && !api.dosage_mode_explicit, "the reproducer's options: dosage fields as constructed");
    require(panvar::execute_genotype_command(api, nullptr) == 0, "v2 API default: predicts");
    audit_ran("sf21_api", "contract_default_sequence_first_v2", "v2 API default");
    v2_inventory("sf21_api", "v2 API default");
    // (2) declared off
    {
        panvar::GenotypeCommandOptions off = fresh("sf21_off", "sequence_first_v2");
        off.dosage_mode_explicit = true;
        require(panvar::execute_genotype_command(off, nullptr) == 0, "v2 API off: predicts");
        const auto locus = sf::kv(path("sf21_off/call.sequence_locus.tsv"));
        require(locus.at("dosage_mode") == "off" && locus.at("dosage_mode_source") == "caller_options" &&
                    !exists(path("sf21_off/call.dosage_status.tsv")) && !exists(path("sf21_off/call.dosage_blocks.tsv")) &&
                    s1d1::locus_row(path("sf21_off/call.locus_genotype.tsv")).at("dosage_mode") == "off",
                "v2 API declared off: off recorded and no audit artefact");
        v2_inventory("sf21_off", "v2 API off");
    }
    // (3) repair refused early, declared or by value; the prefix is untouched
    const auto refused_untouched = [&](panvar::GenotypeCommandOptions o, const panvar::PreparedLocus* locus, const std::string& what) {
        const std::string sentinel = o.out_prefix + ".mosaic_status.tsv";
        { std::ofstream(sentinel) << "sentinel\n"; }
        require(panvar::execute_genotype_command(o, locus) == 2, what + ": exit 2");
        require(slurp(sentinel) == "sentinel\n" && sf::files(fs::path(sentinel).parent_path().string()).size() == 1,
                what + ": refused before the prefix was touched");
    };
    {
        panvar::GenotypeCommandOptions rep = fresh("sf21_repair_declared", "sequence_first_v2");
        rep.dosage_mode = "repair";
        rep.dosage_mode_explicit = true;
        refused_untouched(rep, nullptr, "v2 API declared repair");
        panvar::GenotypeCommandOptions by_value = fresh("sf21_repair_value", "sequence_first_v2");
        by_value.dosage_mode = "repair";   // not flagged explicit: still a declaration, never silently audit
        refused_untouched(by_value, nullptr, "v2 API repair by value");
    }
    // (4) the prepared-locus command route
    {
        const panvar::GenotypeCommandOptions route = fresh("sf21_prepared", "sequence_first_v2");
        const std::shared_ptr<const panvar::PreparedLocus> prepared = panvar::prepare_genotype_locus(route);
        require(panvar::execute_genotype_command(route, prepared.get()) == 0, "v2 prepared default: predicts");
        audit_ran("sf21_prepared", "contract_default_sequence_first_v2", "v2 prepared default");
        v2_inventory("sf21_prepared", "v2 prepared default");
        for (const char* f : {"call.sequence_haplotypes.tsv", "call.sequence_locus.tsv", "call.sequence_blocks.tsv",
                              "call.dosage_blocks.tsv", "call.locus_genotype.tsv"}) {
            require(sf::normalised(path("sf21_api/") + f, path("sf21_api")) == sf::normalised(path("sf21_prepared/") + f, path("sf21_prepared")),
                    std::string("v2 prepared default: ") + f + " identical to the direct API route");
        }
        require(s1d1::fasta_inventory(path("sf21_api")) == s1d1::fasta_inventory(path("sf21_prepared")),
                "v2 prepared default: the same FASTA inventory");
        panvar::GenotypeCommandOptions rep = fresh("sf21_prepared_repair", "sequence_first_v2");
        rep.dosage_mode = "repair";
        rep.dosage_mode_explicit = true;
        refused_untouched(rep, prepared.get(), "v2 prepared repair");
    }
    // (5) a caller-only API run has no dosage stage, and its v2 report says so
    {
        const panvar::GenotypeCommandOptions only = fresh("sf21_caller_only", "sequence_first_v2");
        (void)panvar::run_genotype(only);
        const auto locus = sf::kv(path("sf21_caller_only/call.sequence_locus.tsv"));
        require(locus.at("dosage_mode") == "not_run" && locus.at("dosage_mode_source") == "caller_only_api_without_dosage_stage" &&
                    !exists(path("sf21_caller_only/call.dosage_status.tsv")),
                "caller-only API: the report names no dosage mode it did not run");
    }
    // (6) legacy and v1 keep their FASTA exactly (non-regression), and their links
    for (const std::string schema : {"legacy", "sequence_first_v1"}) {
        const std::string rel = "sf21_" + schema;
        panvar::GenotypeCommandOptions o = fresh(rel, schema);
        o.dosage_mode = "audit";
        o.dosage_mode_explicit = true;
        require(panvar::execute_genotype_command(o, nullptr) == 0, schema + ": calls");
        std::set<std::string> want = {"call.mosaic_haplotypes.fa", "call.dosage_original.fa", "call.dosage_repaired.fa",
                                      "call.dosage_repaired_block_alleles.fa"};
        if (schema == "sequence_first_v1") want.insert("call.sequence_haplotypes.fa");
        require(s1d1::fasta_inventory(path(rel)) == want, schema + ": FASTA inventory unchanged (legacy FASTA and every dosage FASTA kept)");
        const auto row = s1d1::locus_row(path(rel + "/call.locus_genotype.tsv"));
        require(row.at("main_fasta") == "call.mosaic_haplotypes.fa" && row.at("hap1_id") != "NA",
                schema + ": the locus table links the legacy FASTA as before");
        require(s1d1::dangling_references(path(rel)).empty(), schema + ": no dangling output reference");
        require(slurp(path("sf21_api/call.mosaic_scores.tsv")).empty() == false &&
                    sf::normalised(path("sf21_api/call.mosaic_scores.tsv"), path("sf21_api")) ==
                        sf::normalised(path(rel + "/call.mosaic_scores.tsv"), path(rel)),
                schema + " vs v2: the caller's scores are identical");
    }
    // (8) prefix reuse: an off-panel prediction (FASTA written), then a v2 publication failure at the same prefix: exit
    //     2 as before, no FASTA of any kind is left, and nothing references a missing file
    {
        const panvar::GenotypeCommandOptions first = fresh("sf21_failure", "sequence_first_v2");
        require(panvar::execute_genotype_command(first, nullptr) == 0 &&
                    s1d1::fasta_inventory(path("sf21_failure")) == std::set<std::string>{"call.sequence_haplotypes.fa"},
                "v2 prefix reuse: the off-panel prediction writes its FASTA");
        fs::create_directories(path("sf21_failure/call.sequence_alternatives.tsv.writing/blocker"));
        const int code = panvar::execute_genotype_command(first, nullptr);
        require(code == 2 && read_status(path("sf21_failure/call.mosaic_status.tsv")).at("sequence_publication") == "FAILED",
                "v2 publication failure: exit 2 and the status names the failure");
        require(s1d1::fasta_inventory(path("sf21_failure")).empty(), "v2 publication failure: no FASTA is left");
        const auto row = s1d1::locus_row(path("sf21_failure/call.locus_genotype.tsv"));
        require(row.at("main_fasta") == "NA", "v2 publication failure: the locus table links no FASTA");
        require(s1d1::dangling_references(path("sf21_failure")).empty(), "v2 publication failure: no dangling reference");
        fs::remove_all(path("sf21_failure/call.sequence_alternatives.tsv.writing"));   // the test's own blocker
    }
    std::cout << "  sequence_first_v2 contract: API/prepared default audit (artefacts present), declared off, repair "
                 "refused before the prefix; caller-only API says not_run; v2 FASTA inventory = policy, no dangling "
                 "reference (also after a publication failure on a reused prefix); legacy/v1 FASTA unchanged\n";
}

// ---- S2-A: early evaluated incumbent and safe budget-stop reporting (early_incumbent experimental_v1) ----------------
namespace s2a {
// The command in a child (a watchdog or test-hook stop ends the whole process); 99 = an escaped exception.
int execute_in_child(const panvar::GenotypeCommandOptions& options) {
    std::cout.flush();
    const pid_t child = fork();
    require(child >= 0, "fork failed");
    if (child == 0) {
        int code = 99;
        try {
            code = panvar::execute_genotype_command(options, nullptr);
        } catch (...) {
            code = 99;
        }
        std::_Exit(code);
    }
    int wait_status = 0;
    waitpid(child, &wait_status, 0);
    return WIFEXITED(wait_status) ? WEXITSTATUS(wait_status) : -1;
}
std::vector<panvar::ReadPair> fastq_pairs(const std::string& r1, const std::string& r2) {
    const auto load = [](const std::string& file) {
        std::vector<panvar::SequencingRead> out;
        std::ifstream in(file);
        for (std::string head, bases, plus, quals; std::getline(in, head) && std::getline(in, bases) &&
                                                   std::getline(in, plus) && std::getline(in, quals);) {
            panvar::SequencingRead r;
            r.name = head.substr(1, head.find_first_of(" \t") == std::string::npos ? std::string::npos : head.find_first_of(" \t") - 1);
            for (char& c : bases) c = static_cast<char>(std::toupper(c));
            r.bases = bases;
            r.qualities = quals;
            out.push_back(r);
        }
        return out;
    };
    const auto a = load(r1), b = load(r2);
    std::vector<panvar::ReadPair> out;
    for (std::size_t i = 0; i < a.size() && i < b.size(); ++i) out.push_back({a[i], b[i]});
    return out;
}
std::vector<std::size_t> labels(const std::string& text) {
    std::vector<std::size_t> out;
    std::istringstream in(text);
    for (std::string x; std::getline(in, x, ',');) out.push_back(std::stoul(x));
    return out;
}
bool bits(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }
}  // namespace s2a

void test_early_incumbent() {
    namespace fs = std::filesystem;
    const auto fresh = [&](const std::string& rel, bool on) {
        fs::remove_all(path(rel));
        fs::create_directories(path(rel));
        panvar::GenotypeCommandOptions o =
            repair_options(rel + "/call", path("phaseable_r1.fq"), path("phaseable_r2.fq"), "tiny_phaseable.gfa");
        o.output_schema = "sequence_first_v2";
        o.setting_sources["output_schema"] = "test";
        if (on) {
            o.early_incumbent = "experimental_v1";
            o.setting_sources["early_incumbent"] = "test";
        }
        return o;
    };
    static const std::set<std::string> added_status = {"early_incumbent", "checkpoint_commits", "checkpoint_last_stage",
                                                       "checkpoint_evaluated_candidates",
                                                       "checkpoint_seconds_first_evaluated_prediction", "checkpoint_failure"};
    const auto science_status = [&](const std::string& file) {
        std::map<std::string, std::string> out;
        for (const auto& [k, v] : sf::kv(file)) {
            if (sf::volatile_key(k) || k.rfind("early_incumbent", 0) == 0 || added_status.count(k) != 0) continue;
            out[k] = v;
        }
        return out;
    };
    // (1) completed run, off vs on: identical legacy and sequence outputs; only the declared status/settings rows added
    require(panvar::execute_genotype_command(fresh("ei_off", false), nullptr) == 0, "S2-A off: predicts");
    const panvar::GenotypeCommandOptions on = fresh("ei_on", true);
    require(panvar::execute_genotype_command(on, nullptr) == 0, "S2-A on: predicts");
    require(sf::files(path("ei_off")) == sf::files(path("ei_on")), "S2-A off/on: the same files");
    for (const std::string& f : sf::files(path("ei_off"))) {
        if (f == "call.mosaic_status.tsv" || f == "call.settings.tsv" || f == "call.repair_memory.tsv") continue;
        if (f == "call.repair_rounds.tsv") {
            require(rounds_science(path("ei_off/" + f)) == rounds_science(path("ei_on/" + f)), "S2-A off/on: trajectory");
            continue;
        }
        if (f == "call.dosage_status.tsv") continue;   // timing keys only (compared below without them)
        require(sf::normalised(path("ei_off/" + f), path("ei_off")) == sf::normalised(path("ei_on/" + f), path("ei_on")),
                "S2-A off/on: " + f + " identical (completed-run mathematics unchanged)");
    }
    require(science_status(path("ei_off/call.mosaic_status.tsv")) == science_status(path("ei_on/call.mosaic_status.tsv")),
            "S2-A off/on: the status differs only by the declared early-incumbent rows");
    {
        const std::string a = slurp(path("ei_off/call.settings.tsv")), b = slurp(path("ei_on/call.settings.tsv"));
        require(b == a + "early_incumbent\texperimental_v1\ttest\n" || sf::normalised(path("ei_on/call.settings.tsv"), path("ei_on")).find("early_incumbent\texperimental_v1") != std::string::npos,
                "S2-A on: the setting is recorded");
    }
    const auto status = read_status(path("ei_on/call.mosaic_status.tsv"));
    require(status.at("early_incumbent_status") == "evaluated" && status.at("early_incumbent_union_check") == "reproduced" &&
                std::stoul(status.at("checkpoint_commits")) >= 2 && status.at("checkpoint_failure") == "none",
            "S2-A on: early candidate evaluated, reproduced bit for bit by the union pass, checkpoints committed");
    // (2) the early candidate's components, recomputed independently with the ordinary scorer
    {
        const std::shared_ptr<const panvar::PreparedLocus> prepared = panvar::prepare_genotype_locus(on);
        const panvar::PreparedGenotypeModel& model = prepared->model();
        const std::vector<std::size_t> a = s2a::labels(status.at("early_incumbent_hap1"));
        const std::vector<std::size_t> b = s2a::labels(status.at("early_incumbent_hap2"));
        const std::string sa = panvar::spell_mosaic(model, panvar::MosaicPath{a}).bases;
        const std::string sb = panvar::spell_mosaic(model, panvar::MosaicPath{b}).bases;
        const panvar::InternedHaplotypes ih = panvar::intern_haplotypes({panvar::MosaicPath{a}, panvar::MosaicPath{b}}, {sa, sb});
        panvar::FragmentScoringOptions so;
        so.threads = 1;
        so.mismatch_band = on.mismatch_band;
        so.minimum_piece_bp = on.minimum_piece_bp;
        so.max_piece_hits = on.max_piece_hits;
        so.max_proposed_starts_per_haplotype = on.max_proposed_starts_per_haplotype;
        so.max_physical_origins = on.max_physical_origins;
        so.max_state_classes = on.max_state_classes;
        so.max_mass_cells = on.max_mass_cells;
        std::optional<double> background;
        const panvar::LibraryModel library = panvar::read_library(on.library_path, background);
        const panvar::FragmentScoringResult r = panvar::score_diplotypes(
            ih, {{std::min(ih.index_of[0], ih.index_of[1]), std::max(ih.index_of[0], ih.index_of[1])}},
            s2a::fastq_pairs(on.reads1_path, on.reads2_path), library, so);
        const panvar::DiplotypeScore& d = r.scores.front();
        const std::vector<double> sw(model.blocks().size() - 1, on.switch_probability);
        const double prior = panvar::log_haploid_sequence_prior(model, sa, sw) + panvar::log_haploid_sequence_prior(model, sb, sw) +
                             (sa != sb ? std::log(2.0) : 0.0);
        require(s2a::bits(d.log_likelihood, std::stod(status.at("early_incumbent_log_likelihood"))) &&
                    s2a::bits(d.fragment_term, std::stod(status.at("early_incumbent_fragment_term"))) &&
                    s2a::bits(d.exposure_term, std::stod(status.at("early_incumbent_exposure_term"))) &&
                    s2a::bits(d.absolute_term_sum, std::stod(status.at("early_incumbent_absolute_term_sum"))) &&
                    d.fragments_with_mass == std::stoul(status.at("early_incumbent_fragments_with_mass")) &&
                    s2a::bits(prior, std::stod(status.at("early_incumbent_sequence_prior"))),
                "S2-A: the early candidate's components equal an independent ordinary-scorer evaluation, bit for bit");
    }
    // (3) deterministic interruptions
    const auto interrupted = [&](const std::string& rel, const std::string& point, std::size_t threads = 2) {
        panvar::GenotypeCommandOptions o = fresh(rel, true);
        o.test_interrupt_point = point;
        o.threads = threads;
        return s2a::execute_in_child(o);
    };
    {
        const int code = interrupted("ei_before", "before_first_evaluation");
        const auto locus = sf::kv(path("ei_before/call.sequence_locus.tsv"));
        const auto st = read_status(path("ei_before/call.mosaic_status.tsv"));
        require(code == 124 && locus.at("prediction") == "NO_PREDICTION" &&
                    locus.at("completion") == "budget_stop_before_first_evaluated_candidate" &&
                    !exists(path("ei_before/call.sequence_haplotypes.tsv")) && st.at("status") == "TIME_LIMIT" &&
                    st.at("sequence_exit_code") == "1" && st.at("checkpoint_committed") == "no",
                "S2-A stop before the first evaluation: typed NO_PREDICTION, exit 124 (got " + std::to_string(code) + ")");
    }
    std::map<std::string, std::string> early_locus;
    for (const std::string point : {"after_checkpoint", "union_scoring", "adaptive_refinement"}) {
        const std::string rel = "ei_" + point;
        const int code = interrupted(rel, point);
        const auto locus = sf::kv(path(rel + "/call.sequence_locus.tsv"));
        const auto st = read_status(path(rel + "/call.mosaic_status.tsv"));
        const bool early_stage = point != "adaptive_refinement";
        require(code == 3 && locus.at("prediction") == "PROVISIONAL" && locus.at("completion") == "budget_stop" &&
                    locus.at("checkpoint_search_status") == "budget_stop:TIME_LIMIT" &&
                    locus.at("legacy_run_status") == "TIME_LIMIT" &&
                    (early_stage ? locus.at("checkpoint_stage") == "early_panel_incumbent"
                                 : locus.at("checkpoint_stage") != "early_panel_incumbent") &&
                    st.at("status") == "TIME_LIMIT" && st.at("sequence_prediction") == "PROVISIONAL" &&
                    st.at("sequence_exit_code") == "3" && no_ranking_or_call(rel + "/call") &&
                    exists(path(rel + "/call.sequence_haplotypes.tsv")) &&
                    locus.at("checkpoint_semantics").find("search_incomplete") != std::string::npos,
                "S2-A stop " + point + ": the committed checkpoint is kept as PROVISIONAL (exit " + std::to_string(code) + ")");
        const auto rows = read_rows(path(rel + "/call.sequence_haplotypes.tsv"));
        require(rows.size() == 2, "S2-A stop " + point + ": two haplotype rows with references or the v2 FASTA rule");
        require((locus.at("fasta_emitted") == "yes") == exists(path(rel + "/call.sequence_haplotypes.fa")),
                "S2-A stop " + point + ": the v2 FASTA policy holds for the checkpoint");
        if (point == "after_checkpoint") early_locus = locus;
        if (!early_stage) {
            require(std::stod(locus.at("checkpoint_total_score")) >= std::stod(early_locus.at("checkpoint_total_score")),
                    "S2-A: a replacement is never worse than the early incumbent");
        }
    }
    require(early_locus.at("checkpoint_log_likelihood") == status.at("early_incumbent_log_likelihood") &&
                early_locus.at("checkpoint_sequence_prior") == status.at("early_incumbent_sequence_prior") &&
                early_locus.at("checkpoint_evaluated_candidates") == "1" && early_locus.at("checkpoint_tied_groups") == "NA",
            "S2-A: the early checkpoint carries the completed run's early components exactly, one evaluated candidate, no ties claimed");
    // (4) 1 vs 4 threads: the same early checkpoint and the same completed outputs
    {
        require(interrupted("ei_t1", "after_checkpoint", 1) == 3 && interrupted("ei_t4", "after_checkpoint", 4) == 3,
                "S2-A threads: both stop with a checkpoint");
        auto a = sf::kv(path("ei_t1/call.sequence_locus.tsv")), b = sf::kv(path("ei_t4/call.sequence_locus.tsv"));
        for (auto* m : {&a, &b}) {
            m->erase("checkpoint_seconds_first_evaluated_prediction");
            m->erase("checkpoint_seconds_committed");
            m->erase("completion_note");
        }
        require(a == b, "S2-A threads: the early checkpoint is identical at 1 and 4 threads");
        panvar::GenotypeCommandOptions t4 = fresh("ei_on4", true);
        t4.threads = 4;
        require(panvar::execute_genotype_command(t4, nullptr) == 0, "S2-A on, 4 threads: predicts");
        for (const char* f : {"call.sequence_locus.tsv", "call.sequence_haplotypes.tsv", "call.block_genotypes.tsv",
                              "call.mosaic_scores.tsv"}) {
            require(sf::normalised(path("ei_on/") + f, path("ei_on")) == sf::normalised(path("ei_on4/") + f, path("ei_on4")),
                    std::string("S2-A threads: ") + f + " identical at 2 and 4 threads");
        }
    }
    // (5) prefix reuse: a completed run, then a stop after the checkpoint, then a completed run, one prefix
    {
        const panvar::GenotypeCommandOptions first = fresh("ei_reuse", true);
        require(panvar::execute_genotype_command(first, nullptr) == 0, "S2-A reuse: the first run predicts");
        panvar::GenotypeCommandOptions stop = first;
        stop.test_interrupt_point = "after_checkpoint";
        require(s2a::execute_in_child(stop) == 3, "S2-A reuse: the stopped run keeps its checkpoint");
        require(no_ranking_or_call("ei_reuse/call") && !exists(path("ei_reuse/call.locus_genotype.tsv")) &&
                    sf::kv(path("ei_reuse/call.sequence_locus.tsv")).at("checkpoint_stage") == "early_panel_incumbent",
                "S2-A reuse: no earlier final output survives beside the checkpoint");
        require(panvar::execute_genotype_command(first, nullptr) == 0 &&
                    sf::normalised(path("ei_reuse/call.sequence_locus.tsv"), path("ei_reuse")) ==
                        sf::normalised(path("ei_on/call.sequence_locus.tsv"), path("ei_on")),
                "S2-A reuse: a later completed run replaces the checkpoint with the final report");
    }
    // (6) a checkpoint that cannot be published: an explicit failure (exit 2), never a fallback
    {
        panvar::GenotypeCommandOptions o = fresh("ei_publish", true);
        fs::create_directories(path("ei_publish/call.sequence_alternatives.tsv.writing/blocker"));
        const int code = panvar::execute_genotype_command(o, nullptr);
        const auto st = read_status(path("ei_publish/call.mosaic_status.tsv"));
        require(code == 2 && st.at("sequence_publication") == "FAILED" &&
                    st.at("checkpoint_failure").find("checkpoint_publication_failed") == 0 &&
                    !exists(path("ei_publish/call.sequence_haplotypes.tsv")) && no_ranking_or_call("ei_publish/call"),
                "S2-A: a failed checkpoint publication ends the run as an explicit failure (exit " + std::to_string(code) + ")");
        fs::remove_all(path("ei_publish/call.sequence_alternatives.tsv.writing"));
    }
    // (7) the contract: experimental_v1 needs sequence_first_v2 and is refused before the prefix is touched
    {
        panvar::GenotypeCommandOptions o = fresh("ei_refused", true);
        o.output_schema = "sequence_first_v1";
        bool refused = false;
        try {
            (void)panvar::execute_genotype_command(o, nullptr);
        } catch (const std::invalid_argument&) {
            refused = true;
        }
        require(refused && sf::files(path("ei_refused")).empty(), "S2-A: v1 + early_incumbent refused before any output");
    }
    std::cout << "  early incumbent: off/on completed outputs identical; early candidate reproduced by the union and by an "
                 "independent scorer; stops before/after checkpoint, during union scoring and refinement typed correctly; "
                 "1/4 threads, prefix reuse, publication failure, contract refusal\n";
}

// ---- S2-B: maintained regressions for S2-A's declared-cap stops (after Sol's review reproducer) and the explicit
// failure on a scoring inconsistency. Each cap stop after a checkpoint publishes the best fully evaluated candidate as
// PROVISIONAL (exit 3, completion = the cap); the inconsistency is never a provisional fallback.
void test_early_incumbent_caps() {
    namespace fs = std::filesystem;
    const auto fresh = [&](const std::string& rel) {
        fs::remove_all(path(rel));
        fs::create_directories(path(rel));
        panvar::GenotypeCommandOptions o =
            repair_options(rel + "/call", path("phaseable_r1.fq"), path("phaseable_r2.fq"), "tiny_phaseable.gfa");
        o.output_schema = "sequence_first_v2";
        o.setting_sources["output_schema"] = "test";
        o.early_incumbent = "experimental_v1";
        o.setting_sources["early_incumbent"] = "test";
        return o;
    };
    // the reference completed run: the early pass's and the union's physical origins (for the union cap)
    const panvar::GenotypeCommandOptions reference = fresh("eic_reference");
    require(panvar::execute_genotype_command(reference, nullptr) == 0, "S2-B caps: the reference run predicts");
    const auto ref_status = read_status(path("eic_reference/call.mosaic_status.tsv"));
    const std::uint64_t early_origins = std::stoull(ref_status.at("early_incumbent_physical_origins"));
    const std::uint64_t union_origins = std::stoull(ref_status.at("physical_origins"));
    require(early_origins > 0 && union_origins > early_origins, "S2-B caps: the union pass places more origins than the early pass");
    struct Case {
        std::string name, completion, stage;
        std::function<void(panvar::GenotypeCommandOptions&)> configure;
    };
    const std::vector<Case> cases = {
        {"proposal", "proposal_cap", "early_panel_incumbent",
         [](panvar::GenotypeCommandOptions& o) { o.max_mosaic_candidates = 128; o.max_beam_extensions = 1; }},
        {"union", "union_scoring_cap", "early_panel_incumbent",
         [&](panvar::GenotypeCommandOptions& o) { o.max_physical_origins = early_origins; }},
        {"cumulative", "cumulative_cap", "union_best",
         [](panvar::GenotypeCommandOptions& o) { o.repair_max_physical_origins = 1; }},
        {"neighbours", "neighbourhood_bound", "union_best",
         [](panvar::GenotypeCommandOptions& o) { o.repair_max_neighbours = 1; }},
    };
    for (const Case& c : cases) {
        const std::string rel = "eic_" + c.name;
        panvar::GenotypeCommandOptions o = fresh(rel);
        c.configure(o);
        const int code = panvar::execute_genotype_command(o, nullptr);
        const auto locus = sf::kv(path(rel + "/call.sequence_locus.tsv"));
        const auto st = read_status(path(rel + "/call.mosaic_status.tsv"));
        require(code == 3 && locus.at("prediction") == "PROVISIONAL" && locus.at("completion") == c.completion &&
                    locus.at("checkpoint_stage") == c.stage && locus.at("checkpoint_search_status") == c.completion &&
                    st.at("sequence_prediction") == "PROVISIONAL" && st.at("sequence_exit_code") == "3" &&
                    st.at("checkpoint_failure") == "none" && exists(path(rel + "/call.sequence_haplotypes.tsv")),
                "S2-B caps: " + c.name + " stop keeps the best fully evaluated candidate as PROVISIONAL (exit " +
                    std::to_string(code) + ", completion " + (locus.count("completion") ? locus.at("completion") : "NA") +
                    ", stage " + (locus.count("checkpoint_stage") ? locus.at("checkpoint_stage") : "NA") + ")");
        if (c.stage == "early_panel_incumbent") {
            require(locus.at("checkpoint_log_likelihood") == ref_status.at("early_incumbent_log_likelihood"),
                    "S2-B caps: " + c.name + " keeps exactly the early candidate's evaluation");
        }
    }
    // a scoring inconsistency (the early candidate moved one ulp) is an explicit failure, never a provisional fallback
    {
        panvar::GenotypeCommandOptions o = fresh("eic_inconsistent");
        o.test_perturb_early_incumbent_score = true;
        const int code = panvar::execute_genotype_command(o, nullptr);
        const auto locus = sf::kv(path("eic_inconsistent/call.sequence_locus.tsv"));
        const auto st = read_status(path("eic_inconsistent/call.mosaic_status.tsv"));
        require(code != 0 && code != 3 && locus.at("prediction") == "NO_PREDICTION" &&
                    locus.at("legacy_no_call_reason") == "early_incumbent_not_reproduced" &&
                    st.at("early_incumbent_union_check") == "MISMATCH" &&
                    st.at("checkpoint_failure").rfind("scoring_inconsistency", 0) == 0 &&
                    read_rows(path("eic_inconsistent/call.sequence_haplotypes.tsv")).empty() &&
                    no_ranking_or_call("eic_inconsistent/call"),
                "S2-B: a scoring inconsistency ends the run explicitly with no prediction (exit " + std::to_string(code) + ")");
    }
    std::cout << "  early incumbent caps: proposal, union, cumulative and neighbourhood stops keep the best evaluated "
                 "candidate as PROVISIONAL; a scoring inconsistency is an explicit no-prediction failure\n";
}

// ---- S2-D: opt-in progressive scheduling of the indexed union pass (union_progressive_schedule experimental_v1) --------
void test_progressive_schedule() {
    namespace fs = std::filesystem;
    const auto fresh = [&](const std::string& rel, bool progressive) {
        fs::remove_all(path(rel));
        fs::create_directories(path(rel));
        panvar::GenotypeCommandOptions o =
            repair_options(rel + "/call", path("phaseable_r1.fq"), path("phaseable_r2.fq"), "tiny_phaseable.gfa");
        o.output_schema = "sequence_first_v2";
        o.setting_sources["output_schema"] = "test";
        o.early_incumbent = "experimental_v1";
        o.setting_sources["early_incumbent"] = "test";
        o.union_scoring_backend = "indexed_fine_experimental";
        o.indexed_union_store_bytes = std::size_t{256} << 20;
        o.max_union_cached_placements_consumed = 2000000000ULL;
        o.max_union_cached_origins_consumed = 2000000000ULL;
        o.max_union_cached_states_consumed = 2000000000ULL;
        o.max_union_cached_prefilters_consumed = 2000000000ULL;
        o.setting_sources["union_scoring_backend"] = "test";
        o.setting_sources["indexed_union_store_bytes"] = "test";
        if (progressive) {
            o.union_progressive_schedule = "experimental_v1";
            o.setting_sources["union_progressive_schedule"] = "test";
        }
        return o;
    };
    // (1) complete runs: every output identical; only declared status/settings rows (and executed-work rows) differ
    require(panvar::execute_genotype_command(fresh("ps_off", false), nullptr) == 0, "S2-D control predicts");
    require(panvar::execute_genotype_command(fresh("ps_on", true), nullptr) == 0, "S2-D progressive predicts");
    require(sf::files(path("ps_off")) == sf::files(path("ps_on")), "S2-D: the same files");
    for (const std::string& f : sf::files(path("ps_off"))) {
        if (f == "call.mosaic_status.tsv" || f == "call.settings.tsv" || f == "call.repair_memory.tsv" || f == "call.dosage_status.tsv") continue;
        if (f == "call.repair_rounds.tsv") {
            require(rounds_science(path("ps_off/" + f)) == rounds_science(path("ps_on/" + f)), "S2-D: trajectory identical");
            continue;
        }
        require(sf::normalised(path("ps_off/" + f), path("ps_off")) == sf::normalised(path("ps_on/" + f), path("ps_on")),
                "S2-D: " + f + " identical (scheduling only)");
    }
    {
        const auto a = sf::kv(path("ps_off/call.mosaic_status.tsv")), b = sf::kv(path("ps_on/call.mosaic_status.tsv"));
        std::vector<std::string> differing;
        for (const auto& [k, v] : a) {
            if (sf::volatile_key(k)) continue;
            const auto it = b.find(k);
            if (it == b.end() || it->second != v) differing.push_back(k);
        }
        std::string list;
        for (const auto& k : differing) list += k + ",";
        std::cout << "  S2-D control/progressive status keys differing (non-timing): " << (list.empty() ? "none" : list) << '\n';
        // what may differ is executed work, cache accounting and store/aggregate peaks -- never a score, call or filter
        for (const auto& k : differing) {
            const bool resource = k.rfind("checkpoint_", 0) == 0 ||   // declared: the prefix adds a checkpoint commit
                                  k.find("emission") != std::string::npos || k.find("executed") != std::string::npos ||
                                  k.find("indexed_union") != std::string::npos || k.find("cache") != std::string::npos ||
                                  k.find("store") != std::string::npos || k.find("prefilter") != std::string::npos ||
                                  k.find("threads") != std::string::npos || k.find("candidate_fragment_evaluations") != std::string::npos;
            require(resource, "S2-D: status key " + k + " differs and is not a declared resource/work row");
        }
        const auto st = read_status(path("ps_on/call.mosaic_status.tsv"));
        require(st.at("union_progressive_completed") == "yes" && st.at("union_progressive_prefix_checkpoint") == "committed" &&
                    st.at("early_incumbent_union_check") == "reproduced" && st.at("union_progressive_stage1_candidates") == "1",
                "S2-D: schedule completed, prefix checkpoint committed, early candidate reproduced by the union");
    }
    // (2) a stop right after the prefix checkpoint keeps it as PROVISIONAL; a later stop keeps a later checkpoint
    const auto stopped = [&](const std::string& rel, const std::string& point, std::size_t threads = 2) {
        panvar::GenotypeCommandOptions o = fresh(rel, true);
        o.test_interrupt_point = point;
        o.threads = threads;
        return s2a::execute_in_child(o);
    };
    {
        const int code = stopped("ps_prefix", "after_prefix_checkpoint");
        const auto locus = sf::kv(path("ps_prefix/call.sequence_locus.tsv"));
        require(code == 3 && locus.at("prediction") == "PROVISIONAL" && locus.at("checkpoint_stage") == "union_prefix_best" &&
                    locus.at("checkpoint_tied_groups") != "NA" && locus.at("checkpoint_projection_ambiguous") != "NA" &&
                    locus.at("checkpoint_semantics").find("no_confidence_from_the_evaluated_count") != std::string::npos &&
                    no_ranking_or_call("ps_prefix/call"),
                "S2-D: the prefix checkpoint survives a stop as PROVISIONAL, ties and projection recorded (exit " + std::to_string(code) + ")");
        const int later = stopped("ps_later", "adaptive_refinement");
        const auto late = sf::kv(path("ps_later/call.sequence_locus.tsv"));
        require(later == 3 && late.at("checkpoint_stage") != "union_prefix_best" && late.at("checkpoint_stage") != "early_panel_incumbent" &&
                    std::stod(late.at("checkpoint_total_score")) >= std::stod(locus.at("checkpoint_total_score")),
                "S2-D: a later stop keeps the last committed checkpoint, never worse than the prefix's");
        // 1 vs 4 threads: the same prefix checkpoint
        require(stopped("ps_t1", "after_prefix_checkpoint", 1) == 3 && stopped("ps_t4", "after_prefix_checkpoint", 4) == 3,
                "S2-D threads: both stop at the prefix checkpoint");
        auto a = sf::kv(path("ps_t1/call.sequence_locus.tsv")), b = sf::kv(path("ps_t4/call.sequence_locus.tsv"));
        for (auto* m : {&a, &b}) {
            m->erase("checkpoint_seconds_first_evaluated_prediction");
            m->erase("checkpoint_seconds_committed");
            m->erase("completion_note");
        }
        require(a == b, "S2-D threads: the prefix checkpoint is identical at 1 and 4 threads");
    }
    // (3) a union cap during the schedule: PROVISIONAL from the last committed checkpoint, never a full ranking
    {
        const auto ref = read_status(path("ps_on/call.mosaic_status.tsv"));
        panvar::GenotypeCommandOptions o = fresh("ps_cap", true);
        o.max_physical_origins = std::stoull(ref.at("early_incumbent_physical_origins"));
        const int code = panvar::execute_genotype_command(o, nullptr);
        const auto locus = sf::kv(path("ps_cap/call.sequence_locus.tsv"));
        const auto st = read_status(path("ps_cap/call.mosaic_status.tsv"));
        require(code == 3 && locus.at("completion") == "union_scoring_cap" && st.at("union_progressive_completed") == "no" &&
                    (locus.at("checkpoint_stage") == "early_panel_incumbent" || locus.at("checkpoint_stage") == "union_prefix_best"),
                "S2-D: a cap inside the schedule keeps the last committed checkpoint as PROVISIONAL (exit " + std::to_string(code) + ")");
    }
    // (4) the contract: progressive scheduling needs the early incumbent and the indexed backend
    {
        panvar::GenotypeCommandOptions o = fresh("ps_refused", true);
        o.union_scoring_backend = "reference";
        bool refused = false;
        try { (void)panvar::execute_genotype_command(o, nullptr); } catch (const std::invalid_argument&) { refused = true; }
        require(refused && sf::files(path("ps_refused")).empty(), "S2-D: progressive scheduling without the indexed backend is refused");
    }
    std::cout << "  progressive schedule: complete-run outputs identical to the bulk indexed route; prefix checkpoint "
                 "PROVISIONAL with ties/projection; later stop and union cap keep the last checkpoint; 1/4 threads; contract\n";
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: genotype_command_e2e_test <scratch dir> [source dir]\n";
        return 2;
    }
    dir = argv[1];
    // I3-C: a focused run of the prepared-locus gates only (they need nothing but the declared library)
    if (const char* only = std::getenv("PANVAR_E2E_ONLY"); only != nullptr && std::string(only) == "prepared_locus") {
        try {
            std::filesystem::create_directories(dir);
            write_library(path("lib.txt"));
            test_prepared_locus();
            test_prepared_locus_cleanup();
            test_prepared_locus_index();
            std::cout << "genotype command e2e (prepared_locus only): PASS\n";
            return 0;
        } catch (const std::exception& error) {
            std::cerr << "genotype command e2e (prepared_locus only): FAIL: " << error.what() << '\n';
            return 1;
        }
    }
    // S1-D.1: a focused run of the sequence-first gates (their fixtures written here exactly as the full run writes them)
    if (const char* only = std::getenv("PANVAR_E2E_ONLY"); only != nullptr && std::string(only) == "sequence_first") {
        try {
            std::filesystem::create_directories(dir);
            write_library(path("lib.txt"));
            {
                const TinyLocus locus = write_tiny_locus("tiny_phaseable.gfa", 80, 70);
                const auto truths = spell_truths(locus);
                simulate(truths.first, truths.second, path("phaseable_r1.fq"), path("phaseable_r2.fq"), 123);
            }
            {
                const TinyLocus locus = write_tiny_locus();
                const auto truths = spell_truths(locus);
                simulate(truths.first, truths.second, path("unphaseable_r1.fq"), path("unphaseable_r2.fq"), 123);
            }
            {   // test_prepared_locus's sample A
                const TinyLocus locus = write_tiny_locus("prep_locus.gfa", 80, 70);
                const auto truths = spell_truths(locus);
                simulate(truths.first, truths.second, path("prep_A_r1.fq"), path("prep_A_r2.fq"), 123);
            }
            test_sequence_first_schema();
            test_sequence_first_lifecycle();
            test_sequence_first_v2();
            test_sequence_first_v2_contract();
            test_early_incumbent();
            test_early_incumbent_caps();
            test_progressive_schedule();
            std::cout << "genotype command e2e (sequence_first only): PASS\n";
            return 0;
        } catch (const std::exception& error) {
            std::cerr << "genotype command e2e (sequence_first only): FAIL: " << error.what() << '\n';
            return 1;
        }
    }
    try {
        test_end_to_end();
        test_refusals();
        test_fail_closed_statuses();
        test_off_panel_truth_is_recovered();
        test_off_panel_phase_is_not_claimed_from_reads();
        test_proposal_status_rules();
        test_block_genotype_confidence();
        test_prior_only_phase_is_reported();
        test_phase_audit();
        test_candidate_prefixes_match_standalone();
        test_panel_only_path();
        test_block_consensus_rule();
        test_block_consensus_end_to_end();
        test_full_beam_audit();
        test_proposal_payload_ceiling();
        test_curated_panel();
        test_resource_caps_are_thread_invariant();
        test_pooled_decision_rule();
        test_phase_switch_symmetry();
        test_tie_branching_engine();
        test_tie_branch_settle_points();
        test_representation_guard();
        test_adaptive_completion();
        test_representation_guard_outside_pool();
        test_watchdog_and_interruption();
        test_indexed_union_backend();
        test_adaptive_fine_windows();
        test_store_layout_and_timing_detail();
        test_fused_anchor_seeding();
        test_marker_proposal_reuse();
        test_prepared_combination_terms();
        test_gap_recursion_backend();
        test_emission_handoff();
        test_prepared_locus();
        test_prepared_locus_cleanup();
        test_prepared_locus_index();
        test_sequence_first_schema();
        test_sequence_first_lifecycle();
        test_sequence_first_v2();
        test_sequence_first_v2_contract();
        test_early_incumbent();
        test_early_incumbent_caps();
        test_progressive_schedule();
        if (argc >= 3) test_profile(argv[2]);
        test_block_quality_columns();
        std::cout << "genotype command end-to-end: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "genotype command end-to-end: FAIL: " << error.what() << '\n';
        return 1;
    }
}
