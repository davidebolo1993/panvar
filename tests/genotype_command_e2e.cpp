// Phase 7. The command is orchestration only, so its gate is that the orchestration is right: the
// declared execution order runs, every output appears, an incomplete run says so instead of
// presenting a partial ranking, and the model fingerprint identifies the model rather than the run.

#include "panvar/genotype_command.hpp"
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

#include <sys/wait.h>
#include <unistd.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <cstdio>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <optional>
#include <random>
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
    for (const auto& row : read_rows(path(prefix + ".block_genotypes.tsv"))) {
        if (row.empty() || row.front() == "block_id") continue;   // header, if read_rows returns one
        require(row.back() == expected,
                prefix + " block row " + row.front() + " reports scope " + row.back() +
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
                const std::size_t last = line.rfind('\t');
                if (last != std::string::npos) line = line.substr(0, last);
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
        require(s.at("status") == "RESOURCE_REFUSAL", "the status must say RESOURCE_REFUSAL");
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

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: genotype_command_e2e_test <scratch dir> [source dir]\n";
        return 2;
    }
    dir = argv[1];
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
        test_adaptive_completion();
        test_watchdog_and_interruption();
        if (argc >= 3) test_profile(argv[2]);
        std::cout << "genotype command end-to-end: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "genotype command end-to-end: FAIL: " << error.what() << '\n';
        return 1;
    }
}
