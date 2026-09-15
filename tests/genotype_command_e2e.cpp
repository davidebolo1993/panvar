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

#include <cmath>
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
    options.max_candidates = 16;
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
    require(std::stod(reported.at("certified_gap")) > 0.0,
            "a COMPLETE run must carry a positive certified gap");
    require(reported.at("threads_used") == "1",
            "version one is serial; a threads_used above one would contradict the contract");

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
    again.max_candidates = 8;
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
    undeclared.max_candidates = 16;
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
    options.max_candidates = 16;
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
    // two zero-switch phasings and the truth was never offered. Under the same budget of 64 merged states
    // the truth must now be proposed. The locus has 36 distinct diplotypes and 64 merged states hold only
    // four of them, so with K=16 the search is pruned and misses its target: the run must refuse to call,
    // and its diagnostic manifest must still show the truth among what was proposed.
    {
        panvar::GenotypeCommandOptions budget = tiny_options(
            "phaseable_budget", path("phaseable_r1.fq"), path("phaseable_r2.fq"), "tiny_phaseable.gfa");
        budget.proposal_state_limit = 64;
        require(panvar::run_genotype(budget) == panvar::GenotypeRunStatus::CandidateIncomplete,
                "a pruned search that misses K must refuse even when it proposed the truth");
        const auto budget_status = read_status(path("phaseable_budget.mosaic_status.tsv"));
        require(budget_status.at("proposal_search_exhaustive") == "false" &&
                    budget_status.at("proposal_target_reached") == "false",
                "fixture: 64 merged states must prune and miss K=16");
        require(no_ranking_or_call("phaseable_budget"), "a refused proposal must not rank or call");
        const TruthFrame budget_frame = rebuild_truth_frame(
            budget, budget_status.at("model_fingerprint"), kTruth1Nodes, kTruth2Nodes);
        require(truth_phasing_in_candidates(budget_frame, path("phaseable_budget.mosaic_candidates.tsv")),
                "under the budget that used to crowd it out, the merged beam must propose the truth");
    }

    // Recovery, with a budget above every reachable merged state.
    panvar::GenotypeCommandOptions options =
        tiny_options("phaseable", path("phaseable_r1.fq"), path("phaseable_r2.fq"), "tiny_phaseable.gfa");
    options.proposal_state_limit = 1048576;
    const panvar::GenotypeRunStatus status = panvar::run_genotype(options);
    const std::map<std::string, std::string> reported = read_status(path("phaseable.mosaic_status.tsv"));
    require(status == panvar::GenotypeRunStatus::Complete,
            std::string("the phaseable off-panel run must complete, got ") +
                panvar::genotype_run_status_name(status) + ": " +
                (reported.count("note") ? reported.at("note") : ""));
    require(reported.at("proposal_search_exhaustive") == "true" || reported.at("proposal_target_reached") == "true",
            "an accuracy-eligible run must have searched exhaustively or reached its target");
    require(reported.at("phase_support") == "READS",
            "the reads, not the prior, must choose the phaseable off-panel truth");

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
              << " blocks called and phased correctly; truth proposed within 64 merged states; phase support READS, "
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
    options.proposal_state_limit = 1048576;   // far above every reachable merged state: exhaustive
    const panvar::GenotypeRunStatus status = panvar::run_genotype(options);
    const std::map<std::string, std::string> reported = read_status(path("unphaseable.mosaic_status.tsv"));
    require(reported.at("proposal_search_exhaustive") == "true", "fixture: the unphaseable search must be exhaustive");
    require((status == panvar::GenotypeRunStatus::RescoreIncomplete && reported.at("phase_support") == "UNRESOLVED") ||
                (status == panvar::GenotypeRunStatus::Complete && reported.at("phase_support") == "PRIOR_ONLY"),
            std::string("an unphaseable truth must never be reported as read-phased, got ") +
                panvar::genotype_run_status_name(status) + " with phase support " +
                (reported.count("phase_support") ? reported.at("phase_support") : "NA"));
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


// The inverted wrong-call gate. At NA19240 K=8 a pruned proposal that never offered the truth produced a
// certified, fragment-decided, wrong call. A pruned search that misses its target must now end
// CANDIDATE_INCOMPLETE with no ranking and no call; an exhaustive search with fewer candidates than
// requested is complete and unconditional.
void test_proposal_status_rules() {
    const std::string r1 = path("phaseable_r1.fq");
    const std::string r2 = path("phaseable_r2.fq");

    panvar::GenotypeCommandOptions starved = tiny_options("starved", r1, r2, "tiny_phaseable.gfa");
    starved.proposal_state_limit = 1;
    require(panvar::run_genotype(starved) == panvar::GenotypeRunStatus::CandidateIncomplete,
            "a pruned search that misses its target must end CANDIDATE_INCOMPLETE");
    const auto starved_status = read_status(path("starved.mosaic_status.tsv"));
    require(starved_status.at("proposal_search_exhaustive") == "false" &&
                starved_status.at("proposal_target_reached") == "false" &&
                std::stoull(starved_status.at("proposal_states_pruned")) > 0,
            "the status must report the pruned, target-unmet search separately");
    require(no_ranking_or_call("starved"), "a target-unmet pruned search must not rank or call");

    panvar::GenotypeCommandOptions exhaustive = tiny_options("exhaustive", r1, r2, "tiny_phaseable.gfa");
    exhaustive.max_candidates = 64;
    exhaustive.proposal_state_limit = 1048576;
    require(panvar::run_genotype(exhaustive) == panvar::GenotypeRunStatus::Complete,
            "an exhaustive search with fewer candidates than requested must complete");
    const auto exhaustive_status = read_status(path("exhaustive.mosaic_status.tsv"));
    require(exhaustive_status.at("proposal_search_exhaustive") == "true" &&
                exhaustive_status.at("proposal_target_reached") == "false" &&
                exhaustive_status.at("candidate_set_conditional") == "false" &&
                std::stoul(exhaustive_status.at("unique_candidates_found")) < 64,
            "an exhaustive search below K must be reported complete and unconditional");
    require(exhaustive_status.at("candidate_scope") == "COMPLETE" && exhaustive_status.at("global_certified") == "true",
            "an exhaustive, untruncated, complete ranking is the only globally certified result");
    require(exhaustive_status.at("marker_depth_converged") == "true",
            "a completed run must rest on a converged marker-depth fit");
    require(exhaustive_status.at("marker_background_effective") == exhaustive_status.at("marker_background") &&
                std::stod(exhaustive_status.at("marker_background_floor")) == panvar::kMarkerBackgroundFloor,
            "the status must report the model floor and the background proposal scoring actually used");
}

// Prior-only phase. The unphaseable locus with a second carrier of sampleA's haplotype: reads from
// sampleA and sampleB cannot link neighbouring variants, so the zero-switch phasings are fragment
// contenders, and the extra carrier lets the prior pick one. The call is made, but its phase must be
// reported PRIOR_ONLY, while every block genotype -- on which all contenders agree -- stays READS.
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
    const auto slurp = [](const std::string& file, bool status_file) {
        std::ifstream input(file);
        if (!input) return std::string("<absent>");
        std::string line;
        std::string text;
        while (std::getline(input, line)) {
            if (status_file && line.rfind("seconds_", 0) == 0) continue;
            if (!line.empty() && line[0] == '>') line = line.substr(line.find(' '));   // header names the prefix
            text += line + '\n';
        }
        return text;
    };
    const std::vector<std::size_t> prefixes = {2, 4, 16, 64};
    for (const std::size_t state_limit : {std::size_t{64}, std::size_t{1048576}}) {
        const std::string base = "prefixes_" + std::to_string(state_limit);
        panvar::GenotypeCommandOptions ladder = tiny_options(base, r1, r2, "tiny_phaseable.gfa");
        ladder.proposal_state_limit = state_limit;
        ladder.max_candidates = prefixes.back();
        ladder.development_candidate_prefixes = prefixes;
        panvar::run_genotype(ladder);
        std::size_t refused = 0;
        std::size_t complete = 0;
        for (const std::size_t k : prefixes) {
            const std::string name = base + "_standalone_K" + std::to_string(k);
            panvar::GenotypeCommandOptions standalone = tiny_options(name, r1, r2, "tiny_phaseable.gfa");
            standalone.proposal_state_limit = state_limit;
            standalone.max_candidates = k;
            const panvar::GenotypeRunStatus status = panvar::run_genotype(standalone);
            if (status == panvar::GenotypeRunStatus::CandidateIncomplete) ++refused;
            if (status == panvar::GenotypeRunStatus::Complete) ++complete;
            const std::string extracted = path(base + ".K" + std::to_string(k));
            for (const std::string suffix : {".genotypes.tsv", ".block_genotypes.tsv", ".mosaic_candidates.tsv",
                                             ".mosaic_scores.tsv", ".mosaic_status.tsv", ".mosaic_haplotypes.fa"}) {
                const bool status_file = suffix == ".mosaic_status.tsv";
                require(slurp(extracted + suffix, status_file) == slurp(path(name) + suffix, status_file),
                        "prefix K=" + std::to_string(k) + " at state limit " + std::to_string(state_limit) +
                            " must reproduce the standalone " + suffix);
            }
        }
        if (state_limit == 64) {
            require(refused > 0 && complete > 0,
                    "fixture: at 64 merged states the small prefixes must call and the large ones refuse");
        } else {
            require(refused == 0 && complete == prefixes.size(),
                    "fixture: an exhaustive proposal must complete at every prefix");
        }
        std::cout << "  candidate prefixes at state limit " << state_limit << ": " << complete << " complete, "
                  << refused << " refused, each identical to its standalone run\n";
    }
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: genotype_command_e2e_test <scratch dir>\n";
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
        test_prior_only_phase_is_reported();
        test_candidate_prefixes_match_standalone();
        std::cout << "genotype command end-to-end: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "genotype command end-to-end: FAIL: " << error.what() << '\n';
        return 1;
    }
}
