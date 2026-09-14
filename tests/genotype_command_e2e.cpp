// Phase 7. The command is orchestration only, so its gate is that the orchestration is right: the
// declared execution order runs, every output appears, an incomplete run says so instead of
// presenting a partial ranking, and the model fingerprint identifies the model rather than the run.

#include "panvar/genotype_command.hpp"

#include <cstdio>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
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
};

TinyLocus write_tiny_locus() {
    std::mt19937_64 rng(4);
    TinyLocus locus;
    locus.gfa_path = path("tiny.gfa");

    struct Node { std::string id, seq; };
    const std::vector<Node> nodes = {
        {"1", random_bases(rng, 700)},     // left flank
        {"2", random_bases(rng, 300)},     // SV block one, reference allele
        {"3", random_bases(rng, 90)},      // SV block one, alternate allele
        {"4", random_bases(rng, 350)},     // shared context
        {"5", random_bases(rng, 10)},      // sub-threshold variant, allele one
        {"6", random_bases(rng, 10)},      // sub-threshold variant, allele two
        {"7", random_bases(rng, 340)},     // shared context
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
    std::map<std::string, std::string> sequence_of;
    for (const Node& node : nodes) sequence_of[node.id] = node.seq;
    for (const auto& item : templates) {
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
           "substitution_error 0.002\noutlier_weight 0.001\nfragment_start_intensity 0.02\n";
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
    for (const char* suffix : {".genotypes.tsv", ".mosaic_candidates.tsv", ".mosaic_scores.tsv",
                               ".mosaic_status.tsv", ".mosaic_haplotypes.fa"}) {
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
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: genotype_command_e2e_test <scratch dir>\n";
        return 2;
    }
    dir = argv[1];
    try {
        test_end_to_end();
        test_refusals();
        std::cout << "genotype command end-to-end: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "genotype command end-to-end: FAIL: " << error.what() << '\n';
        return 1;
    }
}
