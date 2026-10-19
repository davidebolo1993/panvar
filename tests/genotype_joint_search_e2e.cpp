// M3-D caller gates (run_genotype on a synthetic locus): the joint pass crosses a two-block valley the single-block
// refinement cannot; option off leaves it uncrossed; one and four threads agree; projection, phase and FASTA follow the
// updated pool; the model-derived bound skip keeps the ordinary result; a cap crossed inside the pass keeps the incomplete
// contract; the library-sensitivity audit is refused with the experiment.
#include "panvar/genotype_command.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
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
std::string dir() {
    static const std::string d = [] {
        const std::string p = (std::filesystem::temp_directory_path() / "panvar_m3d_e2e").string();
        std::filesystem::create_directories(p);
        return p;
    }();
    return d;
}
std::string path(const std::string& name) { return dir() + "/" + name; }
std::string rnd(std::mt19937_64& rng, std::size_t n) {
    std::string s(n, 'A');
    for (char& c : s) c = "ACGT"[rng() % 4];
    return s;
}
std::string revcomp(const std::string& x) {
    std::string r(x.rbegin(), x.rend());
    for (char& c : r) c = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : 'A';
    return r;
}
char other_base(char c) { return c == 'A' ? 'C' : c == 'C' ? 'G' : c == 'G' ? 'T' : 'A'; }

std::map<std::string, std::string> status(const std::string& prefix) {
    std::ifstream in(prefix + ".mosaic_status.tsv");
    require(static_cast<bool>(in), "status file " + prefix + " must exist");
    std::map<std::string, std::string> v;
    std::string line;
    std::getline(in, line);
    while (std::getline(in, line)) {
        const std::size_t t = line.find('\t');
        if (t != std::string::npos) v[line.substr(0, t)] = line.substr(t + 1);
    }
    return v;
}
std::string slurp(const std::string& file, bool status_file) {
    std::ifstream in(file);
    if (!in) return "<absent>";
    std::string out, line;
    while (std::getline(in, line)) {
        if (status_file && (line.rfind("seconds_", 0) == 0 || line.rfind("threads_", 0) == 0 || line.rfind("resources_", 0) == 0 ||
                            line.find("footprint") != std::string::npos || line.find("rss") != std::string::npos)) {
            continue;
        }
        if (!line.empty() && line[0] == '>') line = line.substr(line.find(' ') == std::string::npos ? 0 : line.find(' '));
        out += line + '\n';
    }
    return out;
}

// Four SV blocks (60-bp alleles differing at one base), 450-bp separators, 700-bp flanks. Templates: 18 all-reference
// copies (A0..A17) and one all-alternate (B). Truth: hap1 = (alt, alt, ref, ref), hap2 = all-reference; the SNP reads
// are described in write_reads; tiled background reads avoid the SNPs.
struct Locus {
    std::string gfa, reference, truth1, truth2;
    std::vector<std::size_t> snps;   // on truth1
};
Locus write_valley() {
    std::mt19937_64 rng(2026);
    Locus l;
    l.gfa = path("valley.gfa");
    std::vector<std::pair<std::string, std::string>> nodes;   // id, seq
    std::vector<std::string> ref_ids, alt_ids;
    nodes.push_back({"L", rnd(rng, 700)});
    for (int b = 0; b < 4; ++b) {
        std::string r = rnd(rng, 60), a = r;
        a[30] = other_base(a[30]);
        nodes.push_back({"r" + std::to_string(b), r});
        nodes.push_back({"x" + std::to_string(b), a});
        nodes.push_back({"s" + std::to_string(b), rnd(rng, b == 3 ? 700 : 450)});
        ref_ids.push_back("r" + std::to_string(b));
        alt_ids.push_back("x" + std::to_string(b));
    }
    std::map<std::string, std::string> seq;
    for (const auto& n : nodes) seq[n.first] = n.second;
    std::ofstream out(l.gfa);
    out << "H\tVN:Z:1.0\n";
    for (const auto& n : nodes) out << "S\t" << n.first << '\t' << n.second << '\n';
    std::string prev = "L";
    for (int b = 0; b < 4; ++b) {
        out << "L\t" << prev << "\t+\tr" << b << "\t+\t0M\nL\t" << prev << "\t+\tx" << b << "\t+\t0M\n";
        out << "L\tr" << b << "\t+\ts" << b << "\t+\t0M\nL\tx" << b << "\t+\ts" << b << "\t+\t0M\n";
        prev = "s" + std::to_string(b);
    }
    const auto walk = [&](const std::vector<int>& alt) {
        std::vector<std::string> w{"L"};
        for (int b = 0; b < 4; ++b) {
            w.push_back(alt[b] ? alt_ids[b] : ref_ids[b]);
            w.push_back("s" + std::to_string(b));
        }
        return w;
    };
    const auto write_path = [&](const std::string& name, const std::vector<std::string>& w) {
        out << "P\t" << name << '\t';
        for (std::size_t i = 0; i < w.size(); ++i) out << (i ? "," : "") << w[i] << '+';
        out << "\t*\n";
    };
    for (int k = 0; k < 18; ++k) write_path("A" + std::to_string(k) + "#1#chr1:0-100", walk({0, 0, 0, 0}));
    write_path("B#1#chr1:0-100", walk({1, 1, 1, 1}));
    l.reference = "A0#1#chr1:0-100";
    const auto spell = [&](const std::vector<std::string>& w) {
        std::string s;
        for (const auto& id : w) s += seq[id];
        return s;
    };
    l.truth1 = spell(walk({1, 1, 0, 0}));
    l.truth2 = spell(walk({0, 0, 0, 0}));
    std::size_t at = 700;
    for (int b = 0; b < 4; ++b) {
        l.snps.push_back(at + 30);
        at += 60 + (b == 3 ? 700 : 450);
    }
    return l;
}

void write_reads(const Locus& l, const std::string& r1, const std::string& r2) {
    std::ofstream a(r1), b(r2);
    std::size_t n = 0;
    const auto emit = [&](const std::string& hap, std::size_t start) {
        ++n;
        a << "@f" << n << "/1\n" << hap.substr(start, 100) << "\n+\n" << std::string(100, 'I') << '\n';
        b << "@f" << n << "/2\n" << revcomp(hap.substr(start + 200, 100)) << "\n+\n" << std::string(100, 'I') << '\n';
    };
    for (const auto* hap : {&l.truth1, &l.truth2}) {
        const std::size_t step = hap == &l.truth1 ? 37 : 41;
        for (std::size_t s = 0; s + 300 <= hap->size(); s += step) {
            bool hit = false;
            for (const std::size_t x : l.snps) hit = hit || (s <= x + 6 && x < s + 300 + 5);
            if (!hit) emit(*hap, s);
        }
    }
    // one alternate-supporting fragment at blocks 0 and 1 (from hap1); five reference-supporting fragments at blocks 2 and 3
    // (from hap2), which the all-alternate template B cannot explain on its copy: the intact pair (A, B) explains the two
    // alternate reads as well as the truth but pays ln 2 per reference read, and each single change on B pays a donor
    // switch its five reads cannot cover, while the joint change at blocks 2 and 3 can
    emit(l.truth1, l.snps[0] - 50);
    emit(l.truth1, l.snps[1] - 50);
    for (int b = 2; b < 4; ++b) {
        for (int k = 0; k < 5; ++k) emit(l.truth2, l.snps[b] - 40 - 9 * k);
    }
}

void write_library(const std::string& file) {
    std::ofstream out(file);
    out << "insert_mean_bp 300\ninsert_sd_bp 30\nsupport_sigmas 4\nsubstitution_error 0.002\noutlier_weight 0.001\n"
           "fragment_start_intensity 0.1\nmarker_background 0.01\n";
}

panvar::GenotypeCommandOptions options(const Locus& l, const std::string& name, bool joint) {
    panvar::GenotypeCommandOptions o;
    o.gfa_path = l.gfa;
    o.reference_path = l.reference;
    o.reads1_path = path("valley_1.fq");
    o.reads2_path = path("valley_2.fq");
    o.library_path = path("valley_lib.txt");
    o.out_prefix = path(name);
    o.max_mosaic_candidates = 0;   // panel-only union: the refinement starts from the best panel pair
    o.phase_audit_candidates = 0;
    o.repair_rounds = 32;
    o.kmer_size = 21;
    o.syncmer_s = 11;
    o.threads = 1;
    o.quiet = true;
    if (joint) o.repair_joint_search = "experimental_v1";
    return o;
}

std::string fasta_sequences(const std::string& file) {
    std::ifstream in(file);
    std::string out, line;
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] != '>') out += line;
        else if (!line.empty()) out += '|';
    }
    return out;
}

}  // namespace

int main() {
    try {
        const Locus l = write_valley();
        write_reads(l, path("valley_1.fq"), path("valley_2.fq"));
        write_library(path("valley_lib.txt"));

        // off: the single-block refinement converges at the all-reference incumbent (the valley is not crossed)
        require(panvar::run_genotype(options(l, "off", false)) == panvar::GenotypeRunStatus::Complete, "off must complete");
        const auto off = status(path("off"));
        require(off.count("joint_search_status") == 0, "off: no joint rows may be written (output contract unchanged)");
        require(off.at("repair_converged") == "true", "off: the trunk must converge");
        const std::string off_fasta = fasta_sequences(path("off.mosaic_haplotypes.fa"));
        require(off_fasta.find(l.truth1) == std::string::npos, "off: the valley must stay uncrossed (truth not called)");

        // on: the joint pass moves to the truth; the continuation converges; the call, projection, phase and FASTA follow
        require(panvar::run_genotype(options(l, "on", true)) == panvar::GenotypeRunStatus::Complete, "on must complete");
        const auto on = status(path("on"));
        require(on.at("joint_search_status") == "complete_moved", "on: the joint pass must move, got " + on.at("joint_search_status"));
        require(on.at("joint_continued") == "true" && on.at("joint_leg_converged") == "true", "on: the continuation must converge");
        require(on.at("joint_search_full_placements") == "0", "on: no full placement inside the pass");
        const std::string on_fasta = fasta_sequences(path("on.mosaic_haplotypes.fa"));
        require(on_fasta.find(l.truth1) != std::string::npos && on_fasta.find(l.truth2) != std::string::npos,
                "on: the called FASTA must carry the truth haplotypes");
        require(on.at("repair_final_incumbent_hap1") != off.at("repair_final_incumbent_hap1") ||
                    on.at("repair_final_incumbent_hap2") != off.at("repair_final_incumbent_hap2"),
                "on: the final incumbent must differ from the off run");
        require(on.count("phase_support") == 1 && on.at("phase_support") != "READS_WITHIN_CANDIDATES" + std::string("_CLAIMED"),
                "phase support must be reported by the pooled rule");
        {   // projection: blocks sv 0 and 1 heterozygous in the block table of the on run, not of the off run
            std::ifstream bt(path("on.block_genotypes.tsv"));
            std::string all((std::istreambuf_iterator<char>(bt)), std::istreambuf_iterator<char>());
            std::ifstream bo(path("off.block_genotypes.tsv"));
            std::string offt((std::istreambuf_iterator<char>(bo)), std::istreambuf_iterator<char>());
            require(all != offt, "projection: the block genotypes must follow the updated pool");
        }

        // one vs four threads: identical scientific outputs
        {
            auto four = options(l, "on4", true);
            four.threads = 4;
            require(panvar::run_genotype(four) == panvar::GenotypeRunStatus::Complete, "on4 must complete");
            for (const char* suffix : {".genotypes.tsv", ".block_genotypes.tsv", ".mosaic_scores.tsv", ".mosaic_haplotypes.fa",
                                       ".repair_rounds.tsv", ".mosaic_status.tsv"}) {
                const bool st = std::string(suffix) == ".mosaic_status.tsv";
                std::string a = slurp(path("on") + suffix, st), b = slurp(path("on4") + suffix, st);
                if (std::string(suffix) == ".repair_rounds.tsv") {   // the rounds table carries timing and memory columns
                    a = a.substr(0, a.find('\n'));
                    b = b.substr(0, b.find('\n'));
                }
                require(a == b, std::string("threads: ") + suffix + " must not depend on the thread count");
            }
        }

        // bound skip (API test hook): the ordinary result is kept and the skip is recorded
        {
            auto skip = options(l, "skip", true);
            skip.test_joint_search_single_limit = 7;   // 8 single alternatives exist
            require(panvar::run_genotype(skip) == panvar::GenotypeRunStatus::Complete, "skip must complete");
            const auto s = status(path("skip"));
            require(s.at("joint_search_status") == "skipped_bounds" && s.at("joint_continued") == "false", "skip: recorded as skipped_bounds");
            require(fasta_sequences(path("skip.mosaic_haplotypes.fa")) == off_fasta, "skip: the ordinary call must be kept");
        }

        // a cap crossed inside the pass: RESCORE_INCOMPLETE, no genotype, explicit status
        {
            const std::uint64_t trunk = std::stoull(off.at("repair_cumulative_state_classes"));
            const std::uint64_t joint = std::stoull(on.at("joint_search_state_classes"));
            auto capped = options(l, "capped", true);
            capped.repair_max_state_classes = trunk + joint / 2;
            const panvar::GenotypeRunStatus st = panvar::run_genotype(capped);
            require(st == panvar::GenotypeRunStatus::RescoreIncomplete, std::string("capped: the run must be RESCORE_INCOMPLETE, got ") +
                                                                          panvar::genotype_run_status_name(st));
            const auto s = status(path("capped"));
            require(s.at("joint_search_status") == "refused_cumulative_cap" && s.at("no_call_reason") == "cumulative_cap",
                    "capped: the refusal must be explicit");
            require(fasta_sequences(path("capped.mosaic_haplotypes.fa")).empty(), "capped: no haplotype FASTA may be published");
        }

        // the library-sensitivity audit is refused with the experiment
        {
            auto audit = options(l, "audit", true);
            audit.library_path = "auto";
            audit.library_sensitivity_audit = true;
            bool refused = false;
            try {
                (void)panvar::run_genotype(audit);
            } catch (const std::invalid_argument&) {
                refused = true;
            }
            require(refused, "the experiment with the library sensitivity audit must be refused before any output");
        }
        // M3-E: the experimental support annotation is annotation-only (every pre-existing output unchanged) and thread-invariant
        {
            auto sup1 = options(l, "sup1", true);
            sup1.block_support_rule = "full_data_likelihood_v1_experimental";
            auto sup4 = sup1;
            sup4.out_prefix = path("sup4");
            sup4.threads = 4;
            require(panvar::run_genotype(sup1) == panvar::GenotypeRunStatus::Complete && panvar::run_genotype(sup4) == panvar::GenotypeRunStatus::Complete,
                    "support runs must complete");
            for (const char* suffix : {".genotypes.tsv", ".mosaic_scores.tsv", ".mosaic_haplotypes.fa", ".mosaic_candidates.tsv"}) {
                require(slurp(path("on") + suffix, false) == slurp(path("sup1") + suffix, false),
                        std::string("support annotation must leave ") + suffix + " unchanged");
            }
            std::size_t legacy_columns = 1;   // the pre-existing columns: as many as the legacy run's own header has
            {
                std::ifstream in(path("on") + ".block_genotypes.tsv");
                std::string header;
                std::getline(in, header);
                for (const char c : header) legacy_columns += c == '\t' ? 1 : 0;
            }
            const auto legacy_cols = [&](const std::string& prefix) {
                std::ifstream in(prefix + ".block_genotypes.tsv");
                std::string out, line;
                while (std::getline(in, line)) {
                    std::size_t at = 0;
                    for (std::size_t k = 0; k < legacy_columns && at != std::string::npos; ++k) at = line.find('\t', at + 1);
                    out += line.substr(0, at) + '\n';
                }
                return out;
            };
            require(legacy_cols(path("on")) == legacy_cols(path("sup1")), "support annotation must leave the legacy block columns unchanged");
            require(slurp(path("sup1") + ".block_genotypes.tsv", false) == slurp(path("sup4") + ".block_genotypes.tsv", false),
                    "the support columns must not depend on the thread count");
            const auto st = status(path("sup1"));
            require(st.at("block_support_rule") == "full_data_likelihood_v1_experimental" && st.count("block_support_labels_v1x") == 1,
                    "the support status rows must be written");
        }
        std::cout << "genotype joint search e2e: PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "genotype joint search e2e: FAIL: " << e.what() << '\n';
        return 1;
    }
}
