// Gate tool (not part of the suite: it reads sealed campaign outputs). Applies the production headline function to
// copies of sealed run prefixes and writes what it produced, for comparison with an independently derived table.
//
//   genotype_headline_replay <cells.tsv> <sealed replay root> <scratch dir> <out prefix>
//
// cells.tsv: dataset<TAB>job<TAB>dosage_mode<TAB>dosage_audit_status<TAB>comma-separated model block ids. Each prefix's
// status, block tables and dosage blocks are copied into the scratch dir; tables written before the rename carry the
// dosage stage's own filter as `dosage_filter`/`dosage_filter_reason`, which are renamed to `dosage_stage_filter`/
// `dosage_stage_filter_reason` first (the current schema). Writes <out>.blocks.tsv and <out>.loci.tsv.

#include "panvar/genotype_headline.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;
using namespace panvar::headline;

namespace {
std::vector<std::string> split(const std::string& s, char sep) {
    std::vector<std::string> out; std::string cur; std::stringstream in(s);
    while (std::getline(in, cur, sep)) out.push_back(cur);
    return out;
}
void copy_renamed(const fs::path& from, const fs::path& to) {
    if (!fs::exists(from)) return;
    std::ifstream in(from); std::ofstream out(to); std::string line; bool first = true;
    while (std::getline(in, line)) {
        if (first) {
            std::vector<std::string> h = split(line, '\t');
            for (std::string& c : h) {
                if (c == "dosage_filter") c = "dosage_stage_filter";
                else if (c == "dosage_filter_reason") c = "dosage_stage_filter_reason";
            }
            std::string joined;
            for (std::size_t i = 0; i < h.size(); ++i) joined += (i ? "\t" : "") + h[i];
            out << joined << '\n'; first = false;
        } else out << line << '\n';
    }
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 5) { std::cerr << "usage: genotype_headline_replay <cells.tsv> <replay root> <scratch> <out>\n"; return 2; }
    const fs::path root = argv[2], scratch = argv[3];
    fs::create_directories(scratch);
    std::ifstream cells(argv[1]);
    std::ofstream blocks(std::string(argv[4]) + ".blocks.tsv"), loci(std::string(argv[4]) + ".loci.tsv");
    blocks << "dataset\tjob\tblock_id\tgenotype_filter\tgenotype_filter_reasons\tdosage_filter\tdosage_filter_reasons\n";
    loci << "dataset\tjob";
    for (const std::string& c : locus_columns()) loci << '\t' << c;
    loci << '\n';
    std::string line;
    while (std::getline(cells, line)) {
        const std::vector<std::string> f = split(line, '\t');
        if (f.size() < 5) continue;
        const std::string& ds = f[0]; const std::string& job = f[1];
        const fs::path src = root / ds / job, dst = scratch / (ds + "." + job);
        for (const char* suffix : {".mosaic_status.tsv", ".block_genotypes.tsv", ".genotypes.tsv", ".dosage_blocks.tsv"}) {
            copy_renamed(fs::path(src.string() + suffix), fs::path(dst.string() + suffix));
        }
        const LocusHeadline h = apply_headline(dst.string(), split(f[4], ','), f[2], f[3]);
        const std::map<std::string, std::string> v = locus_values(h);
        loci << ds << '\t' << job;
        for (const std::string& c : locus_columns()) loci << '\t' << v.at(c);
        loci << '\n';
        std::ifstream bt(dst.string() + ".block_genotypes.tsv");
        std::string row; bool first = true; std::vector<std::string> head;
        while (std::getline(bt, row)) {
            const std::vector<std::string> cols = split(row, '\t');
            if (first) { head = cols; first = false; continue; }
            const auto at = [&](const std::string& name) {
                for (std::size_t i = 0; i < head.size(); ++i) if (head[i] == name) return cols[i];
                return std::string("NA");
            };
            blocks << ds << '\t' << job << '\t' << at("block_id") << '\t' << at("genotype_filter") << '\t' << at("genotype_filter_reasons")
                   << '\t' << at("dosage_filter") << '\t' << at("dosage_filter_reasons") << '\n';
        }
    }
    return 0;
}
