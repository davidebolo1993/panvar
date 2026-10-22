// The --panel-curation contract: a curated template is removed from the callable panel, not only from the panel-pair
// search. An allele carried only by curated templates is not an allele of the model, so no stage -- panel baseline,
// mosaic proposal, adaptive repair, the dosage census, features and both repair routes, the application guard or any
// provenance -- can reach it. A sequence-identical allele carried by a retained template stays callable under the
// retained source only. The graph's structure is unchanged by curation.

#include "panvar/genotype_command.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::string dir;
std::string path(const std::string& name) { return dir + "/" + name; }
void require(bool ok, const std::string& what) { if (!ok) throw std::runtime_error(what); }
void write(const std::string& file, const std::string& text) { std::ofstream(file) << text; }
std::string slurp(const std::string& file) {
    std::ifstream in(file);
    std::stringstream s;
    s << in.rdbuf();
    return s.str();
}
std::string revcomp(const std::string& bases) {
    std::string out(bases.rbegin(), bases.rend());
    for (char& c : out) c = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : 'A';
    return out;
}

// A locus of shared flanks and two variant sites. Each template is a node path; nodes are random unless given.
struct Locus {
    std::map<std::string, std::string> nodes;
    std::vector<std::pair<std::string, std::vector<std::string>>> templates;
    std::string spell(const std::string& name) const {
        for (const auto& t : templates) {
            if (t.first != name) continue;
            std::string s;
            for (const std::string& n : t.second) s += nodes.at(n);
            return s;
        }
        throw std::runtime_error("no template " + name);
    }
};

Locus base_locus() {
    std::mt19937_64 rng(4);
    const auto bases = [&](std::size_t n) {
        std::string s(n, 'A');
        for (char& c : s) c = "ACGT"[rng() % 4];
        return s;
    };
    Locus l;
    for (const auto& [id, n] : std::vector<std::pair<std::string, std::size_t>>{
             {"1", 700}, {"2", 300}, {"3", 90}, {"4", 350}, {"5", 10}, {"6", 10}, {"7", 340}, {"8", 260}, {"9", 80},
             {"10", 700}, {"11", 1200}}) {
        l.nodes[id] = bases(n);
    }
    l.templates = {{"sampleA#1#chr1:0-100", {"1", "2", "4", "5", "7", "8", "10"}},
                   {"sampleB#1#chr1:0-100", {"1", "3", "4", "6", "7", "9", "10"}},
                   {"sampleC#1#chr1:0-100", {"1", "2", "4", "6", "7", "9", "10"}},
                   {"sampleD#1#chr1:0-100", {"1", "3", "4", "5", "7", "8", "10"}}};
    return l;
}

void write_locus(const Locus& l, const std::string& file) {
    std::set<std::pair<std::string, std::string>> edges;
    for (const auto& t : l.templates) {
        for (std::size_t i = 0; i + 1 < t.second.size(); ++i) edges.insert({t.second[i], t.second[i + 1]});
    }
    std::ofstream out(file);
    out << "H\tVN:Z:1.0\n";
    for (const auto& [id, seq] : l.nodes) out << "S\t" << id << '\t' << seq << '\n';
    for (const auto& [a, b] : edges) out << "L\t" << a << "\t+\t" << b << "\t+\t0M\n";
    for (const auto& t : l.templates) {
        out << "P\t" << t.first << '\t';
        for (std::size_t i = 0; i < t.second.size(); ++i) out << (i ? "," : "") << t.second[i] << '+';
        out << "\t*\n";
    }
}

void simulate(const std::vector<std::string>& haplotypes, const std::string& tag) {
    std::mt19937_64 rng(123);
    std::ofstream r1(path(tag + "_r1.fq")), r2(path(tag + "_r2.fq"));
    std::size_t k = 0;
    for (const std::string& hap : haplotypes) {
        for (std::size_t i = 0; i < 30 * hap.size() / 400; ++i) {
            const std::size_t insert = 260 + rng() % 80;
            const std::size_t start = rng() % (hap.size() - insert);
            r1 << "@p" << k << "/1\n" << hap.substr(start, 100) << "\n+\n" << std::string(100, 'I') << '\n';
            r2 << "@p" << k << "/2\n" << revcomp(hap.substr(start + insert - 100, 100)) << "\n+\n" << std::string(100, 'I') << '\n';
            ++k;
        }
    }
}

const std::string kLibrary = "insert_mean_bp 300\ninsert_sd_bp 30\nsupport_sigmas 4\nsubstitution_error 0.002\n"
                             "outlier_weight 0.001\nfragment_start_intensity 0.085\nmarker_background 0.01\n";

int run(const std::string& tag, const std::string& reads, const std::string& mode, const std::string& curation) {
    std::vector<std::string> args = {"--graph", path(tag + ".gfa"), "--reference-path", "sampleA#1#chr1:0-100",
                                     "--reads1", path(reads + "_r1.fq"), "--reads2", path(reads + "_r2.fq"),
                                     "--library", path("library.txt"), "--threads", "2",
                                     "--out", path(tag + "." + mode + (curation.empty() ? "" : "." + std::filesystem::path(curation).stem().string())),
                                     "--quiet"};
    if (!curation.empty()) {
        args.push_back("--panel-curation");
        args.push_back(path(curation));
    }
    // M15-J: the legacy contract (its FASTA carries the sequences) and the dosage modes, repair included, are internal: set
    // through the API, as the removed --output-schema legacy / --dosage-mode flags did, with the preset's v2-only features
    // (early incumbent, compound annotation) off.
    panvar::GenotypeCommandOptions o;
    std::ostringstream sink;
    panvar::parse_genotype_arguments(args, o, sink);
    o.output_schema = "legacy";
    o.setting_sources["output_schema"] = "test (internal API)";
    o.early_incumbent = "off";
    o.setting_sources["early_incumbent"] = "test (internal API)";
    o.compound_annotation = "off";
    o.setting_sources["compound_annotation"] = "test (internal API)";
    o.dosage_mode = mode;
    o.dosage_mode_explicit = true;
    o.dosage_mode_source = "test (internal API)";
    o.diagnostics_directory = false;   // M15-L: the internal API's flat layout (this test reads the legacy reports by name)
    return panvar::execute_genotype_command(o, nullptr);
}

std::map<std::string, std::string> status_of(const std::string& prefix) {
    std::map<std::string, std::string> kv;
    std::istringstream in(slurp(path(prefix + ".mosaic_status.tsv")));
    for (std::string line; std::getline(in, line);) {
        const std::size_t tab = line.find('\t');
        if (tab != std::string::npos && kv.count(line.substr(0, tab)) == 0) kv[line.substr(0, tab)] = line.substr(tab + 1);
    }
    return kv;
}

// Every output file of a prefix, concatenated with line breaks removed (FASTA wraps at 60 bp): an excluded-only
// sequence or source must appear in none of them.
std::string all_outputs(const std::string& prefix) {
    std::string text;
    for (const auto& entry : std::filesystem::directory_iterator(dir)) {
        const std::string name = entry.path().filename().string();
        if (name.rfind(prefix + ".", 0) == 0) text += slurp(entry.path().string()) + "\t";
    }
    std::string flat;
    for (char c : text) if (c != '\n') flat += c;
    return flat;
}

// The sv:1 row of a genotypes table: (bp, sources) of each reported allele.
std::vector<std::pair<std::string, std::string>> sv1_alleles(const std::string& prefix) {
    std::istringstream in(slurp(path(prefix + ".genotypes.tsv")));
    std::string line;
    std::getline(in, line);
    std::vector<std::string> h;
    for (std::istringstream hs(line); std::getline(hs, line, '\t');) h.push_back(line);
    for (std::string row; std::getline(in, row);) {
        std::vector<std::string> c;
        for (std::istringstream rs(row); std::getline(rs, line, '\t');) c.push_back(line);
        std::map<std::string, std::string> m;
        for (std::size_t i = 0; i < h.size() && i < c.size(); ++i) m[h[i]] = c[i];
        if (m["block_id"] != "sv:1") continue;
        return {{m["allele1_bp"], m["allele1_sources"]}, {m["allele2_bp"], m["allele2_sources"]}};
    }
    throw std::runtime_error("no sv:1 row in " + prefix);
}

// A. The excluded template E carries a unique sv:1 allele (node 2 then the 1,200 bp node 11); reads come from E and A.
void test_unique_excluded_allele() {
    Locus l = base_locus();
    l.templates.push_back({"sampleE#1#chr1:0-100", {"1", "2", "11", "4", "5", "7", "8", "10"}});
    write_locus(l, path("unique.gfa"));
    simulate({l.spell("sampleE#1#chr1:0-100"), l.spell("sampleA#1#chr1:0-100")}, "unique");
    write(path("curate_E.tsv"), "sampleE#1#chr1:0-100\n");
    const std::string marker = l.nodes.at("11").substr(200, 120);   // E-only sequence, in no retained template
    // Without curation the allele is callable and appears: the fixture is sensitive. Its length and source walk (in the
    // model's block-local naming) are read from that run.
    require(run("unique", "unique", "off", "") == 0, "uncurated run must call");
    require(all_outputs("unique.off").find(marker) != std::string::npos, "uncurated: E's allele must be reachable");
    const auto uncurated = sv1_alleles("unique.off");
    require(uncurated[1].first != uncurated[0].first, "uncurated: the call must carry E's allele beside A's");
    const std::string e_bp = "\t" + uncurated[1].first + "\t", e_walk = uncurated[1].second;
    for (const std::string mode : {"off", "audit", "repair"}) {
        const int exit_status = run("unique", "unique", mode, "curate_E.tsv");
        const std::string prefix = "unique." + mode + ".curate_E";
        require(exit_status == 0 || exit_status == 1, mode + ": the curated run must finish (call or no-call)");
        const auto st = status_of(prefix);
        require(st.at("templates_curated_out") == "1" && st.at("curation_inactive_canonical_alleles") == "1" &&
                    st.at("curation_inactive_source_alleles") == "1",
                mode + ": the status must record one inactive allele");
        const std::string outputs = all_outputs(prefix);
        require(outputs.find(marker) == std::string::npos, mode + ": E's sequence must appear in no output file");
        require(outputs.find(e_walk) == std::string::npos, mode + ": E's source walk must appear in no provenance field");
        require(outputs.find(e_bp) == std::string::npos, mode + ": no candidate, repair or call may carry E's allele length");
    }
    std::cout << "  A unique excluded allele: absent from baseline, proposal, repair, dosage candidates, application and "
                 "provenance in off/audit/repair (reachable without curation)\n";
}

// B. Retained R carries a sequence-identical allele through a different walk (one node with the same bytes).
void test_shared_allele() {
    Locus l = base_locus();
    l.nodes["12"] = l.nodes.at("2") + l.nodes.at("11");
    l.templates.push_back({"sampleE#1#chr1:0-100", {"1", "2", "11", "4", "5", "7", "8", "10"}});
    l.templates.push_back({"sampleR#1#chr1:0-100", {"1", "12", "4", "6", "7", "9", "10"}});
    write_locus(l, path("shared.gfa"));
    simulate({l.spell("sampleR#1#chr1:0-100"), l.spell("sampleA#1#chr1:0-100")}, "shared");
    write(path("curate_R.tsv"), "sampleR#1#chr1:0-100\n");
    require(run("shared", "shared", "off", "") == 0, "shared: the uncurated run must call");
    require(run("shared", "shared", "off", "curate_R.tsv") == 0, "shared: the R-curated run must call");
    require(run("shared", "shared", "off", "curate_E.tsv") == 0, "shared: the E-curated run must call");
    const auto st = status_of("shared.off.curate_E");
    require(st.at("curation_inactive_source_alleles") == "1" && st.at("curation_inactive_canonical_alleles") == "0",
            "shared: E's walk is inactive but its sequence stays active through R");
    // The long allele's sources: both walks uncurated, E's alone with R curated, R's alone with E curated.
    const auto long_sources = [](const std::string& prefix) {
        const auto a = sv1_alleles(prefix);
        return std::stoul(a[0].first) > std::stoul(a[1].first) ? a[0] : a[1];
    };
    const auto all = long_sources("shared.off"), e_only = long_sources("shared.off.curate_R"), r_only = long_sources("shared.off.curate_E");
    require(all.first == r_only.first && std::stoul(r_only.first) > 1000, "shared: the long allele must stay callable");
    require(all.second != r_only.second && all.second.find(r_only.second) != std::string::npos,
            "shared: the curated run's provenance must be a strict part of the uncurated one");
    require(r_only.second.find(e_only.second) == std::string::npos && all.second.find(e_only.second) != std::string::npos,
            "shared: the allele's provenance must name R's source, not E's (" + r_only.second + " vs E " + e_only.second + ")");
    std::cout << "  B shared allele: callable, provenance names the retained source only\n";
}

// C. Two curated templates are the only carriers of the allele: it becomes inactive.
void test_all_carriers_excluded() {
    Locus l = base_locus();
    l.templates.push_back({"sampleE#1#chr1:0-100", {"1", "2", "11", "4", "5", "7", "8", "10"}});
    l.templates.push_back({"sampleF#1#chr1:0-100", {"1", "2", "11", "4", "6", "7", "9", "10"}});
    write_locus(l, path("multi.gfa"));
    simulate({l.spell("sampleE#1#chr1:0-100"), l.spell("sampleA#1#chr1:0-100")}, "multi");
    write(path("curate_EF.tsv"), "sampleE#1#chr1:0-100\nsampleF#1#chr1:0-100\n");
    write(path("curate_FE.tsv"), "# reordered\nsampleF#1#chr1:0-100\n\nsampleE#1#chr1:0-100\n");
    write(path("curate_E_only.tsv"), "sampleE#1#chr1:0-100\n");
    require(run("multi", "multi", "repair", "curate_E_only.tsv") <= 1, "one of two carriers curated must run");
    require(status_of("multi.repair.curate_E_only").at("curation_inactive_canonical_alleles") == "0",
            "with F retained the allele stays active");
    require(run("multi", "multi", "repair", "curate_EF.tsv") <= 1, "every carrier curated must run");
    const auto st = status_of("multi.repair.curate_EF");
    require(st.at("curation_inactive_canonical_alleles") == "1", "every carrier curated: the allele is inactive");
    require(all_outputs("multi.repair.curate_EF").find(l.nodes.at("11").substr(200, 120)) == std::string::npos,
            "every carrier curated: the allele's sequence appears in no output");
    // E. The manifest's order and comments do not matter: identical scientific outputs and digest.
    require(run("multi", "multi", "repair", "curate_FE.tsv") <= 1, "reordered manifest must run");
    for (const char* suffix : {".genotypes.tsv", ".block_genotypes.tsv", ".mosaic_scores.tsv", ".mosaic_candidates.tsv",
                               ".equivalence_set.tsv", ".dosage_blocks.tsv", ".dosage_repair_candidates.tsv"}) {
        require(slurp(path(std::string("multi.repair.curate_EF") + suffix)) == slurp(path(std::string("multi.repair.curate_FE") + suffix)),
                std::string("reordered manifest must give an identical ") + suffix);
    }
    require(status_of("multi.repair.curate_FE").at("curation_digest") == st.at("curation_digest"), "same curation digest");
    std::cout << "  C all carriers excluded: inactive; E manifest order: identical outputs and digest\n";
}

// G. An active panel too small to call refuses clearly and never falls back to curated sequence.
void test_insufficient_active_panel() {
    write(path("curate_all.tsv"), "sampleA#1#chr1:0-100\nsampleB#1#chr1:0-100\nsampleC#1#chr1:0-100\nsampleD#1#chr1:0-100\nsampleE#1#chr1:0-100\n");
    write(path("curate_four.tsv"), "sampleB#1#chr1:0-100\nsampleC#1#chr1:0-100\nsampleD#1#chr1:0-100\nsampleE#1#chr1:0-100\n");
    for (const std::string manifest : {"curate_all.tsv", "curate_four.tsv"}) {
        const int exit_status = run("unique", "unique", "repair", manifest);
        const std::string prefix = "unique.repair." + std::filesystem::path(manifest).stem().string();
        require(exit_status == 1, manifest + ": an insufficient active panel must no-call, got " + std::to_string(exit_status));
        const auto st = status_of(prefix);
        require(st.at("status") == "MODEL_INCOMPLETE", manifest + ": the refusal must be a model refusal: " + st.at("status"));
        require(!std::filesystem::exists(path(prefix + ".genotypes.tsv")), manifest + ": no genotype may be written");
        require(st.count("note") && !st.at("note").empty(), manifest + ": the refusal must say why");
    }
    std::cout << "  G insufficient active panel: MODEL_INCOMPLETE with a reason, no genotype, no curated sequence\n";
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: genotype_curation_contract_test <scratch dir>\n";
        return 2;
    }
    dir = argv[1];
    std::filesystem::create_directories(dir);
    try {
        write(path("library.txt"), kLibrary);
        test_unique_excluded_allele();
        test_shared_allele();
        test_all_carriers_excluded();
        test_insufficient_active_panel();
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
    std::cout << "genotype curation contract: all gates passed\n";
    return 0;
}
