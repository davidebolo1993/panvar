// M4-B caller gates (run_genotype on synthetic loci): constant_v1 set explicitly changes nothing but its own record rows;
// distance_ref_core_v1_experimental completes and records its geometry and vector; one and four threads agree; every
// conflicting or unsupported combination is refused explicitly (never silently constant); same-graph curation preserves the
// geometry while the curated alleles leave the callable model; a reverse-written graph and an insertion/deletion locus
// report geometry equal to an independent recomputation from the prepared model and the reference template.
#include "panvar/copying_transitions.hpp"
#include "panvar/genotype_command.hpp"
#include "panvar/genotype_locus.hpp"
#include "panvar/prepared_locus.hpp"

#include <algorithm>
#include <cmath>
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

std::size_t g_failed = 0;
void check(bool ok, const std::string& what) {
    std::cout << (ok ? "PASS\t" : "FAIL\t") << what << '\n';
    g_failed += ok ? 0 : 1;
}
std::string dir() {
    static const std::string d = [] {
        const std::string p = (std::filesystem::temp_directory_path() / "panvar_m4b_e2e").string();
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
std::map<std::string, std::string> status(const std::string& prefix) {
    std::ifstream in(prefix + ".mosaic_status.tsv");
    std::map<std::string, std::string> v;
    std::string line;
    while (std::getline(in, line)) {
        const std::size_t t = line.find('\t');
        if (t != std::string::npos) v[line.substr(0, t)] = line.substr(t + 1);
    }
    return v;
}
bool record_row(const std::string& line) {
    for (const char* k : {"copying_", "locus_geometry_", "transition_"}) if (line.rfind(k, 0) == 0) return true;
    return false;
}
// Every output file of a run, normalised: run prefix removed; volatile status rows dropped; optionally the M4-B record rows.
std::map<std::string, std::string> outputs(const std::string& prefix, bool drop_record_rows) {
    std::map<std::string, std::string> out;
    const std::string name = std::filesystem::path(prefix).filename().string();
    for (const auto& e : std::filesystem::directory_iterator(dir())) {
        const std::string f = e.path().filename().string();
        if (f.rfind(name + ".", 0) != 0) continue;
        const std::string suffix = f.substr(name.size());
        if (suffix.find("partial") != std::string::npos || suffix.find("sorted.gfa") != std::string::npos || suffix == ".repair_memory.tsv") continue;
        std::ifstream in(e.path());
        std::string text, line;
        while (std::getline(in, line)) {
            std::size_t at;
            while ((at = line.find(prefix)) != std::string::npos) line.replace(at, prefix.size(), "<prefix>");
            if (suffix == ".mosaic_status.tsv" || suffix == ".settings.tsv") {
                if (line.rfind("seconds_", 0) == 0 || line.rfind("threads_", 0) == 0 || line.rfind("resources_", 0) == 0 ||
                    line.find("footprint") != std::string::npos || line.find("rss") != std::string::npos || line.rfind("threads\t", 0) == 0) continue;
                if (drop_record_rows && record_row(line)) continue;
            }
            if (suffix == ".repair_rounds.tsv") {   // timing columns from `seconds` on
                const std::size_t s = line.find("\tseconds");
                if (s != std::string::npos && line.rfind("round", 0) == 0) line = line.substr(0, s);
                else {
                    std::size_t tabs = 0, cut = line.size();
                    for (std::size_t i = 0; i < line.size(); ++i) if (line[i] == '\t' && ++tabs == 28) { cut = i; break; }
                    line = line.substr(0, cut);
                }
            }
            text += line + '\n';
        }
        out[suffix] = text;
    }
    return out;
}
bool same(const std::map<std::string, std::string>& a, const std::map<std::string, std::string>& b, std::string& first) {
    if (a.size() != b.size()) { first = "file sets differ (" + std::to_string(a.size()) + " vs " + std::to_string(b.size()) + ")"; return false; }
    for (const auto& [k, v] : a) {
        const auto f = b.find(k);
        if (f == b.end()) { first = "missing " + k; return false; }
        if (f->second != v) { first = "differs: " + k; return false; }
    }
    return true;
}

// The M3-D valley locus (four 60-bp SNP blocks, 18 reference copies + one all-alternate template), optionally written with
// every path reversed (node sequences unchanged, steps reverse-ordered with flipped orientation).
struct Locus { std::string gfa, reference, reads1, reads2, library; };
Locus write_valley(const std::string& tag, bool reversed_paths, bool indels) {
    std::mt19937_64 rng(2026);
    Locus l;
    l.gfa = path(tag + ".gfa");
    std::vector<std::pair<std::string, std::string>> nodes;
    nodes.push_back({"L", rnd(rng, 700)});
    for (int b = 0; b < 4; ++b) {
        std::string r = rnd(rng, 60), a = r;
        a[30] = a[30] == 'A' ? 'C' : 'A';
        if (indels && b == 1) r = "";                      // insertion relative to the reference: the reference skips x1
        if (indels && b == 2) a = "";                      // deletion relative to the reference: the alternate skips r2
        nodes.push_back({"r" + std::to_string(b), r});
        nodes.push_back({"x" + std::to_string(b), a});
        nodes.push_back({"s" + std::to_string(b), rnd(rng, b == 3 ? 700 : 450)});
    }
    std::map<std::string, std::string> seq(nodes.begin(), nodes.end());
    std::ofstream out(l.gfa);
    out << "H\tVN:Z:1.0\n";
    for (const auto& n : nodes) if (!n.second.empty()) out << "S\t" << n.first << '\t' << n.second << '\n';
    std::string prev = "L";
    std::set<std::string> links;
    const auto link = [&](const std::string& a, const std::string& b) {
        if (links.insert(a + ">" + b).second) out << "L\t" << a << "\t+\t" << b << "\t+\t0M\n";
    };
    const auto walk = [&](const std::vector<int>& alt) {
        std::vector<std::string> w{"L"};
        for (int b = 0; b < 4; ++b) {
            const std::string node = (alt[b] ? "x" : "r") + std::to_string(b);
            if (!seq[node].empty()) w.push_back(node);
            w.push_back("s" + std::to_string(b));
        }
        return w;
    };
    std::vector<std::vector<std::string>> walks;
    for (int k = 0; k < 18; ++k) walks.push_back(walk({0, 0, 0, 0}));
    walks.push_back(walk({1, 1, 1, 1}));
    for (const auto& w : walks) for (std::size_t i = 1; i < w.size(); ++i) link(w[i - 1], w[i]);
    const auto write_path = [&](const std::string& name, const std::vector<std::string>& w) {
        out << "P\t" << name << '\t';
        if (!reversed_paths) {
            for (std::size_t i = 0; i < w.size(); ++i) out << (i ? "," : "") << w[i] << '+';
        } else {
            for (std::size_t i = 0; i < w.size(); ++i) out << (i ? "," : "") << w[w.size() - 1 - i] << '-';
        }
        out << "\t*\n";
    };
    for (int k = 0; k < 18; ++k) write_path("A" + std::to_string(k) + "#1#chr1:0-100", walks[k]);
    write_path("B#1#chr1:0-100", walks[18]);
    l.reference = "A0#1#chr1:0-100";
    // reads: tiled from one reference and one alternate haplotype
    std::string h1, h2;
    for (const auto& id : walks[0]) h1 += seq[id];
    for (const auto& id : walks[18]) h2 += seq[id];
    l.reads1 = path(tag + "_1.fq");
    l.reads2 = path(tag + "_2.fq");
    std::ofstream a(l.reads1), bq(l.reads2);
    std::size_t n = 0;
    for (const std::string* hap : {&h1, &h2}) {
        for (std::size_t s = 0; s + 300 <= hap->size(); s += 23) {
            ++n;
            a << "@f" << n << "/1\n" << hap->substr(s, 100) << "\n+\n" << std::string(100, 'I') << '\n';
            bq << "@f" << n << "/2\n" << revcomp(hap->substr(s + 200, 100)) << "\n+\n" << std::string(100, 'I') << '\n';
        }
    }
    l.library = path(tag + "_lib.txt");
    std::ofstream lib(l.library);
    lib << "insert_mean_bp 300\ninsert_sd_bp 30\nsupport_sigmas 4\nsubstitution_error 0.002\noutlier_weight 0.001\n"
           "fragment_start_intensity 0.1\nmarker_background 0.01\n";
    return l;
}

panvar::GenotypeCommandOptions options(const Locus& l, const std::string& name) {
    panvar::GenotypeCommandOptions o;
    o.gfa_path = l.gfa;
    o.reference_path = l.reference;
    o.reads1_path = l.reads1;
    o.reads2_path = l.reads2;
    o.library_path = l.library;
    o.out_prefix = path(name);
    o.max_mosaic_candidates = 16;
    o.phase_audit_candidates = 0;
    o.repair_rounds = 32;
    o.kmer_size = 21;
    o.syncmer_s = 11;
    o.threads = 1;
    o.quiet = true;
    return o;
}
panvar::GenotypeCommandOptions distance(panvar::GenotypeCommandOptions o, double lambda) {
    o.copying_transition = "distance_ref_core_v1_experimental";
    o.copying_length_bp = lambda;
    o.setting_sources["copying_transition"] = "test";
    o.setting_sources["copying_length_bp"] = "test";
    return o;
}
bool refused(const panvar::GenotypeCommandOptions& o, const std::string& must_contain) {
    try {
        (void)panvar::run_genotype(o);
    } catch (const std::invalid_argument& e) {
        return std::string(e.what()).find(must_contain) != std::string::npos;
    }
    return false;
}
std::vector<double> parse_list(const std::string& s) {
    std::vector<double> v;
    std::stringstream in(s);
    for (std::string x; std::getline(in, x, ',');) v.push_back(std::stod(x));
    return v;
}
// Independent recomputation of the geometry from the prepared (uncurated) model and the reference template's alleles.
std::vector<double> recomputed_distances(const panvar::GenotypeCommandOptions& o) {
    const auto locus = panvar::prepare_genotype_locus(o);
    const auto& m = locus->model();
    std::size_t ref = m.panel().template_count();
    for (std::size_t t = 0; t < m.panel().template_count(); ++t) if (m.panel().template_names[t] == locus->resolved_reference()) ref = t;
    if (ref == m.panel().template_count()) throw std::runtime_error("reference is not a template");
    std::vector<std::string> ids, ref_alleles;
    std::vector<std::vector<std::string>> frame;
    for (std::size_t b = 0; b < m.blocks().size(); ++b) {
        ids.push_back(m.blocks()[b].block_id);
        std::vector<std::string> v;
        for (const auto& a : m.blocks()[b].alleles) v.push_back(a.sequence);
        frame.push_back(v);
        ref_alleles.push_back(m.blocks()[b].alleles[m.panel().template_alleles[ref][b]].sequence);
    }
    return panvar::compute_locus_geometry(ids, m.invariant_segments(), frame, &ref_alleles, locus->resolved_reference()).distances();
}

}  // namespace

int main() {
    try {
        const Locus v = write_valley("valley", false, false);
        // E1: constant_v1 set explicitly == default, apart from its own record rows
        const auto d0 = options(v, "const_default");
        auto d1 = options(v, "const_explicit");
        d1.setting_sources["copying_transition"] = "test";
        check(panvar::run_genotype(d0) == panvar::GenotypeRunStatus::Complete && panvar::run_genotype(d1) == panvar::GenotypeRunStatus::Complete,
              "E1 constant runs complete");
        std::string why;
        check(same(outputs(d0.out_prefix, true), outputs(d1.out_prefix, true), why), "E1 explicit constant_v1 equals the default (record rows excluded) " + why);
        check(status(d0.out_prefix).count("copying_transition") == 0, "E1 the default run writes no M4-B record row");
        check(status(d1.out_prefix)["copying_transition"] == "constant_v1" && status(d1.out_prefix)["transition_vector_sha256"].size() == 64,
              "E1 explicit constant_v1 records its mode and vector digest");
        // E2/E3: distance mode completes, records, 1 vs 4 threads identical
        const auto x1 = distance(options(v, "dist_t1"), 500.0);
        auto x4 = distance(options(v, "dist_t4"), 500.0);
        x4.threads = 4;
        const auto s1 = panvar::run_genotype(x1), s4 = panvar::run_genotype(x4);
        check(s1 == panvar::GenotypeRunStatus::Complete && s4 == panvar::GenotypeRunStatus::Complete, "E2 distance runs complete");
        auto st = status(x1.out_prefix);
        check(st["copying_transition"] == "distance_ref_core_v1_experimental" && st["locus_geometry_status"] == "available" &&
                  st["transition_vector_sha256"].size() == 64 && !st["locus_geometry_distances_bp"].empty(),
              "E2 distance mode records geometry status, digest and distances");
        check(same(outputs(x1.out_prefix, false), outputs(x4.out_prefix, false), why), "E3 distance mode: one and four threads agree " + why);
        check(st["transition_vector_sha256"] != status(d1.out_prefix)["transition_vector_sha256"], "E2 distance and constant vectors differ");
        // E4: refusals (never silently constant)
        auto c1 = distance(options(v, "ref_switch"), 500.0);
        c1.setting_sources["switch_probability"] = "profile p.tsv";
        check(refused(c1, "conflicts with switch_probability"), "E4 distance + explicit switch_probability refused");
        auto c2 = distance(options(v, "ok_preset"), 500.0);
        c2.setting_sources["switch_probability"] = "preset short-read (default)";
        check(panvar::run_genotype(c2) == panvar::GenotypeRunStatus::Complete, "E4 distance + preset switch_probability allowed");
        auto c3 = distance(options(v, "ref_joint"), 500.0);
        c3.repair_joint_search = "experimental_v1";
        check(refused(c3, "joint search stays off"), "E4 distance + joint search refused");
        auto c4 = options(v, "ref_length");
        c4.copying_length_bp = 10.0;
        c4.setting_sources["copying_length_bp"] = "test";
        check(refused(c4, "refused with constant_v1"), "E4 copying_length_bp with constant_v1 refused");
        auto c5 = distance(options(v, "ref_nolength"), 500.0);
        c5.setting_sources.erase("copying_length_bp");
        check(refused(c5, "needs copying_length_bp"), "E4 distance without copying_length_bp refused");
        {
            const auto prepared = panvar::prepare_genotype_locus(options(v, "prep"));
            bool ok = false;
            try { (void)panvar::run_genotype_prepared(distance(options(v, "ref_prepared"), 500.0), *prepared, nullptr); }
            catch (const std::invalid_argument& e) { ok = std::string(e.what()).find("prepared locus") != std::string::npos; }
            check(ok, "E4 distance with a prepared locus (no geometry) refused");
            check(panvar::run_genotype_prepared(options(v, "const_prepared"), *prepared, nullptr) == panvar::GenotypeRunStatus::Complete,
                  "E4 constant mode with a prepared locus unchanged");
        }
        // E5: same-graph curation preserves geometry; the curated template's alleles leave the callable model
        {
            std::ofstream cur(path("curate_B.txt"));
            cur << "B#1#chr1:0-100\n";
        }
        auto cu = distance(options(v, "dist_curated"), 500.0);
        cu.curation_manifest_path = path("curate_B.txt");
        check(panvar::run_genotype(cu) == panvar::GenotypeRunStatus::Complete, "E5 curated distance run completes");
        auto sc = status(cu.out_prefix);
        check(sc["locus_geometry_digest"] == st["locus_geometry_digest"] && sc["locus_geometry_frame_digest"] == st["locus_geometry_frame_digest"] &&
                  sc["locus_geometry_distances_bp"] == st["locus_geometry_distances_bp"],
              "E5 curation leaves geometry, frame digest and distances unchanged");
        check(sc["model_fingerprint"] != st["model_fingerprint"] && sc["curation_inactive_canonical_alleles"] != "0" &&
                  !sc["curation_inactive_canonical_alleles"].empty(),
              "E5 the curated alleles are inactive (model fingerprint changes, inactive canonical alleles " + sc["curation_inactive_canonical_alleles"] + ")");
        // E6: geometry equals an independent recomputation (valley), and a reverse-written graph reports the same geometry
        check(parse_list(st["locus_geometry_distances_bp"]) == recomputed_distances(options(v, "recompute")), "E6 adapter geometry = independent recomputation");
        const Locus r = write_valley("valley_rev", true, false);
        const auto xr = distance(options(r, "dist_rev"), 500.0);
        check(panvar::run_genotype(xr) == panvar::GenotypeRunStatus::Complete, "E6 reverse-written graph completes");
        {   // Amendment 2: a reverse-complemented reference reverses the canonical chain, so the distances are the MIRROR (the
            // M4-B freeze expected equality, a test defect; M4-A CONTRACT section 2: reversal mirrors positions)
            auto mirrored = parse_list(status(xr.out_prefix)["locus_geometry_distances_bp"]);
            std::reverse(mirrored.begin(), mirrored.end());
            check(mirrored == parse_list(st["locus_geometry_distances_bp"]),
                  "E6 A2: reverse-written graph: distances mirror the forward ones (" + status(xr.out_prefix)["locus_geometry_distances_bp"] + " vs " +
                      st["locus_geometry_distances_bp"] + ")");
        }
        // E7: insertion and deletion blocks
        const Locus id = write_valley("indel", false, true);
        const auto xi = distance(options(id, "dist_indel"), 500.0);
        check(panvar::run_genotype(xi) == panvar::GenotypeRunStatus::Complete, "E7 insertion/deletion locus completes");
        auto si = status(xi.out_prefix);
        const auto di = parse_list(si["locus_geometry_distances_bp"]);
        bool nonneg = !di.empty();
        for (double x : di) nonneg = nonneg && x >= 0.0;
        check(si["locus_geometry_status"] == "available" && nonneg, "E7 geometry available, distances non-negative (" + si["locus_geometry_distances_bp"] + ")");
        check(di == recomputed_distances(options(id, "recompute_indel")), "E7 adapter geometry = independent recomputation");
    } catch (const std::exception& e) {
        std::cout << "ERROR\t" << e.what() << '\n';
        return 2;
    }
    std::cout << "genotype_copying_transitions_e2e\t" << (g_failed ? "FAIL" : "PASS") << '\n';
    return g_failed ? 1 : 0;
}
