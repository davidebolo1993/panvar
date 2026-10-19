// M6-E diagnostic (not a ctest; panel sequences only): global alignment events of declared component pairs.
//
// For each declared pair <label> <fasta q> <record q> <fasta r> <record r>, the query component is aligned to the
// reference component with the same unbanded global (NW) unit-cost edlib alignment that the evaluator's distances use,
// in path mode, and every difference is written relative to the reference coordinate:
//   SNV  one mismatched column (ref_pos, ref base, query base)
//   INS  a maximal run of query bases absent from the reference, inserted before ref_pos
//   DEL  a maximal run of reference bases absent from the query, starting at ref_pos
// Bases are written for runs of at most 50 bp (longer runs carry their length only). One SUMMARY row per pair.
// A unit-cost optimal alignment is one of possibly several co-optimal alignments; indel placement inside repeats is
// therefore not unique, and the analysis treats indel neighbourhoods as uncertain correspondence.
//
//   unit_correspondence --pairs <tsv> --out <tsv> [--threads 4]
#include <edlib.h>

#include <atomic>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

std::map<std::string, std::string> read_fasta(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot read " + path);
    std::map<std::string, std::string> out;
    std::string line, name;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            name = line.substr(1, line.find_first_of(" \t") == std::string::npos ? std::string::npos : line.find_first_of(" \t") - 1);
            if (out.count(name)) throw std::runtime_error("duplicate record " + name + " in " + path);
            out[name];
        } else {
            out[name] += line;
        }
    }
    return out;
}

struct Pair {
    std::string label, fasta_q, record_q, fasta_r, record_r;
};

std::string align_events(const Pair& p, const std::string& q, const std::string& r) {
    const EdlibAlignConfig cfg = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0);
    EdlibAlignResult res = edlibAlign(q.data(), static_cast<int>(q.size()), r.data(), static_cast<int>(r.size()), cfg);
    if (res.status != EDLIB_STATUS_OK || res.editDistance < 0 || res.alignment == nullptr) {
        edlibFreeAlignResult(res);
        throw std::runtime_error("edlib failed on " + p.label);
    }
    std::ostringstream o;
    const std::string head = p.label + '\t' + p.record_q + '\t' + p.record_r + '\t' + std::to_string(r.size()) + '\t' + std::to_string(q.size()) + '\t' +
                             std::to_string(res.editDistance) + '\t';
    std::size_t qi = 0, ri = 0, snv = 0, ins = 0, del = 0, ins_bp = 0, del_bp = 0, matched = 0;
    const int n = res.alignmentLength;
    for (int k = 0; k < n;) {
        const unsigned char op = res.alignment[k];
        if (op == EDLIB_EDOP_MATCH) { ++qi; ++ri; ++matched; ++k; continue; }
        if (op == EDLIB_EDOP_MISMATCH) {
            o << head << "SNV\t" << ri << "\t1\t" << qi << "\t1\t" << r[ri] << '\t' << q[qi] << '\n';
            ++qi; ++ri; ++snv; ++k;
            continue;
        }
        int e = k;
        while (e < n && res.alignment[e] == op) ++e;
        const std::size_t len = static_cast<std::size_t>(e - k);
        if (op == EDLIB_EDOP_INSERT) {   // query bases absent from the reference
            o << head << "INS\t" << ri << "\t0\t" << qi << '\t' << len << "\t-\t" << (len <= 50 ? q.substr(qi, len) : "*") << '\n';
            qi += len; ++ins; ins_bp += len;
        } else {                         // reference bases absent from the query
            o << head << "DEL\t" << ri << '\t' << len << '\t' << qi << "\t0\t" << (len <= 50 ? r.substr(ri, len) : "*") << "\t-\n";
            ri += len; ++del; del_bp += len;
        }
        k = e;
    }
    edlibFreeAlignResult(res);
    if (qi != q.size() || ri != r.size()) throw std::runtime_error("alignment does not span both sequences: " + p.label);
    o << head << "SUMMARY\t" << matched << '\t' << snv << '\t' << ins << '\t' << ins_bp << '\t' << del << '\t' << del_bp << '\n';
    return o.str();
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::map<std::string, std::string> a;
        for (int i = 1; i + 1 < argc; i += 2) a[argv[i]] = argv[i + 1];
        if (!a.count("--pairs") || !a.count("--out")) {
            std::cerr << "usage: unit_correspondence --pairs <tsv> --out <tsv> [--threads 4]\n";
            return 2;
        }
        const std::size_t threads = a.count("--threads") ? std::stoul(a["--threads"]) : 4;
        std::vector<Pair> pairs;
        std::map<std::string, std::map<std::string, std::string>> fastas;
        std::ifstream in(a["--pairs"]);
        if (!in) throw std::runtime_error("cannot read " + a["--pairs"]);
        for (std::string line; std::getline(in, line);) {
            if (line.empty() || line[0] == '#') continue;
            std::vector<std::string> f;
            std::stringstream ss(line);
            for (std::string x; std::getline(ss, x, '\t');) f.push_back(x);
            if (f.size() != 5) throw std::runtime_error("each pair needs <label> <fasta q> <record q> <fasta r> <record r>");
            for (const auto& fa : {f[1], f[3]}) if (!fastas.count(fa)) fastas[fa] = read_fasta(fa);
            if (!fastas[f[1]].count(f[2]) || !fastas[f[3]].count(f[4])) throw std::runtime_error("missing record in " + line);
            pairs.push_back({f[0], f[1], f[2], f[3], f[4]});
        }
        std::vector<std::string> out(pairs.size());
        std::atomic<std::size_t> next{0};
        std::vector<std::thread> pool;
        std::vector<std::string> errors(threads);
        for (std::size_t t = 0; t < threads; ++t)
            pool.emplace_back([&, t]() {
                try {
                    for (std::size_t i = next++; i < pairs.size(); i = next++)
                        out[i] = align_events(pairs[i], fastas.at(pairs[i].fasta_q).at(pairs[i].record_q), fastas.at(pairs[i].fasta_r).at(pairs[i].record_r));
                } catch (const std::exception& e) { errors[t] = e.what(); }
            });
        for (auto& th : pool) th.join();
        for (const auto& e : errors) if (!e.empty()) throw std::runtime_error(e);
        std::ofstream o(a["--out"]);
        o << "label\tquery\treference\tref_bp\tquery_bp\tdistance\tevent\tref_pos\tref_len\tquery_pos\tquery_len\tref_bases\tquery_bases\n";
        for (const auto& s : out) o << s;
        std::cout << "unit_correspondence: " << pairs.size() << " pairs on " << threads << " threads\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "unit_correspondence: ERROR " << e.what() << '\n';
        return 2;
    }
}
