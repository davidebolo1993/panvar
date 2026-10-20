// M10-C diagnostic (evaluation only; no caller reads it): whole-locus global alignment measurements between a predicted
// haplotype A (query) and a truth haplotype B (target), via the existing panvar::nw_edit_distance (edlib NW + path; the
// benchmark round-trip QV routine) for the headline numbers and one edlib path walk for event detail.
//
// Input manifest (TSV, no header): id  path_to_A  path_to_B  sites  -- A and B are single-line sequence files; sites is a
// comma list of 0-based positions in B (located class-defining sites) or '-'.
// Output (TSV, header): one row per id and direction. 'reverse' aligns reverse(A) to reverse(B) (plain reversal), which
// flips edlib's tie-breaking and so exposes placement-dependent interpretations; its positions are mapped back to B.
//   edits, aln_len, sub_threshold_bp, over_threshold_bp   from nw_edit_distance(A, B, 50) (contiguous non-match blocks)
//   mismatches, ins_bp (bases of A absent from B), del_bp (bases of B absent from A)
//   indel_runs_ge50: contiguous same-op indel runs of >= 50 bp as op:len@b_position (I = extra in A, D = missing from A)
//   max_indel_run, terminal_left / terminal_right: the indel run touching each alignment end (op:len or 0)
//   site_results: per located site pos:truth_base:aligned_A_base (or '-' for a deletion)
#include "panvar/align.hpp"

#include "edlib.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {
std::string slurp(const std::string& path) {
    std::ifstream in(path);
    std::string s, line;
    while (std::getline(in, line)) s += line;
    return s;
}

struct Walk {
    std::size_t mism = 0, ins = 0, del = 0, max_run = 0;
    std::vector<std::string> runs;
    std::string left = "0", right = "0";
    std::vector<std::string> sites;
};

Walk walk(const std::string& a, const std::string& b, const std::vector<std::size_t>& sites, bool reversed) {
    Walk w;
    const EdlibAlignConfig cfg = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0);
    EdlibAlignResult r = edlibAlign(a.data(), static_cast<int>(a.size()), b.data(), static_cast<int>(b.size()), cfg);
    if (r.status != EDLIB_STATUS_OK || r.alignment == nullptr) {
        edlibFreeAlignResult(r);
        throw std::runtime_error("edlib failed");
    }
    std::vector<char> at_b(b.size(), '?');   // aligned A base for every B position
    std::size_t ia = 0, ib = 0;
    char run_op = 0;
    std::size_t run_len = 0, run_b = 0;
    auto close = [&]() {
        if (!run_op) return;
        const std::string tag = std::string(1, run_op) + ":" + std::to_string(run_len);
        w.max_run = std::max(w.max_run, run_len);
        if (run_len >= 50) {
            std::size_t pos = run_b;
            if (reversed) pos = b.size() - run_b - (run_op == 'D' ? run_len : 0);
            w.runs.push_back(tag + "@" + std::to_string(pos));
        }
        run_op = 0;
        run_len = 0;
    };
    for (int k = 0; k < r.alignmentLength; ++k) {
        const unsigned char op = r.alignment[k];
        const char cur = op == EDLIB_EDOP_INSERT ? 'I' : op == EDLIB_EDOP_DELETE ? 'D' : 0;
        if (cur != run_op) {
            close();
            if (cur) {
                run_op = cur;
                run_b = ib;
            }
        }
        if (cur) ++run_len;
        if (op == EDLIB_EDOP_MATCH || op == EDLIB_EDOP_MISMATCH) {
            if (op == EDLIB_EDOP_MISMATCH) ++w.mism;
            at_b[ib] = a[ia];
            ++ia;
            ++ib;
        } else if (op == EDLIB_EDOP_INSERT) {
            ++w.ins;
            ++ia;
        } else {
            ++w.del;
            at_b[ib] = '-';
            ++ib;
        }
    }
    // terminal runs: recompute from the op string ends
    auto end_run = [&](bool from_start) {
        int k = from_start ? 0 : r.alignmentLength - 1;
        const int step = from_start ? 1 : -1;
        const unsigned char op = r.alignment[k];
        if (op != EDLIB_EDOP_INSERT && op != EDLIB_EDOP_DELETE) return std::string("0");
        std::size_t n = 0;
        while (k >= 0 && k < r.alignmentLength && r.alignment[k] == op) { ++n; k += step; }
        return std::string(op == EDLIB_EDOP_INSERT ? "I:" : "D:") + std::to_string(n);
    };
    close();
    const std::string s = end_run(true), e = end_run(false);
    w.left = reversed ? e : s;
    w.right = reversed ? s : e;
    for (const std::size_t p : sites) {
        const std::size_t q = reversed ? b.size() - 1 - p : p;
        const char tb = b[q];
        w.sites.push_back(std::to_string(p) + ":" + tb + ":" + at_b[q]);
    }
    edlibFreeAlignResult(r);
    return w;
}
}  // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "usage: m10c_align_eval <manifest.tsv> <out.tsv>\n";
        return 2;
    }
    std::ifstream in(argv[1]);
    std::ofstream out(argv[2]);
    out << "id\tdirection\tlen_a\tlen_b\tedits\taln_len\tsub_threshold_bp\tover_threshold_bp\tmismatches\tins_bp\tdel_bp\t"
           "indel_runs_ge50\tmax_indel_run\tterminal_left\tterminal_right\tsite_results\n";
    for (std::string line; std::getline(in, line);) {
        if (line.empty()) continue;
        std::vector<std::string> f;
        std::stringstream ss(line);
        for (std::string x; std::getline(ss, x, '\t');) f.push_back(x);
        const std::string a = slurp(f[1]), b = slurp(f[2]);
        std::vector<std::size_t> sites;
        if (f.size() > 3 && f[3] != "-") {
            std::stringstream sp(f[3]);
            for (std::string x; std::getline(sp, x, ',');) sites.push_back(std::stoul(x));
        }
        for (const bool rev : {false, true}) {
            std::string aa = a, bb = b;
            if (rev) {
                std::reverse(aa.begin(), aa.end());
                std::reverse(bb.begin(), bb.end());
            }
            const panvar::NwAlign nw = panvar::nw_edit_distance(aa, bb, 50);
            std::vector<std::size_t> s2 = sites;   // positions are always given in forward B coordinates
            const Walk w = walk(aa, bb, s2, rev);
            std::string runs, sres;
            for (const auto& x : w.runs) runs += (runs.empty() ? "" : ";") + x;
            for (const auto& x : w.sites) sres += (sres.empty() ? "" : ";") + x;
            out << f[0] << '\t' << (rev ? "reverse" : "forward") << '\t' << a.size() << '\t' << b.size() << '\t' << nw.edits << '\t'
                << nw.aln_len << '\t' << nw.sub_threshold_bp << '\t' << nw.over_threshold_bp << '\t' << w.mism << '\t' << w.ins << '\t'
                << w.del << '\t' << (runs.empty() ? "-" : runs) << '\t' << w.max_run << '\t' << w.left << '\t' << w.right << '\t'
                << (sres.empty() ? "-" : sres) << '\n';
        }
    }
    return 0;
}
