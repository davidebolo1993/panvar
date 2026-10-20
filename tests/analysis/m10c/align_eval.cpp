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
// M10-C Amendment 1: `m10c_align_eval --minimap2 <manifest> <out>` runs the existing panvar::minimap2_hits (asm20, all
// regions: primary, secondary and supplementary, unflagged) in both roles (A on B, and B on A) and writes every hit:
//   strand, query start/end, target start/end, matches, block length, and CIGAR I/D ops >= 50 bp as op:len@target_pos
// M10-D: `m10c_align_eval --distance <manifest> <out>` calls panvar::nw_edit_distance(A, B) once per pair (forward only) and
// writes id, len_a, len_b, edits, aln_len -- the exact global edit distance and the benchmark's alignment-column count S.
// M10-C Amendment 2: `m10c_align_eval --poa <manifest> <out>` aligns A and B with the existing panvar::poa_msa (abPOA,
// global, affine gaps) as a two-row alignment, forward and on plain-reversed sequences, and writes the same measurements as
// the default mode from the gapped rows (gap runs >= 50 bp, terminal gap runs, mismatches, located-site bases).
#include "panvar/align.hpp"
#include "panvar/minimap2_align.hpp"
#include "panvar/poa.hpp"

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

int distance_mode(const std::string& manifest, const std::string& out_path) {
    std::ifstream in(manifest);
    std::ofstream out(out_path);
    out << "id\tlen_a\tlen_b\tedits\taln_len\n";
    for (std::string line; std::getline(in, line);) {
        if (line.empty()) continue;
        std::vector<std::string> f;
        std::stringstream ss(line);
        for (std::string x; std::getline(ss, x, '\t');) f.push_back(x);
        const std::string a = slurp(f[1]), b = slurp(f[2]);
        const panvar::NwAlign nw = panvar::nw_edit_distance(a, b);
        out << f[0] << '\t' << a.size() << '\t' << b.size() << '\t' << nw.edits << '\t' << nw.aln_len << '\n';
        out.flush();
    }
    return 0;
}

int poa_mode(const std::string& manifest, const std::string& out_path) {
    std::ifstream in(manifest);
    std::ofstream out(out_path);
    out << "id\tdirection\tlen_a\tlen_b\tcolumns\tmismatches\tins_bp\tdel_bp\tgap_runs_ge50\tmax_gap_run\tterminal_left\tterminal_right\tsite_results\n";
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
            const std::vector<std::string> rows = panvar::poa_msa({aa, bb});
            const std::string& ra = rows.at(0);
            const std::string& rb = rows.at(1);
            std::size_t mism = 0, ins = 0, del = 0, maxrun = 0, ib = 0, run = 0, run_b = 0;
            char op = 0;
            std::vector<std::string> runs;
            std::vector<char> at_b(bb.size(), '?');
            auto close = [&]() {
                if (!op) return;
                maxrun = std::max(maxrun, run);
                if (run >= 50) {
                    std::size_t pos = rev ? bb.size() - run_b - (op == 'D' ? run : 0) : run_b;
                    runs.push_back(std::string(1, op) + ":" + std::to_string(run) + "@" + std::to_string(pos));
                }
                op = 0;
                run = 0;
            };
            for (std::size_t c = 0; c < ra.size(); ++c) {
                const bool ga = ra[c] == '-', gb = rb[c] == '-';
                if (ga && gb) continue;
                const char cur = gb ? 'I' : ga ? 'D' : 0;
                if (cur != op) {
                    close();
                    if (cur) { op = cur; run_b = ib; }
                }
                if (cur) ++run;
                if (!ga && !gb) {
                    if (ra[c] != rb[c]) ++mism;
                    at_b[ib++] = ra[c];
                } else if (gb) {
                    ++ins;
                } else {
                    ++del;
                    at_b[ib++] = '-';
                }
            }
            close();
            auto terminal = [&](bool from_start) {
                std::size_t n = 0;
                char t = 0;
                for (std::size_t k = 0; k < ra.size(); ++k) {
                    const std::size_t c = from_start ? k : ra.size() - 1 - k;
                    const bool ga = ra[c] == '-', gb = rb[c] == '-';
                    if (ga && gb) continue;
                    const char cur = gb ? 'I' : ga ? 'D' : 0;
                    if (!cur || (t && cur != t)) break;
                    t = cur;
                    ++n;
                }
                return t ? std::string(1, t) + ":" + std::to_string(n) : std::string("0");
            };
            const std::string s1 = terminal(true), s2 = terminal(false);
            std::string rs, sr;
            for (const auto& x : runs) rs += (rs.empty() ? "" : ";") + x;
            for (const std::size_t p : sites) {
                const std::size_t q = rev ? bb.size() - 1 - p : p;
                sr += (sr.empty() ? "" : ";") + std::to_string(p) + ":" + bb[q] + ":" + at_b[q];
            }
            out << f[0] << '\t' << (rev ? "reverse" : "forward") << '\t' << a.size() << '\t' << b.size() << '\t' << ra.size() << '\t' << mism << '\t'
                << ins << '\t' << del << '\t' << (rs.empty() ? "-" : rs) << '\t' << maxrun << '\t' << (rev ? s2 : s1) << '\t' << (rev ? s1 : s2)
                << '\t' << (sr.empty() ? "-" : sr) << '\n';
        }
    }
    return 0;
}

int minimap2_mode(const std::string& manifest, const std::string& out_path) {
    std::ifstream in(manifest);
    std::ofstream out(out_path);
    out << "id\trole\thit\tstrand\tqs\tqe\tts\tte\tmatches\tblock_len\tbig_ops\n";
    for (std::string line; std::getline(in, line);) {
        if (line.empty()) continue;
        std::vector<std::string> f;
        std::stringstream ss(line);
        for (std::string x; std::getline(ss, x, '\t');) f.push_back(x);
        const std::string a = slurp(f[1]), b = slurp(f[2]);
        for (const bool swap : {false, true}) {
            const std::string& q = swap ? b : a;
            const std::string& t = swap ? a : b;
            const auto hits = panvar::minimap2_hits("q", q, "t", t, "asm20", 64);
            if (hits.empty()) out << f[0] << '\t' << (swap ? "truth_on_pred" : "pred_on_truth") << "\tnone\t-\t-\t-\t-\t-\t-\t-\t-\n";
            for (std::size_t h = 0; h < hits.size(); ++h) {
                const auto& x = hits[h];
                std::string big;
                std::size_t tp = x.target_start_bp;
                for (const auto& [op, len] : x.cigar) {
                    if ((op == 1 || op == 2) && len >= 50) big += (big.empty() ? "" : ";") + std::string(op == 1 ? "I" : "D") + ":" + std::to_string(len) + "@" + std::to_string(tp);
                    if (op == 7 || op == 8 || op == 2) tp += static_cast<std::size_t>(len);
                }
                out << f[0] << '\t' << (swap ? "truth_on_pred" : "pred_on_truth") << '\t' << h << '\t' << (x.reverse ? '-' : '+') << '\t'
                    << x.query_start_bp << '\t' << x.query_end_bp << '\t' << x.target_start_bp << '\t' << x.target_end_bp << '\t'
                    << x.n_matches << '\t' << x.aln_block_len << '\t' << (big.empty() ? "-" : big) << '\n';
            }
        }
    }
    return 0;
}

int main(int argc, char** argv) {
    if (argc == 4 && std::string(argv[1]) == "--minimap2") return minimap2_mode(argv[2], argv[3]);
    if (argc == 4 && std::string(argv[1]) == "--poa") return poa_mode(argv[2], argv[3]);
    if (argc == 4 && std::string(argv[1]) == "--distance") return distance_mode(argv[2], argv[3]);
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
