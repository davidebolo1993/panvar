#include "panvar/benchmark_command.hpp"

#include "panvar/align.hpp"
#include "panvar/bubble_path.hpp"
#include "panvar/bubbles.hpp"
#include "panvar/cli_utils.hpp"
#include "panvar/gfa.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/output.hpp"
#include "panvar/parallel.hpp"
#include "panvar/ref_path.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace panvar {
namespace {

void print_benchmark_help() {
    std::cout
        << "Usage:\n"
        << "  panvar benchmark -i <graph.gfa> (-b <prefix> | -c <bubbles.csv>) \\\n"
        << "                   --variant-nodes <call.variant_nodes.tsv> -r <name> -o <prefix> [options]\n\n"
        << "Two round-trip scores of the caller's own output, over the same bubbles, the same truth and\n"
        << "the same QV scale, differing only in what the reconstruction is allowed to use.\n\n"
        << "  graph     Reference walk with the haplotype's OWN steps substituted at blocks the calls\n"
        << "            explain (nodes in variant_nodes.tsv), uncalled variation left at reference. It\n"
        << "            asks whether the graph can represent this haplotype and whether the caller\n"
        << "            flagged the divergent blocks. No genotype is read, so it cannot be wrong about\n"
        << "            WHICH haplotype carries what -- it is an upper bound, not a genotyping score.\n"
        << "  genotype  (--vcf) Reference sequence with only the edits this haplotype's GT names,\n"
        << "            applied from the VCF alone. A missed carrier keeps reference and a spurious one\n"
        << "            edits sequence that was already correct, so both error directions cost bases.\n"
        << "            This is what a consumer reconstructing a sample from the VCF actually gets.\n\n"
        << "Both align the reconstruction to the haplotype's true walk (edlib NW) and score\n"
        << "QV = -10*log10(max(0.5, delta)/S), combined per haplotype length-weighted (sum delta / sum S)\n"
        << "and binned into the cosigt bands (<17, 17-23, 23-33, >33). Uncalled variation keeps delta > 0.\n"
        << "A large graph-minus-genotype gap is a genotyping loss; both low is a calling or graph loss.\n\n"
        << "Options:\n"
        << "  -i, --gfa <path>                 Passed graph the calls were made on (required)\n"
        << "  -b, --bubble-prefix-in <prefix>  panphorte output prefix (auto-uses <prefix>.bubbles.csv)\n"
        << "  -c, --bubbles-csv-in <path>      panphorte bubbles CSV (required if no prefix)\n"
        << "      --variant-nodes <path>       call's <prefix>.variant_nodes.tsv (required)\n"
        << "      --vcf <path>                 call's <prefix>.region.vcf. Adds the genotype-aware\n"
        << "                                   reconstruction above; without it only `graph` is scored\n"
        << "  -r, --reference-path <name>      Reference path name/substring, the diff baseline (required)\n"
        << "  -o, --out-prefix <prefix>        Output prefix for the QV tables (required)\n"
        << "      --all-bubbles                Score every bubble in the CSV, not just called ones. A\n"
        << "                                   haplotype diverging from reference at an uncalled bubble\n"
        << "                                   then lands in a low band -- a missed call. Default scores\n"
        << "                                   only bubbles that carry >=1 call.\n"
        << "      --dup-model <cn|cnbp>        How a copy-number record is reconstructed: lay down CN\n"
        << "                                   copies of RU_LEN in place of REF_CN (cn), or apply\n"
        << "                                   the per-sample CNBP bp delta (cnbp, default)\n"
        << "      --min-sv-bp <N>              Threshold for the residual split (match the call run;\n"
        << "                                   default 50): residual blocks < N bp are sub-threshold\n"
        << "                                   variation, >= N bp are callable-size misses.\n"
        << "      --threads <N>                Worker threads over haplotypes (0 = auto)\n"
        << "  -q, --quiet                      Disable progress logs\n"
        << "  -h, --help                       Show this help\n";
}

// cosigt/locityper QV bands.
const char* qv_band(double qv) {
    if (qv < 17.0) return "<17";
    if (qv < 23.0) return "17-23";
    if (qv < 33.0) return "23-33";
    return ">33";
}

// QV = -10*log10(max(0.5, delta)/S). delta<=0.5 (a perfect/near-perfect match) is floored so the score
// stays finite and grows with region length, matching cosigt's convention.
double qv_from(std::size_t delta, std::size_t aln_len) {
    if (aln_len == 0) return 60.0;  // nothing to align; treat as perfect
    const double num = std::max(0.5, static_cast<double>(delta));
    return -10.0 * std::log10(num / static_cast<double>(aln_len));
}

// Ceiling QV for a perfect (delta<=0.5) reconstruction of total aligned length S: QV_max = 10*log10(2S).
// A short region can never reach the high cosigt bands, so we also report qv / qv_max -- length-fair,
// 1.0 for a perfect reconstruction of any size -- and bin that ratio into quintiles.
double qv_max_from(std::size_t aln_len) {
    if (aln_len == 0) return 60.0;
    return 10.0 * std::log10(2.0 * static_cast<double>(aln_len));
}

const char* qv_ratio_quintile(double ratio) {
    if (ratio < 0.2) return "0.0-0.2";
    if (ratio < 0.4) return "0.2-0.4";
    if (ratio < 0.6) return "0.4-0.6";
    if (ratio < 0.8) return "0.6-0.8";
    return "0.8-1.0";
}

const PathRecord* resolve_reference(const Graph& graph, const std::string& wanted) {
    for (const PathRecord& p : graph.paths) {
        if (p.name == wanted) return &p;
    }
    auto lower = [](std::string s) { for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c))); return s; };
    const std::string needle = lower(wanted);
    std::vector<const PathRecord*> hits;
    for (const PathRecord& p : graph.paths) {
        if (lower(p.name).find(needle) != std::string::npos) hits.push_back(&p);
    }
    if (hits.size() == 1) return hits.front();
    if (hits.empty()) throw std::runtime_error("Reference path not found in GFA: " + wanted);
    throw std::runtime_error("Reference path '" + wanted + "' is ambiguous (" +
                             std::to_string(hits.size()) + " matches)");
}

struct CalledVariants {
    std::unordered_set<std::size_t> bubble_ids;                                 // bubbles with >=1 call
    std::unordered_map<std::size_t, std::unordered_set<std::string>> nodes;     // bubble_id -> called nodes
    std::map<std::size_t, std::set<std::string>> svtypes;                       // bubble_id -> svtypes
};

CalledVariants load_variant_nodes(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Failed to open variant nodes: " + path);
    CalledVariants cv;
    std::string line;
    bool header = true;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (header) { header = false; continue; }  // variant_id  bubble_id  svtype  node_ids
        std::vector<std::string> f;
        std::string cur;
        for (char c : line) { if (c == '\t') { f.push_back(cur); cur.clear(); } else cur += c; }
        f.push_back(cur);
        if (f.size() < 4) continue;
        const std::size_t bid = static_cast<std::size_t>(std::stoull(f[1]));
        cv.bubble_ids.insert(bid);
        cv.svtypes[bid].insert(f[2]);
        auto& ns = cv.nodes[bid];
        std::string tok;
        for (char c : f[3]) { if (c == ',') { if (!tok.empty()) ns.insert(tok); tok.clear(); } else tok += c; }
        if (!tok.empty()) ns.insert(tok);
    }
    return cv;
}

std::vector<std::string> split_on(const std::string& s, char sep) {
    std::vector<std::string> f;
    std::string cur;
    for (char c : s) { if (c == sep) { f.push_back(cur); cur.clear(); } else cur += c; }
    f.push_back(cur);
    return f;
}

// One VCF record, reduced to what a reconstruction needs. Genotypes are haploid because each VCF
// sample column IS one haplotype path -- the join to the graph is by name, with no phasing to infer.
struct VcfRecord {
    std::string id;
    std::size_t bubble_id = 0;
    std::size_t pos = 0;               // 1-based, genomic
    std::size_t end = 0;
    long long svlen = 0;
    std::string svtype;
    std::string insseq;
    std::string ref_allele;                 // REF column, when the record carries explicit sequence
    std::vector<std::string> alts;          // ALT alleles, explicit sequence only (empty if symbolic)
    std::size_t ref_cn = 0;            // copy-number records: copies the reference carries
    std::size_t ru_len = 0;            // copy-number records: length of one repeat unit
    // Wide enough for a bubble with hundreds of distinct alleles: the allele VCF indexes every one of
    // them, so an 8-bit slot silently wraps past 127 and drops the edit.
    std::vector<std::int32_t> gt;      // per sample: -1 missing, 0 reference, >=1 alt
    std::vector<long long> cn;         // per sample: copies this haplotype carries (-1 absent)
    std::vector<long long> cnbp;       // per sample: bp this haplotype gains/loses
};

struct VcfData {
    std::vector<std::string> samples;
    std::vector<VcfRecord> records;
};

VcfData load_vcf(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Failed to open VCF: " + path);
    VcfData vd;
    std::string line;
    bool saw_header = false;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line.rfind("##", 0) == 0) continue;
        if (line[0] == '#') {
            const std::vector<std::string> h = split_on(line, '\t');
            for (std::size_t i = 9; i < h.size(); ++i) vd.samples.push_back(h[i]);
            saw_header = true;
            continue;
        }
        const std::vector<std::string> f = split_on(line, '\t');
        if (f.size() < 10) continue;                       // no sample columns: nothing to genotype
        VcfRecord r;
        r.id = f[2];
        r.pos = static_cast<std::size_t>(std::stoull(f[1]));
        // A record whose ALT is literal sequence is applied as a plain REF -> ALT substitution, which
        // is what the allele VCF emits: one record per bubble, every distinct allele spelled out.
        if (!f[4].empty() && f[4] != "." && f[4].find('<') == std::string::npos) {
            r.ref_allele = f[3];
            r.alts = split_on(f[4], ',');
        }
        for (const std::string& kv : split_on(f[7], ';')) {
            const std::size_t eq = kv.find('=');
            if (eq == std::string::npos) continue;
            const std::string k = kv.substr(0, eq), v = kv.substr(eq + 1);
            if (k == "BUBBLE_ID") r.bubble_id = static_cast<std::size_t>(std::stoull(v));
            else if (k == "END") r.end = static_cast<std::size_t>(std::stoull(v));
            else if (k == "SVLEN") r.svlen = std::stoll(v);
            else if (k == "SVTYPE") r.svtype = v;
            else if (k == "INSSEQ") r.insseq = v;
            else if (k == "REF_CN") r.ref_cn = static_cast<std::size_t>(std::stoull(v));
            else if (k == "RU_LEN") r.ru_len = static_cast<std::size_t>(std::stoull(v));
        }
        const std::vector<std::string> fmt = split_on(f[8], ':');
        std::size_t gt_i = fmt.size(), cnbp_i = fmt.size(), cn_i = fmt.size();
        for (std::size_t i = 0; i < fmt.size(); ++i) {
            if (fmt[i] == "GT") gt_i = i;
            else if (fmt[i] == "CNBP") cnbp_i = i;
            else if (fmt[i] == "CN") cn_i = i;
        }
        r.gt.assign(vd.samples.size(), -1);
        r.cn.assign(vd.samples.size(), -1);
        r.cnbp.assign(vd.samples.size(), 0);
        for (std::size_t s = 0; s + 9 < f.size() && s < vd.samples.size(); ++s) {
            const std::vector<std::string> v = split_on(f[s + 9], ':');
            if (gt_i < v.size() && !v[gt_i].empty() && v[gt_i] != ".")
                r.gt[s] = static_cast<std::int32_t>(std::stoi(v[gt_i]));
            if (cnbp_i < v.size() && !v[cnbp_i].empty() && v[cnbp_i] != ".")
                r.cnbp[s] = std::stoll(v[cnbp_i]);
            if (cn_i < v.size() && !v[cn_i].empty() && v[cn_i] != ".")
                r.cn[s] = std::stoll(v[cn_i]);
        }
        vd.records.push_back(std::move(r));
    }
    if (!saw_header) throw std::runtime_error("VCF has no #CHROM header line: " + path);
    return vd;
}

// How many records of each kind could actually be laid down, so a genotype QV is never read without
// knowing what it was computed from.
struct GtStats {
    std::size_t applied = 0;
    std::size_t unplaceable = 0;   // POS outside the bubble's reference span (a merged or shifted bubble)
    std::size_t clamped = 0;       // edit ran past the span end and was truncated
    std::size_t unhandled = 0;     // svtype with no reconstruction rule
    std::size_t ref_mismatch = 0;  // record REF did not match the reference at POS (placement is wrong)
};

// Apply the edits `hap` is genotyped as carrying to the bubble's reference sequence. Nothing from the
// haplotype's own walk enters here -- that is the whole point of the second metric. Edits are laid down
// right to left so earlier ones do not shift the coordinates of later ones.
std::string apply_genotype(const std::string& ref_bubble,
                           std::size_t bubble_start_0,          // 0-based bp offset of ref_bubble in the reference path
                           std::size_t region_start_1based,
                           const std::vector<const VcfRecord*>& recs,
                           std::size_t hap,
                           bool cn_model,
                           GtStats& stats) {
    std::vector<const VcfRecord*> carried;
    for (const VcfRecord* r : recs)
        if (hap < r->gt.size() && r->gt[hap] >= 1) carried.push_back(r);
    std::sort(carried.begin(), carried.end(),
              [](const VcfRecord* a, const VcfRecord* b) { return a->pos > b->pos; });

    std::string out = ref_bubble;
    for (const VcfRecord* r : carried) {
        // POS is the anchor base; the event occupies POS+1 onwards. Everything is expressed relative to
        // the bubble's own reference sequence so the two metrics share a denominator.
        if (r->pos < region_start_1based) { ++stats.unplaceable; continue; }
        const std::size_t g = r->pos - region_start_1based;   // 0-based in the reference path
        if (g < bubble_start_0 || g >= bubble_start_0 + out.size()) { ++stats.unplaceable; continue; }
        const std::size_t at = g - bubble_start_0 + 1;        // first edited base, 0-based in `out`
        if (at > out.size()) { ++stats.unplaceable; continue; }

        if (!r->alts.empty()) {
            // Explicit-sequence allele: replace the record's REF span with the allele the GT names.
            // `at` is one past the anchor, and REF starts AT the anchor, so back up one base.
            const std::size_t start = at - 1;
            const std::int32_t gi = r->gt[hap];
            if (gi < 1 || static_cast<std::size_t>(gi) > r->alts.size()) { ++stats.unhandled; continue; }
            std::size_t n = r->ref_allele.size();
            if (start + n > out.size()) { n = out.size() - start; ++stats.clamped; }
            // The record's REF must be the reference sequence at POS. If it is not, the record is
            // being laid down in the wrong place and every score downstream is fiction.
            if (out.compare(start, n, r->ref_allele, 0, n) != 0) ++stats.ref_mismatch;
            out.erase(start, n);
            out.insert(start, r->alts[static_cast<std::size_t>(gi) - 1]);
            ++stats.applied;
        } else if (r->svtype == "DEL") {
            std::size_t n = static_cast<std::size_t>(r->svlen < 0 ? -r->svlen : r->svlen);
            if (at + n > out.size()) { n = out.size() - at; ++stats.clamped; }
            out.erase(at, n);
            ++stats.applied;
        } else if (r->svtype == "INS") {
            if (r->insseq.empty()) { ++stats.unhandled; continue; }
            out.insert(at, r->insseq);
            ++stats.applied;
        } else if (r->svtype == "DUP") {
            // A copy-number record. The reference carries REF_CN copies of a RU_LEN unit here and this
            // haplotype carries CN of them, so lay down CN copies in place of REF_CN -- the array's
            // LENGTH is what this record is about, and it is what drives the score.
            const std::size_t avail = out.size() - at;
            const long long cn = hap < r->cn.size() ? r->cn[hap] : -1;
            if (cn_model && cn >= 0 && r->ru_len > 0 && r->ref_cn > 0) {
                const std::size_t unit_len = std::min(r->ru_len, avail);
                if (unit_len == 0) { ++stats.unhandled; continue; }
                const std::string unit = out.substr(at, unit_len);
                std::size_t drop = r->ref_cn * r->ru_len;
                if (drop > avail) { drop = avail; ++stats.clamped; }
                out.erase(at, drop);
                const std::size_t want = static_cast<std::size_t>(cn) * r->ru_len;
                std::string add;
                add.reserve(want);
                while (add.size() < want) add += unit.substr(0, std::min(unit.size(), want - add.size()));
                out.insert(at, add);
                ++stats.applied;
                continue;
            }
            // No CN or no unit length: fall back to the record's own bp delta over the annotated span.
            const long long d = hap < r->cnbp.size() ? r->cnbp[hap] : 0;
            if (d == 0) { ++stats.unhandled; continue; }
            const std::size_t unit = std::min(r->end > r->pos ? r->end - r->pos : 0, avail);
            if (unit == 0) { ++stats.unhandled; continue; }
            if (d < 0) {
                std::size_t n = static_cast<std::size_t>(-d);
                if (n > avail) { n = avail; ++stats.clamped; }
                out.erase(at, n);
            } else {
                std::string add;
                add.reserve(static_cast<std::size_t>(d));
                while (add.size() < static_cast<std::size_t>(d))
                    add += out.substr(at, std::min(unit, static_cast<std::size_t>(d) - add.size()));
                out.insert(at, add);
            }
            ++stats.applied;
        } else {
            ++stats.unhandled;
        }
    }
    return out;
}

// Reconstruct a haplotype's walk over a bubble as "reference, with only the called divergences applied".
// Node-align H against the reference walk on shared step tokens (LCS); between consecutive anchors take
// H's block when it is a called event (any node of the ref- or hap-side block is in `called`), otherwise
// revert to the reference block. Coordinate-free, so bubble orientation is irrelevant. `applied` reports
// whether at least one called block was substituted (this haplotype carries a call at the bubble).
std::vector<PathStep> reconstruct(const std::vector<PathStep>& ref_steps,
                                  const std::vector<PathStep>& hap_steps,
                                  const std::unordered_set<std::string>& called,
                                  bool& applied) {
    applied = false;
    const std::vector<std::uint64_t> rt = build_walk_tokens(ref_steps);
    const std::vector<std::uint64_t> ht = build_walk_tokens(hap_steps);
    const std::size_t n = rt.size(), m = ht.size();

    // LCS anchors on step tokens. Guard the O(nm) table on pathological walks (a huge DUP): fall back
    // to treating the whole bubble as one block.
    std::vector<std::pair<std::size_t, std::size_t>> anchors;
    if (n != 0 && m != 0 && static_cast<std::uint64_t>(n) * m <= 8000000ULL) {
        std::vector<std::vector<std::uint32_t>> dp(n + 1, std::vector<std::uint32_t>(m + 1, 0));
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < m; ++j)
                dp[i + 1][j + 1] = (rt[i] == ht[j]) ? dp[i][j] + 1
                                                    : std::max(dp[i][j + 1], dp[i + 1][j]);
        std::size_t i = n, j = m;
        while (i > 0 && j > 0) {
            if (rt[i - 1] == ht[j - 1]) { anchors.emplace_back(i - 1, j - 1); --i; --j; }
            else if (dp[i - 1][j] >= dp[i][j - 1]) --i;
            else --j;
        }
        std::reverse(anchors.begin(), anchors.end());
    }

    auto block_called = [&](std::size_t ri0, std::size_t ri1, std::size_t hi0, std::size_t hi1) {
        for (std::size_t k = ri0; k < ri1; ++k) if (called.count(ref_steps[k].node_id)) return true;
        for (std::size_t k = hi0; k < hi1; ++k) if (called.count(hap_steps[k].node_id)) return true;
        return false;
    };

    std::vector<PathStep> out;
    std::size_t pi = 0, pj = 0;
    auto emit_block = [&](std::size_t ri1, std::size_t hi1) {
        if (block_called(pi, ri1, pj, hi1)) {
            for (std::size_t k = pj; k < hi1; ++k) out.push_back(hap_steps[k]);
            if (ri1 != pi || hi1 != pj) applied = true;
        } else {
            for (std::size_t k = pi; k < ri1; ++k) out.push_back(ref_steps[k]);
        }
    };
    for (const auto& a : anchors) {
        emit_block(a.first, a.second);
        out.push_back(ref_steps[a.first]);  // shared anchor step
        pi = a.first + 1;
        pj = a.second + 1;
    }
    emit_block(n, m);
    return out;
}

struct BubbleObs {
    std::size_t bubble_id = 0;
    bool carrier = false;
    std::size_t delta = 0;
    std::size_t aln_len = 0;
    std::size_t sub_bp = 0;    // residual in blocks < min-sv-bp (variation that couldn't be called)
    std::size_t over_bp = 0;   // residual in blocks >= min-sv-bp (a callable-size event missed/mis-called)
    // Genotype-aware reconstruction over the same bubble against the same truth, plus the do-nothing
    // baseline (plain reference) it has to beat. A VCF that closes none of the gap between the baseline
    // and the graph bound is worth nothing here, and one that scores WORSE than the baseline is applying
    // edits in the wrong place -- which is the invariant that makes this number readable.
    bool gt_scored = false;
    std::size_t gt_delta = 0;
    std::size_t gt_aln_len = 0;
    std::size_t ref_delta = 0;
    std::size_t ref_aln_len = 0;
    bool gt_called_carrier = false;   // the VCF genotypes this haplotype as carrying an alt here
    bool gt_true_carrier = false;     // its walk really does diverge from the reference walk here
};

struct HapResult {
    std::string sample;
    bool scored = false;
    std::size_t sum_delta = 0;
    std::size_t sum_aln = 0;
    std::size_t sum_sub = 0;
    std::size_t sum_over = 0;
    bool gt_scored = false;
    std::size_t gt_sum_delta = 0;
    std::size_t gt_sum_aln = 0;
    std::size_t ref_sum_delta = 0;
    std::size_t ref_sum_aln = 0;
    GtStats gt_stats;
    std::vector<BubbleObs> obs;
};

} // namespace

int run_benchmark_command(const std::vector<std::string>& args) {
    if (args.empty()) { print_benchmark_help(); return 0; }

    std::string gfa_path, bubble_prefix_in, bubbles_csv_in, variant_nodes_in, vcf_in, reference_path, out_prefix;
    std::size_t threads = 0;
    std::size_t min_sv_bp = 50;   // threshold for the residual sub/over split (should match the call run)
    bool quiet = false;
    bool all_bubbles = false;
    // How a copy-number record is laid down. `cnbp` -- the per-sample bp delta -- is the default
    // because it is measured: on acot it closes 95.7% of the gap against 51.5% for `cn`, since
    // CN x RU_LEN understates the true bp change by roughly half at a cyclic array.
    bool dup_cn_model = false;

    for (std::size_t i = 0; i < args.size(); ++i) {
        const std::string& arg = args[i];
        auto require_value = [&](const std::string& flag) -> const std::string& {
            if (i + 1 >= args.size()) throw std::runtime_error("Missing value after " + flag);
            return args[++i];
        };
        if (arg == "-h" || arg == "--help") { print_benchmark_help(); return 0; }
        else if (arg == "-i" || arg == "--gfa") gfa_path = require_value(arg);
        else if (arg == "-b" || arg == "--bubble-prefix-in") bubble_prefix_in = require_value(arg);
        else if (arg == "-c" || arg == "--bubbles-csv-in") bubbles_csv_in = require_value(arg);
        else if (arg == "--variant-nodes") variant_nodes_in = require_value(arg);
        else if (arg == "--vcf") vcf_in = require_value(arg);
        else if (arg == "--dup-model") {
            const std::string v = require_value(arg);
            if (v == "cn") dup_cn_model = true;
            else if (v == "cnbp") dup_cn_model = false;
            else throw std::runtime_error("--dup-model must be cn or cnbp");
        }
        else if (arg == "-r" || arg == "--reference-path") reference_path = require_value(arg);
        else if (arg == "-o" || arg == "--out-prefix") out_prefix = require_value(arg);
        else if (arg == "--all-bubbles") all_bubbles = true;
        else if (arg == "--min-sv-bp") min_sv_bp = cli::parse_size_arg(arg, require_value(arg));
        else if (arg == "--threads") threads = cli::parse_size_arg(arg, require_value(arg));
        else if (arg == "-q" || arg == "--quiet") quiet = true;
        else throw std::runtime_error("Unknown option for benchmark: " + arg);
    }

    if (gfa_path.empty()) throw std::runtime_error("Missing required input: --gfa <path>");
    if (!bubble_prefix_in.empty()) {
        const std::string derived = bubble_prefix_in + ".bubbles.csv";
        if (bubbles_csv_in.empty()) bubbles_csv_in = derived;
        else if (bubbles_csv_in != derived)
            throw std::runtime_error("Conflicting bubble inputs: --bubble-prefix-in resolves to '" +
                                     derived + "' but --bubbles-csv-in is '" + bubbles_csv_in + "'");
    }
    if (bubbles_csv_in.empty())
        throw std::runtime_error("benchmark requires --bubble-prefix-in <prefix> or --bubbles-csv-in <path>");
    if (variant_nodes_in.empty()) throw std::runtime_error("benchmark requires --variant-nodes <path>");
    if (reference_path.empty()) throw std::runtime_error("--reference-path is required for module 'benchmark'");
    if (out_prefix.empty()) throw std::runtime_error("benchmark requires -o/--out-prefix <prefix>");

    ParseGfaOptions parse_options;
    parse_options.include_paths = true;
    parse_options.include_sequences = true;
    const Graph graph = parse_gfa(gfa_path, parse_options);
    if (graph.paths.empty())
        throw std::runtime_error("Input GFA has no P/W paths; benchmarking requires paths");

    cli::RunLog log("benchmark", quiet);
    log.info("input " + gfa_path + " (" + std::to_string(graph.nodes.size()) + " nodes, " +
             std::to_string(graph.paths.size()) + " paths); reference " + reference_path);

    const std::vector<Bubble> all_bubbles_csv = read_bubbles_csv(bubbles_csv_in);
    const CalledVariants cv = load_variant_nodes(variant_nodes_in);
    const PathRecord* ref_path = resolve_reference(graph, reference_path);

    std::vector<Bubble> bubbles;
    for (const Bubble& b : all_bubbles_csv) if (all_bubbles || cv.bubble_ids.count(b.id)) bubbles.push_back(b);
    log.info("scoring " + std::to_string(bubbles.size()) + (all_bubbles ? " bubbles (all) across " : " called bubbles across ") +
             std::to_string(graph.paths.size()) + " haplotypes");

    // Reference walk per scored bubble (the reconstruction backbone). A bubble the reference doesn't
    // traverse has no baseline, so it is dropped from scoring.
    const BubblePathIndex ref_index = build_bubble_path_index(*ref_path);
    std::unordered_map<std::size_t, std::vector<PathStep>> ref_steps;
    for (const Bubble& b : bubbles) {
        const auto s = bubble_steps(*ref_path, ref_index, b);
        if (s) ref_steps[b.id] = *s;
    }

    // Genotype-aware scoring, when a VCF is supplied. VCF POS lives in the reference path's own forward
    // coordinates, so each bubble's reference span is taken as a substring there and the reconstruction
    // is flipped back for bubbles the reference crosses sink->source.
    const bool do_gt = !vcf_in.empty();
    VcfData vcf;
    std::unordered_map<std::size_t, std::vector<const VcfRecord*>> recs_by_bubble;
    std::vector<int> hap_vcf(graph.paths.size(), -1);
    std::unordered_map<std::size_t, std::pair<std::size_t, std::size_t>> ref_span;  // bubble -> [bp start, bp len)
    std::unordered_map<std::size_t, bool> ref_flip;
    std::string ref_seq;
    std::size_t region_start = 1;
    if (do_gt) {
        vcf = load_vcf(vcf_in);
        for (const VcfRecord& r : vcf.records) recs_by_bubble[r.bubble_id].push_back(&r);
        std::unordered_map<std::string, std::size_t> sample_of;
        for (std::size_t i = 0; i < vcf.samples.size(); ++i) sample_of[vcf.samples[i]] = i;
        std::size_t joined = 0;
        for (std::size_t i = 0; i < graph.paths.size(); ++i) {
            const auto it = sample_of.find(graph.paths[i].name);
            if (it != sample_of.end()) { hap_vcf[i] = static_cast<int>(it->second); ++joined; }
        }
        if (joined == 0)
            throw std::runtime_error("No VCF sample column matches a graph path name, so no haplotype "
                                     "can be genotype-scored: " + vcf_in);
        region_start = parse_reference_path_label(ref_path->name).region_start_1based;

        std::vector<std::size_t> prefix(ref_path->steps.size() + 1, 0);
        for (std::size_t i = 0; i < ref_path->steps.size(); ++i) {
            const auto nit = graph.nodes.find(ref_path->steps[i].node_id);
            prefix[i + 1] = prefix[i] + (nit != graph.nodes.end() ? nit->second.sequence.size() : 0);
        }
        ref_seq = spell_path_steps_sequence(graph, ref_path->steps);
        for (const Bubble& b : bubbles) {
            if (!ref_steps.count(b.id)) continue;
            std::size_t lo = 0, hi = 0;
            bool flip = false, ok = false;
            const auto iv = find_best_bubble_path_interval(ref_index, b);
            if (iv.has_value()) {
                lo = std::min(iv->left, iv->right);
                hi = std::max(iv->left, iv->right);
                flip = !iv->source_to_sink;
                ok = true;
            } else {
                // The same fallback bubble_steps uses: an adjacent source/sink pair with no interior.
                const auto si = ref_index.positions.find(b.source);
                const auto ki = ref_index.positions.find(b.sink);
                if (si != ref_index.positions.end() && ki != ref_index.positions.end()) {
                    const std::unordered_set<std::size_t> sink_pos(ki->second.begin(), ki->second.end());
                    for (const std::size_t p : si->second)
                        if (sink_pos.count(p + 1)) { lo = p; hi = p + 1; ok = true; break; }
                    if (!ok) {
                        const std::unordered_set<std::size_t> src_pos(si->second.begin(), si->second.end());
                        for (const std::size_t p : ki->second)
                            if (src_pos.count(p + 1)) { lo = p; hi = p + 1; flip = true; ok = true; break; }
                    }
                }
            }
            if (!ok || hi + 1 >= prefix.size()) continue;
            ref_span[b.id] = { prefix[lo], prefix[hi + 1] - prefix[lo] };
            ref_flip[b.id] = flip;
        }
        log.info("genotype scoring from " + vcf_in + ": " + std::to_string(vcf.records.size()) +
                 " records, " + std::to_string(joined) + "/" + std::to_string(graph.paths.size()) +
                 " haplotypes joined to VCF samples, " + std::to_string(ref_span.size()) +
                 "/" + std::to_string(bubbles.size()) + " bubbles placed on the reference");
    }

    std::vector<HapResult> results(graph.paths.size());
    static const std::unordered_set<std::string> kEmpty;
    run_parallel(graph.paths.size(), threads, [&](std::size_t pi) {
        const PathRecord& path = graph.paths[pi];
        HapResult& hr = results[pi];
        hr.sample = path.name;
        if (&path == ref_path) return;
        const BubblePathIndex idx = build_bubble_path_index(path);
        for (const Bubble& b : bubbles) {
            const auto rit = ref_steps.find(b.id);
            if (rit == ref_steps.end()) continue;                // reference doesn't traverse -> no baseline
            const auto steps = bubble_steps(path, idx, b);
            if (!steps) continue;                                // haplotype doesn't traverse -> skip
            const std::string truth = spell_path_steps_sequence(graph, *steps);
            const auto nit = cv.nodes.find(b.id);
            const std::unordered_set<std::string>& called = nit != cv.nodes.end() ? nit->second : kEmpty;
            bool applied = false;
            const std::vector<PathStep> recon_steps = reconstruct(rit->second, *steps, called, applied);
            const std::string recon = spell_path_steps_sequence(graph, recon_steps);
            const NwAlign nw = nw_edit_distance(recon, truth, min_sv_bp);
            hr.scored = true;
            hr.sum_delta += nw.edits;
            hr.sum_aln += nw.aln_len;
            hr.sum_sub += nw.sub_threshold_bp;
            hr.sum_over += nw.over_threshold_bp;

            BubbleObs o;
            o.bubble_id = b.id;
            o.carrier = applied;
            o.delta = nw.edits;
            o.aln_len = nw.aln_len;
            o.sub_bp = nw.sub_threshold_bp;
            o.over_bp = nw.over_threshold_bp;

            const auto span = do_gt && hap_vcf[pi] >= 0 ? ref_span.find(b.id) : ref_span.end();
            if (span != ref_span.end()) {
                const std::size_t hv = static_cast<std::size_t>(hap_vcf[pi]);
                static const std::vector<const VcfRecord*> kNoRecs;
                const auto brit = recs_by_bubble.find(b.id);
                const std::vector<const VcfRecord*>& brecs =
                    brit != recs_by_bubble.end() ? brit->second : kNoRecs;
                const std::string ref_bubble = ref_seq.substr(span->second.first, span->second.second);
                std::string gt_recon = apply_genotype(ref_bubble, span->second.first, region_start,
                                                      brecs, hv, dup_cn_model, hr.gt_stats);
                std::string ref_recon = ref_bubble;
                if (ref_flip.at(b.id)) {
                    gt_recon = reverse_complement(gt_recon);
                    ref_recon = reverse_complement(ref_recon);
                }
                const NwAlign gnw = nw_edit_distance(gt_recon, truth, min_sv_bp);
                const NwAlign rnw = nw_edit_distance(ref_recon, truth, min_sv_bp);
                hr.gt_scored = true;
                hr.gt_sum_delta += gnw.edits;
                hr.gt_sum_aln += gnw.aln_len;
                hr.ref_sum_delta += rnw.edits;
                hr.ref_sum_aln += rnw.aln_len;
                o.gt_scored = true;
                o.gt_delta = gnw.edits;
                o.gt_aln_len = gnw.aln_len;
                o.ref_delta = rnw.edits;
                o.ref_aln_len = rnw.aln_len;
                for (const VcfRecord* r : brecs)
                    if (hv < r->gt.size() && r->gt[hv] >= 1) { o.gt_called_carrier = true; break; }
                // Truth for the carrier call: is there a callable-SIZE difference between this haplotype
                // and the reference here. Plain walk inequality would count a single SNP as a missed
                // carrier, which is not something the caller ever set out to emit.
                o.gt_true_carrier = rnw.over_threshold_bp > 0;
            }
            hr.obs.push_back(o);
        }
    });

    std::map<std::string, std::size_t> band_count;
    std::map<std::string, std::size_t> quintile_count;
    std::map<std::string, std::size_t> gt_band_count;
    std::map<std::string, std::size_t> gt_quintile_count;
    std::size_t scored_haps = 0, gt_scored_haps = 0;
    double gt_closed_sum = 0.0;
    std::size_t gt_worse_than_ref = 0;
    std::size_t gt_tot_ref_delta = 0, gt_tot_gt_delta = 0, gt_tot_graph_delta = 0;
    std::size_t tot_sub = 0, tot_over = 0;
    GtStats tot_gt;
    {
        std::ofstream by_hap(out_prefix + ".qv_by_haplotype.tsv");
        if (!by_hap) throw std::runtime_error("Failed to write " + out_prefix + ".qv_by_haplotype.tsv");
        by_hap << "sample\tn_bubbles\tsum_delta\tsum_aln_len\tqv\tband\tqv_max\tqv_ratio\tquintile\tidentity\tsub_threshold_bp\tover_threshold_bp"
                  "\tgt_sum_delta\tgt_sum_aln_len\tgt_qv\tgt_band\tgt_quintile\tgt_identity"
                  "\tref_sum_delta\tref_qv\tgap_closed\n";
        for (const HapResult& hr : results) {
            if (!hr.scored) continue;
            const double qv = qv_from(hr.sum_delta, hr.sum_aln);
            const double qmax = qv_max_from(hr.sum_aln);
            const double ratio = std::min(1.0, std::max(0.0, qmax > 0.0 ? qv / qmax : 1.0));
            // Linear reconstruction identity: 1 - error rate. Unlike QV it does not overweight a few
            // uncalled SNPs over a long region (a handful of edits in tens of kb stays ~1.0).
            const double identity = hr.sum_aln ? std::max(0.0, 1.0 - static_cast<double>(hr.sum_delta) / static_cast<double>(hr.sum_aln)) : 1.0;
            const std::string band = qv_band(qv);
            const std::string quintile = qv_ratio_quintile(ratio);
            by_hap << hr.sample << '\t' << hr.obs.size() << '\t' << hr.sum_delta << '\t'
                   << hr.sum_aln << '\t' << qv << '\t' << band << '\t' << qmax << '\t'
                   << ratio << '\t' << quintile << '\t' << identity << '\t'
                   << hr.sum_sub << '\t' << hr.sum_over;
            if (hr.gt_scored) {
                const double gqv = qv_from(hr.gt_sum_delta, hr.gt_sum_aln);
                const double gmax = qv_max_from(hr.gt_sum_aln);
                const double gratio = std::min(1.0, std::max(0.0, gmax > 0.0 ? gqv / gmax : 1.0));
                const double gident = hr.gt_sum_aln ? std::max(0.0, 1.0 - static_cast<double>(hr.gt_sum_delta) / static_cast<double>(hr.gt_sum_aln)) : 1.0;
                const std::string gband = qv_band(gqv);
                const std::string gquint = qv_ratio_quintile(gratio);
                // How much of the distance between doing nothing and the graph bound the VCF closes.
                // 1 = everything the graph could offer, 0 = no better than plain reference, < 0 = the
                // genotypes put edits where they do not belong.
                const double denom = static_cast<double>(hr.ref_sum_delta) - static_cast<double>(hr.sum_delta);
                const double closed = std::fabs(denom) > 1e-9
                    ? (static_cast<double>(hr.ref_sum_delta) - static_cast<double>(hr.gt_sum_delta)) / denom
                    : 1.0;
                by_hap << '\t' << hr.gt_sum_delta << '\t' << hr.gt_sum_aln << '\t' << gqv << '\t'
                       << gband << '\t' << gquint << '\t' << gident << '\t'
                       << hr.ref_sum_delta << '\t' << qv_from(hr.ref_sum_delta, hr.ref_sum_aln) << '\t'
                       << closed;
                gt_closed_sum += closed;
                if (hr.gt_sum_delta > hr.ref_sum_delta) ++gt_worse_than_ref;
                gt_tot_ref_delta += hr.ref_sum_delta;
                gt_tot_gt_delta += hr.gt_sum_delta;
                gt_tot_graph_delta += hr.sum_delta;
                ++gt_band_count[gband];
                ++gt_quintile_count[gquint];
                ++gt_scored_haps;
                tot_gt.applied += hr.gt_stats.applied;
                tot_gt.unplaceable += hr.gt_stats.unplaceable;
                tot_gt.clamped += hr.gt_stats.clamped;
                tot_gt.unhandled += hr.gt_stats.unhandled;
                tot_gt.ref_mismatch += hr.gt_stats.ref_mismatch;
            } else {
                by_hap << "\t.\t.\t.\t.\t.\t.\t.\t.\t.";
            }
            by_hap << '\n';
            ++band_count[band];
            ++quintile_count[quintile];
            ++scored_haps;
            tot_sub += hr.sum_sub;
            tot_over += hr.sum_over;
        }
    }

    std::map<std::string, std::map<std::string, std::size_t>> class_bands;
    // Carrier confusion, per svtype: the genotype call (does the VCF say this haplotype carries an alt
    // at this bubble) against the truth (does its walk really diverge from the reference walk).
    std::map<std::string, std::array<std::size_t, 4>> carrier;   // svtype -> {TP, FP, FN, TN}
    {
        std::ofstream qv_out(out_prefix + ".qv.tsv");
        if (!qv_out) throw std::runtime_error("Failed to write " + out_prefix + ".qv.tsv");
        qv_out << "sample\tbubble_id\tsvtypes\tis_carrier\tdelta\taln_len\tqv\tsub_threshold_bp\tover_threshold_bp"
                  "\tgt_delta\tgt_aln_len\tgt_qv\tref_delta\tgt_called_carrier\tgt_true_carrier\n";
        for (const HapResult& hr : results) {
            if (!hr.scored) continue;
            for (const BubbleObs& o : hr.obs) {
                std::string svt;
                const auto sit = cv.svtypes.find(o.bubble_id);
                if (sit != cv.svtypes.end())
                    for (const std::string& s : sit->second) { if (!svt.empty()) svt += ','; svt += s; }
                const double qv = qv_from(o.delta, o.aln_len);
                const std::string band = qv_band(qv);
                qv_out << hr.sample << '\t' << o.bubble_id << '\t' << (svt.empty() ? "." : svt) << '\t'
                       << (o.carrier ? 1 : 0) << '\t' << o.delta << '\t' << o.aln_len << '\t' << qv << '\t'
                       << o.sub_bp << '\t' << o.over_bp;
                if (o.gt_scored) {
                    qv_out << '\t' << o.gt_delta << '\t' << o.gt_aln_len << '\t'
                           << qv_from(o.gt_delta, o.gt_aln_len) << '\t' << o.ref_delta << '\t'
                           << (o.gt_called_carrier ? 1 : 0) << '\t' << (o.gt_true_carrier ? 1 : 0);
                    const std::size_t cell = o.gt_called_carrier ? (o.gt_true_carrier ? 0 : 1)
                                                                 : (o.gt_true_carrier ? 2 : 3);
                    if (sit != cv.svtypes.end())
                        for (const std::string& s : sit->second) ++carrier[s][cell];
                    ++carrier["ALL"][cell];
                } else {
                    qv_out << "\t.\t.\t.\t.\t.\t.";
                }
                qv_out << '\n';
                if (sit != cv.svtypes.end())
                    for (const std::string& s : sit->second) ++class_bands[s][band];
            }
        }
    }

    static const char* kBands[] = {"<17", "17-23", "23-33", ">33"};
    static const char* kQuintiles[] = {"0.0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0"};
    {
        std::ofstream sum(out_prefix + ".qv_summary.tsv");
        if (!sum) throw std::runtime_error("Failed to write " + out_prefix + ".qv_summary.tsv");
        sum << "scope\tkey\tband\tn\tpct\n";
        // Length-fair headline: % of haplotypes per qv/qv_max quintile (Q5 = 0.8-1.0 is best).
        for (const char* q : kQuintiles) {
            const std::size_t n = quintile_count.count(q) ? quintile_count[q] : 0;
            const double pct = scored_haps ? 100.0 * static_cast<double>(n) / static_cast<double>(scored_haps) : 0.0;
            sum << "quintile\tALL\t" << q << '\t' << n << '\t' << pct << '\n';
        }
        for (const char* band : kBands) {
            const std::size_t n = band_count.count(band) ? band_count[band] : 0;
            const double pct = scored_haps ? 100.0 * static_cast<double>(n) / static_cast<double>(scored_haps) : 0.0;
            sum << "haplotype\tALL\t" << band << '\t' << n << '\t' << pct << '\n';
        }
        // Residual split: of all the sequence we could not reconstruct, how much is sub-threshold
        // variation (couldn't be called) vs callable-size events missed/mis-called.
        {
            const std::size_t tot = tot_sub + tot_over;
            const double sub_pct = tot ? 100.0 * static_cast<double>(tot_sub) / static_cast<double>(tot) : 0.0;
            const double over_pct = tot ? 100.0 * static_cast<double>(tot_over) / static_cast<double>(tot) : 0.0;
            sum << "residual\t<" << min_sv_bp << "bp\tsub_threshold\t" << tot_sub << '\t' << sub_pct << '\n';
            sum << "residual\t>=" << min_sv_bp << "bp\tover_threshold\t" << tot_over << '\t' << over_pct << '\n';
        }
        // Genotype-aware reconstruction, reported on the same scale so the two are read side by side.
        if (gt_scored_haps) {
            for (const char* q : kQuintiles) {
                const std::size_t n = gt_quintile_count.count(q) ? gt_quintile_count[q] : 0;
                const double pct = 100.0 * static_cast<double>(n) / static_cast<double>(gt_scored_haps);
                sum << "gt_quintile\tALL\t" << q << '\t' << n << '\t' << pct << '\n';
            }
            for (const char* band : kBands) {
                const std::size_t n = gt_band_count.count(band) ? gt_band_count[band] : 0;
                const double pct = 100.0 * static_cast<double>(n) / static_cast<double>(gt_scored_haps);
                sum << "gt_haplotype\tALL\t" << band << '\t' << n << '\t' << pct << '\n';
            }
            for (const auto& [svclass, c] : carrier) {
                const std::size_t tp = c[0], fp = c[1], fn = c[2], tn = c[3];
                const double prec = (tp + fp) ? 100.0 * static_cast<double>(tp) / static_cast<double>(tp + fp) : 0.0;
                const double rec = (tp + fn) ? 100.0 * static_cast<double>(tp) / static_cast<double>(tp + fn) : 0.0;
                sum << "gt_carrier\t" << svclass << "\tTP\t" << tp << '\t' << prec << '\n';
                sum << "gt_carrier\t" << svclass << "\tFP\t" << fp << '\t' << prec << '\n';
                sum << "gt_carrier\t" << svclass << "\tFN\t" << fn << '\t' << rec << '\n';
                sum << "gt_carrier\t" << svclass << "\tTN\t" << tn << '\t' << rec << '\n';
            }
            {
                const double denom = static_cast<double>(gt_tot_ref_delta) - static_cast<double>(gt_tot_graph_delta);
                const double closed = std::fabs(denom) > 1e-9
                    ? (static_cast<double>(gt_tot_ref_delta) - static_cast<double>(gt_tot_gt_delta)) / denom : 1.0;
                sum << "gt_gap\tALL\tbaseline_delta\t" << gt_tot_ref_delta << "\t0\n";
                sum << "gt_gap\tALL\tgenotype_delta\t" << gt_tot_gt_delta << "\t0\n";
                sum << "gt_gap\tALL\tgraph_delta\t" << gt_tot_graph_delta << "\t0\n";
                sum << "gt_gap\tALL\tgap_closed_pooled\t0\t" << 100.0 * closed << "\n";
                sum << "gt_gap\tALL\tgap_closed_mean\t0\t"
                    << (gt_scored_haps ? 100.0 * gt_closed_sum / static_cast<double>(gt_scored_haps) : 0.0) << "\n";
                sum << "gt_gap\tALL\tworse_than_baseline\t" << gt_worse_than_ref << "\t"
                    << (gt_scored_haps ? 100.0 * static_cast<double>(gt_worse_than_ref) / static_cast<double>(gt_scored_haps) : 0.0) << "\n";
            }
            // What the genotype QV was actually computed from, so it is never read without that context.
            sum << "gt_records\tALL\tapplied\t" << tot_gt.applied << "\t0\n";
            sum << "gt_records\tALL\tunplaceable\t" << tot_gt.unplaceable << "\t0\n";
            sum << "gt_records\tALL\tclamped\t" << tot_gt.clamped << "\t0\n";
            sum << "gt_records\tALL\tunhandled\t" << tot_gt.unhandled << "\t0\n";
            sum << "gt_records\tALL\tref_mismatch\t" << tot_gt.ref_mismatch << "\t0\n";
        }
        for (const auto& [svclass, bands] : class_bands) {
            std::size_t total = 0;
            for (const auto& [b, n] : bands) total += n;
            for (const char* band : kBands) {
                const auto bit = bands.find(band);
                const std::size_t n = bit != bands.end() ? bit->second : 0;
                const double pct = total ? 100.0 * static_cast<double>(n) / static_cast<double>(total) : 0.0;
                sum << "svclass\t" << svclass << '\t' << band << '\t' << n << '\t' << pct << '\n';
            }
        }
    }

    auto quintile_recap = [&](const std::map<std::string, std::size_t>& counts, std::size_t denom) {
        std::string s;
        for (const char* q : kQuintiles) {
            const auto it = counts.find(q);
            const std::size_t n = it != counts.end() ? it->second : 0;
            const double pct = denom ? 100.0 * static_cast<double>(n) / static_cast<double>(denom) : 0.0;
            char buf[64];
            std::snprintf(buf, sizeof(buf), " %s=%.1f%%", q, pct);
            s += buf;
        }
        return s;
    };
    log.info("graph (self-consistency) " + std::to_string(scored_haps) + " haplotypes, qv/qv_max quintiles:" +
             quintile_recap(quintile_count, scored_haps));
    if (gt_scored_haps) {
        log.info("genotype (from VCF) " + std::to_string(gt_scored_haps) + " haplotypes, qv/qv_max quintiles:" +
                 quintile_recap(gt_quintile_count, gt_scored_haps));
        {
            const double denom = static_cast<double>(gt_tot_ref_delta) - static_cast<double>(gt_tot_graph_delta);
            const double closed = std::fabs(denom) > 1e-9
                ? (static_cast<double>(gt_tot_ref_delta) - static_cast<double>(gt_tot_gt_delta)) / denom : 1.0;
            char gb[224];
            std::snprintf(gb, sizeof(gb),
                          "residual bases: baseline(no edits)=%zu genotype=%zu graph=%zu -> gap closed %.1f%%"
                          " (%zu haplotypes worse than baseline)",
                          gt_tot_ref_delta, gt_tot_gt_delta, gt_tot_graph_delta, 100.0 * closed, gt_worse_than_ref);
            log.info(gb);
        }
        const auto& c = carrier["ALL"];
        const double prec = (c[0] + c[1]) ? 100.0 * static_cast<double>(c[0]) / static_cast<double>(c[0] + c[1]) : 0.0;
        const double rec = (c[0] + c[2]) ? 100.0 * static_cast<double>(c[0]) / static_cast<double>(c[0] + c[2]) : 0.0;
        char buf[192];
        std::snprintf(buf, sizeof(buf),
                      "carrier calls: TP=%zu FP=%zu FN=%zu TN=%zu (precision %.1f%%, recall %.1f%%)",
                      c[0], c[1], c[2], c[3], prec, rec);
        log.info(buf);
    } else if (do_gt) {
        log.info("genotype scoring produced no scored haplotype (no bubble placed on the reference)");
    }
    log.wrote({out_prefix + ".qv.tsv", out_prefix + ".qv_by_haplotype.tsv", out_prefix + ".qv_summary.tsv"});
    log.done();
    return 0;
}

} // namespace panvar
