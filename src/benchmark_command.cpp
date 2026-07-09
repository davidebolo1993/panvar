#include "panvar/benchmark_command.hpp"

#include "panvar/align.hpp"
#include "panvar/bubble_path.hpp"
#include "panvar/bubbles.hpp"
#include "panvar/cli_utils.hpp"
#include "panvar/gfa.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/output.hpp"
#include "panvar/parallel.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
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
        << "Round-trip QV of the caller's own output. For each called bubble and each haplotype that\n"
        << "traverses it, reconstruct the haplotype on the passed graph -- take the reference walk and\n"
        << "substitute only the divergent blocks the calls explain (nodes in variant_nodes.tsv), leaving\n"
        << "uncalled variation (SNPs, sub-threshold indels) at reference -- then align the reconstruction\n"
        << "to the haplotype's true walk (edlib NW) and score QV = -10*log10(max(0.5, delta)/S).\n"
        << "Per-haplotype QV is the length-weighted combine over its bubbles (sum delta / sum S), binned\n"
        << "into the cosigt bands (<17, 17-23, 23-33, >33). Uncalled variation is what keeps delta > 0,\n"
        << "so well-reconstructed haplotypes land high; a missed event drops the band.\n\n"
        << "Options:\n"
        << "  -i, --gfa <path>                 Passed graph the calls were made on (required)\n"
        << "  -b, --bubble-prefix-in <prefix>  panphorte output prefix (auto-uses <prefix>.bubbles.csv)\n"
        << "  -c, --bubbles-csv-in <path>      panphorte bubbles CSV (required if no prefix)\n"
        << "      --variant-nodes <path>       call's <prefix>.variant_nodes.tsv (required)\n"
        << "  -r, --reference-path <name>      Reference path name/substring, the diff baseline (required)\n"
        << "  -o, --out-prefix <prefix>        Output prefix for the QV tables (required)\n"
        << "      --all-bubbles                Score every bubble in the CSV, not just called ones. A\n"
        << "                                   haplotype diverging from reference at an uncalled bubble\n"
        << "                                   then lands in a low band -- a missed call. Default scores\n"
        << "                                   only bubbles that carry >=1 call.\n"
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

// Steps of `path` across `bubble` (canonical source->sink), mirroring call's bubble_steps: fall back
// to a bare [source, sink] pair for haplotypes that cross with no inside node (the short side of a
// deletion), which the inside-node interval finder would otherwise drop.
std::optional<std::vector<PathStep>> bubble_steps(
    const PathRecord& path, const BubblePathIndex& index, const Bubble& bubble) {
    const auto iv = find_best_bubble_path_interval(index, bubble);
    if (iv.has_value()) {
        std::vector<PathStep> s = canonical_bubble_path_steps(path, bubble, *iv);
        if (!s.empty()) return s;
    }
    const auto si = index.positions.find(bubble.source);
    const auto ki = index.positions.find(bubble.sink);
    if (si == index.positions.end() || ki == index.positions.end()) return std::nullopt;
    const std::unordered_set<std::size_t> sink_pos(ki->second.begin(), ki->second.end());
    for (const std::size_t p : si->second) {
        if (sink_pos.count(p + 1)) return std::vector<PathStep>{ path.steps[p], path.steps[p + 1] };
    }
    const std::unordered_set<std::size_t> src_pos(si->second.begin(), si->second.end());
    for (const std::size_t p : ki->second) {
        if (src_pos.count(p + 1)) {
            return std::vector<PathStep>{
                PathStep{ path.steps[p + 1].node_id, !path.steps[p + 1].reverse },
                PathStep{ path.steps[p].node_id, !path.steps[p].reverse } };
        }
    }
    return std::nullopt;
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
};

struct HapResult {
    std::string sample;
    bool scored = false;
    std::size_t sum_delta = 0;
    std::size_t sum_aln = 0;
    std::vector<BubbleObs> obs;
};

} // namespace

int run_benchmark_command(const std::vector<std::string>& args) {
    if (args.empty()) { print_benchmark_help(); return 0; }

    std::string gfa_path, bubble_prefix_in, bubbles_csv_in, variant_nodes_in, reference_path, out_prefix;
    std::size_t threads = 0;
    bool quiet = false;
    bool all_bubbles = false;

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
        else if (arg == "-r" || arg == "--reference-path") reference_path = require_value(arg);
        else if (arg == "-o" || arg == "--out-prefix") out_prefix = require_value(arg);
        else if (arg == "--all-bubbles") all_bubbles = true;
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
            const NwAlign nw = nw_edit_distance(recon, truth);
            hr.scored = true;
            hr.sum_delta += nw.edits;
            hr.sum_aln += nw.aln_len;
            hr.obs.push_back({b.id, applied, nw.edits, nw.aln_len});
        }
    });

    std::map<std::string, std::size_t> band_count;
    std::map<std::string, std::size_t> quintile_count;
    std::size_t scored_haps = 0;
    {
        std::ofstream by_hap(out_prefix + ".qv_by_haplotype.tsv");
        if (!by_hap) throw std::runtime_error("Failed to write " + out_prefix + ".qv_by_haplotype.tsv");
        by_hap << "sample\tn_bubbles\tsum_delta\tsum_aln_len\tqv\tband\tqv_max\tqv_ratio\tquintile\tidentity\n";
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
                   << ratio << '\t' << quintile << '\t' << identity << '\n';
            ++band_count[band];
            ++quintile_count[quintile];
            ++scored_haps;
        }
    }

    std::map<std::string, std::map<std::string, std::size_t>> class_bands;
    {
        std::ofstream qv_out(out_prefix + ".qv.tsv");
        if (!qv_out) throw std::runtime_error("Failed to write " + out_prefix + ".qv.tsv");
        qv_out << "sample\tbubble_id\tsvtypes\tis_carrier\tdelta\taln_len\tqv\n";
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
                       << (o.carrier ? 1 : 0) << '\t' << o.delta << '\t' << o.aln_len << '\t' << qv << '\n';
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

    std::string recap = std::to_string(scored_haps) + " haplotypes (qv/qv_max quintiles):";
    for (const char* q : kQuintiles) {
        const std::size_t n = quintile_count.count(q) ? quintile_count[q] : 0;
        const double pct = scored_haps ? 100.0 * static_cast<double>(n) / static_cast<double>(scored_haps) : 0.0;
        char buf[64];
        std::snprintf(buf, sizeof(buf), " %s=%.1f%%", q, pct);
        recap += buf;
    }
    log.info(recap);
    log.wrote({out_prefix + ".qv.tsv", out_prefix + ".qv_by_haplotype.tsv", out_prefix + ".qv_summary.tsv"});
    log.done();
    return 0;
}

} // namespace panvar
