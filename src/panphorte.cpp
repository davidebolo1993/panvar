#include "panvar/panphorte.hpp"

#include "panvar/align.hpp"
#include "panvar/bubble_path.hpp"
#include "panvar/cli_utils.hpp"
#include "panvar/gfa.hpp"
#include "panvar/gfa_io.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/output.hpp"

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace panvar {
namespace {

// Run fn(i) for i in [0, n) across worker threads (work-stealing on an atomic
// counter). threads==0 -> hardware concurrency. Results must be written to
// distinct indices per i (no shared mutation) for thread safety + determinism.
void run_parallel(std::size_t n, std::size_t threads, const std::function<void(std::size_t)>& fn) {
    std::size_t nthreads = threads != 0 ? threads : std::thread::hardware_concurrency();
    if (nthreads == 0) nthreads = 1;
    if (n > 0) nthreads = std::min(nthreads, n);
    if (nthreads <= 1) {
        for (std::size_t i = 0; i < n; ++i) fn(i);
        return;
    }
    std::atomic<std::size_t> next{0};
    std::vector<std::thread> pool;
    pool.reserve(nthreads);
    for (std::size_t t = 0; t < nthreads; ++t) {
        pool.emplace_back([&] {
            for (;;) {
                const std::size_t i = next.fetch_add(1);
                if (i >= n) break;
                fn(i);
            }
        });
    }
    for (std::thread& th : pool) th.join();
}

constexpr std::size_t kMaxGapSteps = 32; // bound interruption length (steps) between copies

std::string step_sequence(const Graph& graph, const PathStep& step) {
    const auto it = graph.nodes.find(step.node_id);
    if (it == graph.nodes.end()) {
        return {};
    }
    return step.reverse ? reverse_complement(it->second.sequence) : it->second.sequence;
}

std::uint64_t hash_sequence(const std::string& s) {
    std::uint64_t h = 1469598103934665603ULL;
    for (const unsigned char c : s) {
        h ^= static_cast<std::uint64_t>(c);
        h *= 1099511628211ULL;
    }
    return h;
}

char orient_char(bool reverse) { return reverse ? '-' : '+'; }

// Per-node precomputed tokens (forward + reverse-complement sequence hash) and
// length, so step tokenization is O(1) lookups instead of re-hashing sequences.
struct NodeTok {
    std::uint64_t fwd = 0;
    std::uint64_t rc = 0;
    std::size_t len = 0;
};

std::unordered_map<std::string, NodeTok> build_node_tokens(const Graph& graph) {
    std::unordered_map<std::string, NodeTok> out;
    out.reserve(graph.nodes.size() * 2);
    for (const auto& [id, node] : graph.nodes) {
        NodeTok nt;
        nt.fwd = hash_sequence(node.sequence);
        nt.rc = hash_sequence(reverse_complement(node.sequence));
        nt.len = node.sequence.size();
        out.emplace(id, nt);
    }
    return out;
}

// One element of a (possibly interrupted) tandem array, in emission order: a unit
// copy (replaced by the REP node) or an interruption kept as literal path steps.
struct ArrayElement {
    bool is_copy = false;
    std::vector<PathStep> steps; // literal steps for interruptions
};

struct TandemArray {
    std::size_t start = 0;       // absolute index into path.steps (array start)
    std::size_t length = 0;      // number of steps replaced
    std::size_t copies = 0;      // unit copies (REP traversals)
    std::size_t interruption_bp = 0;
    std::string unit_seq;        // one unit, native orientation
    std::vector<ArrayElement> elements;                        // ordered emission plan
    std::vector<std::pair<std::size_t, std::size_t>> copy_ranges; // absolute [lo,hi) of copies (for removal)
};

// Detect non-overlapping tandem arrays within the bubble interval [left, right],
// comparing units by spelled-sequence tokens so identical-sequence copies match
// even as distinct node ids. Tolerates short interruptions between copies (kept
// as literal nested steps), e.g. CGG (A) CGG CGG CGG.
std::vector<TandemArray> detect_tandems(
    const Graph& graph,
    const std::unordered_map<std::string, NodeTok>& node_tok,
    const PathRecord& path,
    std::size_t left,
    std::size_t right,
    std::size_t min_unit_bp,
    std::size_t min_copies,
    double max_interruption_frac) {

    std::vector<TandemArray> arrays;
    if (right < left || right >= path.steps.size()) {
        return arrays;
    }
    const std::size_t m = right - left + 1;

    std::vector<std::uint64_t> tok(m);
    std::vector<std::size_t> bp(m);
    for (std::size_t k = 0; k < m; ++k) {
        const PathStep& st = path.steps[left + k];
        const auto it = node_tok.find(st.node_id);
        const NodeTok nt = (it != node_tok.end()) ? it->second : NodeTok{};
        tok[k] = st.reverse ? nt.rc : nt.fwd;
        bp[k] = nt.len;
    }

    // Materialize an oriented step sequence on demand (only for accepted arrays).
    auto seq_at = [&](std::size_t rel) { return step_sequence(graph, path.steps[left + rel]); };

    std::unordered_map<std::uint64_t, std::vector<std::size_t>> token_positions;
    for (std::size_t k = 0; k < m; ++k) {
        token_positions[tok[k]].push_back(k);
    }

    // Fast equality by 64-bit sequence token during detection; accepted arrays
    // are verified by actual sequence below, and the sequence-preserving invariant
    // is the final backstop against the (~1/2^64) token collision.
    auto block_equal = [&](std::size_t a, std::size_t b, std::size_t p) {
        for (std::size_t t = 0; t < p; ++t) {
            if (tok[a + t] != tok[b + t]) {
                return false;
            }
        }
        return true;
    };
    auto block_equal_seq = [&](std::size_t a, std::size_t b, std::size_t p) {
        for (std::size_t t = 0; t < p; ++t) {
            if (seq_at(a + t) != seq_at(b + t)) {
                return false;
            }
        }
        return true;
    };

    std::size_t lower_bound = 0; // do not extend an array left of the previous array
    std::size_t i = 0;
    while (i < m) {
        // Establish a unit: smallest period p with a clean adjacent pair at i.
        std::size_t period = 0;
        const auto pit = token_positions.find(tok[i]);
        if (pit != token_positions.end()) {
            std::size_t tried = 0;
            for (const std::size_t j : pit->second) {
                if (j <= i) {
                    continue;
                }
                const std::size_t p = j - i;
                if (i + 2 * p > m) {
                    break;
                }
                if (++tried > 8) {
                    break;
                }
                if (block_equal(i, i + p, p)) {
                    period = p;
                    break;
                }
            }
        }
        if (period == 0) {
            ++i;
            continue;
        }
        const std::size_t p = period;

        std::size_t unit_bp = 0;
        for (std::size_t t = 0; t < p; ++t) {
            unit_bp += bp[i + t];
        }

        // Collect copy start positions around the anchor, tolerating gaps.
        std::vector<std::size_t> copy_starts = {i};

        // Extend right.
        std::size_t last = i;
        while (true) {
            const std::size_t adj = last + p;
            if (adj + p <= m && block_equal(i, adj, p)) {
                copy_starts.push_back(adj);
                last = adj;
                continue;
            }
            bool found = false;
            std::size_t gap_bp = 0;
            for (std::size_t g = 1; g <= kMaxGapSteps; ++g) {
                const std::size_t cand = last + p + g;
                if (cand + p > m) {
                    break;
                }
                gap_bp += bp[last + p + (g - 1)];
                if (gap_bp > unit_bp) {
                    break;
                }
                if (block_equal(i, cand, p)) {
                    copy_starts.push_back(cand);
                    last = cand;
                    found = true;
                    break;
                }
            }
            if (!found) {
                break;
            }
        }

        // Extend left (not below lower_bound).
        std::size_t first = i;
        while (first >= p && first - p >= lower_bound) {
            const std::size_t adj = first - p;
            if (block_equal(i, adj, p)) {
                copy_starts.insert(copy_starts.begin(), adj);
                first = adj;
                continue;
            }
            bool found = false;
            std::size_t gap_bp = 0;
            for (std::size_t g = 1; g <= kMaxGapSteps; ++g) {
                if (first < p + g || first - p - g < lower_bound) {
                    break;
                }
                const std::size_t cand = first - p - g;
                gap_bp += bp[first - g];
                if (gap_bp > unit_bp) {
                    break;
                }
                if (block_equal(i, cand, p)) {
                    copy_starts.insert(copy_starts.begin(), cand);
                    first = cand;
                    found = true;
                    break;
                }
            }
            if (!found) {
                break;
            }
        }

        const std::size_t array_start = copy_starts.front();
        const std::size_t array_end = copy_starts.back() + p; // exclusive
        const std::size_t copies = copy_starts.size();

        // Build emission plan + measure interruptions.
        std::vector<ArrayElement> elements;
        std::vector<std::pair<std::size_t, std::size_t>> copy_ranges;
        std::size_t interruption_bp = 0;
        for (std::size_t k = 0; k < copy_starts.size(); ++k) {
            ArrayElement copy_el;
            copy_el.is_copy = true;
            elements.push_back(std::move(copy_el));
            copy_ranges.emplace_back(left + copy_starts[k], left + copy_starts[k] + p);
            if (k + 1 < copy_starts.size()) {
                const std::size_t gap_lo = copy_starts[k] + p;
                const std::size_t gap_hi = copy_starts[k + 1];
                if (gap_hi > gap_lo) {
                    ArrayElement gap_el;
                    gap_el.is_copy = false;
                    for (std::size_t s = gap_lo; s < gap_hi; ++s) {
                        gap_el.steps.push_back(path.steps[left + s]);
                        interruption_bp += bp[s];
                    }
                    elements.push_back(std::move(gap_el));
                }
            }
        }

        std::size_t array_bp = 0;
        for (std::size_t s = array_start; s < array_end; ++s) {
            array_bp += bp[s];
        }

        bool seq_ok = true;
        for (std::size_t k = 0; seq_ok && k < copy_starts.size(); ++k) {
            if (!block_equal_seq(i, copy_starts[k], p)) {
                seq_ok = false; // token collision: not a real identical copy
            }
        }

        const bool ok = seq_ok && copies >= min_copies && unit_bp >= min_unit_bp &&
                        (array_bp == 0 ||
                         static_cast<double>(interruption_bp) <= max_interruption_frac * static_cast<double>(array_bp));
        if (ok) {
            TandemArray arr;
            arr.start = left + array_start;
            arr.length = array_end - array_start;
            arr.copies = copies;
            arr.interruption_bp = interruption_bp;
            arr.elements = std::move(elements);
            arr.copy_ranges = std::move(copy_ranges);
            for (std::size_t t = 0; t < p; ++t) {
                arr.unit_seq += seq_at(i + t);
            }
            arrays.push_back(std::move(arr));
            lower_bound = array_end;
            i = array_end;
        } else {
            ++i;
        }
    }
    return arrays;
}

std::string spell_model_path(const GfaModel& model, const GfaPath& path) {
    std::string out;
    for (const PathStep& s : path.steps) {
        const auto it = model.seq.find(s.node_id);
        if (it == model.seq.end()) {
            throw std::runtime_error("Path references missing node: " + s.node_id);
        }
        out += s.reverse ? reverse_complement(it->second) : it->second;
    }
    return out;
}

// Ensure every consecutive step pair in every path has a backing L edge.
struct EdgeSet {
    std::unordered_set<std::string> keys;
    GfaModel& model;
    std::size_t added = 0;

    explicit EdgeSet(GfaModel& m) : model(m) {
        for (const GfaEdge& e : m.edges) {
            keys.insert(key(e.from, e.from_orient, e.to, e.to_orient));
        }
    }
    static std::string key(const std::string& f, char fo, const std::string& t, char to) {
        return f + fo + ">" + t + to;
    }
    void ensure(const PathStep& a, const PathStep& b) {
        const std::string k = key(a.node_id, orient_char(a.reverse), b.node_id, orient_char(b.reverse));
        if (keys.insert(k).second) {
            model.edges.push_back(GfaEdge{a.node_id, orient_char(a.reverse), b.node_id, orient_char(b.reverse), "0M"});
            ++added;
        }
    }
};

struct PathArray {
    std::size_t start = 0;
    std::size_t length = 0;
    std::vector<PathStep> replacement;
};

// ---- Approximate (similarity-based) single-block collapse ----
//
// One representative repeat unit per bubble. Near-identical copies are found by
// banded alignment of the unit against the path's spelled sequence (so a large
// indel between divergent copies, e.g. the C4 long/short HERV, is bridged when
// --min-similarity is low enough), then collapsed (lossy) to one REP node.
// Per-path detection is multithreaded; the graph mutation is serial.

constexpr std::size_t kAnchorSeeds = 8;     // seeds spread across the unit for anchoring

std::string sample_of(const std::string& path_name) {
    const std::size_t h = path_name.find('#');
    return h == std::string::npos ? path_name : path_name.substr(0, h);
}

int base2(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1;
    }
}

// Start positions where `kmer` (<=31 bp, ACGT) occurs in s, via a 2-bit roll.
std::vector<std::size_t> find_kmer_starts(const std::string& s, const std::string& kmer) {
    std::vector<std::size_t> out;
    const std::size_t k = kmer.size();
    if (k == 0 || k > 31 || s.size() < k) return out;
    std::uint64_t target = 0;
    for (const char c : kmer) {
        const int b = base2(c);
        if (b < 0) return out;
        target = (target << 2) | static_cast<std::uint64_t>(b);
    }
    const std::uint64_t mask = (k < 32) ? ((1ULL << (2 * k)) - 1) : ~0ULL;
    std::uint64_t code = 0;
    std::size_t valid = 0;
    for (std::size_t i = 0; i < s.size(); ++i) {
        const int b = base2(s[i]);
        if (b < 0) { valid = 0; code = 0; continue; }
        code = ((code << 2) | static_cast<std::uint64_t>(b)) & mask;
        if (++valid >= k && code == target) out.push_back(i - k + 1);
    }
    return out;
}

struct ApproxCopy {
    std::size_t bp_lo = 0;
    std::size_t bp_hi = 0;
    bool rev = false;
    double identity = 0.0;
};

// Near-identical copies of unit r (and its RC) in spelled sequence s: multi-seed
// anchors propose copy starts, then a banded global alignment of r into the
// window decides each copy and its extent/orientation. Band is (1-min_sim)*|r|
// (uncapped), so divergent copies with a large internal indel still align.
std::vector<ApproxCopy> detect_copies_native(const std::string& s, const std::string& r, double min_sim) {
    std::vector<ApproxCopy> copies;
    const std::size_t k = 16;
    if (r.size() < k || s.size() < r.size()) return copies;
    const std::string rrc = reverse_complement(r);

    std::vector<std::size_t> anchors;
    auto add_seeds = [&](const std::string& ref) {
        const std::size_t span = ref.size() - k;
        for (std::size_t i = 0; i < kAnchorSeeds; ++i) {
            const std::size_t off = (kAnchorSeeds <= 1) ? 0 : span * i / (kAnchorSeeds - 1);
            for (const std::size_t pos : find_kmer_starts(s, ref.substr(off, k))) {
                if (pos >= off) anchors.push_back(pos - off);
            }
        }
    };
    add_seeds(r);
    add_seeds(rrc);
    std::sort(anchors.begin(), anchors.end());
    anchors.erase(std::unique(anchors.begin(), anchors.end()), anchors.end());

    const std::size_t band =
        std::max<std::size_t>(8, static_cast<std::size_t>((1.0 - min_sim) * static_cast<double>(r.size())) + 8);
    std::size_t last_end = 0;
    bool have_last = false;
    for (const std::size_t st : anchors) {
        if (have_last && st < last_end) continue;
        const std::size_t winlen = std::min(s.size() - st, r.size() + band);
        const std::string win = s.substr(st, winlen);
        const FitAlignResult ff = fit_align(r, win, band);
        const FitAlignResult fr = fit_align(rrc, win, band);
        const bool use_rev = fr.identity > ff.identity;
        const FitAlignResult& best = use_rev ? fr : ff;
        if (best.ok && best.identity >= min_sim && best.target_consumed >= r.size() / 2) {
            copies.push_back({st, st + best.target_consumed, use_rev, best.identity});
            last_end = st + best.target_consumed;
            have_last = true;
        }
    }
    return copies;
}

// Nearest step-boundary offset to a bp position, given a prefix-bp table.
std::size_t nearest_step_off(const std::vector<std::size_t>& prefix, std::size_t bp) {
    const auto it = std::lower_bound(prefix.begin(), prefix.end(), bp);
    if (it == prefix.begin()) return 0;
    if (it == prefix.end()) return prefix.size() - 1;
    const std::size_t hi = static_cast<std::size_t>(it - prefix.begin());
    const std::size_t lo = hi - 1;
    return (bp - prefix[lo] <= prefix[hi] - bp) ? lo : hi;
}

// Seed the bubble's repeat unit from the exact tandem detector: a unit must come
// from a clean ADJACENT identical pair (so we get the real repeat period, e.g.
// the whole ~32 kb C4 module, not a small sub-segment that merely recurs far
// apart). Across sampled paths, tally each detected unit and return the most
// supported one (tie-break to the longer unit), or "" if none qualifies.
std::string pick_reference_unit(
    const Graph& graph,
    const std::unordered_map<std::string, NodeTok>& node_tok,
    const std::vector<BubblePathIndex>& path_indexes,
    const Bubble& bubble,
    const PanphorteOptions& options) {

    std::unordered_map<std::string, std::size_t> tally;
    for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
        const auto interval = find_best_bubble_path_interval(path_indexes[pi], bubble);
        if (!interval.has_value()) continue;
        // Exact tandem arrays only need an adjacent identical pair (min_copies 2),
        // independent of the approximate --min-copies the caller will enforce.
        const auto arrays = detect_tandems(
            graph, node_tok, graph.paths[pi], interval->left, interval->right,
            options.min_unit_bp, 2, options.max_interruption_frac);
        for (const TandemArray& arr : arrays) {
            if (arr.unit_seq.size() >= options.min_unit_bp) tally[arr.unit_seq] += arr.copies;
        }
    }
    std::string best;
    std::size_t best_n = 0;
    for (const auto& [seq, cnt] : tally) {
        if (cnt > best_n || (cnt == best_n && seq.size() > best.size())) {
            best_n = cnt;
            best = seq;
        }
    }
    return best;
}

} // namespace

void panphorte_normalize(const PanphorteOptions& options, PanphorteSummary* summary_out) {
    if (options.gfa_path.empty()) {
        throw std::runtime_error("panphorte requires -i/--gfa <graph.gfa>");
    }
    if (options.out_prefix.empty()) {
        throw std::runtime_error("panphorte requires -o/--out-prefix <prefix>");
    }

    ParseGfaOptions parse_opts;
    parse_opts.include_paths = true;
    parse_opts.include_sequences = true;
    const Graph graph = parse_gfa(options.gfa_path, parse_opts);

    GfaModel model = read_gfa_model(options.gfa_path);
    if (model.paths.size() != graph.paths.size()) {
        throw std::runtime_error("Internal: path count mismatch between Graph and GfaModel");
    }

    const std::vector<Bubble> bubbles = read_bubbles_csv(options.bubbles_csv_in);
    std::unordered_set<std::size_t> bubble_filter(options.bubble_ids.begin(), options.bubble_ids.end());

    const std::unordered_map<std::string, NodeTok> node_tok = build_node_tokens(graph);

    std::vector<BubblePathIndex> path_indexes;
    path_indexes.reserve(graph.paths.size());
    for (const auto& p : graph.paths) {
        path_indexes.push_back(build_bubble_path_index(p));
    }

    // Approximate mode (--min-similarity < 1.0): collapse near-identical copies
    // (lossy) to a single representative unit looped N times, preserving per-copy
    // orientation. Exact mode (== 1.0) keeps the sequence-preserving collapse.
    const bool approximate = options.min_similarity < 1.0;
    std::ofstream copies_out;
    if (approximate) {
        copies_out.open(options.out_prefix + ".panphorte.copies.tsv");
        if (!copies_out) {
            throw std::runtime_error("Failed to write panphorte copies TSV");
        }
        copies_out << "path_name\tsample\tbubble_id\tcopies\tunit_bp\torientations\t"
                   << "mean_identity\tregion_bp\tfrom_node\tto_node\tn_long\tn_short\n";
    }

    std::unordered_map<std::string, std::string> rep_by_unit; // canonical unit seq -> REP node id
    std::unordered_set<std::string> removal_candidates;
    std::vector<std::vector<PathArray>> per_path_arrays(graph.paths.size());

    PanphorteSummary summary;
    std::ofstream report(options.out_prefix + ".panphorte.report.tsv");
    if (!report) {
        throw std::runtime_error("Failed to write panphorte report: " + options.out_prefix + ".panphorte.report.tsv");
    }
    report << "bubble_id\tnormalized\tunit_bp\tpaths_normalized\tmin_copies\tmax_copies\t"
           << "interruptions_bp\tnodes_collapsed\n";

    for (const Bubble& bubble : bubbles) {
        if (!bubble_filter.empty() && bubble_filter.find(bubble.id) == bubble_filter.end()) {
            continue;
        }
        ++summary.bubbles_seen;

        std::size_t paths_norm = 0;
        std::size_t min_copies_seen = 0;
        std::size_t max_copies_seen = 0;
        std::size_t unit_bp_report = 0;
        std::size_t interruptions_bp = 0;
        std::unordered_set<std::string> collapsed_nodes;

        if (approximate) {
            // Single-block: seed ONE representative repeat unit for the bubble, then
            // find its near-identical copies (any orientation, not necessarily
            // adjacent) per path by banded sequence alignment, and collapse them to
            // one REP node looped per-copy (lossy). Sequence between copies is kept
            // as literal steps, so only the repeat copies lose within-copy detail.
            const std::string ref_unit =
                pick_reference_unit(graph, node_tok, path_indexes, bubble, options);
            if (ref_unit.size() >= options.min_unit_bp) {
                struct Mapped { std::size_t off_lo, off_hi; bool rev; double id; };
                struct Res { bool has = false; std::size_t left = 0; std::vector<Mapped> copies; };
                std::vector<Res> results(graph.paths.size());

                auto detect_one = [&](std::size_t pi) {
                    const auto interval = find_best_bubble_path_interval(path_indexes[pi], bubble);
                    if (!interval.has_value()) return;
                    const std::size_t left = interval->left;
                    const std::size_t n = interval->right - left + 1;
                    const std::vector<PathStep> nat(
                        graph.paths[pi].steps.begin() + static_cast<std::ptrdiff_t>(left),
                        graph.paths[pi].steps.begin() + static_cast<std::ptrdiff_t>(left + n));
                    const std::string snat = spell_path_steps_sequence(graph, nat);
                    std::vector<std::size_t> prefix(n + 1, 0); // cumulative bp at step boundaries
                    for (std::size_t off = 0; off < n; ++off) {
                        const auto it = node_tok.find(nat[off].node_id);
                        prefix[off + 1] = prefix[off] + (it != node_tok.end() ? it->second.len : 0);
                    }
                    const auto bp_copies = detect_copies_native(snat, ref_unit, options.min_similarity);
                    std::vector<Mapped> copies;
                    std::size_t prev_hi = 0;
                    for (const ApproxCopy& c : bp_copies) {
                        std::size_t lo = nearest_step_off(prefix, c.bp_lo);
                        std::size_t hi = nearest_step_off(prefix, c.bp_hi);
                        if (lo < prev_hi) lo = prev_hi;       // keep monotonic / non-overlapping
                        if (hi <= lo) continue;
                        copies.push_back({lo, hi, c.rev, c.identity});
                        prev_hi = hi;
                    }
                    if (copies.size() >= options.min_copies) {
                        results[pi] = Res{true, left, std::move(copies)};
                    }
                };

                run_parallel(graph.paths.size(), options.threads, detect_one);

                // Serial: one REP node per bubble; build PathArrays + copies.tsv rows.
                std::string rep_id;
                for (std::size_t pi = 0; pi < results.size(); ++pi) {
                    if (!results[pi].has) continue;
                    const std::size_t left = results[pi].left;
                    const std::vector<Mapped>& copies = results[pi].copies;
                    if (rep_id.empty()) {
                        const auto rit = rep_by_unit.find(ref_unit);
                        if (rit == rep_by_unit.end()) {
                            rep_id = add_new_node(model, ref_unit);
                            rep_by_unit.emplace(ref_unit, rep_id);
                            ++summary.nodes_added;
                        } else {
                            rep_id = rit->second;
                        }
                    }

                    PathArray pa;
                    pa.start = left + copies.front().off_lo;
                    std::size_t prev_hi = copies.front().off_lo;
                    std::string orients;
                    double id_sum = 0.0;
                    // Per-copy length classes. The reference unit is the long form; a copy
                    // whose traversed sequence is substantially shorter carries a large
                    // internal deletion (e.g. the C4 HERV-K ~6.4 kb) and is the SHORT form.
                    // This recovers the long/short composition the collapse would otherwise
                    // hide (the short copy is folded onto the same REP node).
                    constexpr double kShortMaxFrac = 0.90;
                    std::size_t n_short = 0;
                    for (const Mapped& m : copies) {
                        for (std::size_t off = prev_hi; off < m.off_lo; ++off) {
                            pa.replacement.push_back(graph.paths[pi].steps[left + off]);
                        }
                        pa.replacement.push_back(PathStep{rep_id, m.rev});
                        std::size_t copy_bp = 0;
                        for (std::size_t off = m.off_lo; off < m.off_hi; ++off) {
                            removal_candidates.insert(graph.paths[pi].steps[left + off].node_id);
                            collapsed_nodes.insert(graph.paths[pi].steps[left + off].node_id);
                            const auto cit = node_tok.find(graph.paths[pi].steps[left + off].node_id);
                            if (cit != node_tok.end()) copy_bp += cit->second.len;
                        }
                        if (ref_unit.size() > 0 &&
                            static_cast<double>(copy_bp) < kShortMaxFrac * static_cast<double>(ref_unit.size())) {
                            ++n_short;
                        }
                        orients += (m.rev ? '-' : '+');
                        id_sum += m.id;
                        prev_hi = m.off_hi;
                    }
                    pa.length = (left + prev_hi) - pa.start;
                    per_path_arrays[pi].push_back(std::move(pa));
                    ++paths_norm;

                    const std::size_t occ = copies.size();
                    if (min_copies_seen == 0 || occ < min_copies_seen) min_copies_seen = occ;
                    if (occ > max_copies_seen) max_copies_seen = occ;
                    if (unit_bp_report == 0) unit_bp_report = ref_unit.size();
                    std::size_t region_bp = 0;
                    for (std::size_t off = copies.front().off_lo; off < copies.back().off_hi; ++off) {
                        const auto it = node_tok.find(graph.paths[pi].steps[left + off].node_id);
                        if (it != node_tok.end()) region_bp += it->second.len;
                    }
                    copies_out << graph.paths[pi].name << '\t' << sample_of(graph.paths[pi].name) << '\t'
                               << bubble.id << '\t' << occ << '\t' << ref_unit.size() << '\t'
                               << orients << '\t' << (id_sum / static_cast<double>(occ)) << '\t'
                               << region_bp << '\t'
                               << graph.paths[pi].steps[left + copies.front().off_lo].node_id << '\t'
                               << graph.paths[pi].steps[left + copies.back().off_hi - 1].node_id << '\t'
                               << (occ - n_short) << '\t' << n_short << '\n';
                }
            }

            const bool norm_ok = paths_norm > 0;
            if (norm_ok) ++summary.bubbles_normalized;
            report << bubble.id << '\t' << (norm_ok ? "yes" : "no") << '\t'
                   << unit_bp_report << '\t' << paths_norm << '\t'
                   << min_copies_seen << '\t' << max_copies_seen << '\t'
                   << interruptions_bp << '\t' << collapsed_nodes.size() << '\n';
            continue; // bubble fully handled in approximate mode
        }

        for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
            const auto interval = find_best_bubble_path_interval(path_indexes[pi], bubble);
            if (!interval.has_value()) {
                continue;
            }

            const auto arrays = detect_tandems(
                graph, node_tok, graph.paths[pi], interval->left, interval->right,
                options.min_unit_bp, options.min_copies, options.max_interruption_frac);
            if (arrays.empty()) {
                continue;
            }

            for (const TandemArray& arr : arrays) {
                const std::string& unit = arr.unit_seq;
                const std::string rc = reverse_complement(unit);
                const bool forward = unit <= rc;
                const std::string canonical = forward ? unit : rc;
                const bool rev = !forward;

                auto rit = rep_by_unit.find(canonical);
                std::string rep_id;
                if (rit == rep_by_unit.end()) {
                    rep_id = add_new_node(model, canonical);
                    rep_by_unit.emplace(canonical, rep_id);
                    ++summary.nodes_added;
                } else {
                    rep_id = rit->second;
                }

                PathArray pa;
                pa.start = arr.start;
                pa.length = arr.length;
                for (const ArrayElement& el : arr.elements) {
                    if (el.is_copy) {
                        pa.replacement.push_back(PathStep{rep_id, rev});
                    } else {
                        pa.replacement.insert(pa.replacement.end(), el.steps.begin(), el.steps.end());
                    }
                }
                per_path_arrays[pi].push_back(std::move(pa));

                for (const auto& r : arr.copy_ranges) {
                    for (std::size_t s = r.first; s < r.second; ++s) {
                        removal_candidates.insert(graph.paths[pi].steps[s].node_id);
                        collapsed_nodes.insert(graph.paths[pi].steps[s].node_id);
                    }
                }
                if (min_copies_seen == 0 || arr.copies < min_copies_seen) min_copies_seen = arr.copies;
                if (arr.copies > max_copies_seen) max_copies_seen = arr.copies;
                if (unit_bp_report == 0) unit_bp_report = unit.size();
                interruptions_bp += arr.interruption_bp;
            }
            ++paths_norm;
        }

        const bool normalized = paths_norm > 0;
        if (normalized) {
            ++summary.bubbles_normalized;
        }
        report << bubble.id << '\t' << (normalized ? "yes" : "no") << '\t'
               << unit_bp_report << '\t' << paths_norm << '\t'
               << min_copies_seen << '\t' << max_copies_seen << '\t'
               << interruptions_bp << '\t' << collapsed_nodes.size() << '\n';
    }

    // Splice replacements into each path's steps.
    for (std::size_t pi = 0; pi < model.paths.size(); ++pi) {
        auto& arrays = per_path_arrays[pi];
        if (arrays.empty()) {
            continue;
        }
        std::sort(arrays.begin(), arrays.end(),
                  [](const PathArray& a, const PathArray& b) { return a.start < b.start; });

        const std::vector<PathStep>& old_steps = model.paths[pi].steps;
        std::vector<PathStep> new_steps;
        new_steps.reserve(old_steps.size());
        std::size_t cursor = 0;
        for (const PathArray& a : arrays) {
            for (; cursor < a.start; ++cursor) {
                new_steps.push_back(old_steps[cursor]);
            }
            new_steps.insert(new_steps.end(), a.replacement.begin(), a.replacement.end());
            cursor = a.start + a.length;
        }
        for (; cursor < old_steps.size(); ++cursor) {
            new_steps.push_back(old_steps[cursor]);
        }
        model.paths[pi].steps = std::move(new_steps);
        summary.paths_rewritten += 1;
    }

    // Ensure every consecutive step pair in every (rewritten) path has an L edge.
    EdgeSet edges(model);
    for (const GfaPath& p : model.paths) {
        for (std::size_t k = 1; k < p.steps.size(); ++k) {
            edges.ensure(p.steps[k - 1], p.steps[k]);
        }
    }
    summary.edges_added = edges.added;

    // Exact mode is sequence-preserving: every rewritten path must spell the same
    // sequence. Approximate mode is intentionally lossy (copies are canonicalized
    // to the representative unit), so this invariant does not apply there.
    if (!approximate) {
        for (std::size_t i = 0; i < model.paths.size(); ++i) {
            if (per_path_arrays[i].empty()) {
                continue;
            }
            const std::string original = spell_path_steps_sequence(graph, graph.paths[i].steps);
            if (spell_model_path(model, model.paths[i]) != original) {
                throw std::runtime_error(
                    "panphorte sequence-preservation violated for path index " + std::to_string(i));
            }
        }
    }

    // Remove collapsed nodes no longer referenced by any path; drop their edges.
    std::unordered_set<std::string> referenced;
    for (const GfaPath& p : model.paths) {
        for (const PathStep& s : p.steps) {
            referenced.insert(s.node_id);
        }
    }
    std::unordered_set<std::string> to_remove;
    for (const std::string& id : removal_candidates) {
        if (referenced.find(id) == referenced.end()) {
            to_remove.insert(id);
        }
    }
    if (!to_remove.empty()) {
        std::vector<std::string> kept_order;
        kept_order.reserve(model.node_order.size());
        for (const std::string& id : model.node_order) {
            if (to_remove.find(id) == to_remove.end()) {
                kept_order.push_back(id);
            } else {
                model.seq.erase(id);
                ++summary.nodes_removed;
            }
        }
        model.node_order = std::move(kept_order);
        std::vector<GfaEdge> kept_edges;
        kept_edges.reserve(model.edges.size());
        for (const GfaEdge& e : model.edges) {
            if (to_remove.find(e.from) == to_remove.end() && to_remove.find(e.to) == to_remove.end()) {
                kept_edges.push_back(e);
            }
        }
        model.edges = std::move(kept_edges);
    }

    const std::string normalized_gfa = options.out_prefix + ".normalized.gfa";
    cli::ensure_parent_dir_for_file(normalized_gfa);
    write_gfa_model(normalized_gfa, model);

    if (summary_out != nullptr) {
        *summary_out = summary;
    }
}

} // namespace panvar
