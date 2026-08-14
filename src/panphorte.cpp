#include "panvar/panphorte.hpp"

#include "panvar/align.hpp"
#include "panvar/bubble_path.hpp"
#include "panvar/bubbles.hpp"
#include "panvar/cli_utils.hpp"
#include "panvar/gfa.hpp"
#include "panvar/gfa_io.hpp"
#include "panvar/graph_sort.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/gtf.hpp"
#include "panvar/integrated_snarls.hpp"
#include "panvar/output.hpp"
#include "panvar/parallel.hpp"

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace panvar {
namespace {

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

// The identity of a repeat unit, independent of where the array happens to be cut, which strand it is
// stored on, and whether the caller handed us one copy or three.
//
// Prevalence used to ask "does this haplotype contain SOME tandem", which is a different question from
// "do these haplotypes share a repeat". One bubble carrying a poly-C array on one haplotype and a
// poly-T array on another satisfied the first and was folded into two unrelated REP nodes; only the
// second question distinguishes a population VNTR from two coincidental arrays.
std::string canonical_motif_key(const std::string& unit) {
    if (unit.empty()) return unit;
    // Primitive root: if the unit is itself w repeated, the motif is w.
    std::string prim = unit;
    for (std::size_t p = 1; p <= unit.size() / 2; ++p) {
        if (unit.size() % p != 0) continue;
        bool ok = true;
        for (std::size_t i = p; i < unit.size() && ok; ++i) ok = (unit[i] == unit[i - p]);
        if (ok) { prim = unit.substr(0, p); break; }
    }
    // Least rotation, on both strands.
    const auto least_rotation = [](const std::string& t) {
        std::string best = t;
        const std::string dbl = t + t;
        for (std::size_t i = 1; i < t.size(); ++i) {
            const std::string rot = dbl.substr(i, t.size());
            if (rot < best) best = rot;
        }
        return best;
    };
    const std::string f = least_rotation(prim);
    const std::string r = least_rotation(reverse_complement(prim));
    return f <= r ? f : r;
}


// Copies of a KNOWN motif that sit exactly on node boundaries, including a lone copy.
//
// detect_tandems finds a PERIOD, which needs at least two copies to exist -- so a haplotype carrying
// one copy has no array and, in exact mode, was left on the original nodes while its neighbours folded.
// Once the site is confirmed the motif is known, so its copies can be matched directly instead of
// rediscovered. A copy spanning several nodes is not matched here; that is the multi-node case
// detect_tandems already covers whenever two or more copies are present.
// Copies of a KNOWN motif that sit exactly on node boundaries, including a lone copy and a copy split
// across several nodes.
//
// detect_tandems finds a PERIOD, which needs two copies to exist, so a haplotype carrying one copy has
// no array. Matching only a single node then left a copy split across nodes literal while its
// neighbours folded -- and a mixed representation is what makes `call` report CN 0 for a haplotype that
// carries one copy. The unit is already known once the site is confirmed, so this is a direct match
// rather than a rediscovery: spans of whole consecutive steps whose spelling IS the motif.
std::vector<TandemArray> find_motif_copies_by_node(
    const Graph& graph,
    const PathRecord& path,
    std::size_t left,
    std::size_t right,
    const std::string& motif_key,
    std::size_t min_unit_bp) {

    std::vector<TandemArray> out;
    if (right <= left || right >= path.steps.size()) return out;

    const auto spell_span = [&](std::size_t i, std::size_t j) {
        std::string acc;
        for (std::size_t k = i; k < j; ++k) {
            const auto nit = graph.nodes.find(path.steps[k].node_id);
            if (nit == graph.nodes.end()) return std::string();
            acc += path.steps[k].reverse ? reverse_complement(nit->second.sequence)
                                         : nit->second.sequence;
        }
        return acc;
    };

    // Accepted spans, left to right. Both ends are step boundaries by construction, so no copy can be
    // accepted whose edge falls inside a node.
    constexpr std::size_t kMaxSpanSteps = 64;
    std::vector<std::pair<std::size_t, std::size_t>> spans;
    std::vector<std::string> span_seq;
    for (std::size_t i = left + 1; i < right;) {
        std::size_t end = 0;
        std::string seq_at_end;
        for (std::size_t j = i + 1; j <= right && j - i <= kMaxSpanSteps; ++j) {
            const std::string seq = spell_span(i, j);
            if (seq.empty()) break;
            if (seq.size() >= min_unit_bp && canonical_motif_key(seq) == motif_key) {
                end = j;
                seq_at_end = seq;
                break;
            }
        }
        if (end == 0) { ++i; continue; }
        spans.emplace_back(i, end);
        span_seq.push_back(std::move(seq_at_end));
        i = end;
    }

    // Contiguous spans form one array; a gap starts a new one.
    for (std::size_t k = 0; k < spans.size();) {
        std::size_t j = k;
        while (j + 1 < spans.size() && spans[j + 1].first == spans[j].second) ++j;
        TandemArray a;
        a.start = spans[k].first;
        a.length = spans[j].second - spans[k].first;
        a.copies = j - k + 1;
        a.unit_seq = span_seq[k];
        a.elements.assign(a.copies, ArrayElement{true, {}});
        a.copy_ranges.emplace_back(spans[k].first, spans[j].second);
        out.push_back(std::move(a));
        k = j + 1;
    }
    return out;
}

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
    // ONE arc, read from either end: a+ -> b+ and b- -> a- are the same link, and a GFA stores whichever
    // of the two it happens to store. Keying on the stored direction alone means an edge recorded from
    // one side is not found from the other -- which left reverse-traversed obsolete arcs in the graph
    // and let ensure() add a second L line for a link that was already there.
    static std::string key(const std::string& f, char fo, const std::string& t, char to) {
        const std::string fwd = f + fo + ">" + t + to;
        const char rf = fo == '+' ? '-' : '+';
        const char rt = to == '+' ? '-' : '+';
        const std::string rev = t + rt + ">" + f + rf;
        return fwd <= rev ? fwd : rev;
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
    std::size_t bubble_id = 0;   // the site that produced this edit, for the acceptance check
    std::vector<PathStep> replacement;
};

// Approximate (similarity-based) collapse: pick one representative repeat unit per bubble, find
// near-identical copies by banded alignment of the unit against each path, then collapse
// them lossily to one REP node. Detection is per-path/multithreaded; the graph mutation is serial.

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


// Base-level tandem seeding: the smallest lag (>= min_unit_bp) at which the spelled interval is
// self-similar. Fallback for graphs whose node boundaries do not follow repeat-unit boundaries, where
// `detect_tandems` -- which measures its period in node STEPS -- sees no period at all. Cost is bounded
// and independent of unit size: 2-bit k-mer scans propose candidate lags, each settled by a sampled
// identity check rather than a full O(len^2) scan.
std::string seed_unit_by_period(const std::string& s, std::size_t min_unit_bp, double min_identity) {
    constexpr std::size_t k = 16;
    constexpr std::size_t kProbes = 16;        // seed offsets spread ACROSS the interval: the array is
                                               // usually interior (the interval carries the source/sink
                                               // flanks, often kb of backbone), and probing only the
                                               // start would seed inside the flank and find nothing
    constexpr std::size_t kMaxCandidates = 8;  // candidate lags examined per offset
    constexpr std::size_t kSamples = 2048;     // sampled identity keeps each candidate ~O(1)
    if (min_unit_bp == 0 || s.size() < 2 * min_unit_bp || s.size() < 2 * k) {
        return {};
    }
    const std::size_t stride = std::max<std::size_t>(k, s.size() / kProbes);
    for (std::size_t off = 0, tried_off = 0; tried_off < kProbes && off + k <= s.size();
         off += stride, ++tried_off) {
        const std::vector<std::size_t> occ = find_kmer_starts(s, s.substr(off, k));
        std::size_t tried = 0;
        for (const std::size_t j : occ) {
            if (j <= off) {
                continue;
            }
            const std::size_t p = j - off;
            if (p < min_unit_bp) {
                continue;
            }
            if (off + 2 * p > s.size()) {
                break; // need two full copies to call it a tandem
            }
            if (++tried > kMaxCandidates) {
                break;
            }
            // Verify LOCALLY, on adjacent copy pairs walking right from the seed. Measuring identity
            // across the whole interval would drown the signal: the interval carries the source/sink
            // flanks (often kb of non-periodic backbone), so a real array's lag scores near-random there.
            const std::size_t sample_step = std::max<std::size_t>(1, p / kSamples);
            std::size_t copies = 1;
            for (std::size_t a = off; a + 2 * p <= s.size(); a += p) {
                std::size_t seen = 0;
                std::size_t hit = 0;
                for (std::size_t i = 0; i < p; i += sample_step) {
                    ++seen;
                    if (s[a + i] == s[a + p + i]) {
                        ++hit;
                    }
                }
                if (seen == 0 || static_cast<double>(hit) / static_cast<double>(seen) < min_identity) {
                    break;
                }
                ++copies;
            }
            if (copies >= 2) {
                return s.substr(off, p); // a rotation of the unit; any rotation folds equivalently
            }
        }
    }
    return {};
}

// Seed the repeat unit from the exact tandem detector: a unit must come from a clean adjacent
// identical pair, so we get the true period (the whole ~32 kb C4 module, not a recurring sub-segment).
// Per-path scan is parallel; the tally merge is serial in path order and the winner is the
// most-supported unit (tie-break: longer, then lexicographic) so it's deterministic across threads.
std::string pick_reference_unit(
    const Graph& graph,
    const std::unordered_map<std::string, NodeTok>& node_tok,
    const std::vector<BubblePathIndex>& path_indexes,
    const Bubble& bubble,
    const PanphorteOptions& options) {

    using Cand = std::pair<std::string, std::size_t>; // (unit_seq, copies)
    std::vector<std::vector<Cand>> per_path(graph.paths.size());
    run_parallel(graph.paths.size(), options.threads, [&](std::size_t pi) {
        // The same INTERIOR window detection uses. Seeding from the anchor-inclusive interval let a
        // periodic backbone anchor win unit selection, so a site could be seeded on sequence that is
        // then never rewritten -- a false negative produced by the anchors alone.
        BubblePathInterval used{};
        if (!bubble_steps(graph.paths[pi], path_indexes[pi], bubble, &used).has_value()) return;
        if (used.inside_count == 0) return;
        const std::optional<BubblePathInterval> interval = used;
        // Exact tandem arrays only need an adjacent identical pair (min_copies 2),
        // independent of the approximate --min-copies the caller will enforce.
        const auto arrays = detect_tandems(
            graph, node_tok, graph.paths[pi], interval->left + 1, interval->right - 1,
            options.min_unit_bp, 2, options.max_interruption_frac);
        std::vector<Cand> cands;
        for (const TandemArray& arr : arrays) {
            if (arr.unit_seq.size() >= options.min_unit_bp) {
                cands.emplace_back(arr.unit_seq, arr.copies);
            }
        }
        per_path[pi] = std::move(cands);
    });
    std::unordered_map<std::string, std::size_t> tally;
    for (const auto& v : per_path) {
        for (const auto& [seq, w] : v) tally[seq] += w;
    }
    // Fallback: no node-step period anywhere in this bubble. Retry at base level, where a unit whose
    // boundaries straddle node boundaries is still visible. Gated on the token detector coming up empty,
    // so graphs that split at repeat boundaries never reach this and are unchanged.
    if (tally.empty()) {
        std::vector<std::vector<Cand>> per_path_bp(graph.paths.size());
        run_parallel(graph.paths.size(), options.threads, [&](std::size_t pi) {
            BubblePathInterval bused{};
            if (!bubble_steps(graph.paths[pi], path_indexes[pi], bubble, &bused).has_value()) return;
            if (bused.inside_count == 0) return;
            const std::size_t left = bused.left + 1;          // interior only, as detection uses
            if (bused.right <= left) return;
            const std::size_t n = bused.right - left;
            const std::vector<PathStep> nat(
                graph.paths[pi].steps.begin() + static_cast<std::ptrdiff_t>(left),
                graph.paths[pi].steps.begin() + static_cast<std::ptrdiff_t>(left + n));
            const std::string snat = spell_path_steps_sequence(graph, nat);
            const std::string unit =
                seed_unit_by_period(snat, options.min_unit_bp, options.min_similarity);
            if (!unit.empty()) {
                per_path_bp[pi].emplace_back(unit, snat.size() / unit.size());
            }
        });
        for (const auto& v : per_path_bp) {
            for (const auto& [seq, w] : v) tally[seq] += w;
        }
    }
    std::string best;
    std::size_t best_n = 0;
    for (const auto& [seq, cnt] : tally) {
        const bool better = cnt > best_n ||
                            (cnt == best_n && seq.size() > best.size()) ||
                            (cnt == best_n && seq.size() == best.size() && seq < best);
        if (best.empty() || better) {
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
    // Create the output directory up front (like the other modules) so the report /
    // copies / GFA writers below don't fail when the prefix points at a missing dir.
    cli::ensure_parent_dir_for_file(options.out_prefix);

    ParseGfaOptions parse_opts;
    parse_opts.include_paths = true;
    parse_opts.include_sequences = true;
    const Graph graph = parse_gfa(options.gfa_path, parse_opts);
    // The same contract every other module applies. panphorte REWRITES the graph, so accepting a
    // malformed one means emitting a repaired-looking graph that was never validated: a path step with
    // no link behind it was silently given the edge it lacked.
    validate_graph_paths(graph, "panphorte", true, true);

    GfaModel model = read_gfa_model(options.gfa_path);
    if (model.paths.size() != graph.paths.size()) {
        throw std::runtime_error("Internal: path count mismatch between Graph and GfaModel");
    }

    // No output may name an input: panphorte reads the graph and the CSV throughout the run, and the
    // report and copies files used to be opened long before detection, sorting, re-snarl or GTF work
    // could fail.
    {
        const std::string outs[] = {
            options.out_prefix + ".panphorte.report.tsv",
            options.out_prefix + ".panphorte.copies.tsv",
            options.out_prefix + ".normalized.gfa",
            options.out_prefix + ".normalized.sorted.gfa",
            options.out_prefix + ".bubbles.csv",
            options.out_prefix + ".bandage_nodes.csv",
            options.out_prefix + ".bandage_genes.csv",
            options.out_prefix + ".panphorte.rep_provenance.tsv",
        };
        // The GTF is an input like any other. Checked only against the gene CSV, it could name the
        // report or the provenance table and be overwritten by the commit, after it had been read.
        for (const std::string& in : {options.gfa_path, options.bubbles_csv_in, options.gtf_path}) {
            if (in.empty()) continue;
            std::error_code ec;
            const auto ip = std::filesystem::weakly_canonical(in, ec);
            if (ec || ip.empty()) continue;
            for (const std::string& o : outs) {
                std::error_code e2;
                const auto op = std::filesystem::weakly_canonical(o, e2);
                if (!e2 && ip == op)
                    throw std::runtime_error("panphorte: output '" + o + "' is the same file as input '" +
                                             in + "'");
            }
        }
    }

    const std::vector<Bubble> bubbles = read_bubbles_csv(options.bubbles_csv_in);
    // The CSV and the graph are separate inputs with nothing tying them together, so a CSV built from
    // a different graph -- or from this one before a rewrite renumbered it -- would otherwise run to
    // completion against nodes that do not exist.
    {
        std::vector<std::string> missing;
        std::unordered_set<std::size_t> seen_ids;
        for (const Bubble& b : bubbles) {
            if (!seen_ids.insert(b.id).second)
                throw std::runtime_error("panphorte: duplicate bubble id in the CSV: " +
                                         std::to_string(b.id));
            const auto check = [&](const std::string& n) {
                if (graph.nodes.find(n) == graph.nodes.end() && missing.size() < 8)
                    missing.push_back("bubble " + std::to_string(b.id) + " node " + n);
            };
            check(b.source);
            check(b.sink);
            for (const std::string& n : b.inside) check(n);
        }
        if (!missing.empty())
            throw std::runtime_error(
                "panphorte: the bubbles CSV describes nodes the GFA does not contain (" +
                cli::join_with_comma(missing) + "); the CSV and the graph are not the same graph");
    }
    std::unordered_set<std::size_t> bubble_filter(options.bubble_ids.begin(), options.bubble_ids.end());
    // A --bubble-id that names nothing used to run to completion over an empty selection.
    if (!bubble_filter.empty()) {
        std::unordered_set<std::size_t> present;
        for (const Bubble& b : bubbles) present.insert(b.id);
        std::vector<std::string> absent;
        for (const std::size_t want : bubble_filter)
            if (present.find(want) == present.end()) absent.push_back(std::to_string(want));
        if (!absent.empty())
            throw std::runtime_error("panphorte: --bubble-id not present in the bubbles CSV: " +
                                     cli::join_with_comma(absent));
    }
    // Sites are meant to be disjoint. Two that share interior nodes describe the same sequence twice,
    // and folding both rewrites one span inside another -- the same array folded at two scales, which
    // the acceptance check later catches as overlapping edits, but only after every haplotype has been
    // aligned. Caught here instead, naming the pair, since the input is what has to be fixed.
    //
    // Only among the SELECTED bubbles: what is not folded cannot overlap anything, and refusing on a
    // pair the run was never going to touch blocks the one safe way to work at a locus whose CSV
    // carries an overlap -- naming one member of it.
    {
        std::unordered_map<std::string, std::size_t> owner;
        for (const Bubble& b : bubbles) {
            if (!bubble_filter.empty() && bubble_filter.find(b.id) == bubble_filter.end()) continue;
            for (const std::string& n : b.inside) {
                const auto it = owner.find(n);
                if (it != owner.end()) {
                    throw std::runtime_error(
                        "panphorte: bubbles " + std::to_string(it->second) + " and " +
                        std::to_string(b.id) + " both claim interior node " + n +
                        "; overlapping sites describe the same sequence twice and folding both would "
                        "rewrite one span inside the other");
                }
                owner.emplace(n, b.id);
            }
        }
    }

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
    // Nothing lands until the run succeeds: the report and copies table were opened before even the
    // --bubble-id check ran, and the GFA landed before re-snarl or GTF work could fail.
    cli::StagedOutputs staged("panphorte");

    // What each REP node stands for. Two phase-rotated units at one site cannot share an unsplit node
    // while exact spelling is preserved, so they become separate REP nodes -- and without this table a
    // consumer counting REP occurrences sees two independent DUPs instead of one site's copy number.
    //
    // Buffered rather than streamed, because the id a REP is created with is not the id it is delivered
    // under: --reference-path renumbers the graph, so a table written during the loop named nodes that
    // do not exist in the GFA beside it and could not be joined to anything. Both ids are emitted --
    // created_rep_node for the input-side numbering, output_rep_node for the delivered graph.
    struct ProvRecord {
        std::string created_node;
        std::size_t bubble_id = 0;
        std::string canonical_motif;
        std::string phase_unit;
    };
    std::vector<ProvRecord> provenance;

    std::ofstream copies_out;
    if (approximate) {
        copies_out.open(staged.stage(options.out_prefix + ".panphorte.copies.tsv"));
        if (!copies_out) {
            throw std::runtime_error("Failed to write panphorte copies TSV");
        }
        // These ids describe the INPUT graph. Under --reference-path the normalized graph is sorted,
        // which renumbers nodes and reassigns bubble ids -- so the same numbers still exist afterwards
        // and refer to different sequence. Naming them input_* is the difference between a join key and
        // a trap.
        copies_out << "path_name\tsample\tinput_bubble_id\tcopies\tunit_bp\torientations\t"
                   << "mean_identity\tregion_bp\tinput_from_node\tinput_to_node\n";
    }

    // REP nodes are SITE-LOCAL: the key is (bubble id, canonical unit), not the unit alone. Keyed by
    // sequence only, two distinct loci that happen to carry the same repeat unit were folded onto the
    // SAME node -- which joins them in the graph, so a walk can leave one locus and arrive in the
    // other. That corrupts sorting, snarl decomposition and any copy number read from node
    // multiplicity, and identical repeat units at different loci are the normal case, not a corner
    // one.
    std::unordered_map<std::string, std::string> rep_by_unit; // "<bubble>\t<canonical unit>" -> node
    std::unordered_set<std::string> removal_candidates;
    // Oriented links a path traversed across a span this run replaced, and the site that replaced it. A
    // replaced node that survives because some other site still uses it would otherwise keep its old
    // arcs, so the branch the rewrite was supposed to normalize away stays reachable in the delivered
    // graph. The bubble id is carried so a link that cannot be resolved can name where it came from.
    std::unordered_map<std::string, std::size_t> obsolete_edge_keys;
    std::vector<std::vector<PathArray>> per_path_arrays(graph.paths.size());

    // Every consecutive pair the ORIGINAL path walks across [lo, hi), plus the arcs into and out of the
    // span: those are the ones the replacement is meant to supersede.
    const auto record_replaced_span = [&](const PathRecord& p, std::size_t bubble_id, std::size_t lo,
                                          std::size_t hi) {
        if (p.steps.empty()) return;
        const std::size_t first = lo == 0 ? 1 : lo;
        const std::size_t last = std::min(hi + 1, p.steps.size());
        for (std::size_t k = first; k < last; ++k) {
            obsolete_edge_keys.emplace(EdgeSet::key(p.steps[k - 1].node_id,
                                                    orient_char(p.steps[k - 1].reverse),
                                                    p.steps[k].node_id,
                                                    orient_char(p.steps[k].reverse)),
                                       bubble_id);
        }
    };

    PanphorteSummary summary;
    std::ofstream report(staged.stage(options.out_prefix + ".panphorte.report.tsv"));
    if (!report) {
        throw std::runtime_error("Failed to write panphorte report: " + options.out_prefix + ".panphorte.report.tsv");
    }
    // Enough to explain a decision without re-running: how many haplotypes were at the site, how many
    // carried the winning motif, what that fraction was, and the reason the bubble was or was not
    // folded. "normalized=no" on its own gives a reader nothing to act on.
    report << "bubble_id\tnormalized\tunit_bp\tpaths_normalized\tmin_copies\tmax_copies\t"
           << "interruptions_bp\tnodes_collapsed\tn_traversing\tn_motif_carriers\tprevalence\t"
           << "n_motifs\tcopies_declined_partial_boundary\tpaths_with_partial_boundary\tstatus\n";

    std::size_t total_bubbles = 0;
    for (const Bubble& b : bubbles) {
        if (bubble_filter.empty() || bubble_filter.find(b.id) != bubble_filter.end()) ++total_bubbles;
    }
    // Exact mode: a coarse per-bubble bar. Approximate mode shows a richer per-bubble,
    // per-haplotype line during alignment (below), so the coarse bar is disabled there.
    const bool bubble_bar = !options.quiet && !approximate;
    cli::ProgressBar progress(bubble_bar ? "Normalizing bubbles" : std::string(),
                              bubble_bar ? total_bubbles : 0);

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
            // Live, per-bubble progress on stderr (interactive only): makes it obvious
            // whether time goes into seeding or the per-haplotype alignment, and how far
            // the alignment has progressed. Suppressed by --quiet / non-TTY.
            const bool show_detail = !options.quiet && cli::stderr_is_tty();
            if (show_detail) {
                std::cerr << "\r[bubble " << bubble.id << "] seeding...                    "
                          << std::flush;
            }

            // Single-block: seed one representative unit, find its near-identical copies per path (any
            // orientation, banded alignment), collapse to one REP node looped per copy. Lossy --
            // within-copy detail is dropped; inter-copy sequence is kept as literal steps.
            const std::string ref_unit =
                pick_reference_unit(graph, node_tok, path_indexes, bubble, options);
            if (!options.quiet && ref_unit.size() >= options.min_unit_bp) {
                // One concise line per SEEDED bubble (works in non-interactive logs too): what was
                // seeded and how much alignment work it implies, printed BEFORE the expensive
                // per-haplotype alignment so it is visible even if that phase is slow. Bubbles with
                // no foldable unit (sub-threshold seed) stay silent to keep the log concise.
                std::cerr << "[bubble " << bubble.id << "] seeded unit_bp=" << ref_unit.size()
                          << "; aligning across " << graph.paths.size() << " haplotypes\n"
                          << std::flush;
            }
            // Hoisted so the report can state them whether or not the inner block runs.
            std::size_t approx_traversing = 0, approx_carriers = 0;
            std::size_t approx_partial_boundary = 0;   // copies refused for an in-node boundary
            std::size_t approx_partial_paths = 0;      // haplotypes carrying at least one such copy
            const bool seeded = ref_unit.size() >= options.min_unit_bp;
            double approx_prevalence = 0.0;
            if (ref_unit.size() >= options.min_unit_bp) {
                // One replaced step range holding one or more copies. Several copies can share a range
                // -- two copies inside a single node have no step boundary between them -- so the range
                // is emitted as fragment, REP, fragment, REP, ... fragment: `frag` always holds one
                // more entry than there are copies, and the outer two are the bases before the first
                // copy and after the last.
                struct Mapped {
                    std::size_t off_lo = 0, off_hi = 0;
                    std::vector<bool> rev;
                    std::vector<double> id;
                    std::vector<std::string> frag;
                    std::size_t last_bp_hi = 0;   // end of the most recent copy, in interval bp
                    std::size_t n_copies() const { return rev.size(); }
                };
                struct Res { bool has = false; std::size_t left = 0; std::vector<Mapped> copies; };
                std::vector<Res> results(graph.paths.size());

                std::atomic<std::size_t> detect_done{0};
                std::atomic<std::size_t> n_traverse{0};   // haplotypes that traverse this bubble
                std::size_t partial_boundary_declined = 0;  // copies refused for a mid-node boundary
                std::size_t partial_boundary_paths = 0;     // haplotypes carrying at least one
                std::size_t zero_after_decline = 0;         // haplotypes left with NO foldable copy
                std::size_t approx_interruptions_bp = 0;    // real interrupting bases within arrays
                std::size_t group_carriers = 0;             // paths with an array reaching min_copies
                std::mutex detail_mtx;
                const std::size_t n_paths = graph.paths.size();

                auto detect_one = [&](std::size_t pi) {
                    // bubble_steps, not the inside-node-only interval finder: a haplotype crossing the
                    // site with no interior carries ZERO copies and is still a haplotype at this site.
                    // The exact branch was fixed for this; the approximate branch kept the old finder,
                    // so one array carrier among nine direct alleles still read prevalence 1/1.
                    BubblePathInterval used{};
                    const auto steps_opt = bubble_steps(graph.paths[pi], path_indexes[pi], bubble, &used);
                    if (!steps_opt.has_value()) return;
                    n_traverse.fetch_add(1);
                    if (used.inside_count == 0) return;      // crosses with no copies
                    // Search the INTERIOR only. Including the source and sink steps let the anchors --
                    // which are backbone, not bubble content -- contribute copies, and a long periodic
                    // backbone is exactly where that goes wrong.
                    const std::size_t left = used.left + 1;
                    if (used.right <= left) return;
                    const std::size_t n = used.right - left;
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
                    // A copy boundary rarely lands on a node boundary. Rounding to the nearest one
                    // and replacing whole steps DELETED the bases between the node edge and the copy
                    // edge -- sequence outside the copy, which approximate mode may not discard (it may
                    // drop differences WITHIN a copy; that is a different claim).
                    //
                    // Declining those copies instead was worse: a declined copy at a tandem array is a
                    // whole repeat unit missing from the copy number, and the bases at stake are a
                    // handful. So the containing step range is taken -- floor for the start, ceiling
                    // for the end -- and the bases inside it that fall outside the copy are kept as
                    // literal fragments, emitted around the REP step at rewrite time. Nothing is lost
                    // and nothing is undercounted.
                    const auto floor_off = [&](std::size_t bp) {
                        std::size_t lo = 0;
                        while (lo + 1 < prefix.size() && prefix[lo + 1] <= bp) ++lo;
                        return lo;
                    };
                    const auto ceil_off = [&](std::size_t bp) {
                        std::size_t hi = 0;
                        while (hi < prefix.size() && prefix[hi] < bp) ++hi;
                        return hi < prefix.size() ? hi : prefix.size() - 1;
                    };
                    // Copies are grouped into BLOCKS by the step range that contains them. Two copies
                    // inside one node share a range, and mapping each copy to a range independently
                    // declined the second one (its floor collided with the first one's ceiling) --
                    // which at an array is a whole repeat unit missing from the copy number.
                    std::vector<Mapped> copies;
                    std::size_t declined = 0;
                    for (const ApproxCopy& c : bp_copies) {
                        const std::size_t lo = floor_off(c.bp_lo);
                        const std::size_t hi = ceil_off(c.bp_hi);
                        if (hi <= lo) { ++declined; continue; }
                        if (!copies.empty() && lo < copies.back().off_hi) {
                            // Shares the open block: extend it and record the gap since the last copy.
                            Mapped& blk = copies.back();
                            if (hi > blk.off_hi) blk.off_hi = hi;
                            const std::size_t gap_lo = blk.last_bp_hi;
                            const std::size_t gap_hi = std::max(c.bp_lo, gap_lo);
                            blk.frag.push_back(snat.substr(gap_lo, gap_hi - gap_lo));
                            blk.rev.push_back(c.rev);
                            blk.id.push_back(c.identity);
                            blk.last_bp_hi = std::max(c.bp_hi, gap_hi);
                            continue;
                        }
                        Mapped blk;
                        blk.off_lo = lo;
                        blk.off_hi = hi;
                        const std::size_t range_lo = prefix[lo];
                        const std::size_t copy_lo = std::max(c.bp_lo, range_lo);
                        blk.frag.push_back(snat.substr(range_lo, copy_lo - range_lo));
                        blk.rev.push_back(c.rev);
                        blk.id.push_back(c.identity);
                        blk.last_bp_hi = std::min(c.bp_hi, prefix[hi]);
                        copies.push_back(std::move(blk));
                    }
                    // Close each block with the bases between its last copy and the end of its range.
                    for (Mapped& blk : copies) {
                        const std::size_t range_hi = prefix[blk.off_hi];
                        const std::size_t tail_lo = std::min(blk.last_bp_hi, range_hi);
                        blk.frag.push_back(snat.substr(tail_lo, range_hi - tail_lo));
                    }
                    if (declined > 0) {
                        std::lock_guard<std::mutex> lk(detail_mtx);
                        partial_boundary_declined += declined;
                        ++partial_boundary_paths;
                        // The condition that actually breaks `call`: this haplotype carries the motif
                        // but ends with NOTHING foldable, so it reaches the site literally while its
                        // neighbours reach it through the REP node -- and REP occurrences are what copy
                        // number is counted from, so it reads 0 while carrying copies. A haplotype that
                        // folds some copies and declines others is merely undercounted by those copies,
                        // which is a different and much smaller harm.
                        if (copies.empty()) ++zero_after_decline;
                    }
                    // Copies are grouped into tandem ARRAYS by the gap between them, and
                    // --max-interruption-frac decides where a group ends. The option was applied when
                    // seeding and nowhere afterwards, so two copies with a 64 bp gap between them were
                    // accepted as one array at every setting and interruptions_bp reported 0.
                    std::size_t path_interruptions = 0;
                    std::size_t best_group = 0;
                    {
                        // Group sizes count COPIES, not blocks: a block sharing one node can hold
                        // several, and they are adjacent within it by construction.
                        std::size_t gstart = 0;
                        std::size_t ginterrupt = 0;
                        std::size_t group_copies = 0;
                        const auto blk_copies = [&](std::size_t a, std::size_t b) {
                            std::size_t n = 0;
                            for (std::size_t q = a; q < b; ++q) n += copies[q].n_copies();
                            return n;
                        };
                        (void)group_copies;
                        for (std::size_t k = 0; k < copies.size(); ++k) {
                            if (k > 0) {
                                const std::size_t gap =
                                    prefix[copies[k].off_lo] - prefix[copies[k - 1].off_hi];
                                const std::size_t span =
                                    prefix[copies[k].off_hi] - prefix[copies[gstart].off_lo];
                                const double frac = span > 0
                                    ? static_cast<double>(ginterrupt + gap) / static_cast<double>(span)
                                    : 0.0;
                                if (frac > options.max_interruption_frac) {
                                    best_group = std::max(best_group, blk_copies(gstart, k));
                                    path_interruptions += ginterrupt;
                                    gstart = k;
                                    ginterrupt = 0;
                                    continue;
                                }
                                ginterrupt += gap;
                            }
                        }
                        best_group = std::max(best_group, blk_copies(gstart, copies.size()));
                        path_interruptions += ginterrupt;
                    }
                    {
                        std::lock_guard<std::mutex> lk(detail_mtx);
                        approx_interruptions_bp += path_interruptions;
                        if (best_group >= options.min_copies) ++group_carriers;
                    }
                    // Record any path with >=1 detected copy. The bubble-level min_copies gate is
                    // applied below: once an array is confirmed (some path reaches min_copies), even
                    // single-copy haplotypes are folded through the REP node so their copy number is
                    // 1 rather than 0 (they would otherwise bypass the REP node entirely).
                    if (!copies.empty()) {
                        results[pi] = Res{true, left, std::move(copies)};
                    }
                    if (show_detail) {
                        const std::size_t d = detect_done.fetch_add(1) + 1;
                        if ((d & 0x1F) == 0 || d == n_paths) { // throttle redraws
                            std::lock_guard<std::mutex> lk(detail_mtx);
                            std::cerr << "\r[bubble " << bubble.id << "] unit_bp=" << ref_unit.size()
                                      << " aligning " << d << '/' << n_paths << " haplotypes      "
                                      << std::flush;
                        }
                    }
                };

                if (show_detail) {
                    std::cerr << "\r[bubble " << bubble.id << "] unit_bp=" << ref_unit.size()
                              << " aligning 0/" << n_paths << " haplotypes      " << std::flush;
                }
                run_parallel(graph.paths.size(), options.threads, detect_one);
                if (show_detail) {
                    std::cerr << '\n';
                }

                // Confirm a real tandem array before folding: some path must carry >= min_copies copies.
                // Then a cohort-prevalence gate -- require >= min_array_prevalence of traversing
                // haplotypes to carry one -- so a rare private duplication (a CYP2D6x2 in a few
                // haplotypes) isn't folded, collapsing genes `call` resolves downstream.
                // Carriers are counted by GROUP, not by total copy count: two copies separated by more
                // than --max-interruption-frac are two arrays of one, not one array of two.
                std::size_t max_copies_path = 0;
                for (const Res& r : results) {
                    if (!r.has) continue;
                    std::size_t nc = 0;
                    for (const auto& blk : r.copies) nc += blk.n_copies();
                    max_copies_path = std::max(max_copies_path, nc);
                }
                const std::size_t traversing = n_traverse.load();
                const double prevalence =
                    traversing > 0 ? static_cast<double>(group_carriers) / static_cast<double>(traversing) : 0.0;
                // SITE-WIDE refusal, not per copy. Declining one copy while folding the rest leaves a
                // mixed representation: some haplotypes reach the motif through the REP node and one
                // still spells it literally. `call` counts REP occurrences, so that haplotype is
                // reported CN=0 with the literal node as an insertion, when it carries one copy. The
                // sequence was safe and the CALL was not, so the whole bubble is left alone.
                // Refuse the site only where the mixed representation would produce a FALSE ZERO, not
                // wherever any boundary happens to fall mid-node. Blanket refusal cost the principal
                // case for nothing: at a real tandem array every haplotype still folds most of its
                // copies (measured at KIV-2: 466 of 466 haplotypes fold at least one, median 20), so
                // refusing the whole site to prevent an undercount of one copy in twenty threw away
                // the fold for all 466.
                const bool partial_boundary = zero_after_decline > 0 && !options.allow_partial_boundary;
                const bool array_confirmed = !partial_boundary &&
                                             max_copies_path >= options.min_copies &&
                                             prevalence >= options.min_array_prevalence;
                approx_traversing = traversing;
                approx_carriers = group_carriers;
                approx_partial_boundary = partial_boundary_declined;
                approx_partial_paths = partial_boundary_paths;
                if (partial_boundary_declined > 0 && zero_after_decline == 0 && !options.quiet) {
                    std::cerr << "[bubble " << bubble.id << "] " << partial_boundary_declined
                              << " copy/copies declined (boundary inside a node) across "
                              << partial_boundary_paths
                              << " haplotype(s); each still folds other copies, so copy number is "
                                 "undercounted by those copies rather than zeroed\n";
                }
                if (zero_after_decline > 0 && options.allow_partial_boundary && !options.quiet) {
                    std::cerr << "[bubble " << bubble.id << "] --allow-partial-boundary: folding with "
                              << zero_after_decline << " haplotype(s) left with no foldable copy. They "
                                 "keep the motif literally while the rest reach it through the REP "
                                 "node, so `call` will report them CN 0 despite carrying copies\n";
                }
                if (partial_boundary && !options.quiet) {
                    std::cerr << "[bubble " << bubble.id << "] not normalized: " << zero_after_decline
                              << " haplotype(s) carry the motif but have NO copy that can be folded "
                                 "(every boundary falls inside a node). Folding the rest would leave "
                                 "them literal while their neighbours use the REP node, and `call` "
                                 "counts REP occurrences -- so they would read CN 0 while carrying "
                                 "copies. --allow-partial-boundary folds anyway\n";
                }
                // The report's interruptions column was never fed by the approximate branch, so it
                // read 0 whatever the arrays actually contained.
                interruptions_bp = approx_interruptions_bp;
                approx_prevalence = prevalence;

                // Serial: one REP node per bubble; build PathArrays + copies.tsv rows.
                std::string rep_id;
                for (std::size_t pi = 0; array_confirmed && pi < results.size(); ++pi) {
                    if (!results[pi].has) continue;
                    const std::size_t left = results[pi].left;
                    const std::vector<Mapped>& copies = results[pi].copies;
                    if (rep_id.empty()) {
                        const std::string rep_key = std::to_string(bubble.id) + "\t" + ref_unit;
                        const auto rit = rep_by_unit.find(rep_key);
                        if (rit == rep_by_unit.end()) {
                            rep_id = add_new_node(model, ref_unit);
                            rep_by_unit.emplace(rep_key, rep_id);
                            ++summary.nodes_added;
                            provenance.push_back(
                                ProvRecord{rep_id, bubble.id, canonical_motif_key(ref_unit), ref_unit});
                        } else {
                            rep_id = rit->second;
                        }
                    }

                    PathArray pa;
                    pa.bubble_id = bubble.id;
                    pa.start = left + copies.front().off_lo;
                    std::size_t prev_hi = copies.front().off_lo;
                    std::string orients;
                    double id_sum = 0.0;
                    for (const Mapped& m : copies) {
                        for (std::size_t off = prev_hi; off < m.off_lo; ++off) {
                            pa.replacement.push_back(graph.paths[pi].steps[left + off]);
                        }
                        // fragment, REP, fragment, REP, ... fragment. The fragments carry every base
                        // in the replaced range that is not inside a copy, verbatim, as
                        // per-occurrence nodes: nothing outside an accepted copy is lost, and copies
                        // that share a node are all folded rather than the second one declined.
                        for (std::size_t ci = 0; ci < m.n_copies(); ++ci) {
                            if (!m.frag[ci].empty()) {
                                pa.replacement.push_back(PathStep{add_new_node(model, m.frag[ci]), false});
                                ++summary.nodes_added;
                                ++summary.fragment_nodes_added;
                            }
                            pa.replacement.push_back(PathStep{rep_id, m.rev[ci]});
                        }
                        if (!m.frag.back().empty()) {
                            pa.replacement.push_back(PathStep{add_new_node(model, m.frag.back()), false});
                            ++summary.nodes_added;
                            ++summary.fragment_nodes_added;
                        }
                        // The steps this block replaces are gone from the path. Marking them (and the
                        // arcs across them) makes them candidates for removal; whether they actually
                        // go is decided globally, since another site may still use them. Without this
                        // the approximate branch left every replaced node and link in the graph beside
                        // its own normalized route -- a dead branch the re-snarl then picked up -- and
                        // reported nodes_collapsed=0 while doing it.
                        for (std::size_t off = m.off_lo; off < m.off_hi; ++off) {
                            const std::string& nid = graph.paths[pi].steps[left + off].node_id;
                            removal_candidates.insert(nid);
                            collapsed_nodes.insert(nid);
                        }
                        record_replaced_span(graph.paths[pi], bubble.id, left + m.off_lo, left + m.off_hi);
                        for (std::size_t ci = 0; ci < m.n_copies(); ++ci) {
                            orients += (m.rev[ci] ? '-' : '+');
                            id_sum += m.id[ci];
                        }
                        prev_hi = m.off_hi;
                    }
                    pa.length = (left + prev_hi) - pa.start;
                    per_path_arrays[pi].push_back(std::move(pa));
                    ++paths_norm;

                    std::size_t occ = 0;
                    for (const Mapped& m : copies) occ += m.n_copies();
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
                               << graph.paths[pi].steps[left + copies.back().off_hi - 1].node_id << '\n';
                }
            }

            const bool norm_ok = paths_norm > 0;
            if (norm_ok) ++summary.bubbles_normalized;
            report << bubble.id << '\t' << (norm_ok ? "yes" : "no") << '\t'
                   << unit_bp_report << '\t' << paths_norm << '\t'
                   << min_copies_seen << '\t' << max_copies_seen << '\t'
                   << interruptions_bp << '\t' << collapsed_nodes.size() << '\t'
                   // Without a seed nothing was counted, so these are unknown rather than zero. A
                   // bubble reported "traversing=0" that 466 haplotypes cross is a misleading number,
                   // and I introduced it with these columns.
                   << (seeded ? std::to_string(approx_traversing) : ".") << '\t'
                   << (seeded ? std::to_string(approx_carriers) : ".") << '\t'
                   << (seeded ? std::to_string(approx_prevalence) : ".") << '\t'
                   << (ref_unit.empty() ? 0 : 1) << '\t' << approx_partial_boundary << '\t'
                   << approx_partial_paths << '\t'
                   << (norm_ok ? "normalized"
                       : approx_partial_boundary > 0 ? "partial_boundary"
                       : ref_unit.size() < options.min_unit_bp ? "no_seed" : "below_prevalence")
                   << '\n';
            progress.tick();
            continue; // bubble fully handled in approximate mode
        }

        // Detect per-path tandem arrays first, then apply the SAME cohort-prevalence gate as the
        // approximate branch: fold the bubble only if a >= min_copies array is carried by
        // >= --min-array-prevalence of the traversing haplotypes. This keeps a rare private
        // duplication of a gene/segmental module (e.g. a CYP2D6x2 in a handful of haplotypes) from
        // being collapsed in exact mode too, so the gate behaves identically at every similarity.
        std::vector<std::vector<TandemArray>> path_arrays(graph.paths.size());
        std::size_t n_traverse = 0;
        // Carriers of each candidate motif, and the motif's own copy-number range, so prevalence can be
        // asked per motif rather than of "any tandem at all".
        std::unordered_map<std::string, std::size_t> carriers_by_motif;
        for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
            // bubble_steps, not the inside-node-only interval finder: a haplotype crossing the site
            // with no interior node carries ZERO copies, and it is still a haplotype at this site.
            // Counting only paths with an interior made the denominator the set of carriers, so one
            // array carrier among nine deletion alleles read as prevalence 1/1 instead of 1/10 and was
            // folded straight through a 50% gate.
            BubblePathInterval used{};
            const auto steps_opt = bubble_steps(graph.paths[pi], path_indexes[pi], bubble, &used);
            if (!steps_opt.has_value()) continue;
            ++n_traverse;
            if (used.inside_count == 0) continue;             // crosses with 0 copies
            auto arrays = detect_tandems(
                graph, node_tok, graph.paths[pi], used.left, used.right,
                options.min_unit_bp, options.min_copies, options.max_interruption_frac);
            if (arrays.empty()) continue;
            std::unordered_set<std::string> motifs_here;
            for (const TandemArray& a : arrays) motifs_here.insert(canonical_motif_key(a.unit_seq));
            for (const std::string& m : motifs_here) ++carriers_by_motif[m];
            path_arrays[pi] = std::move(arrays);
        }
        // Prevalence is asked of each MOTIF on its own. Two unrelated repeats in one bubble are two
        // separate questions -- neither lends its support to the other -- and a motif that clears the
        // gate is folded into its own site-local REP.
        std::unordered_set<std::string> confirmed_motifs;
        double prevalence = 0.0;
        for (const auto& [m, c] : carriers_by_motif) {
            const double p = n_traverse > 0 ? static_cast<double>(c) / static_cast<double>(n_traverse) : 0.0;
            if (p >= options.min_array_prevalence) confirmed_motifs.insert(m);
            if (p > prevalence) prevalence = p;
        }
        const bool array_confirmed = !confirmed_motifs.empty();
        // Once the site is confirmed, every haplotype carrying at least ONE copy of the chosen motif
        // folds, not only those reaching --min-copies. --min-copies confirms the SITE; it is not a
        // per-haplotype threshold. Leaving one-copy haplotypes on the original nodes contradicted the
        // documented contract and left them outside the REP node, so a copy number read from REP
        // multiplicity called them 0 when they carry 1.
        if (array_confirmed) {
            for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
                if (!path_arrays[pi].empty()) continue;       // already has >= min_copies
                BubblePathInterval used{};
                const auto steps_opt = bubble_steps(graph.paths[pi], path_indexes[pi], bubble, &used);
                if (!steps_opt.has_value() || used.inside_count == 0) continue;
                for (const std::string& m : confirmed_motifs) {
                    auto singles = find_motif_copies_by_node(
                        graph, graph.paths[pi], used.left, used.right, m, options.min_unit_bp);
                    for (TandemArray& a : singles) path_arrays[pi].push_back(std::move(a));
                }
            }
        }
        for (std::size_t pi = 0; array_confirmed && pi < graph.paths.size(); ++pi) {
            if (path_arrays[pi].empty()) continue;
            for (const TandemArray& arr : path_arrays[pi]) {
                // Only the site's motif is folded. Without this a second, unrelated array in the same
                // bubble got its own REP node on the strength of the first one's prevalence.
                if (confirmed_motifs.find(canonical_motif_key(arr.unit_seq)) == confirmed_motifs.end())
                    continue;
                const std::string& unit = arr.unit_seq;
                const std::string rc = reverse_complement(unit);
                const bool forward = unit <= rc;
                const std::string canonical = forward ? unit : rc;
                const bool rev = !forward;

                const std::string rep_key = std::to_string(bubble.id) + "\t" + canonical;
                auto rit = rep_by_unit.find(rep_key);
                std::string rep_id;
                if (rit == rep_by_unit.end()) {
                    rep_id = add_new_node(model, canonical);
                    rep_by_unit.emplace(rep_key, rep_id);
                    ++summary.nodes_added;
                    provenance.push_back(
                        ProvRecord{rep_id, bubble.id, canonical_motif_key(unit), canonical});
                } else {
                    rep_id = rit->second;
                }

                PathArray pa;
                pa.bubble_id = bubble.id;
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
                    record_replaced_span(graph.paths[pi], bubble.id, r.first, r.second);
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
        std::size_t top_carriers = 0;
        for (const auto& [_m, c] : carriers_by_motif) top_carriers = std::max(top_carriers, c);
        const char* status = normalized ? "normalized"
                             : carriers_by_motif.empty() ? "no_tandem_detected"
                             : confirmed_motifs.empty() ? "below_prevalence"
                                                        : "no_copies_rewritten";
        report << bubble.id << '\t' << (normalized ? "yes" : "no") << '\t'
               << unit_bp_report << '\t' << paths_norm << '\t'
               << min_copies_seen << '\t' << max_copies_seen << '\t'
               << interruptions_bp << '\t' << collapsed_nodes.size() << '\t'
               << n_traverse << '\t' << top_carriers << '\t' << prevalence << '\t'
               << carriers_by_motif.size() << '\t' << 0 << '\t' << 0 << '\t' << status << '\n';
        progress.tick();
    }
    progress.done();

    // Splice replacements into each path's steps.
    for (std::size_t pi = 0; pi < model.paths.size(); ++pi) {
        auto& arrays = per_path_arrays[pi];
        if (arrays.empty()) {
            continue;
        }
        std::sort(arrays.begin(), arrays.end(),
                  [](const PathArray& a, const PathArray& b) { return a.start < b.start; });

        // ACCEPTANCE. Exact mode proves whole-sequence equality afterwards; approximate mode
        // deliberately changes sequence INSIDE accepted copies, so what it can prove is structural:
        // the edits are ordered, non-overlapping and in bounds, and the spliced result is exactly the
        // original with those spans -- and only those spans -- substituted. A general backstop: the
        // boundary-deletion bug would have failed here rather than in the one fixture covering it.
        std::string expected_after;
        {
            const std::vector<PathStep>& old_steps = model.paths[pi].steps;
            const auto spell_steps = [&](const std::vector<PathStep>& st, std::size_t lo, std::size_t hi) {
                std::string acc;
                for (std::size_t k = lo; k < hi; ++k) {
                    const auto it = model.seq.find(st[k].node_id);
                    if (it == model.seq.end()) continue;
                    acc += st[k].reverse ? reverse_complement(it->second) : it->second;
                }
                return acc;
            };
            std::size_t prev_end = 0;
            for (const PathArray& a : arrays) {
                if (a.start < prev_end)
                    throw std::runtime_error("panphorte: overlapping edits on path " + graph.paths[pi].name);
                if (a.start + a.length > old_steps.size())
                    throw std::runtime_error("panphorte: edit runs past the end of path " +
                                             graph.paths[pi].name);
                if (a.length == 0 || a.replacement.empty())
                    throw std::runtime_error("panphorte: empty edit on path " + graph.paths[pi].name);
                expected_after += spell_steps(old_steps, prev_end, a.start);          // outside, verbatim
                expected_after += spell_steps(a.replacement, 0, a.replacement.size()); // the substitution
                prev_end = a.start + a.length;
            }
            expected_after += spell_steps(old_steps, prev_end, old_steps.size());
        }

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
        {
            std::string got;
            for (const PathStep& st : model.paths[pi].steps) {
                const auto it = model.seq.find(st.node_id);
                if (it == model.seq.end()) continue;
                got += st.reverse ? reverse_complement(it->second) : it->second;
            }
            if (got != expected_after)
                throw std::runtime_error(
                    "panphorte: the rewrite of path " + graph.paths[pi].name +
                    " is not the original with the accepted copies substituted (" +
                    std::to_string(expected_after.size()) + " bp expected, " +
                    std::to_string(got.size()) + " bp produced)");
        }
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
    // The arcs the FINAL paths walk, under the same canonical key as everything else, so a step pair
    // b-,a- is recognised as the stored a+>b+ rather than treated as a different link.
    std::unordered_set<std::string> traversed_edges;
    for (const GfaPath& p : model.paths) {
        for (const PathStep& s : p.steps) {
            referenced.insert(s.node_id);
        }
        for (std::size_t k = 1; k < p.steps.size(); ++k) {
            const PathStep& a = p.steps[k - 1];
            const PathStep& b = p.steps[k];
            traversed_edges.insert(
                EdgeSet::key(a.node_id, orient_char(a.reverse), b.node_id, orient_char(b.reverse)));
        }
    }
    std::unordered_set<std::string> to_remove;
    for (const std::string& id : removal_candidates) {
        if (referenced.find(id) == referenced.end()) {
            to_remove.insert(id);
        }
    }
    {
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
            if (to_remove.find(e.from) != to_remove.end() || to_remove.find(e.to) != to_remove.end()) {
                continue;
            }
            // An arc only a pre-rewrite path used, between two nodes that both survive for other
            // reasons. Removing the nodes is not enough on its own: the link is what keeps the old
            // route walkable, and the re-snarl reads links, not paths.
            // A link is dropped only when NO final path walks it. One a rewritten span used that some
            // path still walks stays: a GFA link is global and deleting it would break that path. At an
            // array that is ordinary rather than exceptional -- LPA keeps 2009 of them, because an arc
            // inside the repeat unit is crossed once per copy and folding replaces all but the one
            // crossing that falls outside the folded span. Whether that matters is not a property of
            // any single link, so it is not judged here; the acceptance check below asks the question
            // that does matter.
            const std::string k = EdgeSet::key(e.from, e.from_orient, e.to, e.to_orient);
            if (obsolete_edge_keys.find(k) != obsolete_edge_keys.end() &&
                traversed_edges.find(k) == traversed_edges.end()) {
                ++summary.edges_removed;
                continue;
            }
            kept_edges.push_back(e);
        }
        model.edges = std::move(kept_edges);
    }
    // ACCEPTANCE, at the level the question is actually asked: is the route a path used BEFORE the
    // rewrite still walkable AFTER it? A link that survives because another path needs it is not on its
    // own a defect -- at a real array that is the normal case, 2009 such links at LPA -- and refusing on
    // the count would block the locus for nothing. What matters is whether the arcs and nodes that
    // survive still spell the old allele end to end.
    //
    // Asked PER EDIT, not per path. A haplotype normalized at two sites has two independent old routes,
    // and testing the whole original walk conflates them: the first site's route being properly gone
    // makes the walk unwalkable and hides the second site's route surviving intact. Each replaced span
    // is therefore checked on its own -- the arc into it, its interior, and the arc out of it, which
    // together are what make the old branch a detour a walk can actually take.
    {
        std::unordered_set<std::string> live_nodes(model.node_order.begin(), model.node_order.end());
        std::unordered_set<std::string> live_edges;
        live_edges.reserve(model.edges.size() * 2);
        for (const GfaEdge& e : model.edges)
            live_edges.insert(EdgeSet::key(e.from, e.from_orient, e.to, e.to_orient));
        std::vector<std::string> survivors;
        for (std::size_t i = 0; i < model.paths.size(); ++i) {
            const std::vector<PathStep>& was = graph.paths[i].steps;
            const auto step_live = [&](std::size_t k) {
                return live_nodes.find(was[k].node_id) != live_nodes.end();
            };
            const auto arc_live = [&](std::size_t k) {   // the arc from step k-1 to step k
                return live_edges.find(EdgeSet::key(was[k - 1].node_id, orient_char(was[k - 1].reverse),
                                                    was[k].node_id, orient_char(was[k].reverse))) !=
                       live_edges.end();
            };
            for (const PathArray& a : per_path_arrays[i]) {
                const std::size_t lo = a.start, hi = a.start + a.length;   // [lo, hi) was replaced
                bool walkable = true;
                for (std::size_t k = lo; walkable && k < hi; ++k) {
                    if (!step_live(k)) walkable = false;
                    else if (k > lo && !arc_live(k)) walkable = false;
                }
                // The anchors: without an arc in and an arc out the old interior is stranded, not a
                // route anything can walk.
                if (walkable && lo > 0 && !(step_live(lo - 1) && arc_live(lo))) walkable = false;
                if (walkable && hi < was.size() && !(step_live(hi) && arc_live(hi))) walkable = false;
                if (walkable) {
                    ++summary.routes_surviving;
                    if (survivors.size() < 4)
                        survivors.push_back(graph.paths[i].name + " (bubble " +
                                            std::to_string(a.bubble_id) + ")");
                }
            }
        }
        // No override. This is an internal invariant, not a choice about how to fold: a graph carrying
        // the same site both folded and unfolded is wrong for every consumer, and the answer is to fix
        // whatever produced it, not to accept it. Nothing in the six reference loci reaches it.
        if (summary.routes_surviving > 0) {
            throw std::runtime_error(
                "panphorte: " + std::to_string(summary.routes_surviving) +
                " rewritten path(s) can still be walked along their ORIGINAL route in the normalized "
                "graph (" + cli::join_with_comma(survivors) +
                "), so the site is represented twice -- once folded and once as the branch it replaced, "
                "which the re-snarl will pick up as a second allele. This happens when every node and "
                "link of the replaced span is also used somewhere else, and separating the two would "
                "need the reused context cloned before the rewrite, which is not implemented");
        }
    }

    // old REP id -> the id it is delivered under. Empty means the ids did not move.
    std::unordered_map<std::string, std::string> rep_id_remap;

    // When a reference is given we make the output call-ready with no external tools: internally
    // sort+flip along the reference (placing the appended REP nodes into reference order) and
    // re-snarl with the cactus finder, writing only the sorted GFA + a bubbles CSV + Bandage colors.
    // The unsorted normalized GFA is written only when no reference is supplied (nothing to sort by).
    if (options.reference_path.empty()) {
        write_gfa_model(staged.stage(options.out_prefix + ".normalized.gfa"), model);
    } else {
        GraphSortOptions sort_opts;
        sort_opts.reference_path = options.reference_path;
        sort_opts.flip = !options.no_flip;
        // The sorter resolves an exact name, a unique case-insensitive match or a unique substring, and
        // everything after this compares path names EXACTLY. Discarding the resolved name meant
        // `-r FULL` for a path named `full` sorted correctly and then oriented nothing: the re-snarled
        // bubbles came back reversed with ref_allele_support 0.
        const GraphSortResult sort_result = sort_graph_reference(model, sort_opts);
        const std::string resolved_reference = sort_result.resolved_reference;
        rep_id_remap = sort_result.id_remap;

        const std::string sorted_gfa = staged.stage(options.out_prefix + ".normalized.sorted.gfa");
        write_gfa_model(sorted_gfa, model);

        BubbleCallOptions bopts;
        bopts.reference_path = resolved_reference;
        // The re-snarl used BubbleCallOptions' own default (50 bp), silently re-filtering a CSV whose
        // bubbles may have been produced under a different threshold -- so a small bubble handed to
        // panphorte could disappear from the call-ready CSV with nothing said. The threshold is an
        // option now, and any input bubble missing afterwards is reported.
        bopts.min_variant_bp = options.resnarl_min_variant_bp;
        bopts.snarl_pairs_override = find_top_level_snarls_cactus(snarl_input_from_model(model));
    bopts.snarl_source_supplied = true;
        bopts.quiet = options.quiet;

        ParseGfaOptions parse_options;
        parse_options.include_paths = true;
        parse_options.include_sequences = true;
        const Graph sorted_graph = parse_gfa(sorted_gfa, parse_options);
        const BubbleCallReport report = call_bubbles_report(sorted_graph, bopts);
        write_bubbles_csv(staged.stage(options.out_prefix + ".bubbles.csv"), report.bubbles);
        // Bubble ids are REASSIGNED by the re-snarl, so comparing input ids against output ids says
        // nothing: id 3 afterwards is not id 3 before. Only the count is comparable without remapping
        // anchors through the sort, so only the count is reported.
        if (!options.quiet) {
            // Total in against total out. In targeted mode the output still contains the WHOLE graph,
            // so comparing the selected count against every re-snarled bubble was comparing a subset
            // with a superset and always looked like loss.
            const std::size_t before = bubbles.size();
            const std::size_t selected = bubble_filter.empty()
                ? bubbles.size()
                : std::count_if(bubbles.begin(), bubbles.end(), [&](const Bubble& b) {
                      return bubble_filter.find(b.id) != bubble_filter.end();
                  });
            if (report.bubbles.size() != before) {
                std::cerr << "[panphorte] bubble count " << before << " -> " << report.bubbles.size()
                          << " after re-snarling (normalization merges sites, and "
                             "--resnarl-min-variant-bp " << options.resnarl_min_variant_bp
                          << " filters them); ids are reassigned, so they do not correspond"
                          << (bubble_filter.empty() ? std::string()
                                                    : " [" + std::to_string(selected) + " selected]")
                          << "\n";
            }
        }
        write_bandage_node_colors_csv(staged.stage(options.out_prefix + ".bandage_nodes.csv"),
                                      report.bubbles, report.non_snp_bubbles);

        // Gene annotation belongs inside the transaction. Run by the caller after this function
        // returned, a malformed GTF failed with the whole normalized family already on disk -- so the
        // one output family the staging exists to keep whole was not whole.
        if (!options.gtf_path.empty()) {
            summary.genes_written = emit_gene_annotation(
                sorted_graph, resolved_reference, options.gtf_path,
                staged.stage(options.out_prefix + ".bandage_genes.csv"));
        }

        summary.sorted = true;
        summary.resnarled_bubbles = report.bubbles.size();
    }

    // The REP ids, mapped through the sort so they name nodes in the graph delivered beside this file.
    {
        std::ofstream prov(staged.stage(options.out_prefix + ".panphorte.rep_provenance.tsv"));
        if (!prov) throw std::runtime_error("panphorte: cannot write the REP provenance table");
        // Two unit sequences, because the sort may FLIP a node: created_phase_unit is what panphorte
        // built, output_phase_unit is what the delivered node actually spells. One column called
        // phase_unit could not be both, and a consumer comparing it against the graph would have been
        // wrong for every flipped REP. canonical_motif is the min over both strands, so it is the same
        // either way and is the safe key for grouping phase-rotated units at one site.
        prov << "created_rep_node\toutput_rep_node\tinput_bubble_id\tcanonical_motif\t"
                "created_phase_unit\toutput_phase_unit\tunit_bp\tcopy_quantum\n";
        for (const ProvRecord& r : provenance) {
            const auto it = rep_id_remap.find(r.created_node);
            const std::string out_id = it == rep_id_remap.end() ? r.created_node : it->second;
            const auto sit = model.seq.find(out_id);
            const std::string out_unit = sit == model.seq.end() ? r.phase_unit : sit->second;
            prov << r.created_node << '\t' << out_id << '\t' << r.bubble_id << '\t'
                 << r.canonical_motif << '\t' << r.phase_unit << '\t' << out_unit << '\t'
                 << r.phase_unit.size() << '\t' << 1 << '\n';
        }
        prov.close();
        if (!prov) throw std::runtime_error("panphorte: failed to write the REP provenance table");
    }

    // Everything succeeded: move the family into place. Closed and checked, not just flushed -- a
    // failure that surfaces on close would otherwise be committed as a complete file.
    report.close();
    if (!report) throw std::runtime_error("panphorte: failed to write the report table");
    if (copies_out.is_open()) {
        copies_out.close();
        if (!copies_out) throw std::runtime_error("panphorte: failed to write the copies table");
    }
    staged.commit();

    if (summary_out != nullptr) {
        *summary_out = summary;
    }
}

} // namespace panvar
