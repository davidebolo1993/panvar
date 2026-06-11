#include "panvar/variant_call.hpp"

#include "panvar/align.hpp"
#include "panvar/bubble_path.hpp"
#include "panvar/bubbles.hpp"
#include "panvar/cli_utils.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/minimap2_align.hpp"
#include "panvar/output.hpp"
#include "panvar/ref_path.hpp"

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace panvar {
namespace {

constexpr std::size_t kAlignCellCap = 4000000; // skip walk-diff if |ref|*|hap| exceeds this

enum class EvType { Del, Ins, Inv, Dup };

const char* ev_svtype(EvType t) {
    switch (t) {
        case EvType::Del: return "DEL";
        case EvType::Ins: return "INS";
        case EvType::Inv: return "INV";
        case EvType::Dup: return "DUP";
    }
    return "DEL";
}

std::size_t node_len(const Graph& graph, const std::string& id) {
    const auto it = graph.nodes.find(id);
    return it == graph.nodes.end() ? 0 : it->second.sequence.size();
}

// Nodes carrying a self-loop edge (a panphorte REP, or any tandem-unit node):
// these are copy-number loci even when a haplotype traverses them only once.
std::unordered_set<std::string> self_loop_nodes(const Graph& graph) {
    std::unordered_set<std::string> out;
    for (const auto& [id, node] : graph.nodes) {
        bool loop = false;
        for (const Neighbor& nb : node.start) if (nb.node_id == id) { loop = true; break; }
        if (!loop) for (const Neighbor& nb : node.end) if (nb.node_id == id) { loop = true; break; }
        if (loop) out.insert(id);
    }
    return out;
}

// One node-token position in a walk (oriented), for the DEL/INS/INV alignment.
struct Tok {
    std::uint64_t token = 0;
    std::string node_id;
    bool reverse = false;
};

// Token walk for the DEL/INS/INV alignment. Copy-number nodes (REP self-loops)
// are dropped entirely — they are handled separately as count-based DUP/CN events,
// so a haplotype's extra REP copies are never mistyped as INS.
std::vector<Tok> collapse_walk(
    const std::vector<PathStep>& steps,
    const std::unordered_set<std::string>& cn_nodes) {

    std::vector<Tok> out;
    for (const PathStep& s : steps) {
        if (cn_nodes.count(s.node_id) != 0) {
            continue;
        }
        Tok t;
        t.token = hash_step_token(s);
        t.node_id = s.node_id;
        t.reverse = s.reverse;
        out.push_back(std::move(t));
    }
    return out;
}

// A typed event derived from one haplotype walk vs the reference walk.
struct Event {
    EvType type = EvType::Del;
    std::vector<std::string> nodes;   // variant node set (length-weighted Jaccard key)
    std::string start_node;           // graph anchors (first/last involved node)
    std::string end_node;
    std::string seq;                  // event sequence (DELSEQ/INSSEQ/INVSEQ)
    std::size_t ref_cn = 0;           // DUP only
    std::size_t alt_cn = 0;           // DUP only
    std::string ins_subtype;          // INS only: "", "NOVEL", "DUP"
    std::string link_id;              // shared id for a co-located DEL+INS substitution (EVENTID)
    // Reference anchoring for VCF coordinates.
    std::string anchor_node;          // ref node POS is taken from
    bool anchor_after = false;        // true: POS = last base of anchor (INS); false: first base
    std::size_t size_bp = 0;          // |event| for min_sv filtering / SVLEN magnitude
    std::size_t ref_pos = 0;          // reference genomic position of the anchor (merge window)
};

// Spell a token run into sequence.
std::string spell_toks(const Graph& graph, const std::vector<const Tok*>& toks) {
    std::vector<PathStep> steps;
    for (const Tok* t : toks) steps.push_back(PathStep{t->node_id, t->reverse});
    return spell_path_steps_sequence(graph, steps);
}

std::size_t toks_bp(const Graph& graph, const std::vector<const Tok*>& toks) {
    std::size_t bp = 0;
    for (const Tok* t : toks) bp += node_len(graph, t->node_id);
    return bp;
}

// Is rc(hap block) == ref block (an inversion)? Compares node-id + orientation
// of the reference block against the reverse-complemented haplotype block.
bool is_inversion(const std::vector<const Tok*>& ref_blk, const std::vector<const Tok*>& hap_blk) {
    if (ref_blk.size() != hap_blk.size() || ref_blk.empty()) return false;
    for (std::size_t i = 0; i < ref_blk.size(); ++i) {
        const Tok* r = ref_blk[i];
        const Tok* h = hap_blk[hap_blk.size() - 1 - i]; // reverse order
        if (r->node_id != h->node_id || r->reverse == h->reverse) return false; // orientation must flip
    }
    return true;
}

// Global node-token alignment (diagonal only on equal tokens; else gaps), then
// read off DEL/INS/INV events of the haplotype relative to the reference. DUP/CN
// is handled separately (count-based), so CN nodes are already excluded from R/H.
std::vector<Event> diff_walks(
    const Graph& graph,
    const std::vector<Tok>& R,
    const std::vector<Tok>& H,
    const std::string& bubble_source) {

    std::vector<Event> events;
    const std::size_t m = R.size();
    const std::size_t n = H.size();
    if (static_cast<std::size_t>(m) * static_cast<std::size_t>(n) > kAlignCellCap) {
        return events; // pathologically large bubble; skip (intended input is normalized)
    }

    // DP: match (+2) only when tokens equal, gaps (-1). Mismatch diagonals forbidden
    // so substitutions surface as separate ref-only / hap-only blocks (for INV check).
    const int kMatch = 2;
    const int kGap = -1;
    const int kNeg = -1000000;
    std::vector<std::vector<int>> dp(m + 1, std::vector<int>(n + 1, 0));
    for (std::size_t i = 1; i <= m; ++i) dp[i][0] = static_cast<int>(i) * kGap;
    for (std::size_t j = 1; j <= n; ++j) dp[0][j] = static_cast<int>(j) * kGap;
    for (std::size_t i = 1; i <= m; ++i) {
        for (std::size_t j = 1; j <= n; ++j) {
            const int diag = (R[i - 1].token == H[j - 1].token) ? dp[i - 1][j - 1] + kMatch : kNeg;
            const int up = dp[i - 1][j] + kGap;
            const int left = dp[i][j - 1] + kGap;
            dp[i][j] = std::max(diag, std::max(up, left));
        }
    }

    // Traceback into columns (ri, hi); -1 means a gap on that side.
    struct Col { long long ri; long long hi; };
    std::vector<Col> cols;
    std::size_t i = m, j = n;
    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 && R[i - 1].token == H[j - 1].token && dp[i][j] == dp[i - 1][j - 1] + kMatch) {
            cols.push_back({static_cast<long long>(i - 1), static_cast<long long>(j - 1)});
            --i; --j;
        } else if (i > 0 && dp[i][j] == dp[i - 1][j] + kGap) {
            cols.push_back({static_cast<long long>(i - 1), -1});
            --i;
        } else {
            cols.push_back({-1, static_cast<long long>(j - 1)});
            --j;
        }
    }
    std::reverse(cols.begin(), cols.end());

    // Emit events from maximal gap blocks -> DEL / INS / INV. Track the last matched
    // ref node as the anchor for an INS that follows it.
    std::string last_ref_node = bubble_source;
    auto flush_block = [&](std::vector<const Tok*>& ref_blk, std::vector<const Tok*>& hap_blk) {
        if (ref_blk.empty() && hap_blk.empty()) return;
        if (!ref_blk.empty() && !hap_blk.empty() && is_inversion(ref_blk, hap_blk)) {
            Event e;
            e.type = EvType::Inv;
            for (const Tok* t : ref_blk) e.nodes.push_back(t->node_id);
            e.start_node = ref_blk.front()->node_id;
            e.end_node = ref_blk.back()->node_id;
            e.seq = spell_toks(graph, ref_blk);
            e.size_bp = toks_bp(graph, ref_blk);
            e.anchor_node = ref_blk.front()->node_id;
            events.push_back(std::move(e));
        } else {
            // A gap-block with both ref and hap content is a substitution: emit DEL+INS
            // but link them with a shared EVENTID so downstream can pair them.
            const bool substitution = !ref_blk.empty() && !hap_blk.empty();
            const std::string link =
                substitution ? ("sub_" + ref_blk.front()->node_id + "_" + hap_blk.front()->node_id) : "";
            if (!ref_blk.empty()) {
                Event e;
                e.type = EvType::Del;
                for (const Tok* t : ref_blk) e.nodes.push_back(t->node_id);
                e.start_node = ref_blk.front()->node_id;
                e.end_node = ref_blk.back()->node_id;
                e.seq = spell_toks(graph, ref_blk);
                e.size_bp = toks_bp(graph, ref_blk);
                e.anchor_node = ref_blk.front()->node_id;
                e.link_id = link;
                events.push_back(std::move(e));
            }
            if (!hap_blk.empty()) {
                Event e;
                e.type = EvType::Ins;
                for (const Tok* t : hap_blk) e.nodes.push_back(t->node_id);
                e.start_node = hap_blk.front()->node_id;
                e.end_node = hap_blk.back()->node_id;
                e.seq = spell_toks(graph, hap_blk);
                e.size_bp = toks_bp(graph, hap_blk);
                e.anchor_node = last_ref_node;
                e.anchor_after = true;
                e.link_id = link;
                events.push_back(std::move(e));
            }
        }
        ref_blk.clear();
        hap_blk.clear();
    };

    std::vector<const Tok*> ref_blk;
    std::vector<const Tok*> hap_blk;
    for (const Col& c : cols) {
        if (c.ri >= 0 && c.hi >= 0) {
            flush_block(ref_blk, hap_blk);
            last_ref_node = R[static_cast<std::size_t>(c.ri)].node_id;
        } else if (c.ri >= 0) {
            ref_blk.push_back(&R[static_cast<std::size_t>(c.ri)]);
        } else {
            hap_blk.push_back(&H[static_cast<std::size_t>(c.hi)]);
        }
    }
    flush_block(ref_blk, hap_blk);
    return events;
}

// Within a haplotype, coalesce consecutive same-type DEL/INS/INV events whose
// reference gap is <= merge_distance_bp into one event (union node sets). The gap
// is measured in reference bp via the bubble's node->genomic-start map.
void coalesce_events(
    const Graph& graph,
    std::vector<Event>& events,
    const std::unordered_map<std::string, std::size_t>& ref_node_pos,
    std::size_t merge_distance_bp) {

    if (events.size() < 2) return;
    auto pos = [&](const std::string& node) -> long long {
        const auto it = ref_node_pos.find(node);
        return it == ref_node_pos.end() ? -1 : static_cast<long long>(it->second);
    };
    // Reference bp interval [lo, hi) that an event occupies on the reference.
    auto ref_span = [&](const Event& e, long long& lo, long long& hi) {
        if (e.type == EvType::Ins) {
            const long long p = pos(e.anchor_node);
            const long long end = (p < 0) ? p : p + static_cast<long long>(node_len(graph, e.anchor_node));
            lo = end; hi = end;
            return;
        }
        lo = -1; hi = -1;
        for (const std::string& n : e.nodes) {
            const long long p = pos(n);
            if (p < 0) continue;
            const long long q = p + static_cast<long long>(node_len(graph, n));
            if (lo < 0 || p < lo) lo = p;
            if (hi < 0 || q > hi) hi = q;
        }
    };
    std::vector<Event> out;
    long long prev_lo = -1, prev_hi = -1;
    for (Event& e : events) {
        long long lo = -1, hi = -1;
        ref_span(e, lo, hi);
        if (e.type != EvType::Dup && !out.empty() && out.back().type == e.type &&
            prev_hi >= 0 && lo >= 0 && lo - prev_hi <= static_cast<long long>(merge_distance_bp) &&
            lo - prev_hi >= -static_cast<long long>(merge_distance_bp)) {
            Event& prev = out.back();
            for (const std::string& nd : e.nodes) prev.nodes.push_back(nd);
            prev.end_node = e.end_node;
            if (!e.seq.empty()) prev.seq += e.seq;
            prev.size_bp += e.size_bp;
            prev_hi = std::max(prev_hi, hi);
            if (lo < prev_lo) prev_lo = lo;
        } else {
            out.push_back(std::move(e));
            prev_lo = lo;
            prev_hi = hi;
        }
    }
    events = std::move(out);
}

double weighted_jaccard(
    const Graph& graph,
    const std::vector<std::string>& a,
    const std::vector<std::string>& b) {

    std::unordered_set<std::string> sa(a.begin(), a.end());
    std::unordered_set<std::string> sb(b.begin(), b.end());
    std::size_t inter = 0;
    std::size_t uni = 0;
    std::unordered_set<std::string> seen;
    for (const std::string& n : sa) {
        seen.insert(n);
        const std::size_t l = node_len(graph, n);
        uni += l;
        if (sb.count(n)) inter += l;
    }
    for (const std::string& n : sb) {
        if (seen.count(n)) continue;
        uni += node_len(graph, n);
    }
    return uni == 0 ? 0.0 : static_cast<double>(inter) / static_cast<double>(uni);
}

// Banded alignment identity between two event sequences, gated on a length ratio
// (so wildly different sizes are not compared). 0 when either is empty.
double seq_identity(const std::string& a, const std::string& b, double min_id) {
    if (a.empty() || b.empty()) return 0.0;
    const std::size_t lo = std::min(a.size(), b.size());
    const std::size_t hi = std::max(a.size(), b.size());
    if (static_cast<double>(lo) < min_id * static_cast<double>(hi)) return 0.0;
    const std::string& shorter = a.size() <= b.size() ? a : b;
    const std::string& longer = a.size() <= b.size() ? b : a;
    const std::size_t band = std::max<std::size_t>(
        8, static_cast<std::size_t>((1.0 - min_id) * static_cast<double>(shorter.size())) + 8);
    const FitAlignResult fa = fit_align(shorter, longer, band);
    return fa.ok ? fa.identity : 0.0;
}

// A merged event across haplotypes: a representative event + its carriers.
struct MergedRecord {
    Event seed;                                             // representative (largest member)
    std::vector<std::string> carriers;                      // haplotype path names with GT=1
    std::unordered_set<std::size_t> member_alleles;         // allele indices already merged in
    std::unordered_map<std::string, std::size_t> sample_cn; // DUP per-sample copy number
};

std::string upper_base(char c) {
    const char u = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    if (u == 'A' || u == 'C' || u == 'G' || u == 'T' || u == 'N') return std::string(1, u);
    return "N";
}

} // namespace

void call_variants(
    const Graph& graph,
    const VariantCallOptions& options,
    VariantCallSummary* summary_out) {

    if (options.reference_path.empty()) {
        throw std::runtime_error("call requires --reference-path");
    }
    if (options.out_prefix.empty()) {
        throw std::runtime_error("call requires -o/--out-prefix");
    }

    // Locate the reference path record.
    const PathRecord* ref_path = nullptr;
    for (const PathRecord& p : graph.paths) {
        if (p.name == options.reference_path) { ref_path = &p; break; }
    }
    if (ref_path == nullptr) {
        throw std::runtime_error("Reference path not found in GFA: " + options.reference_path);
    }

    const std::vector<Bubble> bubbles = read_bubbles_csv(options.bubbles_csv_in);
    std::unordered_set<std::size_t> bubble_filter(options.bubble_ids.begin(), options.bubble_ids.end());

    const ParsedReferencePath ref_meta = parse_reference_path_label(options.reference_path);
    const std::vector<std::size_t> ref_prefix = path_prefix_bp(*ref_path, graph.nodes);

    // Fixed sample list: every haplotype path (the reference included).
    std::vector<std::string> sample_names;
    sample_names.reserve(graph.paths.size());
    for (const PathRecord& p : graph.paths) sample_names.push_back(p.name);

    // Per-path bubble index, reused across bubbles.
    std::vector<BubblePathIndex> path_indexes;
    path_indexes.reserve(graph.paths.size());
    for (const PathRecord& p : graph.paths) path_indexes.push_back(build_bubble_path_index(p));

    const std::unordered_set<std::string> selfloops = self_loop_nodes(graph);

    cli::ensure_parent_dir_for_file(options.out_prefix + ".region.vcf");

    auto write_vcf_header = [&](std::ostream& out) {
        out << "##fileformat=VCFv4.2\n";
        out << "##source=panvar call\n";
        out << "##reference=" << options.reference_path << "\n";
        out << "##contig=<ID=" << ref_meta.chrom << ">\n";
        out << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant\">\n";
        out << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Structural variant type\">\n";
        out << "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length difference ALT-REF\">\n";
        out << "##INFO=<ID=BUBBLE_ID,Number=1,Type=Integer,Description=\"panvar bubble identifier\">\n";
        out << "##INFO=<ID=START_NODE,Number=1,Type=String,Description=\"First graph node of the event\">\n";
        out << "##INFO=<ID=END_NODE,Number=1,Type=String,Description=\"Last graph node of the event\">\n";
        out << "##INFO=<ID=EVENT_NODES,Number=.,Type=String,Description=\"Variant node set\">\n";
        out << "##INFO=<ID=INS_SUBTYPE,Number=1,Type=String,Description=\"INS subtype: NOVEL or DUP (minimap2 refined)\">\n";
        out << "##INFO=<ID=REF_CN,Number=1,Type=Integer,Description=\"Reference copy number of the repeat unit (DUP)\">\n";
        out << "##INFO=<ID=NMERGED,Number=1,Type=Integer,Description=\"Haplotype carriers merged into this record\">\n";
        out << "##INFO=<ID=EVENTID,Number=1,Type=String,Description=\"Shared id linking a co-located DEL+INS substitution\">\n";
        out << "##INFO=<ID=INSSEQ,Number=1,Type=String,Description=\"Inserted sequence\">\n";
        out << "##INFO=<ID=DELSEQ,Number=1,Type=String,Description=\"Deleted reference sequence\">\n";
        out << "##INFO=<ID=INVSEQ,Number=1,Type=String,Description=\"Inverted reference sequence\">\n";
        out << "##ALT=<ID=DEL,Description=\"Deletion\">\n";
        out << "##ALT=<ID=INS,Description=\"Insertion\">\n";
        out << "##ALT=<ID=INV,Description=\"Inversion\">\n";
        out << "##ALT=<ID=DUP,Description=\"Copy-number / tandem duplication locus\">\n";
        out << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"1=carrier, 0=reference-like, .=bubble not traversed\">\n";
        out << "##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number of the repeat unit (DUP records)\">\n";
        out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
        for (const std::string& s : sample_names) out << '\t' << s;
        out << '\n';
    };

    std::ofstream region_out(options.out_prefix + ".region.vcf");
    if (!region_out) {
        throw std::runtime_error("Failed to write region VCF: " + options.out_prefix + ".region.vcf");
    }
    write_vcf_header(region_out);

    VariantCallSummary summary;

    // Genomic 1-based start of a bubble node from the reference path (first occurrence).
    auto build_ref_node_pos = [&](const Bubble& bubble) {
        std::unordered_map<std::string, std::size_t> pos; // node -> 1-based genomic start
        std::unordered_set<std::string> bnodes(bubble.inside.begin(), bubble.inside.end());
        bnodes.insert(bubble.source);
        bnodes.insert(bubble.sink);
        for (std::size_t k = 0; k < ref_path->steps.size(); ++k) {
            const std::string& id = ref_path->steps[k].node_id;
            if (bnodes.count(id) && !pos.count(id)) {
                pos[id] = ref_meta.region_start_1based + ref_prefix[k];
            }
        }
        return pos;
    };

    for (const Bubble& bubble : bubbles) {
        if (!bubble_filter.empty() && bubble_filter.find(bubble.id) == bubble_filter.end()) continue;
        ++summary.bubbles_seen;

        // Reference walk through this bubble. path_indexes is parallel to graph.paths.
        std::size_t ref_idx = graph.paths.size();
        for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
            if (graph.paths[pi].name == options.reference_path) { ref_idx = pi; break; }
        }
        const auto ref_iv = (ref_idx < graph.paths.size())
            ? find_best_bubble_path_interval(path_indexes[ref_idx], bubble)
            : std::optional<BubblePathInterval>{};
        if (!ref_iv.has_value()) {
            continue; // reference does not traverse this bubble; cannot type events
        }
        ++summary.bubbles_with_reference;

        const std::vector<PathStep> ref_steps =
            canonical_bubble_path_steps(graph.paths[ref_idx], bubble, *ref_iv);
        if (ref_steps.empty()) continue;

        // Distinct alleles: group paths by canonical-walk signature.
        struct Allele { std::vector<PathStep> steps; std::vector<std::string> members; };
        std::unordered_map<std::string, std::size_t> sig_to_allele;
        std::vector<Allele> alleles;
        std::unordered_set<std::string> traverses; // path names that cross the bubble
        const std::string ref_sig = build_walk_signature(ref_steps);
        for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
            const auto iv = find_best_bubble_path_interval(path_indexes[pi], bubble);
            if (!iv.has_value()) continue;
            const std::vector<PathStep> steps = canonical_bubble_path_steps(graph.paths[pi], bubble, *iv);
            if (steps.empty()) continue;
            traverses.insert(graph.paths[pi].name);
            const std::string sig = build_walk_signature(steps);
            auto it = sig_to_allele.find(sig);
            if (it == sig_to_allele.end()) {
                sig_to_allele.emplace(sig, alleles.size());
                Allele a; a.steps = steps; a.members.push_back(graph.paths[pi].name);
                alleles.push_back(std::move(a));
            } else {
                alleles[it->second].members.push_back(graph.paths[pi].name);
            }
        }

        // CN nodes: self-loop nodes present in this bubble (a panphorte REP / genuine
        // tandem unit). They are handled as count-based DUP/CN events and excluded from
        // the DEL/INS/INV alignment. Ordinary nodes that merely recur (e.g. shared module
        // content) are NOT copy-number loci and stay in the alignment.
        std::unordered_map<std::string, std::size_t> ref_count;
        for (const PathStep& s : ref_steps) ++ref_count[s.node_id];
        std::unordered_set<std::string> cn_nodes;
        for (const auto& [id, c] : ref_count) { (void)c; if (selfloops.count(id)) cn_nodes.insert(id); }
        for (const Allele& a : alleles) {
            for (const PathStep& s : a.steps) if (selfloops.count(s.node_id)) cn_nodes.insert(s.node_id);
        }

        const std::vector<Tok> Rtok = collapse_walk(ref_steps, cn_nodes);

        std::vector<MergedRecord> merged;
        const std::unordered_map<std::string, std::size_t> ref_node_pos = build_ref_node_pos(bubble);
        const std::size_t rescue_floor =
            options.rescue_min_bp != 0 ? options.rescue_min_bp : std::max<std::size_t>(1, options.min_sv_bp / 2);

        auto ev_ref_pos = [&](const Event& e) -> long long {
            const auto it = ref_node_pos.find(e.anchor_node);
            if (it == ref_node_pos.end()) return -1;
            const std::size_t glen = node_len(graph, e.anchor_node);
            return static_cast<long long>(it->second + (e.anchor_after && glen > 0 ? glen - 1 : 0));
        };
        // Two non-DUP events are the same site: same type, anchors within the window,
        // and either node sets overlap (Jaccard) OR sequences are similar.
        auto events_match = [&](const Event& a, const Event& b) {
            if (a.type != b.type) return false;
            if (a.ref_pos == 0 || b.ref_pos == 0) {
                // fall back to node/seq only when no coordinate
            } else {
                const long long d = static_cast<long long>(a.ref_pos) - static_cast<long long>(b.ref_pos);
                if (d > static_cast<long long>(options.merge_distance_bp) ||
                    d < -static_cast<long long>(options.merge_distance_bp)) return false;
            }
            if (weighted_jaccard(graph, a.nodes, b.nodes) >= options.merge_jaccard) return true;
            return seq_identity(a.seq, b.seq, options.merge_seq_identity) >= options.merge_seq_identity;
        };

        // ---- DUP/CN events: count self-loop traversals per allele vs reference. Merge
        // on shared REP node; per-sample CN. Independent of the walk-diff alignment.
        for (std::size_t ai = 0; ai < alleles.size(); ++ai) {
            std::unordered_map<std::string, std::size_t> alt_count;
            for (const PathStep& s : alleles[ai].steps) ++alt_count[s.node_id];
            for (const std::string& cn : cn_nodes) {
                const std::size_t rc = ref_count.count(cn) ? ref_count.at(cn) : 0;
                const std::size_t ac = alt_count.count(cn) ? alt_count.at(cn) : 0;
                if (rc == ac) continue;
                const std::size_t unit = node_len(graph, cn);
                const std::size_t delta = rc > ac ? rc - ac : ac - rc;
                if (unit * delta < options.min_sv_bp) continue;
                MergedRecord* grp = nullptr;
                for (MergedRecord& mr : merged) {
                    if (mr.seed.type == EvType::Dup && mr.seed.nodes.front() == cn) { grp = &mr; break; }
                }
                if (grp == nullptr) {
                    MergedRecord mr;
                    mr.seed.type = EvType::Dup;
                    mr.seed.nodes.push_back(cn);
                    mr.seed.start_node = cn; mr.seed.end_node = cn;
                    mr.seed.ref_cn = rc; mr.seed.alt_cn = ac;
                    mr.seed.anchor_node = bubble.source;
                    mr.seed.size_bp = unit * delta;
                    merged.push_back(std::move(mr));
                    grp = &merged.back();
                }
                for (const std::string& m : alleles[ai].members) {
                    grp->carriers.push_back(m);
                    grp->sample_cn[m] = ac;
                }
            }
        }

        // ---- DEL/INS/INV events: derive per allele, keep ALL (down to the rescue floor)
        // for the re-scan, then merge events >= floor by position + sequence/node match.
        std::vector<std::vector<Event>> allele_events(alleles.size());
        for (std::size_t ai = 0; ai < alleles.size(); ++ai) {
            if (build_walk_signature(alleles[ai].steps) == ref_sig) continue;
            const std::vector<Tok> Htok = collapse_walk(alleles[ai].steps, cn_nodes);
            std::vector<Event> events = diff_walks(graph, Rtok, Htok, bubble.source);
            coalesce_events(graph, events, ref_node_pos, options.merge_distance_bp);
            for (Event& e : events) {
                const long long p = ev_ref_pos(e);
                e.ref_pos = p < 0 ? 0 : static_cast<std::size_t>(p);
            }
            allele_events[ai] = std::move(events);
        }

        for (std::size_t ai = 0; ai < alleles.size(); ++ai) {
            for (const Event& e : allele_events[ai]) {
                if (e.size_bp < rescue_floor) continue;
                MergedRecord* grp = nullptr;
                for (MergedRecord& mr : merged) {
                    if (mr.seed.type == EvType::Dup) continue;
                    if (events_match(mr.seed, e)) { grp = &mr; break; }
                }
                if (grp == nullptr) {
                    MergedRecord mr;
                    mr.seed = e;
                    mr.member_alleles.insert(ai);
                    for (const std::string& m : alleles[ai].members) mr.carriers.push_back(m);
                    merged.push_back(std::move(mr));
                } else {
                    if (e.size_bp > grp->seed.size_bp) {
                        const std::string keep_link = grp->seed.link_id;
                        grp->seed = e; // largest member represents the record
                        if (grp->seed.link_id.empty()) grp->seed.link_id = keep_link;
                    }
                    if (grp->member_alleles.insert(ai).second) {
                        for (const std::string& m : alleles[ai].members) grp->carriers.push_back(m);
                    }
                }
            }
        }

        // ---- Joint re-scan: rescue any haplotype whose walk carries a comparable
        // signature for a called record, even if its own event was sub-threshold.
        for (MergedRecord& mr : merged) {
            if (mr.seed.type == EvType::Dup) continue;
            for (std::size_t ai = 0; ai < alleles.size(); ++ai) {
                if (mr.member_alleles.count(ai)) continue;
                for (const Event& e : allele_events[ai]) {
                    if (events_match(mr.seed, e)) {
                        mr.member_alleles.insert(ai);
                        for (const std::string& m : alleles[ai].members) mr.carriers.push_back(m);
                        break;
                    }
                }
            }
        }

        // ---- INS subtype refinement on the representative only (bounded minimap2 calls).
        if (options.classify_ins) {
            std::string ref_window;
            for (MergedRecord& mr : merged) {
                if (mr.seed.type != EvType::Ins || mr.seed.seq.empty()) continue;
                if (ref_window.empty()) ref_window = spell_path_steps_sequence(graph, ref_steps);
                const Minimap2Hit hit = minimap2_best_hit(
                    "ins", mr.seed.seq, "refwin", ref_window, options.minimap_preset, options.minimap_best_n);
                mr.seed.ins_subtype = (hit.ok && hit.identity() >= options.ins_dup_min_identity &&
                                       hit.query_end_bp - hit.query_start_bp >= mr.seed.seq.size() / 2)
                                          ? "DUP" : "NOVEL";
            }
        }

        // ---- Keep records reaching min_sv_bp (by representative) and --min-haplotypes.
        {
            std::vector<MergedRecord> kept;
            for (MergedRecord& mr : merged) {
                std::unordered_set<std::string> uniq(mr.carriers.begin(), mr.carriers.end());
                std::size_t support = uniq.size();
                if (mr.seed.type == EvType::Dup) {
                    support = 0;
                    for (const auto& [s, cn] : mr.sample_cn) { (void)s; if (cn != mr.seed.ref_cn) ++support; }
                } else if (mr.seed.size_bp < options.min_sv_bp) {
                    continue;
                }
                if (support < options.min_haplotypes) continue;
                kept.push_back(std::move(mr));
            }
            merged = std::move(kept);
        }
        if (merged.empty()) continue;
        ++summary.bubbles_with_calls;

        // Emit VCF rows for this bubble (and append to the region VCF).
        std::ofstream bubble_out;
        if (options.write_per_bubble_vcf) {
            const std::string path = options.out_prefix + ".bubble_" + std::to_string(bubble.id) + ".vcf";
            bubble_out.open(path);
            if (!bubble_out) throw std::runtime_error("Failed to write bubble VCF: " + path);
            write_vcf_header(bubble_out);
        }

        for (const MergedRecord& mr : merged) {
            const Event& e = mr.seed;
            // POS / REF base from the reference node map.
            std::size_t pos = 1;
            std::string ref_base = "N";
            const auto ait = ref_node_pos.find(e.anchor_node);
            if (ait != ref_node_pos.end()) {
                const std::size_t glen = node_len(graph, e.anchor_node);
                if (e.anchor_after && glen > 0) {
                    pos = ait->second + glen - 1; // last base of preceding ref node
                } else {
                    pos = ait->second;            // first base of first event node
                }
                const auto nit = graph.nodes.find(e.anchor_node);
                if (nit != graph.nodes.end() && !nit->second.sequence.empty()) {
                    const std::size_t bi = (e.anchor_after && glen > 0) ? glen - 1 : 0;
                    ref_base = upper_base(nit->second.sequence[bi]);
                }
            }
            std::size_t end = pos;
            long long svlen = 0;
            const char* svt = ev_svtype(e.type);
            if (e.type == EvType::Del) { svlen = -static_cast<long long>(e.size_bp); end = pos + e.size_bp; }
            else if (e.type == EvType::Ins) { svlen = static_cast<long long>(e.size_bp); end = pos; }
            else if (e.type == EvType::Inv) { svlen = static_cast<long long>(e.size_bp); end = pos + e.size_bp; }
            else { // DUP
                svlen = (static_cast<long long>(e.alt_cn) - static_cast<long long>(e.ref_cn)) *
                        static_cast<long long>(node_len(graph, e.nodes.empty() ? std::string() : e.nodes.front()));
                end = pos + node_len(graph, e.nodes.empty() ? std::string() : e.nodes.front());
            }

            std::unordered_set<std::string> carrier_set(mr.carriers.begin(), mr.carriers.end());

            std::ostringstream info;
            info << "END=" << end << ";SVTYPE=" << svt << ";SVLEN=" << svlen
                 << ";BUBBLE_ID=" << bubble.id
                 << ";START_NODE=" << e.start_node << ";END_NODE=" << e.end_node
                 << ";NMERGED=" << carrier_set.size();
            info << ";EVENT_NODES=";
            for (std::size_t k = 0; k < e.nodes.size(); ++k) { if (k) info << ','; info << e.nodes[k]; }
            if (!e.link_id.empty()) info << ";EVENTID=bubble" << bubble.id << "_" << e.link_id;
            if (e.type == EvType::Dup) info << ";REF_CN=" << e.ref_cn;
            if (e.type == EvType::Ins && !e.ins_subtype.empty()) info << ";INS_SUBTYPE=" << e.ins_subtype;
            if (!e.seq.empty() && e.seq.size() <= 20000) {
                if (e.type == EvType::Ins) info << ";INSSEQ=" << e.seq;
                else if (e.type == EvType::Del) info << ";DELSEQ=" << e.seq;
                else if (e.type == EvType::Inv) info << ";INVSEQ=" << e.seq;
            }

            const std::string id = "bubble" + std::to_string(bubble.id) + "_" + svt + "_" + e.start_node;

            std::ostringstream row;
            row << ref_meta.chrom << '\t' << pos << '\t' << id << '\t' << ref_base
                << "\t<" << svt << ">\t.\t.\t" << info.str() << "\tGT:CN";
            for (const std::string& s : sample_names) {
                row << '\t';
                if (!traverses.count(s)) { row << ".:."; continue; }
                const bool carrier = carrier_set.count(s) != 0;
                row << (carrier ? "1" : "0") << ':';
                if (e.type == EvType::Dup) {
                    const auto cit = mr.sample_cn.find(s);
                    row << (cit != mr.sample_cn.end() ? std::to_string(cit->second) : std::to_string(e.ref_cn));
                } else {
                    row << '.';
                }
            }
            row << '\n';
            region_out << row.str();
            if (options.write_per_bubble_vcf) bubble_out << row.str();

            ++summary.records_written;
            if (e.type == EvType::Del) ++summary.del;
            else if (e.type == EvType::Ins) ++summary.ins;
            else if (e.type == EvType::Inv) ++summary.inv;
            else ++summary.dup;
        }
    }

    if (summary_out) *summary_out = summary;
    if (!options.quiet) {
        std::cerr << "[call] bubbles=" << summary.bubbles_seen
                  << " with_ref=" << summary.bubbles_with_reference
                  << " with_calls=" << summary.bubbles_with_calls
                  << " records=" << summary.records_written
                  << " (DEL=" << summary.del << " INS=" << summary.ins
                  << " INV=" << summary.inv << " DUP=" << summary.dup << ")\n";
    }
}

} // namespace panvar
