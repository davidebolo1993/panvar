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
#include <numeric>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <system_error>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace panvar {
namespace {

void print_benchmark_help() {
    std::cout
        << "Usage:\n"
        << "  panvar benchmark -i <graph.gfa> (-b <prefix> | -c <bubbles.csv>) -r <name> --calls <prefix>\n"
        << "                   -o <prefix> [options]\n\n"
        << "Check how much of every haplotype's large variation call's output accounts for. Each haplotype\n"
        << "is rebuilt from the reference twice, with call's own events for it and with the VCF records it\n"
        << "carries, and each rebuild is aligned to the haplotype's own sequence: what each leaves is listed\n"
        << "with why, a VCF leftover as one of call's misses again or as lost in writing the VCF.\n\n"
        << "Options:\n"
        << "  -i, --gfa <path>                  Graph the calls were made on (required)\n"
        << "  -b, --bubble-prefix <prefix>      Sites from <prefix>.bubbles.csv (required unless -c)\n"
        << "  -c, --bubbles-csv <path>          Site table, by path\n"
        << "  -r, --reference-path <name>       Reference path, as given to call (required)\n"
        << "      --calls <prefix>              call's output: <prefix>.region.vcf and\n"
        << "                                    <prefix>.carrier_events.tsv, without --diploid\n"
        << "                                    (required)\n"
        << "  -o, --out-prefix <prefix>         Output prefix (required)\n"
        << "      --min-sv-bp <N>               Smallest large event (default 50); below call's own,\n"
        << "                                    what lies under call's threshold is listed too\n"
        << "      --merge-gap <N>               Differences N or fewer aligned bases apart are one event\n"
        << "                                    (default 10)\n"
        << "  -t, --threads <N>                 Worker threads (default 0 = all cores)\n"
        << "  -q, --quiet                       No progress logs\n"
        << "  -h, --help                        Show this help\n";
}

std::vector<std::string> split_on(const std::string& s, char sep) {
    std::vector<std::string> out;
    std::size_t start = 0;
    for (std::size_t i = 0; i <= s.size(); ++i)
        if (i == s.size() || s[i] == sep) { out.push_back(s.substr(start, i - start)); start = i + 1; }
    return out;
}

std::size_t parse_count(const std::string& t, const std::string& what) {
    char* end = nullptr;
    const unsigned long long v = std::strtoull(t.c_str(), &end, 10);
    if (t.empty() || end == t.c_str() || *end != '\0' || t[0] == '-')
        throw std::runtime_error(what + ": '" + t + "' is not a count");
    return static_cast<std::size_t>(v);
}

// One row of call's carrier_events.tsv: the event the listed carriers were called for, as the reference
// interval [ref_start, ref_end] (1-based; ref_end = ref_start - 1 for an insertion before ref_start)
// and the forward-strand bases that replace it.
struct CarrierEvent {
    std::string variant_id;
    std::size_t bubble_id = 0;
    std::string kind;
    std::size_t ref_start = 0, ref_end = 0;
    std::string cn;
    std::string walk_cn;
    long long event_bp = 0;        // the event as call measured it: net change, an inversion's length
    std::string alt;
    std::vector<std::string> haplotypes;
    std::size_t order = 0;         // --vcf: the record's place in the file
    bool insertion() const { return ref_end + 1 == ref_start; }
    std::size_t ref_len() const { return insertion() ? 0 : ref_end - ref_start + 1; }
};

std::vector<CarrierEvent> load_carrier_events(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path + " (written by call; rerun call to produce it)");
    std::string line;
    if (!std::getline(in, line) ||
        line != "variant_id\tbubble_id\tkind\tref_start\tref_end\tcn\twalk_cn\tevent_bp\talt_bp\talt\thaplotypes")
        throw std::runtime_error(path + ": not a carrier_events.tsv header");
    std::vector<CarrierEvent> out;
    std::size_t lineno = 1;
    while (std::getline(in, line)) {
        ++lineno;
        if (line.empty()) continue;
        const std::vector<std::string> f = split_on(line, '\t');
        const std::string where = path + " line " + std::to_string(lineno);
        if (f.size() != 11) throw std::runtime_error(where + ": expected 11 columns");
        CarrierEvent e;
        e.variant_id = f[0];
        e.bubble_id = parse_count(f[1], where);
        e.kind = f[2];
        e.ref_start = parse_count(f[3], where);
        e.ref_end = parse_count(f[4], where);
        e.cn = f[5];
        e.walk_cn = f[6];
        if (f[7] != ".") {
            char* end = nullptr;
            e.event_bp = std::strtoll(f[7].c_str(), &end, 10);
            if (end == f[7].c_str() || *end != '\0') throw std::runtime_error(where + ": event_bp is not a number");
        }
        e.alt = f[9] == "." ? std::string() : f[9];
        if (parse_count(f[8], where) != e.alt.size()) throw std::runtime_error(where + ": alt_bp does not match alt");
        if (e.kind != "unresolved" && (e.ref_start == 0 || e.ref_end + 1 < e.ref_start))
            throw std::runtime_error(where + ": not a reference interval");
        e.haplotypes = split_on(f[10], ',');
        out.push_back(std::move(e));
    }
    return out;
}

// What the region VCF says about each record: its size and, per carrier, its copy number.
struct VcfRecordInfo {
    std::string id;
    std::size_t pos = 0, end = 0;                         // POS and INFO/END (END = POS without one)
    std::string svtype;
    bool has_svlen = false;
    long long svlen = 0;
    long long ref_cn = -1;
    long long unit_bp = 0;                                // RU_LEN or CN_UNIT_BP: one copy
    long long cn_ref_bp = 0;                              // a CNV's span, POS+1..END: the reference's copies
    std::string ref, info;
    std::vector<std::string> alts;
    std::unordered_set<std::string> carriers;             // GT names an ALT
    std::unordered_map<std::string, int> sample_allele;   // GT's allele index; -1 for '.'
    std::unordered_map<std::string, long long> sample_cn;
    std::unordered_map<std::string, long long> sample_cnbp;   // FORMAT:CNBP, the carrier's own change in bp
    // MODULE: the copies carriers gained or lost (INFO/CN_FORMS), reference [start, end] replaced by seq
    // (end = start - 1 for an insertion; spelled false where call left the sequence out), and each carrier's
    // (FORMAT/CNFORM, 0-based)
    struct Form { std::size_t start = 0, end = 0; std::string seq; bool spelled = true; };
    std::vector<Form> cn_forms;
    std::unordered_map<std::string, std::vector<std::size_t>> sample_cnform;
};

// Whether a record's change lies inside genomic [lo, hi]: an insertion goes in before POS+1, so it is inside
// only when POS+1 is past lo and no further than hi; any other change runs from POS+1 to END.
bool change_inside(const VcfRecordInfo& o, std::size_t lo, std::size_t hi) {
    const std::size_t a = o.pos + 1, b = std::max(o.end, o.pos);
    if (b < a) return a > lo && a <= hi;
    return a <= hi && b >= lo;
}

// A record's length change as far as it falls inside [lo, hi]: a deletion reaching past an edge counts
// only its bases inside; anything else counts whole.
long long change_within(const VcfRecordInfo& o, long long change, std::size_t lo, std::size_t hi) {
    const std::size_t a = o.pos + 1, b = std::max(o.end, o.pos);
    if (change >= 0 || b < a || (a >= lo && b <= hi) || b < lo || a > hi) return change;
    const std::size_t in = std::min(b, hi) + 1 - std::max(a, lo);
    return std::max(change, -static_cast<long long>(in));
}

std::string info_field(const std::string& info, const std::string& key) {
    for (const std::string& kv : split_on(info, ';')) {
        const std::size_t eq = kv.find('=');
        if (eq != std::string::npos && kv.compare(0, eq, key) == 0) return kv.substr(eq + 1);
    }
    return std::string();
}

bool has_flag(const std::string& info, const std::string& key) {
    for (const std::string& kv : split_on(info, ';')) if (kv == key) return true;
    return false;
}

std::vector<VcfRecordInfo> load_region_vcf(const std::string& path, std::vector<std::string>& samples) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    std::vector<VcfRecordInfo> out;
    std::string line;
    bool header = false;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        // A --diploid VCF has a column per sample, not per path, even where its GTs show one haplotype.
        if (line.rfind("##panvar_columns=diploid", 0) == 0)
            throw std::runtime_error(path + " was written by call --diploid (one column per sample); benchmark scores "
                                     "haplotypes (one column per path): rerun call without --diploid");
        if (line.empty() || line.rfind("##", 0) == 0) continue;
        const std::vector<std::string> f = split_on(line, '\t');
        if (line.rfind("#CHROM", 0) == 0) { samples.assign(f.begin() + std::min<std::size_t>(9, f.size()), f.end()); header = true; continue; }
        if (!header) throw std::runtime_error(path + ": a record before the #CHROM line; not a VCF call wrote");
        if (f.size() < 9 + samples.size()) throw std::runtime_error(path + ": a record has fewer columns than samples");
        VcfRecordInfo r;
        r.id = f[2];
        r.pos = std::stoull(f[1]);
        const std::string en = info_field(f[7], "END");
        r.end = en.empty() ? r.pos + (f[3].size() > 0 ? f[3].size() - 1 : 0) : std::stoull(en);
        r.svtype = info_field(f[7], "SVTYPE");
        r.ref = f[3];
        r.alts = split_on(f[4], ',');
        r.info = f[7];
        const std::string sl = info_field(f[7], "SVLEN");
        if (!sl.empty() && sl.find(',') == std::string::npos) { r.has_svlen = true; r.svlen = std::stoll(sl); }
        const std::string rc = info_field(f[7], "REF_CN");
        if (!rc.empty()) r.ref_cn = std::stoll(rc);
        std::string unit = info_field(f[7], "RU_LEN");
        if (unit.empty()) unit = info_field(f[7], "CN_UNIT_BP");
        if (!unit.empty()) r.unit_bp = std::stoll(unit);
        if (r.svtype == "CNV" && r.end > r.pos) r.cn_ref_bp = static_cast<long long>(r.end - r.pos);
        const std::string forms = info_field(f[7], "CN_FORMS");
        if (!forms.empty())
            for (const std::string& x : split_on(forms, ',')) {
                const std::size_t colon = x.find(':'), dash = x.find('-');
                if (colon == std::string::npos || dash == std::string::npos || dash > colon)
                    throw std::runtime_error(path + ": record " + r.id + " has a CN_FORMS entry that is not START-END:SEQ ('" + x + "')");
                VcfRecordInfo::Form fm;
                fm.start = parse_count(x.substr(0, dash), path + " CN_FORMS");
                fm.end = parse_count(x.substr(dash + 1, colon - dash - 1), path + " CN_FORMS");
                fm.seq = x.substr(colon + 1);
                if (fm.seq == "*") { fm.seq.clear(); fm.spelled = false; }
                if (fm.end + 1 < fm.start) throw std::runtime_error(path + ": record " + r.id + " has a CN_FORMS entry ending before it starts");
                r.cn_forms.push_back(std::move(fm));
            }
        const std::vector<std::string> fmt = split_on(f[8], ':');
        const auto gt_i = std::find(fmt.begin(), fmt.end(), "GT") - fmt.begin();
        const auto cn_i = std::find(fmt.begin(), fmt.end(), "CN") - fmt.begin();
        const auto cnbp_i = std::find(fmt.begin(), fmt.end(), "CNBP") - fmt.begin();
        const auto cnform_i = std::find(fmt.begin(), fmt.end(), "CNFORM") - fmt.begin();
        for (std::size_t s = 0; s < samples.size(); ++s) {
            const std::vector<std::string> v = split_on(f[9 + s], ':');
            if (static_cast<std::size_t>(gt_i) < v.size()) {
                const std::string& gt = v[static_cast<std::size_t>(gt_i)];
                if (gt.find_first_of("/|") != std::string::npos)
                    throw std::runtime_error(path + ": record " + r.id + " has a diploid GT ('" + gt +
                                             "'); benchmark scores haplotypes (one column per path)");
                if (gt != "." && gt != "0") r.carriers.insert(samples[s]);
                r.sample_allele[samples[s]] = gt == "." ? -1 : static_cast<int>(parse_count(gt, path + " GT"));
            }
            if (static_cast<std::size_t>(cn_i) < v.size() && v[static_cast<std::size_t>(cn_i)] != ".")
                r.sample_cn[samples[s]] = std::stoll(v[static_cast<std::size_t>(cn_i)]);
            if (static_cast<std::size_t>(cnbp_i) < v.size() && v[static_cast<std::size_t>(cnbp_i)] != ".")
                r.sample_cnbp[samples[s]] = std::stoll(v[static_cast<std::size_t>(cnbp_i)]);
            if (static_cast<std::size_t>(cnform_i) < v.size() && v[static_cast<std::size_t>(cnform_i)] != ".")
                for (const std::string& k : split_on(v[static_cast<std::size_t>(cnform_i)], ',')) {
                    const std::size_t i = parse_count(k, path + " CNFORM");
                    if (i == 0 || i > r.cn_forms.size())
                        throw std::runtime_error(path + ": record " + r.id + " names CN_FORMS entry " + k + ", which it does not have");
                    r.sample_cnform[samples[s]].push_back(i - 1);
                }
        }
        out.push_back(std::move(r));
    }
    return out;
}

// A difference between two sequences, from their alignment. Gaps (I/D columns) no more than `merge_gap`
// aligned columns apart are one event, with the mismatches within `merge_gap` of either end; a run of
// mismatches with no gap is an event where it is dense, one column in five or more over `min_event`
// columns or more (an inversion or an equal-length replacement can align with few gaps). An event is
// sized as call sizes a replaced run, by the larger of the bases each side holds over it: its gaps and
// mismatches, and the chance matches between them.
struct Diff {
    std::size_t t0 = 0, t1 = 0;     // target interval [t0, t1)
    std::size_t q0 = 0, q1 = 0;     // query interval [q0, q1)
    std::size_t del_bp = 0, ins_bp = 0;
    bool has_del = false, has_ins = false;
    std::size_t size() const { return std::max(del_bp, ins_bp); }
    const char* type() const { return !has_ins && has_del ? "DEL" : !has_del && has_ins ? "INS" : "SUB"; }
};

// Matches inside an event in runs shorter than this are chance matches; a longer run is shared sequence,
// as a repeat lining up between scattered changes.
constexpr std::size_t kChanceRun = 5;
// A dense window is searched over at most this many mismatches from its start.
constexpr std::size_t kMaxWindowMismatches = 2000;

struct DiffSet {
    std::vector<Diff> events;       // every event, any size
    std::size_t mismatches = 0;     // mismatches outside every event (SNVs)
};

DiffSet diff_of(const std::string& ops, std::size_t merge_gap, std::size_t min_event) {
    DiffSet out;
    const std::size_t n = ops.size();
    const auto gap = [&](std::size_t c) { return ops[c] == 'I' || ops[c] == 'D'; };
    std::vector<std::pair<std::size_t, std::size_t>> spans;   // events as column intervals [c0, c1]
    std::vector<char> dense;                                   // per span: a dense mismatch window
    std::vector<std::pair<std::size_t, std::size_t>> cores;    // per span: its first and last gap (or the span)
    std::vector<char> in_event(n, 0);
    // Gap clusters, widened by the mismatches within merge_gap of either end.
    for (std::size_t c = 0; c < n;) {
        if (!gap(c)) { ++c; continue; }
        std::size_t c0 = c, c1 = c;
        for (std::size_t k = c + 1; k < n && k <= c1 + merge_gap + 1; ++k) if (gap(k)) c1 = k;
        const std::size_t g0 = c0, g1 = c1;
        for (std::size_t k = g0; k-- > 0 && g0 - k <= merge_gap && !gap(k);) if (ops[k] == 'X') c0 = k;
        for (std::size_t k = g1 + 1; k < n && k - g1 <= merge_gap && !gap(k); ++k) if (ops[k] == 'X') c1 = k;
        spans.emplace_back(c0, c1);
        dense.push_back(0);
        cores.emplace_back(g0, g1);
        for (std::size_t k = c0; k <= c1; ++k) in_event[k] = 1;
        c = c1 + 1;
    }
    // Mismatch runs outside the gap clusters: each within merge_gap of the next. A dense window of one is an
    // event; a window is grown from each mismatch as far as it stays dense, so sparse mismatches beside a
    // dense stretch do not dilute it.
    std::vector<std::size_t> lone;
    for (std::size_t c = 0; c < n; ++c) if (ops[c] == 'X' && !in_event[c]) lone.push_back(c);
    for (std::size_t a = 0; a < lone.size();) {
        std::size_t b = a;
        while (b + 1 < lone.size() && lone[b + 1] - lone[b] <= merge_gap + 1) ++b;
        for (std::size_t i = a; i <= b;) {
            // The furthest mismatch keeping the window dense, within a bounded number of mismatches.
            std::size_t best = i;
            const std::size_t last = std::min(b, i + kMaxWindowMismatches);
            for (std::size_t j = last; j > i; --j)
                if (5 * (j - i + 1) >= lone[j] - lone[i] + 1) { best = j; break; }
            const std::size_t span = lone[best] - lone[i] + 1;
            if (best > i && span >= min_event) {
                spans.emplace_back(lone[i], lone[best]);
                dense.push_back(1);
                cores.emplace_back(lone[i], lone[best]);
                i = best + 1;
            }
            else ++i;
        }
        a = b + 1;
    }
    // Spans within merge_gap of each other: a dense window joins the spans beside it; two gap clusters are
    // one event only when the stretch between their gaps, any dense window there included, is one mismatch
    // in five or more (an inverted or replaced stretch). Gap clusters with mostly matches between their
    // gaps stay apart, however short a dense window between them.
    std::vector<std::size_t> order(spans.size());
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::sort(order.begin(), order.end(), [&](std::size_t x, std::size_t y) { return spans[x] < spans[y]; });
    std::vector<std::pair<std::size_t, std::size_t>> merged;
    constexpr std::size_t kNoGap = static_cast<std::size_t>(-1);
    std::size_t gap_core = kNoGap;   // the last gap in the event being built, kNoGap while it has none
    for (const std::size_t i : order) {
        const auto& sp = spans[i];
        bool join = !merged.empty() && sp.first <= merged.back().second + merge_gap + 1;
        if (join && !dense[i] && gap_core != kNoGap) {
            std::size_t x = 0;
            const std::size_t from = gap_core + 1, to = cores[i].first;   // columns strictly between the gaps
            for (std::size_t k = from; k < to; ++k) x += ops[k] == 'X';
            join = to <= from || 5 * x >= to - from;
        }
        if (join) merged.back().second = std::max(merged.back().second, sp.second);
        else { merged.push_back(sp); gap_core = kNoGap; }
        if (!dense[i]) gap_core = gap_core == kNoGap ? cores[i].second : std::max(gap_core, cores[i].second);
    }
    // Mismatches outside every event are SNVs.
    {
        std::size_t e = 0;
        for (std::size_t c = 0; c < n; ++c) {
            while (e < merged.size() && merged[e].second < c) ++e;
            if (ops[c] == 'X' && !(e < merged.size() && merged[e].first <= c)) ++out.mismatches;
        }
    }
    // Target and query positions per column, then each event's extent on both sides.
    std::size_t t = 0, q = 0, m = 0;
    for (std::size_t c = 0; c < n; ++c) {
        while (m < merged.size() && merged[m].first == c) {
            Diff d;
            d.t0 = t;
            d.q0 = q;
            std::size_t tt = t, qq = q, run = 0, chance = 0, nd = 0, ni = 0, nx = 0;
            for (std::size_t k = merged[m].first; k <= merged[m].second; ++k) {
                const char o = ops[k];
                if (o != 'I') ++tt;
                if (o != 'D') ++qq;
                if (o == 'M') { ++run; continue; }
                if (run < kChanceRun) chance += run;
                run = 0;
                if (o == 'X') ++nx;
                if (o == 'D') ++nd;
                if (o == 'I') ++ni;
            }
            d.has_del = nd > 0;
            d.has_ins = ni > 0;
            // A deletion or insertion alone: its gap bases (the aligner scatters chance matches through
            // it too). A replacement or inversion -- mismatches or gaps on both sides, in 40% of its
            // columns or more -- the bases each side holds, chance matches included. Scattered changes
            // in a repeat, with fewer edits between longer matches, count as their edits.
            const std::size_t cols = merged[m].second - merged[m].first + 1;
            const bool replaced = 4 * (std::min(nd, ni) + nx) >= std::max(nd, ni) && 5 * (nd + ni + nx) >= 2 * cols;
            d.del_bp = nd + nx + (replaced ? chance : 0);
            d.ins_bp = ni + nx + (replaced ? chance : 0);
            d.t1 = tt;
            d.q1 = qq;
            out.events.push_back(d);
            ++m;
        }
        if (ops[c] != 'I') ++t;
        if (ops[c] != 'D') ++q;
    }
    return out;
}

// Where on the target an insertion or deletion alone can lie with the same result: inside a repeat its
// bases can be written at any copy, and an aligner picks one end. [s0, s1) spans every placement; false
// for any other event, which keeps its own place.
bool slide_of(const std::string& target, const std::string& query, const Diff& d, std::size_t& s0, std::size_t& s1) {
    if (d.has_del == d.has_ins || d.t1 > target.size() || d.q1 > query.size()) return false;
    // The bases the two sides share at either end are not part of it (an aligner can thread matches
    // through a gap); what remains must be bases on one side only.
    std::size_t t0 = d.t0, t1 = d.t1, q0 = d.q0, q1 = d.q1;
    while (t0 < t1 && q0 < q1 && target[t0] == query[q0]) { ++t0; ++q0; }
    while (t0 < t1 && q0 < q1 && target[t1 - 1] == query[q1 - 1]) { --t1; --q1; }
    if ((t0 < t1) == (q0 < q1)) return false;
    if (t0 < t1) {   // target[a, b) removed
        const std::size_t a = t0, b = t1, len = b - a;
        std::size_t lo = a, hi = b;
        while (lo > 0 && target[lo - 1] == target[lo - 1 + len]) --lo;
        while (hi < target.size() && target[hi] == target[hi - len]) ++hi;
        s0 = lo; s1 = hi;
        return true;
    }
    // query[q0, q1) put in before target[t0]: moved right by k, it starts at ins[k % len]; moved left, it
    // ends at ins[len - 1 - k % len]
    const std::string ins = query.substr(q0, q1 - q0);
    const std::size_t len = ins.size();
    std::size_t lo = t0, hi = t0;
    while (hi < target.size() && target[hi] == ins[(hi - t0) % len]) ++hi;
    while (lo > 0 && target[lo - 1] == ins[len - 1 - (t0 - lo) % len]) --lo;
    s0 = lo; s1 = hi;
    return true;
}

// Global alignment of `target` and `query` cut at blocks both hold exactly: (target pos, query pos,
// length), increasing in both. Each stretch between blocks is aligned on its own and the blocks are
// matches, so a large copy-number difference costs only its own stretch, not the whole region.
std::string aligned_ops(const std::string& target, const std::string& query,
                        const std::vector<std::array<std::size_t, 3>>& blocks) {
    std::string ops, part;
    std::size_t t = 0, q = 0;
    for (const auto& [bt, bq, len] : blocks) {
        if (bt < t || bq < q || bt + len > target.size() || bq + len > query.size()) continue;
        if (target.compare(bt, len, query, bq, len) != 0) continue;
        nw_edit_ops(target.substr(t, bt - t), query.substr(q, bq - q), part);
        ops += part;
        ops.append(len, 'M');
        t = bt + len;
        q = bq + len;
    }
    nw_edit_ops(target.substr(t), query.substr(q), part);
    ops += part;
    return ops;
}

// Exact blocks to cut an alignment at: runs of 500 bases or more whose 31-mers occur once in `target` and
// once in `query`, on one diagonal, chained co-linearly to hold the most bases. Sequence anchors, not
// graph nodes, so a copy call placed at another offset in a repeat than the walk does is not forced
// there; a k-mer inside a repeated copy is not unique, so the aligner stays free to place copies.
// `conflict` is set when a run left out of the chain crosses it: the sequences hold a block in another
// order (a moved or swapped segment), and the chain may then not lie on the best alignment.
std::vector<std::array<std::size_t, 3>> kmer_blocks(const std::string& target, const std::string& query,
                                                    bool& conflict) {
    conflict = false;
    constexpr std::size_t k = 31;
    const auto code = [](char c) -> int {
        switch (c) { case 'A': case 'a': return 0; case 'C': case 'c': return 1; case 'G': case 'g': return 2;
                     case 'T': case 't': return 3; default: return -1; }
    };
    const auto scan = [&](const std::string& s, const auto& fn) {
        std::uint64_t v = 0;
        std::size_t run = 0;
        const std::uint64_t mask = (std::uint64_t{1} << (2 * k)) - 1;
        for (std::size_t i = 0; i < s.size(); ++i) {
            const int c = code(s[i]);
            if (c < 0) { run = 0; v = 0; continue; }
            v = ((v << 2) | static_cast<std::uint64_t>(c)) & mask;
            if (++run >= k) fn(v, i + 1 - k);
        }
    };
    std::unordered_map<std::uint64_t, std::size_t> tpos;   // k-mer -> position, SIZE_MAX when repeated
    tpos.reserve(target.size());
    scan(target, [&](std::uint64_t v, std::size_t p) {
        auto [it, fresh] = tpos.try_emplace(v, p);
        if (!fresh) it->second = SIZE_MAX;
    });
    std::unordered_map<std::uint64_t, std::size_t> qcount;
    qcount.reserve(query.size());
    scan(query, [&](std::uint64_t v, std::size_t) { ++qcount[v]; });
    std::vector<std::pair<std::size_t, std::size_t>> hits;   // (query pos, target pos), query order
    scan(query, [&](std::uint64_t v, std::size_t p) {
        const auto it = tpos.find(v);
        if (it != tpos.end() && it->second != SIZE_MAX && qcount[v] == 1) hits.emplace_back(p, it->second);
    });
    // Only long runs anchor: in unique sequence they run for as long as the sequence is unique, while
    // inside a repeat a copy-specific difference makes a run of a few dozen bases.
    constexpr std::size_t kMinRun = 500;
    std::vector<std::array<std::size_t, 3>> runs;   // (target, query, length), query order
    for (const auto& [q, t] : hits) {
        if (!runs.empty()) {
            auto& r = runs.back();
            if (q == r[1] + r[2] - k + 1 && t == r[0] + r[2] - k + 1) { ++r[2]; continue; }
        }
        runs.push_back({t, q, k});
    }
    // Shorter unique runs still show order: a segment that moved can hold no run of kMinRun (a few
    // differences in it break its runs up) and still cross the chain. Runs of kMinCross or more are kept to
    // look for that.
    constexpr std::size_t kMinCross = 200;
    std::vector<std::array<std::size_t, 3>> short_runs;
    for (const auto& r : runs) if (r[2] >= kMinCross && r[2] < kMinRun) short_runs.push_back(r);
    runs.erase(std::remove_if(runs.begin(), runs.end(), [](const auto& r) { return r[2] < kMinRun; }), runs.end());
    // The chain of non-overlapping runs, increasing in both, holding the most bases.
    const std::size_t nr = runs.size();
    std::vector<std::size_t> best(nr, 0), prev(nr, SIZE_MAX);
    for (std::size_t a = 0; a < nr; ++a) {
        best[a] = runs[a][2];
        for (std::size_t b = 0; b < a; ++b)
            if (runs[b][0] + runs[b][2] <= runs[a][0] && runs[b][1] + runs[b][2] <= runs[a][1] &&
                best[b] + runs[a][2] > best[a]) { best[a] = best[b] + runs[a][2]; prev[a] = b; }
    }
    std::vector<std::array<std::size_t, 3>> blocks;
    if (nr == 0) return blocks;
    std::size_t top = 0;
    for (std::size_t a = 1; a < nr; ++a) if (best[a] > best[top]) top = a;
    std::vector<char> chained(nr, 0);
    for (std::size_t a = top; a != SIZE_MAX; a = prev[a]) { blocks.push_back(runs[a]); chained[a] = 1; }
    std::reverse(blocks.begin(), blocks.end());
    // A run left out crosses the chain when its middle lies before a chained run on one sequence and
    // after it on the other; two runs beside a small indel only overlap by a few bases.
    const auto crosses = [&](const std::array<std::size_t, 3>& r) {
        const std::size_t mt = r[0] + r[2] / 2, mq = r[1] + r[2] / 2;
        for (const auto& b : blocks)
            if ((mt < b[0] && mq >= b[1] + b[2]) || (mt >= b[0] + b[2] && mq < b[1])) return true;
        return false;
    };
    for (std::size_t a = 0; a < nr && !conflict; ++a) if (!chained[a] && crosses(runs[a])) conflict = true;
    for (std::size_t a = 0; a < short_runs.size() && !conflict; ++a) if (crosses(short_runs[a])) conflict = true;
    return blocks;
}

// The alignment of `target` and `query`, as edit operations. Cut at exact blocks when they are all in
// order; when one is out of order, the cut alignment is only an upper bound, and the exact one is found
// within it.
std::string align_ops(const std::string& target, const std::string& query) {
    bool conflict = false;
    const std::vector<std::array<std::size_t, 3>> blocks = kmer_blocks(target, query, conflict);
    std::string ops = aligned_ops(target, query, blocks);
    if (!conflict) return ops;
    const std::size_t bound = static_cast<std::size_t>(std::count_if(ops.begin(), ops.end(), [](char o) { return o != 'M'; }));
    if (bound == 0) return ops;
    std::string exact;
    if (nw_edit_ops_bounded(target, query, bound - 1, exact) != SIZE_MAX) return exact;
    return ops;
}

// Where a stretch of the reference repeats itself with a period near `unit`: the offset d at which the
// most 31-mers recur -- within a tenth of `unit` if any recur there, else from unit/2 to 3 unit/2 -- and
// the first position x0 of a k-mer recurring at d. Copies that share segments at more than one offset
// (one copy carrying an insertion the other lacks) give several; `unit`, call's own measure of a copy,
// picks between them. A copy is then [x0, x0 + d): deleting it or adding it there is what a recombination
// between two copies does. `forms` gets every offset with a quarter of the strongest one's support, with
// its x0 and support: the reference's copies when they come in several lengths. False when nothing recurs
// (a reference with one copy).
struct CopyForm { std::size_t start = 0, length = 0, support = 0; };
bool period_of(const std::string& seq, std::size_t unit, std::size_t& d, std::size_t& x0,
               std::vector<CopyForm>& forms) {
    constexpr std::size_t k = 31, kMaxOcc = 64;
    if (seq.size() < 2 * k || unit < 2) return false;
    const std::size_t dlo = unit / 2, dhi = unit + unit / 2;
    std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> at;
    std::uint64_t v = 0;
    std::size_t run = 0;
    const std::uint64_t mask = (std::uint64_t{1} << (2 * k)) - 1;
    for (std::size_t i = 0; i < seq.size(); ++i) {
        int c = -1;
        switch (seq[i]) { case 'A': case 'a': c = 0; break; case 'C': case 'c': c = 1; break;
                          case 'G': case 'g': c = 2; break; case 'T': case 't': c = 3; break; default: break; }
        if (c < 0) { run = 0; v = 0; continue; }
        v = ((v << 2) | static_cast<std::uint64_t>(c)) & mask;
        if (++run >= k) at[v].push_back(static_cast<std::uint32_t>(i + 1 - k));
    }
    std::unordered_map<std::size_t, std::size_t> votes;
    std::unordered_map<std::size_t, std::size_t> first;
    for (const auto& [kmer, pos] : at) {
        if (pos.size() < 2 || pos.size() > kMaxOcc) continue;
        for (std::size_t a = 0; a < pos.size(); ++a)
            for (std::size_t b = a + 1; b < pos.size(); ++b) {
                const std::size_t dd = pos[b] - pos[a];
                if (dd < dlo || dd > dhi) continue;
                ++votes[dd];
                const auto [it, fresh] = first.try_emplace(dd, pos[a]);
                if (!fresh) it->second = std::min<std::size_t>(it->second, pos[a]);
            }
    }
    // An offset near `unit` is taken only with a quarter of the strongest offset's support: a few
    // k-mers recurring there by chance must not stand for the copies.
    std::size_t strongest = 0;
    for (const auto& [dd, n] : votes) strongest = std::max(strongest, n);
    std::size_t best = 0;
    for (const bool near : {true, false}) {
        for (const auto& [dd, n] : votes) {
            if (near && 10 * (dd > unit ? dd - unit : unit - dd) > unit) continue;
            if (n > best || (n == best && dd < d)) { best = n; d = dd; }
        }
        if (near && 4 * best < strongest) best = 0;
        if (best > 0) break;
    }
    if (best == 0) return false;
    x0 = first.at(d);
    forms.clear();
    for (const auto& [dd, n] : votes) if (4 * n >= strongest) forms.push_back({first.at(dd), dd, n});
    std::sort(forms.begin(), forms.end(), [](const CopyForm& a, const CopyForm& b) { return a.length < b.length; });
    return true;
}

// A piece of a reconstruction: target positions [r0, r1) copy the reference from g0 (genomic), or are
// an applied event replacing [g0, g1] (genomic, g1 = g0 - 1 for an insertion).
struct Piece {
    std::size_t r0 = 0, r1 = 0;
    std::size_t g0 = 0, g1 = 0;
    int event = -1;                 // index into the haplotype's applied events, -1 for reference
};

struct LargeRow {
    std::string haplotype;
    std::size_t g0 = 0, g1 = 0;     // genomic interval the difference sits at
    // where it may sit: [g0, g1], or the whole of an applied event it lies in, whose bases are the
    // carrier's own and map to the reference only in proportion
    std::size_t m0 = 0, m1 = 0;
    std::size_t w0 = 0, w1 = 0;     // [g0, g1] widened over every place an indel alone can slide to
    std::string type;
    std::size_t size = 0, del_bp = 0, ins_bp = 0;
    std::string cause;              // incomplete, not_carrier or not_reported
    bool unwritten = false;         // at a record of its own whose event could not be placed or spelled
    std::string bubbles = ".";
    std::string records = ".";
};

struct HapResult {
    std::string name;
    std::string status = "scored";
    std::size_t g_start = 0, g_end = 0;
    std::size_t applied = 0, conflicts = 0, outside_extent = 0, unresolved = 0, absorbed = 0;
    std::size_t present_n = 0, present_bp = 0, present_small_bp = 0, present_snv = 0;
    std::size_t left_n = 0, left_bp = 0, left_small_bp = 0, left_snv = 0;
    struct Span { std::string id; std::size_t lo = 0, hi = 0; };
    std::vector<Span> applied_spans;        // where each applied event sits (genomic)
    std::set<std::string> outside_ids;      // records whose event lies past the scored stretch
    std::vector<LargeRow> rows;
    // the large differences between the plain reference and the haplotype: what there was to find
    struct Present { std::size_t g0 = 0, g1 = 0, size = 0, del_bp = 0, ins_bp = 0; const char* type = "DEL"; std::string bubbles = "."; };
    std::vector<Present> present_rows;
};

struct CheckRow {
    std::string variant_id, haplotype, kind, reported, own;
    long long difference = 0;
};

} // namespace

int run_benchmark_command(const std::vector<std::string>& args) {
    if (args.empty()) { print_benchmark_help(); return 0; }

    std::string gfa_path, bubble_prefix_in, bubbles_csv_in, calls_prefix, reference_path, out_prefix;
    std::size_t threads = 0, min_sv_bp = 50, merge_gap = 10;
    bool quiet = false;
    for (std::size_t i = 0; i < args.size(); ++i) {
        const std::string arg = cli::renamed_input_flag(args[i], "benchmark");
        auto require_value = [&](const std::string& flag) -> const std::string& {
            if (i + 1 >= args.size()) throw std::runtime_error("Missing value after " + flag);
            return args[++i];
        };
        if (arg == "-h" || arg == "--help") { print_benchmark_help(); return 0; }
        else if (arg == "-i" || arg == "--gfa") gfa_path = require_value(arg);
        else if (arg == "-b" || arg == "--bubble-prefix") bubble_prefix_in = require_value(arg);
        else if (arg == "-c" || arg == "--bubbles-csv") bubbles_csv_in = require_value(arg);
        else if (arg == "--calls") calls_prefix = require_value(arg);
        else if (arg == "-r" || arg == "--reference-path") reference_path = require_value(arg);
        else if (arg == "-o" || arg == "--out-prefix") out_prefix = require_value(arg);
        else if (arg == "--min-sv-bp") min_sv_bp = cli::parse_size_arg(arg, require_value(arg));
        else if (arg == "--merge-gap") merge_gap = cli::parse_size_arg(arg, require_value(arg));
        else if (arg == "-t" || arg == "--threads") threads = cli::parse_size_arg(arg, require_value(arg));
        else if (arg == "-q" || arg == "--quiet") quiet = true;
        else throw std::runtime_error("Unknown option for benchmark: " + arg);
    }
    if (gfa_path.empty()) throw std::runtime_error("Missing required input: --gfa <path>");
    if (!bubble_prefix_in.empty()) {
        const std::string derived = bubble_prefix_in + ".bubbles.csv";
        if (bubbles_csv_in.empty()) bubbles_csv_in = derived;
        else if (bubbles_csv_in != derived)
            throw std::runtime_error("Conflicting bubble inputs: --bubble-prefix resolves to '" + derived +
                                     "' but --bubbles-csv is '" + bubbles_csv_in + "'");
    }
    if (bubbles_csv_in.empty()) throw std::runtime_error("benchmark requires --bubble-prefix <prefix> or --bubbles-csv <path>");
    if (calls_prefix.empty()) throw std::runtime_error("benchmark requires --calls <call_prefix>");
    if (reference_path.empty()) throw std::runtime_error("benchmark requires --reference-path <name>");
    if (out_prefix.empty()) throw std::runtime_error("benchmark requires -o/--out-prefix <prefix>");
    if (min_sv_bp == 0) throw std::runtime_error("--min-sv-bp must be at least 1");
    // Two stages from one call run: each carrier's own events (what call found), then the VCF's records
    // (what the file keeps). Both are scored and written together.
    struct ModeSpec {
        bool vcf = false;
        std::string vcf_in, ce_in;
    };
    const std::string vcf_in_path = calls_prefix + ".region.vcf", ce_in_path = calls_prefix + ".carrier_events.tsv";
    const std::vector<ModeSpec> specs = {{false, vcf_in_path, ce_in_path}, {true, vcf_in_path, std::string()}};
    const std::vector<std::string> finals = {out_prefix + ".summary.tsv", out_prefix + ".haplotypes.tsv",
                                             out_prefix + ".events.tsv", out_prefix + ".carrier_checks.tsv",
                                             out_prefix + ".present.tsv"};
    for (const std::string& f : finals)
        for (const std::string& in : {gfa_path, bubbles_csv_in, vcf_in_path, ce_in_path}) {
            std::error_code ec;
            if (f == in || std::filesystem::equivalent(f, in, ec))
                throw std::runtime_error("benchmark: output '" + f + "' is also an input; refusing to "
                                         "overwrite what the run reads");
        }

    ParseGfaOptions parse_options;
    parse_options.include_paths = true;
    parse_options.include_sequences = true;
    const Graph graph = parse_gfa(gfa_path, parse_options);
    validate_graph_paths(graph, "benchmark", true, true);
    const std::string ref_name = resolve_reference_path_name(graph, reference_path, "benchmark");
    const PathRecord* ref_path = nullptr;
    for (const PathRecord& p : graph.paths) if (p.name == ref_name) ref_path = &p;

    cli::RunLog log("benchmark", quiet);
    log.info("reading: " + gfa_path);
    log.info("input: haplotypes: " + std::to_string(graph.paths.size()) + "; nodes: " +
             std::to_string(graph.nodes.size()));
    log.info("reference: " + ref_name);
    log.info("calls: " + calls_prefix);

    // ---- the reference, and each bubble's span on it ----
    const std::string ref_seq = spell_path_steps_sequence(graph, ref_path->steps);
    // The name's interval convention follows from the path's length (see parse_reference_path_label).
    const ParsedReferencePath ref_meta = parse_reference_path_label(ref_name, ref_seq.size(), ref_path->type == 'W');
    if (ref_meta.length_mismatch)
        std::cerr << cli::log_prefix("benchmark") << "warning: the interval in reference path name '" << ref_name
                  << "' does not match the path's " << ref_seq.size() << " bp; reading its start as "
                  << (ref_meta.half_open ? "0-based (half-open)" : "the W line's start, plus one") << "\n";
    if (!ref_meta.has_interval)
        std::cerr << cli::log_prefix("benchmark") << "warning: reference path name '" << ref_name << "' has no :start-end interval; "
                  << "positions are counted from the path's first base (1)\n";
    const std::size_t region_start = ref_meta.region_start_1based;
    std::vector<std::size_t> ref_prefix(ref_path->steps.size() + 1, 0);
    for (std::size_t k = 0; k < ref_path->steps.size(); ++k)
        ref_prefix[k + 1] = ref_prefix[k] + graph.nodes.at(ref_path->steps[k].node_id).sequence.size();
    const std::vector<Bubble> bubbles = read_bubbles_csv(bubbles_csv_in, /*with_inversion_sites=*/true);
    // The sites CSV and the GFA are separate inputs: a table from another graph would leave every
    // difference outside the sites, so one naming nodes this graph lacks, or one id twice, is refused.
    {
        std::unordered_set<std::size_t> seen_ids;
        std::vector<std::string> missing;
        for (const Bubble& b : bubbles) {
            if (!seen_ids.insert(b.id).second)
                throw std::runtime_error("benchmark: duplicate bubble id in " + bubbles_csv_in + ": " + std::to_string(b.id));
            const auto check = [&](const std::string& n) {
                if (graph.nodes.find(n) == graph.nodes.end() && missing.size() < 8)
                    missing.push_back("bubble " + std::to_string(b.id) + " node " + n);
            };
            check(b.source);
            check(b.sink);
            for (const std::string& n : b.inside) check(n);
        }
        if (!missing.empty())
            throw std::runtime_error("benchmark: the bubbles CSV describes nodes the GFA does not contain (" +
                                     cli::join_with_comma(missing) + "); the CSV and the graph are not the same graph");
    }
    std::vector<std::pair<std::size_t, std::pair<std::size_t, std::size_t>>> bubble_span;   // id -> genomic [lo, hi]
    {
        const BubblePathIndex ridx = build_bubble_path_index(*ref_path);
        for (const Bubble& b : bubbles) {
            std::vector<std::size_t> used;
            if (!bubble_steps(*ref_path, ridx, b, nullptr, &used) || used.empty()) continue;
            // The interior, from the last base of the near boundary to the first of the far one: the
            // boundary nodes are shared sequence.
            const auto [lo, hi] = std::minmax_element(used.begin(), used.end());
            bubble_span.push_back({b.id, {region_start + ref_prefix[*lo + 1] - 1, region_start + ref_prefix[*hi]}});
        }
    }

    // ---- each haplotype, once for every mode: its strand, the stretch it shares with the reference, and
    // the alignment of the plain reference over it ----
    std::unordered_map<std::string, std::size_t> ref_count;
    for (const PathStep& s : ref_path->steps) ++ref_count[s.node_id];
    std::unordered_map<std::string, std::size_t> ref_pos_of;   // unique reference node -> step index
    for (std::size_t k = 0; k < ref_path->steps.size(); ++k)
        if (ref_count[ref_path->steps[k].node_id] == 1) ref_pos_of[ref_path->steps[k].node_id] = k;
    struct Prep {
        std::string status = "not_prepared";   // "scored" once its alignment is done
        std::size_t g_start = 0, g_end = 0, o0 = 0, o1 = 0;
        std::string truth;
        DiffSet base;
    };
    std::vector<Prep> preps(graph.paths.size());
    // Only the paths the VCF has a column for are scored.
    std::unordered_set<std::string> in_vcf_header;
    {
        std::ifstream vh(vcf_in_path);
        if (!vh) throw std::runtime_error("cannot open " + vcf_in_path);
        for (std::string l; std::getline(vh, l);) {
            if (!l.empty() && l.back() == '\r') l.pop_back();
            if (l.rfind("#CHROM", 0) == 0) {
                const std::vector<std::string> f = split_on(l, '\t');
                in_vcf_header.insert(f.begin() + std::min<std::size_t>(9, f.size()), f.end());
                break;
            }
        }
    }
    {
        std::size_t to_align = 0;
        for (const PathRecord& p : graph.paths) to_align += &p != ref_path && in_vcf_header.count(p.name);
        log.info("aligning to the reference: " + std::to_string(to_align) + " haplotypes");
    }
    run_parallel(graph.paths.size(), threads, [&](std::size_t pi) {
        const PathRecord& path = graph.paths[pi];
        if (&path == ref_path || !in_vcf_header.count(path.name)) return;
        Prep& pr = preps[pi];
        // Orientation and extent from the nodes the haplotype and the reference each hold once: the
        // haplotype is turned to the reference's strand, and both are cut to the longest co-linear
        // chain of such nodes, so a partial assembly is not charged for the sequence it does not have.
        std::unordered_map<std::string, std::size_t> hcount;
        for (const PathStep& s : path.steps) ++hcount[s.node_id];
        std::size_t same = 0, opposite = 0;
        int first_same = -1, last_same = -1;   // the orientation of the haplotype's first and last anchors
        for (const PathStep& s : path.steps) {
            const auto it = ref_pos_of.find(s.node_id);
            if (it == ref_pos_of.end() || hcount[s.node_id] != 1) continue;
            const bool sm = s.reverse == ref_path->steps[it->second].reverse;
            (sm ? same : opposite) += 1;
            if (first_same < 0) first_same = sm;
            last_same = sm;
        }
        if (same + opposite == 0) { pr.status = "no_shared_anchor"; return; }
        std::vector<PathStep> hsteps = path.steps;
        // The haplotype's direction comes from its two ends when they agree, else from most of its anchors:
        // an inversion longer than both flanks together would outvote them.
        const bool flip = first_same == last_same ? first_same == 0 : opposite > same;
        if (flip) {
            std::reverse(hsteps.begin(), hsteps.end());
            for (PathStep& s : hsteps) s.reverse = !s.reverse;
        }
        std::vector<std::pair<std::size_t, std::size_t>> anchors;   // (hap step, ref step)
        for (std::size_t k = 0; k < hsteps.size(); ++k) {
            const auto it = ref_pos_of.find(hsteps[k].node_id);
            if (it != ref_pos_of.end() && hcount[hsteps[k].node_id] == 1 &&
                hsteps[k].reverse == ref_path->steps[it->second].reverse)
                anchors.emplace_back(k, it->second);
        }
        // longest chain increasing in both (patience sort on the reference index)
        std::vector<std::size_t> tails, tail_idx, prev(anchors.size(), SIZE_MAX);
        for (std::size_t a = 0; a < anchors.size(); ++a) {
            const auto pos = std::lower_bound(tails.begin(), tails.end(), anchors[a].second) - tails.begin();
            if (static_cast<std::size_t>(pos) == tails.size()) { tails.push_back(anchors[a].second); tail_idx.push_back(a); }
            else { tails[static_cast<std::size_t>(pos)] = anchors[a].second; tail_idx[static_cast<std::size_t>(pos)] = a; }
            if (pos > 0) prev[a] = tail_idx[static_cast<std::size_t>(pos) - 1];
        }
        if (tail_idx.empty()) { pr.status = "no_shared_anchor"; return; }
        std::size_t last = tail_idx.back(), first = last;
        while (prev[first] != SIZE_MAX) first = prev[first];
        std::size_t h0 = anchors[first].first, h1 = anchors[last].first;
        std::size_t r0 = anchors[first].second, r1 = anchors[last].second;
        // A haplotype whose walk starts (ends) on the reference's own first (last) node, the same way, is
        // complete at that end: scored to the region's edge, so a deletion right after it counts.
        const auto same_step = [](const PathStep& a, const PathStep& b) {
            return a.node_id == b.node_id && a.reverse == b.reverse;
        };
        if (same_step(hsteps.front(), ref_path->steps.front())) { h0 = 0; r0 = 0; }
        if (same_step(hsteps.back(), ref_path->steps.back())) { h1 = hsteps.size() - 1; r1 = ref_path->steps.size() - 1; }
        pr.truth = spell_path_steps_sequence(
            graph, std::vector<PathStep>(hsteps.begin() + static_cast<std::ptrdiff_t>(h0),
                                         hsteps.begin() + static_cast<std::ptrdiff_t>(h1) + 1));
        pr.o0 = ref_prefix[r0];
        pr.o1 = ref_prefix[r1 + 1];   // reference offsets [o0, o1)
        pr.g_start = region_start + pr.o0;
        pr.g_end = region_start + pr.o1 - 1;
        pr.base = diff_of(align_ops(ref_seq.substr(pr.o0, pr.o1 - pr.o0), pr.truth), merge_gap, min_sv_bp);
        pr.status = "scored";
    });

    // What one stage leaves, per haplotype, and its counts.
    struct Stage {
        std::vector<HapResult> results;
        std::vector<CheckRow> checks;
        std::size_t checked = 0, unchecked = 0, cn_checked = 0, cn_unchecked = 0, cn_exact_checked = 0, unresolved_rows = 0;
        std::size_t not_reconstructable = 0, missing_gt = 0, not_in_vcf = 0;
        std::vector<VcfRecordInfo> vcf;
        // A record that could not be spelled for some carrier, and why: "sequence" (left out above call's
        // --max-seq-bp), "no_copy" (a copy-number record with no reference copy to paste), "other".
        std::unordered_map<std::string, std::string> unspelled;
    };
    Stage stages[2];
    const auto score_mode = [&](const ModeSpec& spec, Stage& out) {
    const bool vcf_mode = spec.vcf;
    const std::string& vcf_in = spec.vcf_in;
    const std::string& ce_in = spec.ce_in;
    // ---- call's output: the VCF and each carrier's own event ----
    std::vector<std::string> vcf_samples;
    const std::vector<VcfRecordInfo> vcf = load_region_vcf(vcf_in, vcf_samples);
    std::unordered_map<std::string, const VcfRecordInfo*> vcf_by_id;
    // Records are matched to carrier events and to what is left by their ID.
    for (const VcfRecordInfo& r : vcf)
        if (r.id.empty() || r.id == "." || !vcf_by_id.emplace(r.id, &r).second)
            throw std::runtime_error(vcf_in + ": record at POS " + std::to_string(r.pos) + " has " +
                                     (r.id.empty() || r.id == "." ? std::string("no ID") : "the repeated ID '" + r.id + "'") +
                                     "; benchmark names each record by a unique ID");
    // --vcf: each record's written allele, for every haplotype whose GT names it. A copy-number record has
    // no sequence: a gain is (CN - REF_CN) tandem copies of one reference copy, a loss deletes that many.
    // One reference copy, genomic (start, length): one period where the reference's copies repeat, else
    // the copies' last RU_LEN / CN_UNIT_BP bases; with one reference copy, its span or those bases. The
    // copies run from POS+1 to END. Real copies differ in length and content, so this approximates.
    // A carrier of a MODULE record with FORMAT:CNBP gets the form of the reference's copies that fits its own
    // change per copy: paralogous copies come in several lengths, and one length cannot stand for all its
    // carriers. CNBP is the net change over the copies, so the carrier's other records there come off it
    // first. A UNIT record's copies are one period.
    // period: the copy is a period at which the reference's copies repeat, so a whole copy away is another
    struct CopyForms { std::pair<std::size_t, std::size_t> fallback; std::vector<CopyForm> forms; bool period = false; };
    std::map<std::string, CopyForms> copy_cache;
    const auto copy_of = [&](const VcfRecordInfo& r, const std::string& h) -> std::pair<std::size_t, std::size_t> {
        auto it = copy_cache.find(r.id);
        if (it == copy_cache.end()) {
            CopyForms cf;
            const std::size_t span = r.cn_ref_bp > 0 ? static_cast<std::size_t>(r.cn_ref_bp) : std::max(r.end, r.pos) - r.pos;
            const std::size_t aend = r.pos + span, unit = static_cast<std::size_t>(std::max(0LL, r.unit_bp));
            cf.fallback = {0, unit};   // start 0: its length known, no reference copy to spell
            if (r.ref_cn == 1 && r.cn_ref_bp > 0 && r.pos >= region_start && aend < region_start + ref_seq.size()) {
                cf.fallback = {r.pos + 1, span};   // the reference holds one copy: the span is it
                // unless the copies are a shorter stretch inside it: the carrier's change then tells
                if (unit > 0 && unit < span) cf.forms = {{r.pos + 1, span, 1}, {aend - unit + 1, unit, 1}};
            } else if (unit > 0 && unit <= span && r.pos >= region_start && aend < region_start + ref_seq.size()) {
                cf.fallback = {aend - unit + 1, unit};
                // The period is looked for near the span over the reference's copies: CN_UNIT_BP is only
                // what the copies share, short of a copy where they differ.
                const std::size_t hint = r.ref_cn >= 2 && r.cn_ref_bp > 0 ? span / static_cast<std::size_t>(r.ref_cn) : unit;
                std::size_t pd = 0, px = 0;
                std::vector<CopyForm> forms;
                if (r.ref_cn >= 2 && period_of(ref_seq.substr(r.pos + 1 - region_start, span), hint, pd, px, forms)) {
                    cf.fallback = {r.pos + 1 + px, pd};
                    cf.period = true;
                    for (CopyForm& f : forms) { f.start += r.pos + 1; cf.forms.push_back(f); }
                }
            }
            it = copy_cache.emplace(r.id, std::move(cf)).first;
        }
        const CopyForms& cf = it->second;
        const auto cn = r.sample_cn.find(h);
        const auto bp = r.sample_cnbp.find(h);
        if (cf.forms.empty() || cn == r.sample_cn.end() || bp == r.sample_cnbp.end() || cn->second == r.ref_cn ||
            info_field(r.info, "CN_METHOD") != "MODULE")
            return cf.fallback;
        long long own = bp->second;
        const std::size_t span_end = r.pos + (r.cn_ref_bp > 0 ? static_cast<std::size_t>(r.cn_ref_bp) : r.end - r.pos);
        for (const VcfRecordInfo& o : vcf) {
            if (&o == &r || !change_inside(o, r.pos + 1, span_end)) continue;   // changing the copies
            const auto a = o.sample_allele.find(h);
            if (a == o.sample_allele.end() || a->second <= 0 || static_cast<std::size_t>(a->second) > o.alts.size()) continue;
            const std::string& alt = o.alts[static_cast<std::size_t>(a->second) - 1];
            if (alt.empty() || alt[0] != '<')
                own -= change_within(o, static_cast<long long>(alt.size()) - static_cast<long long>(o.ref.size()), r.pos + 1, span_end);
            else if (alt.rfind("<CNV", 0) != 0 && alt != "<INV>" && o.has_svlen) own -= change_within(o, o.svlen, r.pos + 1, span_end);   // an INV's SVLEN changes no length
        }
        const long long per = std::llabs(own) / std::llabs(cn->second - r.ref_cn);
        // Lengths within 2% (100 bp at least) of the carrier's are one form, many offsets of one copy: the
        // best supported of them, so a shifted copy is not taken for a closer length; of equal support, the nearer.
        const auto gap = [&](const CopyForm& f) { return std::llabs(static_cast<long long>(f.length) - per); };
        const long long tol = std::max(100LL, per / 50);
        const CopyForm* best = nullptr;
        for (const CopyForm& f : cf.forms)
            if (gap(f) <= tol && (best == nullptr || f.support > best->support || (f.support == best->support && gap(f) < gap(*best)))) best = &f;
        if (best == nullptr)
            best = &*std::min_element(cf.forms.begin(), cf.forms.end(),
                                      [&](const CopyForm& x, const CopyForm& y) { return gap(x) < gap(y); });
        return {best->start, best->length};
    };
    std::size_t vcf_not_reconstructable = 0, vcf_missing_gt = 0;
    std::unordered_map<std::string, std::string> unspelled;
    std::vector<CarrierEvent> cevents;
    if (vcf_mode) {
        const std::size_t ref_lo = region_start, ref_hi = region_start + ref_seq.size() - 1;   // genomic
        auto ref_at = [&](std::size_t g0, std::size_t g1) {   // genomic [g0, g1], g1 >= g0 - 1
            return g1 + 1 <= g0 ? std::string() : ref_seq.substr(g0 - region_start, g1 - g0 + 1);
        };
        std::map<std::tuple<std::string, std::size_t, std::size_t, std::string>, std::size_t> row_of;
        auto add = [&](const VcfRecordInfo& r, const char* kind, std::size_t g0, std::size_t g1,
                       const std::string& alt, const std::string& cn, const std::string& h) {
            const auto key = std::make_tuple(r.id + '\t' + kind + '\t' + cn, g0, g1, alt);
            auto it = row_of.find(key);
            if (it == row_of.end()) {
                CarrierEvent e;
                e.order = static_cast<std::size_t>(&r - vcf.data());
                e.variant_id = r.id; e.kind = kind; e.ref_start = g0; e.ref_end = g1; e.cn = cn;
                e.walk_cn = "."; e.alt = alt;
                e.event_bp = static_cast<long long>(alt.size()) - static_cast<long long>(g1 + 1 - g0);
                it = row_of.emplace(key, cevents.size()).first;
                cevents.push_back(std::move(e));
            }
            cevents[it->second].haplotypes.push_back(h);
        };
        for (const VcfRecordInfo& r : vcf) {
            bool reconstructable = true;
            for (const auto& [h, a] : r.sample_allele) {
                if (a < 0) { ++vcf_missing_gt; continue; }
                if (a == 0) continue;
                if (static_cast<std::size_t>(a) > r.alts.size())
                    throw std::runtime_error(vcf_in + ": record " + r.id + " GT names an allele it does not have");
                const std::string& alt = r.alts[static_cast<std::size_t>(a) - 1];
                const bool subst = has_flag(r.info, "SUBSTITUTION");
                const std::size_t end = std::max(r.end, r.pos);
                if (r.pos < ref_lo || end > ref_hi) { unspelled.emplace(r.id, "other"); reconstructable = false; add(r, "unresolved", 0, 0, std::string(), ".", h); continue; }
                if (alt.empty() || alt[0] != '<') {
                    // sequence-resolved: REF at POS replaced by ALT, their shared first bases trimmed
                    const auto upper = [](std::string x) {
                        for (char& ch : x) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
                        return x;
                    };
                    if (upper(ref_at(r.pos, r.pos + r.ref.size() - 1)) != upper(r.ref))
                        throw std::runtime_error(vcf_in + ": record " + r.id + " REF does not match the reference at POS "
                                                 + std::to_string(r.pos));
                    std::size_t c = 0;
                    while (c < r.ref.size() && c < alt.size() &&
                           std::toupper(static_cast<unsigned char>(r.ref[c])) == std::toupper(static_cast<unsigned char>(alt[c]))) ++c;
                    add(r, "seq", r.pos + c, r.pos + r.ref.size() - 1, alt.substr(c), ".", h);
                } else if (alt == "<DEL>" || alt == "<INS>") {
                    const std::string ins = info_field(r.info, "INSSEQ");
                    // An insertion, or a replacement, spells its new sequence from INSSEQ.
                    if ((alt == "<INS>" || subst) && ins.empty()) { unspelled.emplace(r.id, "sequence"); reconstructable = false; add(r, "unresolved", 0, 0, std::string(), ".", h); continue; }
                    if (alt == "<INS>" && !subst) add(r, "ins", r.pos + 1, r.pos, ins, ".", h);
                    else add(r, alt == "<DEL>" && !subst ? "del" : "sub", r.pos + 1, end, subst ? ins : std::string(), ".", h);
                } else if (alt == "<INV>") {
                    add(r, "inv", r.pos + 1, end, reverse_complement(ref_at(r.pos + 1, end)), ".", h);
                } else if (alt.rfind("<CNV", 0) == 0 && !r.cn_forms.empty()) {
                    // The record names the carrier's copies: each pasted as written. A carrier with none named
                    // holds its copies in no change of its walk.
                    const auto cit = r.sample_cn.find(h);
                    const std::string cn = cit == r.sample_cn.end() ? std::string(".") : std::to_string(cit->second);
                    const auto fit = r.sample_cnform.find(h);
                    if (fit == r.sample_cnform.end()) continue;
                    for (const std::size_t k : fit->second) {
                        const VcfRecordInfo::Form& fm = r.cn_forms[k];
                        if (!fm.spelled || fm.start < ref_lo || fm.end > ref_hi) {
                            unspelled.emplace(r.id, fm.spelled ? "other" : "sequence");
                            reconstructable = false; add(r, "unresolved", 0, 0, std::string(), ".", h); continue;
                        }
                        add(r, "cnv", fm.start, fm.end, fm.seq, cn, h);
                    }
                } else if (alt.rfind("<CNV", 0) == 0) {
                    const auto cit = r.sample_cn.find(h);
                    const auto [c0, clen] = copy_of(r, h);
                    const std::size_t aend = r.pos + (r.cn_ref_bp > 0 ? static_cast<std::size_t>(r.cn_ref_bp) : end - r.pos);
                    if (cit == r.sample_cn.end() || r.ref_cn <= 0 || c0 == 0 || clen == 0) {
                        unspelled.emplace(r.id, cit == r.sample_cn.end() ? "other" : "no_copy");
                        reconstructable = false; add(r, "unresolved", 0, 0, std::string(), ".", h); continue;
                    }
                    const long long d = cit->second - r.ref_cn;
                    const std::string cn = std::to_string(cit->second);
                    if (d == 0) continue;
                    // Copies gained or lost go where the reference's copies repeat: at c0, or a whole copy
                    // either side of it, the first of these clear of the carrier's other records there (a
                    // copy lost over a deletion the carrier also has would undo that record).
                    const std::size_t want = d > 0 ? 0 : static_cast<std::size_t>(-d) * clen;
                    struct Change { std::size_t lo, hi; bool ins; };
                    std::vector<Change> others;   // genomic [lo, hi] the carrier changes; an insertion goes in before lo
                    for (const VcfRecordInfo& o : vcf) {
                        if (&o == &r || !change_inside(o, r.pos + 1, aend)) continue;
                        const auto oa = o.sample_allele.find(h);
                        if (oa != o.sample_allele.end() && oa->second > 0)
                            others.push_back({o.pos + 1, std::max(o.end, o.pos + 1), o.end <= o.pos});
                    }
                    const auto clear = [&](std::size_t s) {
                        const std::size_t e = want ? s + want - 1 : s;
                        for (const Change& o : others)
                            if (want ? (o.ins ? s < o.lo && o.lo <= e : o.lo <= e && s <= o.hi) : (o.lo < s && s <= o.hi)) return false;
                        return true;
                    };
                    // A window must lie inside the copies; of those, the nearest clear one, else the nearest
                    // that fits.
                    const auto fits = [&](std::size_t s) { return s >= r.pos + 1 && s + (want ? want - 1 : 0) <= aend; };
                    // with one reference copy there is no other copy to step to
                    const std::size_t steps = r.ref_cn == 1 ? 0 : (aend + 1 - std::min(aend + 1, r.pos + 1)) / clen + 1;
                    std::size_t at = c0;
                    bool found = false;
                    for (int pass = 0; pass < 2 && !found; ++pass)
                        for (std::size_t step = 0; step <= steps && !found; ++step)
                            for (const int sign : {1, -1}) {
                                if (step == 0 && sign < 0) continue;
                                if (sign < 0 && c0 < step * clen) continue;
                                const std::size_t s0c = sign > 0 ? c0 + step * clen : c0 - step * clen;
                                if (fits(s0c) && (pass == 1 || clear(s0c))) { at = s0c; found = true; break; }
                            }
                    // copies lost that fit at no copy's start (as many as the reference holds): the copies' end
                    if (!found && want) at = aend + 1 >= r.pos + 1 + want ? aend + 1 - want : r.pos + 1;
                    if (d > 0) {
                        // the copy beside the insertion point, so the gain is a tandem copy wherever it goes
                        const std::size_t s0 = at + clen - 1 <= aend ? at : at - clen;
                        const std::string copy = ref_at(s0, s0 + clen - 1);
                        // A gained copy can be a shorter form than the reference's. FORMAT:CNBP says how much
                        // shorter (its change, less the carrier's other records there, less the copies pasted);
                        // where a deletion record inside the copies accounts for that, k copies short by it,
                        // those copies are pasted without it. Without such a record the VCF holds the length
                        // only, and what is left counts as copy_count.
                        std::string short_copy;
                        long long short_k = 0;
                        const auto bpit = r.sample_cnbp.find(h);
                        if (bpit != r.sample_cnbp.end() && info_field(r.info, "CN_METHOD") == "MODULE") {
                            long long own = bpit->second;
                            // deletions inside the copies (a replacement too, by its net loss): the reference they
                            // replace, genomic [lo, lo + span), and what they put there
                            struct Shorter { std::size_t lo = 0, span = 0; std::string with; long long net = 0; };
                            std::vector<Shorter> dels;
                            for (const VcfRecordInfo& o : vcf) {
                                if (&o == &r || !change_inside(o, r.pos + 1, aend)) continue;
                                const bool symbolic = !o.alts.empty() && !o.alts[0].empty() && o.alts[0][0] == '<';
                                if (o.svtype == "DEL" && !o.alts.empty()) {
                                    Shorter sh;
                                    if (symbolic) {
                                        sh.lo = o.pos + 1; sh.span = std::max(o.end, o.pos) - o.pos;
                                        if (has_flag(o.info, "SUBSTITUTION")) sh.with = info_field(o.info, "INSSEQ");
                                    } else if (o.ref.size() > o.alts[0].size()) {
                                        sh.lo = o.pos + 1; sh.span = o.ref.size() - 1; sh.with = o.alts[0].substr(1);
                                    }
                                    sh.net = static_cast<long long>(sh.span) - static_cast<long long>(sh.with.size());
                                    if (sh.net >= static_cast<long long>(min_sv_bp) && sh.span < clen && sh.lo >= r.pos + 1 && sh.lo + sh.span - 1 <= aend &&
                                        !(has_flag(o.info, "SUBSTITUTION") && sh.with.empty()))
                                        dels.push_back(sh);
                                }
                                const auto oa = o.sample_allele.find(h);   // the carrier's own other records come off CNBP
                                if (oa == o.sample_allele.end() || oa->second <= 0 || static_cast<std::size_t>(oa->second) > o.alts.size()) continue;
                                const std::string& alt = o.alts[static_cast<std::size_t>(oa->second) - 1];
                                if (alt.empty() || alt[0] != '<')
                                    own -= change_within(o, static_cast<long long>(alt.size()) - static_cast<long long>(o.ref.size()), r.pos + 1, aend);
                                else if (alt.rfind("<CNV", 0) != 0 && alt != "<INV>" && o.has_svlen) own -= change_within(o, o.svlen, r.pos + 1, aend);
                            }
                            const long long resid = own - d * static_cast<long long>(clen);
                            // Where a deletion sits in the pasted copy: a whole number of copies away when the copy
                            // is a period of the reference's copies, else inside the copy itself; npos if neither.
                            const bool periodic = copy_cache.at(r.id).period;
                            const auto offset = [&](const Shorter& del) -> std::size_t {
                                const std::size_t off = periodic ? (del.lo + clen * (s0 / clen + 2) - s0) % clen
                                                                 : (del.lo >= s0 ? del.lo - s0 : std::string::npos);
                                return off != std::string::npos && off + del.span <= clen ? off : std::string::npos;
                            };
                            // Only a shortfall beyond the form's own tolerance (lengths within 2%, 100 bp at least,
                            // are one form); only a deletion that brings the pasted copies nearer CNBP, within
                            // min_sv plus 1% of a copy per copy; and only when every such deletion is one stretch of the
                            // copy, within min_sv in place and length (two stretches: the VCF does not say which, and
                            // none is applied; a smaller shift leaves less than min_sv).
                            const long long form_tol = std::max(static_cast<long long>(min_sv_bp),
                                                                d * std::max(100LL, static_cast<long long>(clen) / 50));
                            long long best_diff = -1;
                            // every deletion that would explain it, to tell whether they are one stretch
                            std::vector<std::pair<std::size_t, std::size_t>> fits_off_span;
                            const auto apart = [&](std::size_t a, std::size_t b, bool circular) {
                                std::size_t g = a > b ? a - b : b - a;
                                if (circular) g = std::min(g, clen - g);
                                return g >= min_sv_bp;
                            };
                            if (resid < -form_tol)
                                for (const Shorter& del : dels) {
                                    const std::size_t off = offset(del);
                                    if (off == std::string::npos) continue;
                                    const std::string made = copy.substr(0, off) + del.with + copy.substr(off + del.span);
                                    for (long long k = 1; k <= d; ++k) {
                                        const long long diff = std::llabs(resid + k * del.net);
                                        const long long tol = static_cast<long long>(min_sv_bp) + k * static_cast<long long>(clen) / 100;
                                        if (diff > tol || diff >= std::llabs(resid)) continue;
                                        fits_off_span.push_back({off, del.span});
                                        if (best_diff < 0 || diff < best_diff) { best_diff = diff; short_copy = made; short_k = k; }
                                    }
                                }
                            bool ambiguous = false;
                            for (std::size_t i = 0; i < fits_off_span.size() && !ambiguous; ++i)
                                for (std::size_t j = i + 1; j < fits_off_span.size() && !ambiguous; ++j)
                                    ambiguous = apart(fits_off_span[i].first, fits_off_span[j].first, periodic) ||
                                                apart(fits_off_span[i].second, fits_off_span[j].second, false);
                            if (ambiguous || best_diff < 0) { short_k = 0; short_copy.clear(); }
                        }
                        std::string ins;
                        for (long long k = 0; k < d; ++k) ins += k < short_k ? short_copy : copy;
                        add(r, "cnv", at, at - 1, ins, cn, h);
                    } else {
                        const std::size_t del = std::min<std::size_t>(want, aend + 1 - at);
                        add(r, "cnv", at, at + del - 1, std::string(), cn, h);
                    }
                } else { unspelled.emplace(r.id, "other"); reconstructable = false; add(r, "unresolved", 0, 0, std::string(), ".", h); }
            }
            if (!reconstructable) ++vcf_not_reconstructable;
        }
    } else cevents = load_carrier_events(ce_in);
    // The two files must describe one run: every carrier of every record has exactly one event.
    if (!vcf_mode) {
        std::unordered_map<std::string, std::unordered_set<std::string>> ce_carriers;
        for (const CarrierEvent& e : cevents) {
            if (!vcf_by_id.count(e.variant_id))
                throw std::runtime_error(ce_in + " names record '" + e.variant_id + "', which " + vcf_in +
                                         " does not have: the two files are from different runs");
            for (const std::string& h : e.haplotypes) ce_carriers[e.variant_id].insert(h);   // can be several events
        }
        for (const VcfRecordInfo& r : vcf) {
            const auto it = ce_carriers.find(r.id);
            const std::size_t n = it == ce_carriers.end() ? 0 : it->second.size();
            bool same = n == r.carriers.size();
            for (const std::string& h : r.carriers) if (same && !it->second.count(h)) same = false;
            if (!same)
                throw std::runtime_error("record '" + r.id + "' has " + std::to_string(r.carriers.size()) +
                                         " carrier(s) in " + vcf_in + " and " + std::to_string(n) + " in " + ce_in +
                                         ": the two files are from different runs");
        }
    }
    // Each record's reference span, for attributing what is left: its POS..END, widened to wherever any
    // carrier's own event sits.
    std::map<std::string, std::pair<std::size_t, std::size_t>> record_span, event_span;
    for (const VcfRecordInfo& r : vcf) record_span[r.id] = {r.pos, std::max(r.pos, r.end)};
    std::unordered_map<std::string, std::vector<const CarrierEvent*>> events_of;   // haplotype -> events
    std::size_t unresolved_rows = 0;
    for (const CarrierEvent& e : cevents) {
        if (e.kind == "unresolved") { unresolved_rows += e.haplotypes.size(); continue; }
        std::size_t lo = e.insertion() ? e.ref_start - 1 : e.ref_start;
        std::size_t hi = e.insertion() ? e.ref_start : e.ref_end;
        // A --vcf copy gain or loss sits at one copy, but stands for the record's whole run of copies.
        if (e.kind == "cnv") {
            const VcfRecordInfo& r = *vcf_by_id.at(e.variant_id);
            lo = std::min(lo, r.pos + 1);
            hi = std::max(hi, r.pos + (r.cn_ref_bp > 0 ? static_cast<std::size_t>(r.cn_ref_bp) : std::max(r.end, r.pos) - r.pos));
        }
        auto& sp = record_span[e.variant_id];
        sp.first = std::min(sp.first, lo);
        sp.second = std::max(sp.second, hi);
        const auto [es, fresh] = event_span.try_emplace(e.variant_id, lo, hi);
        if (!fresh) { es->second.first = std::min(es->second.first, lo); es->second.second = std::max(es->second.second, hi); }
        for (const std::string& h : e.haplotypes) events_of[h].push_back(&e);
    }
    std::unordered_map<std::string, std::unordered_set<std::string>> unresolved_of;
    for (const CarrierEvent& e : cevents)
        if (e.kind == "unresolved") for (const std::string& h : e.haplotypes) unresolved_of[h].insert(e.variant_id);

    // A carrier's copy-number event is its own stretch: its copies, and for an allele call compared node by
    // node the sequence beside them up to the steps its walk and the reference share, so its other events
    // inside that stretch -- an insertion between or beside its copies, a deletion of a copy's part -- are
    // already in it. An insertion is inside when it goes in after the stretch's first base and before the
    // base after its last; one at either edge may be beside the copies or in the stretch.
    const auto inside_at = [](const CarrierEvent& c, const CarrierEvent& e) {
        if (&c == &e || e.kind == "cn") return false;
        if (c.insertion()) return e.insertion() && e.ref_start == c.ref_start;
        if (e.insertion()) return c.ref_start <= e.ref_start && e.ref_start <= c.ref_end + 1;
        return c.ref_start <= e.ref_start && e.ref_end <= c.ref_end;
    };
    // An insertion that repeats the sequence beside it can go in anywhere over its record's CIPOS. When that
    // range reaches a stretch's edge, its listed place does not say which side it is on: it is inside when
    // counting it brings the stretch's expected length (as the carrier check computes it) nearer its own.
    const auto place_range = [&](const CarrierEvent& e) {
        std::size_t lo = e.ref_start, hi = e.ref_start;
        const VcfRecordInfo& r = *vcf_by_id.at(e.variant_id);
        const std::vector<std::string> ci = split_on(info_field(r.info, "CIPOS"), ',');
        if (ci.size() == 2) {
            char* e0 = nullptr; char* e1 = nullptr;
            const long long a = std::strtoll(ci[0].c_str(), &e0, 10), b = std::strtoll(ci[1].c_str(), &e1, 10);
            const long long p = static_cast<long long>(r.pos) + 1;
            if (*e0 == '\0' && *e1 == '\0' && p + a > 0 && a <= b) {
                lo = std::min(lo, static_cast<std::size_t>(p + a));
                hi = std::max(hi, static_cast<std::size_t>(p + b));
            }
        }
        return std::make_pair(lo, hi);
    };
    std::map<std::pair<const CarrierEvent*, std::string>, std::unordered_set<const CarrierEvent*>> inside_of;
    for (const auto& [h, evs] : events_of) {
        for (const CarrierEvent* c : evs) {
            if (c->kind != "cn") continue;
            auto& in = inside_of[{c, h}];
            std::vector<const CarrierEvent*> edge;
            for (const CarrierEvent* e : evs) {
                if (e != c && e->kind != "cn" && e->insertion() && !c->insertion()) {
                    const auto [lo, hi] = place_range(*e);
                    if (lo <= c->ref_end + 1 && hi >= c->ref_start && (lo <= c->ref_start || hi >= c->ref_end + 1)) {
                        edge.push_back(e);
                        continue;
                    }
                }
                if (inside_at(*c, *e)) in.insert(e);
            }
            if (edge.empty()) continue;
            const VcfRecordInfo& r = *vcf_by_id.at(c->variant_id);
            const auto cit = r.sample_cn.find(h);
            const long long copy_bp = static_cast<long long>(copy_of(r, h).second);
            if (cit == r.sample_cn.end() || r.ref_cn < 0 || copy_bp <= 0) {
                for (const CarrierEvent* e : edge) if (inside_at(*c, *e)) in.insert(e);
                continue;
            }
            long long expect = static_cast<long long>(c->ref_len()) + (cit->second - r.ref_cn) * copy_bp;
            for (const CarrierEvent* e : in) expect += static_cast<long long>(e->alt.size()) - static_cast<long long>(e->ref_len());
            const long long got = static_cast<long long>(c->alt.size());
            for (const CarrierEvent* e : edge) {
                const long long with = expect + static_cast<long long>(e->alt.size());
                if (std::llabs(got - with) < std::llabs(got - expect)) { in.insert(e); expect = with; }
            }
        }
    }
    // Inversions whose alt is the plain reference interval turned round.
    std::unordered_set<const CarrierEvent*> plain_inv;
    for (const CarrierEvent& e : cevents)
        if (e.kind == "inv" && !e.insertion() && e.ref_start >= region_start && e.ref_end < region_start + ref_seq.size() &&
            e.alt == reverse_complement(ref_seq.substr(e.ref_start - region_start, e.ref_len())))
            plain_inv.insert(&e);
    const auto inside_cn = [&](const CarrierEvent& c, const CarrierEvent& e, const std::string& h) {
        const auto it = inside_of.find({&c, h});
        return it != inside_of.end() && it->second.count(&e) > 0;
    };

    // ---- the carrier check: each carrier's own event against what the VCF reports for it ----
    std::vector<CheckRow> checks;
    std::size_t checked = 0, unchecked = 0, cn_stretch_checked = 0, cn_stretch_unchecked = 0, cn_exact_checked = 0;
    // A carrier can hold several events under one record (two pieces call filed together): their net change
    // together is what the record's SVLEN stands for, so they are checked once, summed.
    std::map<std::pair<std::string, std::string>, std::pair<long long, std::size_t>> seq_net;   // (record, carrier)
    if (!vcf_mode)
        for (const CarrierEvent& e : cevents)
            if (e.kind != "unresolved" && e.kind != "cn" && e.kind != "inv")
                for (const std::string& h : e.haplotypes) {
                    auto& v = seq_net[{e.variant_id, h}];
                    v.first += e.event_bp; ++v.second;
                }
    std::set<std::pair<std::string, std::string>> seq_done;
    for (const CarrierEvent& e : cevents) {
        if (vcf_mode) break;   // the events ARE the VCF's
        if (e.kind == "unresolved") continue;
        const VcfRecordInfo& r = *vcf_by_id.at(e.variant_id);
        const long long own_net = e.event_bp;
        for (const std::string& h : e.haplotypes) {
            long long reported = 0, own = own_net, limit = static_cast<long long>(min_sv_bp);
            std::string rep_s, own_s;
            if (e.kind == "cn") {
                // The CN reported against the copies the carrier's own walk holds.
                const auto cit = r.sample_cn.find(h);
                if (cit == r.sample_cn.end() || e.walk_cn == ".") { ++unchecked; continue; }
                reported = cit->second;
                own = static_cast<long long>(parse_count(e.walk_cn, ce_in));
                limit = 1;
                rep_s = "CN=" + std::to_string(cit->second);
                own_s = std::to_string(own) + " copies in its walk";
            } else if (e.kind == "inv" || r.svtype == "INV") {
                if (!r.has_svlen) { ++unchecked; continue; }
                reported = r.svlen < 0 ? -r.svlen : r.svlen;
                own = e.event_bp;
                rep_s = "SVLEN=" + std::to_string(r.svlen);
                own_s = std::to_string(own) + " bp inverted";
            } else {
                if (!seq_done.insert({e.variant_id, h}).second) continue;   // checked with its first event
                if (!r.has_svlen) { ++unchecked; continue; }
                const auto& net = seq_net.at({e.variant_id, h});
                own = net.first;
                reported = r.svlen;
                rep_s = "SVLEN=" + std::to_string(r.svlen);
                own_s = std::to_string(own) + " bp" + (net.second > 1 ? " over " + std::to_string(net.second) + " events" : "");
            }
            ++checked;
            const long long diff = own - reported;
            if ((diff < 0 ? -diff : diff) >= limit) checks.push_back({e.variant_id, h, e.kind, rep_s, own_s, diff});
            if (e.kind != "cn") continue;
            const auto eit = events_of.find(h);
            // A record naming each carrier's copies (CN_FORMS): the reference's stretch with the carrier's copies
            // and its own events inside the stretch applied must make the stretch its walk holds, but for changes
            // under --min-sv-bp, which call does not report: flagged when they overlap, or when aligned the two
            // differ by --min-sv-bp or more somewhere.
            if (!r.cn_forms.empty()) {
                const auto fit = r.sample_cnform.find(h);
                if (fit == r.sample_cnform.end() || e.ref_start < region_start || e.ref_end >= region_start + ref_seq.size()) continue;
                ++cn_exact_checked;
                struct Edit { std::size_t g0, g1; const std::string* alt; };   // [g0, g1] replaced; g1 = g0 - 1 inserts
                std::vector<Edit> edits;
                bool spelled = true;
                for (const std::size_t k : fit->second) {
                    const VcfRecordInfo::Form& fm = r.cn_forms[k];
                    spelled = spelled && fm.spelled;
                    edits.push_back({fm.start, fm.end, &fm.seq});
                }
                if (eit != events_of.end())
                    for (const CarrierEvent* o : eit->second)
                        if (inside_cn(e, *o, h)) edits.push_back({o->ref_start, o->ref_end, &o->alt});
                std::stable_sort(edits.begin(), edits.end(), [](const Edit& a, const Edit& b) {
                    if (a.g0 != b.g0) return a.g0 < b.g0;
                    return a.g1 < b.g1;   // an insertion before a change starting at its base
                });
                const auto upper = [](std::string x) {
                    for (char& c : x) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
                    return x;
                };
                // Two insertions at one point go in either order (the VCF does not say which comes first): the
                // stretch is made both ways, and the one nearer the walk's counts.
                std::string made;
                bool clash = false, large = false;
                for (int pass = 0; pass < 2; ++pass) {
                    if (pass == 1) {
                        bool shared = false;
                        for (std::size_t i = 0; i + 1 < edits.size(); ++i)
                            if (edits[i].g0 == edits[i + 1].g0 && edits[i].g1 + 1 == edits[i].g0 && edits[i + 1].g1 + 1 == edits[i + 1].g0) {
                                std::swap(edits[i], edits[i + 1]); shared = true; ++i;
                            }
                        if (!shared) break;
                    }
                    made.clear();
                    std::size_t at = e.ref_start;
                    clash = false;
                    for (const Edit& x : edits) {
                        if (x.g0 < at || x.g0 < e.ref_start || x.g1 > e.ref_end) { clash = true; break; }
                        made += ref_seq.substr(at - region_start, x.g0 - at);
                        made += *x.alt;
                        at = x.g1 + 1;
                    }
                    if (!clash) made += ref_seq.substr(at - region_start, e.ref_end + 1 - at);
                    large = false;
                    if (!clash && upper(made) != upper(e.alt))
                        for (const Diff& d : diff_of(align_ops(upper(made), upper(e.alt)), merge_gap, min_sv_bp).events)
                            large = large || d.size() >= min_sv_bp;
                    if (!clash && !large) break;
                }
                const long long gap = static_cast<long long>(e.alt.size()) - static_cast<long long>(made.size());
                if (spelled && (clash || large))
                    checks.push_back({e.variant_id, h, "cn_exact",
                                      clash ? std::string("copies and events that overlap") :
                                              std::to_string(made.size()) + " bp from its copies and events",
                                      std::to_string(e.alt.size()) + " bp stretch", gap});
                continue;
            }
            // And the stretch's length against the CN: the reference's stretch, (CN - REF_CN) copies, and
            // the carrier's own events inside it. Sequence the stretch holds beyond that is in no record.
            const auto cit = r.sample_cn.find(h);
            const long long copy_bp = static_cast<long long>(copy_of(r, h).second);
            if (cit == r.sample_cn.end() || r.ref_cn < 0 || copy_bp <= 0) { ++cn_stretch_unchecked; continue; }
            ++cn_stretch_checked;
            long long expect = static_cast<long long>(e.ref_len()) + (cit->second - r.ref_cn) * copy_bp, inner = 0;
            if (eit != events_of.end())
                for (const CarrierEvent* o : eit->second)
                    if (inside_cn(e, *o, h)) inner += static_cast<long long>(o->alt.size()) - static_cast<long long>(o->ref_len());
            expect += inner;
            const long long got = static_cast<long long>(e.alt.size());
            // Copies of one record differ by small indels that add up over many copies: each copy changed
            // may differ by 1% of a copy, on top of --min-sv-bp.
            const long long dcn = cit->second > r.ref_cn ? cit->second - r.ref_cn : r.ref_cn - cit->second;
            const long long off = got - expect < 0 ? expect - got : got - expect;
            if (100 * off >= 100 * static_cast<long long>(min_sv_bp) + dcn * copy_bp)
                checks.push_back({e.variant_id, h, "cn_stretch",
                                  std::to_string(e.ref_len()) + " bp + (" + std::to_string(cit->second) + "-" +
                                      std::to_string(r.ref_cn) + ") x " + std::to_string(copy_bp) + " bp copies + " +
                                      std::to_string(inner) + " bp of its events = " + std::to_string(expect) + " bp",
                                  std::to_string(got) + " bp stretch", got - expect});
        }
    }

    // ---- each haplotype: apply its events, align, and read what is left ----
    std::vector<const PathRecord*> haps;
    std::size_t not_in_vcf = 0;
    {
        std::unordered_set<std::string> in_vcf(vcf_samples.begin(), vcf_samples.end());
        for (const PathRecord& p : graph.paths) {
            if (&p == ref_path) continue;
            if (!in_vcf.count(p.name)) { ++not_in_vcf; continue; }
            haps.push_back(&p);
        }
    }

    if (!vcf_mode) {
        std::size_t unanchored = 0;   // a haplotype sharing no anchor with the reference is not scored
        for (const PathRecord* p : haps)
            unanchored += preps[static_cast<std::size_t>(p - graph.paths.data())].status != "scored";
        if (unanchored > 0)
            std::cerr << cli::log_prefix("benchmark") << "warning: " << unanchored
                      << " haplotype(s) share no anchor node with the reference and are not scored\n";
        log.info("haplotypes scored: " + std::to_string(haps.size() - unanchored) + "; large events: " +
                 std::to_string(min_sv_bp) + " bp or more");
        log.info("call's events: records: " + std::to_string(vcf.size()) + "; carrier events: " +
                 std::to_string(cevents.size()));
    } else {
        log.info("VCF: records: " + std::to_string(vcf.size()) + "; events rebuilt from them: " +
                 std::to_string(cevents.size()));
    }

    std::unordered_map<std::string, const HapResult*> call_stage;   // the VCF stage reads the call stage's
    if (vcf_mode) for (const HapResult& r : stages[0].results) call_stage[r.name] = &r;
    std::vector<HapResult> results(haps.size());
    run_parallel(haps.size(), threads, [&](std::size_t hi) {
        const PathRecord& path = *haps[hi];
        HapResult& hr = results[hi];
        hr.name = path.name;
        const std::set<std::string>* call_outside = nullptr;
        if (const auto it = call_stage.find(path.name); it != call_stage.end()) call_outside = &it->second->outside_ids;
        const Prep& pr = preps[static_cast<std::size_t>(&path - graph.paths.data())];
        if (pr.status != "scored") { hr.status = pr.status; return; }
        hr.g_start = pr.g_start;
        hr.g_end = pr.g_end;
        const std::size_t o0 = pr.o0, o1 = pr.o1;
        const std::string& truth = pr.truth;

        // The events, in reference order; an event outside the extent, or overlapping one already
        // applied, is left out and counted.
        std::vector<const CarrierEvent*> evs;
        if (const auto it = events_of.find(path.name); it != events_of.end()) evs = it->second;
        if (const auto it = unresolved_of.find(path.name); it != unresolved_of.end()) hr.unresolved = it->second.size();
        // Events outside the scored stretch leave first, so they neither absorb nor outweigh the events
        // inside it. The VCF stage leaves out what the call stage left out for this haplotype, so an edge
        // event is in or out the same way in both.
        {
            std::vector<const CarrierEvent*> in;
            for (const CarrierEvent* e : evs) {
                const std::size_t lo = e->ref_start, hi_end = e->ref_end + 1;
                const bool outside = e->insertion() ? (lo <= hr.g_start || lo > hr.g_end) : (lo < hr.g_start || hi_end > hr.g_end + 1);
                const bool out_in_call = vcf_mode && call_outside != nullptr && call_outside->count(e->variant_id);
                if (outside || out_in_call) { ++hr.outside_extent; hr.outside_ids.insert(e->variant_id); continue; }
                in.push_back(e);
            }
            evs.swap(in);
        }
        const auto by_order = [vcf_mode](const CarrierEvent* a, const CarrierEvent* b) {
            if (a->ref_start != b->ref_start) return a->ref_start < b->ref_start;
            if (a->insertion() != b->insertion()) return a->insertion();   // goes in before ref_start
            if (a->ref_end != b->ref_end) return a->ref_end > b->ref_end;  // the enclosing one first
            if (vcf_mode) return a->order < b->order;   // two records at one point: the VCF's order
            return a->alt.size() > b->alt.size();   // one interval: the replacement holding more of the haplotype
        };
        std::sort(evs.begin(), evs.end(), by_order);
        // Its events inside its own copy-number stretch are part of that stretch (see inside_cn).
        {
            std::vector<const CarrierEvent*> cns;
            for (const CarrierEvent* e : evs) if (e->kind == "cn") cns.push_back(e);
            if (!cns.empty()) {
                const auto n0 = evs.size();
                evs.erase(std::remove_if(evs.begin(), evs.end(), [&](const CarrierEvent* e) {
                    for (const CarrierEvent* c : cns) if (inside_cn(*c, *e, path.name)) return true;
                    return false;
                }), evs.end());
                hr.absorbed = n0 - evs.size();
            }
        }
        // Other overlapping events of one haplotype: the one with the longer reference interval is kept,
        // and of two over one interval the first in the order above. Two insertions at one point, or two
        // crossing intervals, are conflicts too; each left out is counted. An event inside an inversion
        // that turns the plain reference is not a conflict: the sites inside an inversion keep their own
        // events, so they are applied to the reference stretch first, and the inversion turns the result.
        std::unordered_map<const CarrierEvent*, std::vector<const CarrierEvent*>> nested;
        // k starts no later than e. An insertion goes in before its ref_start, so it is inside an interval
        // [S, E] only when S < ref_start <= E, and two insertions clash at one point (in --vcf mode two
        // records inserting at one point are both applied, in the file's order).
        const auto overlaps = [&](const CarrierEvent* k, const CarrierEvent* e) {
            if (e->insertion() && k->insertion()) return !vcf_mode && e->ref_start == k->ref_start;
            if (e->insertion()) return k->ref_start < e->ref_start && e->ref_start <= k->ref_end;
            if (k->insertion()) return false;
            return e->ref_start <= k->ref_end;
        };
        // Longest first, each event is kept unless it overlaps one already kept, so of two overlapping events
        // the longer is kept whatever lies around them. `nest` sets events inside a kept plain inversion
        // aside in `nested` instead of weighing them against it.
        std::unordered_map<const CarrierEvent*, std::size_t> rank;   // place in the order above
        for (std::size_t i = 0; i < evs.size(); ++i) rank[evs[i]] = i;
        const auto before = [&](const CarrierEvent* a, const CarrierEvent* b) { return rank.at(a) < rank.at(b); };
        const auto keep_apart = [&](std::vector<const CarrierEvent*> list, bool nest) {
            std::sort(list.begin(), list.end(), [&](const CarrierEvent* a, const CarrierEvent* b) {
                return a->ref_len() != b->ref_len() ? a->ref_len() > b->ref_len() : before(a, b);
            });
            std::vector<const CarrierEvent*> kept;   // apart, in the order above
            for (const CarrierEvent* e : list) {
                const auto at = std::upper_bound(kept.begin(), kept.end(), e, before);
                // Kept events are apart, so only the one before e and those starting within e can overlap it.
                const CarrierEvent* hit = nullptr;
                if (at != kept.begin() && overlaps(*(at - 1), e)) hit = *(at - 1);
                for (auto it = at; !hit && it != kept.end() && (*it)->ref_start <= std::max(e->ref_start, e->ref_end); ++it)
                    if (overlaps(e, *it)) hit = *it;
                if (!hit) kept.insert(at, e);
                else if (nest && plain_inv.count(hit) && before(hit, e) && e->kind != "inv" &&
                         (e->insertion() || e->ref_end <= hit->ref_end))
                    nested[hit].push_back(e);
                else ++hr.conflicts;
            }
            return kept;
        };
        evs = keep_apart(evs, true);
        // inside an inversion, the same rule among its own events
        for (auto& [inv, inner] : nested) inner = keep_apart(inner, false);
        std::string recon;
        recon.reserve(o1 - o0);
        std::vector<Piece> pieces;
        std::vector<const CarrierEvent*> applied;
        std::size_t cursor = hr.g_start;          // next genomic base to copy
        bool last_was_insertion_here = false;
        std::size_t last_insertion_at = 0;
        auto copy_ref = [&](std::size_t g_to) {   // copy [cursor, g_to)
            if (g_to <= cursor) return;
            Piece p; p.r0 = recon.size(); p.g0 = cursor; p.g1 = g_to - 1;
            recon.append(ref_seq, cursor - region_start, g_to - cursor);
            p.r1 = recon.size();
            pieces.push_back(p);
            cursor = g_to;
        };
        for (const CarrierEvent* e : evs) {
            const std::size_t lo = e->ref_start, hi_end = e->ref_end + 1;   // genomic [lo, hi_end)
            if (lo < cursor || (!vcf_mode && e->insertion() && last_was_insertion_here && last_insertion_at == lo)) { ++hr.conflicts; continue; }   // a guard: the events kept above are apart
            copy_ref(lo);
            Piece p; p.r0 = recon.size(); p.g0 = lo; p.g1 = e->ref_end; p.event = static_cast<int>(applied.size());
            applied.push_back(e);
            if (const auto it = nested.find(e); it != nested.end()) {
                std::string inner;
                std::size_t c = lo;
                bool ins_here = false;
                for (const CarrierEvent* n : it->second) {
                    if (n->ref_start < c || (!vcf_mode && n->insertion() && ins_here && n->ref_start == c)) { ++hr.conflicts; continue; }
                    inner.append(ref_seq, c - region_start, n->ref_start - c);
                    inner += n->alt;
                    c = n->ref_end + 1;
                    ins_here = n->insertion();
                    applied.push_back(n);
                }
                inner.append(ref_seq, c - region_start, hi_end - c);
                recon += reverse_complement(inner);
            } else recon += e->alt;
            p.r1 = recon.size();
            pieces.push_back(p);
            cursor = hi_end;
            last_was_insertion_here = e->insertion();
            last_insertion_at = lo;
        }
        copy_ref(hr.g_end + 1);
        hr.applied = applied.size();
        for (const CarrierEvent* e : applied) {
            hr.applied_spans.push_back({e->variant_id, e->insertion() ? e->ref_start - 1 : e->ref_start,
                                        e->insertion() ? e->ref_start : e->ref_end});
        }

        // reconstruction position -> (genomic position, applied event or -1)
        auto locate = [&](std::size_t r) -> std::pair<std::size_t, int> {
            auto it = std::upper_bound(pieces.begin(), pieces.end(), r,
                                       [](std::size_t v, const Piece& p) { return v < p.r0; });
            if (it == pieces.begin()) return {hr.g_start, -1};
            --it;
            if (it->event >= 0) {   // the carrier's own bases: placed in proportion over the reference it replaces
                if (it->g1 < it->g0 || it->r1 <= it->r0) return {it->g0, it->event};
                const std::size_t off = std::min((r - it->r0) * (it->g1 - it->g0 + 1) / (it->r1 - it->r0), it->g1 - it->g0);
                // an inversion's bases run backwards over the reference
                const bool turned = applied[static_cast<std::size_t>(it->event)]->kind == "inv";
                return {turned ? it->g1 - off : it->g0 + off, it->event};
            }
            const std::size_t g = it->g0 + std::min(r - it->r0, it->r1 > it->r0 ? it->r1 - it->r0 : 0);
            return {g, -1};
        };

        const DiffSet& base = pr.base;
        hr.present_snv = base.mismatches;
        const std::size_t stretch = pr.o1 - pr.o0;
        const std::string ref_stretch = ref_seq.substr(pr.o0, stretch);
        for (const Diff& d : base.events) {
            if (d.size() < min_sv_bp) { hr.present_small_bp += d.size(); continue; }
            ++hr.present_n; hr.present_bp += d.size();
            // Placed as the leftovers are: from the base before the event to the base after it.
            HapResult::Present p;
            p.g0 = pr.g_start + (d.t0 > 0 ? d.t0 - 1 : 0);
            p.g1 = pr.g_start + std::min(d.t1, stretch ? stretch - 1 : 0);
            p.size = d.size(); p.del_bp = d.del_bp; p.ins_bp = d.ins_bp; p.type = d.type();
            // the sites it lies in, by the leftovers' rule, anywhere it can slide to
            std::size_t w0 = p.g0, w1 = p.g1, s0 = 0, s1 = 0;
            if (slide_of(ref_stretch, pr.truth, d, s0, s1)) {
                w0 = std::min(w0, pr.g_start + (s0 > 0 ? s0 - 1 : 0));
                w1 = std::max(w1, pr.g_start + std::min(s1, stretch ? stretch - 1 : 0));
            }
            const std::size_t lo = w0 > merge_gap ? w0 - merge_gap : 0, hi = w1 + merge_gap;
            for (const auto& [bid, sp] : bubble_span)
                if (sp.first <= hi && lo <= sp.second) p.bubbles = (p.bubbles == "." ? "" : p.bubbles + ";") + std::to_string(bid);
            hr.present_rows.push_back(std::move(p));
        }
        // With nothing applied the reconstruction is the reference, already aligned.
        const DiffSet left = applied.empty() ? base : diff_of(align_ops(recon, truth), merge_gap, min_sv_bp);
        hr.left_snv = left.mismatches;
        std::unordered_set<std::string> carried;   // every record the haplotype carries, applied or not
        if (const auto it = events_of.find(path.name); it != events_of.end())
            for (const CarrierEvent* e : it->second) carried.insert(e->variant_id);
        for (const Diff& d : left.events) {
            if (d.size() < min_sv_bp) { hr.left_small_bp += d.size(); continue; }
            ++hr.left_n;
            hr.left_bp += d.size();
            LargeRow row;
            row.haplotype = path.name;
            row.type = d.type();
            row.size = d.size(); row.del_bp = d.del_bp; row.ins_bp = d.ins_bp;
            // where it sits, and whether it touches an event applied to this haplotype
            int touched = -1;
            std::size_t g_lo = SIZE_MAX, g_hi = 0, m_lo = SIZE_MAX, m_hi = 0;
            const std::size_t r_hi = d.t1 > d.t0 ? d.t1 - 1 : d.t0;
            // its first and last bases, and the bases before and after it (an insertion has only these)
            for (std::size_t r : {d.t0, r_hi, d.t0 > 0 ? d.t0 - 1 : d.t0, std::min(d.t1, recon.size() ? recon.size() - 1 : 0)}) {
                const auto [g, ev] = locate(r);
                const bool own = d.t1 > d.t0 && (r == d.t0 || r == r_hi);   // its own bases, not its flanks
                if (ev >= 0 && touched < 0 && (own || (d.t1 == d.t0 && r == d.t0))) touched = ev;
                g_lo = std::min(g_lo, g); g_hi = std::max(g_hi, g);
                if (ev >= 0 && own) {
                    const CarrierEvent& ae = *applied[static_cast<std::size_t>(ev)];
                    m_lo = std::min(m_lo, ae.ref_start); m_hi = std::max(m_hi, std::max(ae.ref_end, ae.ref_start));
                }
            }
            // An insertion has no bases of its own: it lies in an applied event when both its flanks do.
            if (d.t1 == d.t0 && d.t0 > 0) {
                const int ea = locate(d.t0 - 1).second, eb = locate(std::min(d.t0, recon.size() ? recon.size() - 1 : 0)).second;
                if (ea >= 0 && ea == eb) {
                    const CarrierEvent& ae = *applied[static_cast<std::size_t>(ea)];
                    m_lo = std::min(m_lo, ae.ref_start); m_hi = std::max(m_hi, std::max(ae.ref_end, ae.ref_start));
                }
            }
            for (auto it = std::upper_bound(pieces.begin(), pieces.end(), d.t0,
                                            [](std::size_t v, const Piece& p) { return v < p.r0; });
                 touched < 0 && it != pieces.end() && it->r0 < d.t1; ++it)
                if (it->event >= 0) touched = it->event;
            row.g0 = g_lo; row.g1 = g_hi;
            row.m0 = std::min(m_lo, g_lo); row.m1 = std::max(m_hi, g_hi);
            // Sites and records are looked for anywhere it can slide to.
            std::size_t w_lo = g_lo, w_hi = g_hi;
            if (std::size_t s0 = 0, s1 = 0; slide_of(recon, truth, d, s0, s1)) {
                const std::size_t last = recon.size() ? recon.size() - 1 : 0;
                const std::size_t r_lo = s0 > 0 ? s0 - 1 : 0, r_hi = std::min(s1, last);
                const auto take = [&](std::size_t r) {
                    const std::size_t g = locate(r).first;
                    w_lo = std::min(w_lo, g); w_hi = std::max(w_hi, g);
                };
                take(r_lo);
                take(r_hi);
                // and the ends of every piece inside, as an applied inversion maps backwards
                for (auto it = std::upper_bound(pieces.begin(), pieces.end(), r_lo,
                                                [](std::size_t v, const Piece& p) { return v < p.r0; });
                     it != pieces.end() && it->r0 <= r_hi; ++it) {
                    take(it->r0);
                    if (it->r1 > it->r0) take(std::min(it->r1 - 1, r_hi));
                }
                if (auto it = std::upper_bound(pieces.begin(), pieces.end(), r_lo,
                                               [](std::size_t v, const Piece& p) { return v < p.r0; });
                    it != pieces.begin() && (--it)->r1 > it->r0)
                    take(std::min(it->r1 - 1, r_hi));
            }
            row.w0 = w_lo; row.w1 = w_hi;
            row.m0 = std::min(row.m0, w_lo); row.m1 = std::max(row.m1, w_hi);
            const std::size_t slo = w_lo > merge_gap ? w_lo - merge_gap : 0, shi = w_hi + merge_gap;
            const auto near = [&](const std::map<std::string, std::pair<std::size_t, std::size_t>>& spans) {
                std::vector<std::string> v;
                for (const auto& [id, sp] : spans) if (sp.first <= shi && slo <= sp.second) v.push_back(id);
                return v;
            };
            std::vector<std::string> recs, bids;
            for (const auto& [bid, sp] : bubble_span)
                if (sp.first <= shi && slo <= sp.second) bids.push_back(std::to_string(bid));
            if (touched >= 0) {
                const CarrierEvent& te = *applied[static_cast<std::size_t>(touched)];
                row.cause = "incomplete";
                recs = {te.variant_id};
            } else {
                // Where the records' events sit first, then their POS..END (a copy-number record's span is its
                // whole module). A record this haplotype carries there: its event is not all of what the
                // haplotype has, or it could not be written: both incomplete. Otherwise a record it does not
                // carry, or none at all: a miss.
                const std::unordered_set<std::string>* unres = nullptr;
                if (const auto it = unresolved_of.find(path.name); it != unresolved_of.end()) unres = &it->second;
                // A record of its own that could not be written, there, comes first: nothing of it was applied.
                std::vector<std::string> mine, gone, others;
                if (unres != nullptr)
                    for (const std::string& id : near(record_span)) if (unres->count(id)) gone.push_back(id);
                for (const auto* spans : {&event_span, &record_span}) {
                    if (!gone.empty()) break;
                    for (const std::string& id : near(*spans)) (carried.count(id) ? mine : others).push_back(id);
                    if (!mine.empty() || !others.empty()) break;
                }
                if (!gone.empty()) { row.cause = "incomplete"; row.unwritten = true; recs = gone; }
                else if (!mine.empty()) { row.cause = "incomplete"; recs = mine; }
                else if (!others.empty()) { row.cause = "not_carrier"; recs = others; }
                else { row.cause = "not_reported"; recs.clear(); }
            }
            auto join = [](const std::vector<std::string>& v) {
                if (v.empty()) return std::string(".");
                std::string s;
                for (std::size_t k = 0; k < v.size(); ++k) { if (k) s += ';'; s += v[k]; }
                return s;
            };
            row.records = join(recs);
            row.bubbles = join(bids);
            hr.rows.push_back(std::move(row));
        }
    });

    out.results = std::move(results);
    out.checks = std::move(checks);
    out.checked = checked; out.unchecked = unchecked;
    out.cn_checked = cn_stretch_checked; out.cn_unchecked = cn_stretch_unchecked; out.cn_exact_checked = cn_exact_checked;
    out.unresolved_rows = unresolved_rows;
    out.not_reconstructable = vcf_not_reconstructable; out.missing_gt = vcf_missing_gt;
    out.unspelled = std::move(unspelled);
    out.not_in_vcf = not_in_vcf;
    out.vcf = vcf;
    };
    // ---- outputs: what call's events leave, then what the VCF's records leave ----
    // A large difference left by the VCF's records that overlaps one left by call's own events (on the
    // reference, within --merge-gap) is the same miss, with its cause. The rest is what the step from
    // call's findings to the file lost, with a VCF reason: sequence the record does not hold
    // (no_sequence), copies it counts but cannot spell (copy_count), a carrier's own allele its
    // representative does not match (representative), or anything else (other).
    const auto write_outputs = [&]() {
        cli::StagedOutputs staged("benchmark");
        const Stage& sc = stages[0];
        const Stage& sv = stages[1];
        std::unordered_map<std::string, const HapResult*> vres;
        for (const HapResult& r : sv.results) vres[r.name] = &r;
        std::unordered_map<std::string, const VcfRecordInfo*> rec;
        for (const VcfRecordInfo& r : sv.vcf) rec[r.id] = &r;
        // Why the VCF rebuild leaves `e` where call's own events do not; `records` may be filled in.
        const auto lost_reason = [&](const LargeRow& e, const HapResult& c, std::string& records) -> std::string {
            if (e.unwritten) return "no_sequence";
            const std::vector<std::string> ids = e.records == "." ? std::vector<std::string>() : split_on(e.records, ';');
            for (const std::string& id : ids)
                if (const auto it = rec.find(id); it != rec.end() && (it->second->svtype == "CNV" || it->second->svtype == "DUP"))
                    return "copy_count";
            // approximate breakpoints: not the representative's fault
            for (const std::string& id : ids)
                if (const auto it = rec.find(id); it != rec.end() && has_flag(it->second->info, "IMPRECISE")) return "other";
            if (e.cause == "incomplete") return "representative";
            // A record whose own event for this haplotype sits there, written elsewhere by the VCF.
            for (const HapResult::Span& sp : c.applied_spans)
                if (e.w0 <= sp.hi + merge_gap && sp.lo <= e.w1 + merge_gap) { records = sp.id; return "representative"; }
            return "other";
        };
        const std::vector<std::string> causes = {"incomplete", "not_carrier", "not_reported"};
        const std::vector<std::string> vcf_reasons = {"no_sequence", "copy_count", "representative", "other"};
        std::map<std::string, std::size_t> call_by_cause, vcf_by_reason;
        for (const std::string& c : causes) call_by_cause[c] = 0;
        for (const std::string& r : vcf_reasons) vcf_by_reason[r] = 0;
        std::map<std::string, std::size_t> unspelled_lost;   // no_sequence rows, by why the record is unspelled
        std::map<std::string, std::size_t> skipped;
        std::size_t scored = 0, present_n = 0, present_bp = 0, call_n = 0, call_bp = 0, vcf_only_n = 0, vcf_only_bp = 0;
        std::size_t vcf_n = 0, vcf_bp = 0, clean_call = 0, clean_vcf = 0, present_small = 0, small_call = 0, small_vcf = 0;
        std::size_t present_snv = 0, snv_call = 0, snv_vcf = 0, applied_call = 0, applied_vcf = 0, conflicts = 0;
        std::size_t outside = 0, absorbed = 0, conflicts_vcf = 0;
        std::set<std::pair<std::string, std::string>> mismatching;
        for (const CheckRow& c : sc.checks) mismatching.emplace(c.variant_id, c.haplotype);
        std::ofstream ho(staged.stage(finals[1]));
        std::ofstream eo(staged.stage(finals[2]));
        if (!ho || !eo) throw std::runtime_error("cannot write benchmark outputs at " + out_prefix);
        ho << "haplotype\tstatus\tref_start\tref_end\tlarge_present\tlarge_present_bp\tlarge_left_call\tlarge_left_call_bp"
              "\tlarge_left_vcf\tlarge_left_vcf_bp\tlarge_left_vcf_only\tlarge_left_vcf_only_bp\tsmall_bp_present"
              "\tsmall_bp_left_call\tsmall_bp_left_vcf\tsnv_present\tsnv_left_call\tsnv_left_vcf\n";
        eo << "haplotype\tstage\tref_start\tref_end\ttype\tsize_bp\tdel_bp\tins_bp\treason\tbubbles\trecords\n";
        const auto row_out = [&](const LargeRow& e, const char* stage, const std::string& reason) {
            eo << e.haplotype << '\t' << stage << '\t' << e.g0 << '\t' << e.g1 << '\t' << e.type << '\t' << e.size
               << '\t' << e.del_bp << '\t' << e.ins_bp << '\t' << reason << '\t' << e.bubbles << '\t' << e.records << '\n';
        };
        for (const HapResult& c : sc.results) {
            const auto vit = vres.find(c.name);
            if (c.status != "scored" || vit == vres.end()) {
                ++skipped[c.status];
                ho << c.name << '\t' << c.status;
                for (int k = 0; k < 16; ++k) ho << "\t.";
                ho << '\n';
                continue;
            }
            const HapResult& v = *vit->second;
            ++scored;
            std::size_t h_only = 0, h_only_bp = 0;
            for (const LargeRow& e : c.rows) {
                ++call_by_cause[e.cause];
                row_out(e, "call", e.cause);
            }
            // Each of the VCF rebuild's leftovers is listed whole. One lying on call's leftovers (within
            // --merge-gap; a leftover inside an applied event may lie anywhere over it) is those misses again,
            // up to their size, and takes the cause of the miss it draws most on, unless what goes beyond them
            // reaches --min-sv-bp: then it was lost in the VCF, with a VCF reason. Each of call's leftovers is
            // drawn on once in all, so a miss is never counted twice.
            std::vector<std::size_t> held(c.rows.size(), 0);
            for (const LargeRow& e0 : v.rows) {
                std::size_t beyond = e0.size, most = 0, most_k = 0;
                for (std::size_t k = 0; k < c.rows.size() && beyond > 0; ++k) {
                    const LargeRow& m = c.rows[k];
                    if (!(e0.m0 <= m.m1 + merge_gap && m.m0 <= e0.m1 + merge_gap) || held[k] >= m.size) continue;
                    // the same direction: a loss is not a gain call missed (a replacement is either)
                    if (e0.type != m.type && e0.type != "SUB" && m.type != "SUB") continue;
                    const std::size_t take = std::min(beyond, m.size - held[k]);
                    held[k] += take;
                    beyond -= take;
                    if (take > most) { most = take; most_k = k; }
                }
                LargeRow row = e0;
                if (beyond < min_sv_bp) {
                    row.records = c.rows[most_k].records;
                    row_out(row, "vcf", c.rows[most_k].cause);
                    continue;
                }
                const std::string why = lost_reason(e0, c, row.records);
                ++vcf_by_reason[why];
                if (why == "no_sequence") {   // by why its record could not be spelled
                    const auto ids = split_on(row.records, ';');
                    const auto it = ids.empty() ? sv.unspelled.end() : sv.unspelled.find(ids.front());
                    ++unspelled_lost[it == sv.unspelled.end() ? "other" : it->second];
                }
                ++h_only; h_only_bp += e0.size;
                row_out(row, "vcf", why);
            }
            present_n += c.present_n; present_bp += c.present_bp;
            call_n += c.left_n; call_bp += c.left_bp;
            vcf_n += v.left_n; vcf_bp += v.left_bp;
            vcf_only_n += h_only; vcf_only_bp += h_only_bp;
            clean_call += c.left_n == 0; clean_vcf += v.left_n == 0;
            present_small += c.present_small_bp; small_call += c.left_small_bp; small_vcf += v.left_small_bp;
            present_snv += c.present_snv; snv_call += c.left_snv; snv_vcf += v.left_snv;
            applied_call += c.applied; applied_vcf += v.applied;
            conflicts += c.conflicts; conflicts_vcf += v.conflicts; outside += c.outside_extent;
            absorbed += c.absorbed;
            ho << c.name << '\t' << c.status << '\t' << c.g_start << '\t' << c.g_end << '\t' << c.present_n << '\t'
               << c.present_bp << '\t' << c.left_n << '\t' << c.left_bp << '\t' << v.left_n << '\t' << v.left_bp << '\t'
               << h_only << '\t' << h_only_bp << '\t' << c.present_small_bp << '\t' << c.left_small_bp << '\t'
               << v.left_small_bp << '\t' << c.present_snv << '\t' << c.left_snv << '\t' << v.left_snv << '\n';
        }
        ho.close(); eo.close();
        if (!ho || !eo) throw std::runtime_error("failed writing benchmark outputs at " + out_prefix);
        {
            // What there was to find, per scored haplotype: the truth set a comparison with other callers uses.
            std::ofstream po(staged.stage(finals[4]));
            if (!po) throw std::runtime_error("cannot write " + finals[4]);
            po << "haplotype\tref_start\tref_end\ttype\tsize_bp\tdel_bp\tins_bp\tbubbles\n";
            for (const HapResult& c : sc.results) {
                if (c.status != "scored" || !vres.count(c.name)) continue;
                for (const HapResult::Present& p : c.present_rows)
                    po << c.name << '\t' << p.g0 << '\t' << p.g1 << '\t' << p.type << '\t' << p.size << '\t'
                       << p.del_bp << '\t' << p.ins_bp << '\t' << p.bubbles << '\n';
            }
            po.close();
            if (!po) throw std::runtime_error("failed writing " + finals[4]);
        }
        {
            std::ofstream co(staged.stage(finals[3]));
            if (!co) throw std::runtime_error("cannot write " + finals[3]);
            co << "variant_id\thaplotype\tkind\treported\town_event\tdifference_bp\n";
            for (const CheckRow& c : sc.checks)
                co << c.variant_id << '\t' << c.haplotype << '\t' << c.kind << '\t' << c.reported << '\t' << c.own
                   << '\t' << c.difference << '\n';
            co.close();
            if (!co) throw std::runtime_error("failed writing " + finals[3]);
        }
        // 1 - left / present; negative when more is left than was there
        const auto recovered = [&](std::size_t left) {
            if (present_bp == 0) return std::string("NA");
            char b[32];
            std::snprintf(b, sizeof(b), "%.4f", 1.0 - static_cast<double>(left) / static_cast<double>(present_bp));
            return std::string(b);
        };
        std::ofstream so(staged.stage(finals[0]));
        if (!so) throw std::runtime_error("cannot write " + finals[0]);
        so << "key\tvalue\n"
           << "min_sv_bp\t" << min_sv_bp << '\n'
           << "merge_gap\t" << merge_gap << '\n'
           << "haplotypes_scored\t" << scored << '\n';
        for (const auto& [why, n] : skipped) so << "haplotypes_skipped_" << why << '\t' << n << '\n';
        so << "large_present\t" << present_n << '\n'
           << "large_present_bp\t" << present_bp << '\n'
           << "large_left_call\t" << call_n << '\n'
           << "large_left_call_bp\t" << call_bp << '\n';
        for (const std::string& c : causes) so << "large_left_call_" << c << '\t' << call_by_cause[c] << '\n';
        so << "large_left_vcf\t" << vcf_n << '\n'
           << "large_left_vcf_bp\t" << vcf_bp << '\n'
           << "large_left_vcf_only\t" << vcf_only_n << '\n'
           << "large_left_vcf_only_bp\t" << vcf_only_bp << '\n';
        for (const std::string& r : vcf_reasons) so << "large_left_vcf_" << r << '\t' << vcf_by_reason[r] << '\n';
        so << "recovered_by_call_fraction\t" << recovered(call_bp) << '\n'
           << "recovered_by_vcf_fraction\t" << recovered(vcf_bp) << '\n'
           << "haplotypes_clean_call\t" << clean_call << '\n'
           << "haplotypes_clean_vcf\t" << clean_vcf << '\n'
           << "small_bp_present\t" << present_small << '\n'
           << "small_bp_left_call\t" << small_call << '\n'
           << "small_bp_left_vcf\t" << small_vcf << '\n'
           << "snv_present\t" << present_snv << '\n'
           << "snv_left_call\t" << snv_call << '\n'
           << "snv_left_vcf\t" << snv_vcf << '\n'
           << "carriers_checked\t" << sc.checked << '\n'
           << "carrier_mismatches\t" << mismatching.size() << '\n'
           << "cn_exact_checked\t" << sc.cn_exact_checked << '\n';
        so.close();
        if (!so) throw std::runtime_error("failed writing " + finals[0]);
        log.info("present: " + std::to_string(present_bp) + " bp of large events");
        log.info("left by call: " + std::to_string(call_bp) + " bp");
        log.info("left by the VCF: " + std::to_string(vcf_bp) + " bp; by the VCF only: " + std::to_string(vcf_only_bp) + " bp");
        log.info("carriers disagreeing with their record: " + std::to_string(mismatching.size()));
        // How each rebuild was made: counts to read when a result looks wrong, warned about where they can
        // mean a loss.
        log.info("events applied: call: " + std::to_string(applied_call) + "; VCF: " + std::to_string(applied_vcf) +
                 "; left out for overlapping another: call: " + std::to_string(conflicts) + "; VCF: " +
                 std::to_string(conflicts_vcf));
        log.info("events not applied: inside a carrier's own copy stretch: " + std::to_string(absorbed) +
                 "; past the scored stretch: " + std::to_string(outside));
        log.info("VCF genotypes . (no allele there): " + std::to_string(sv.missing_gt));
        const auto warn = [&](std::size_t n, const std::string& what) {
            if (n > 0) std::cerr << cli::log_prefix("benchmark") << "warning: " << n << ' ' << what << '\n';
        };
        warn(sc.unresolved_rows, "carrier event(s) call could not place on the reference");
        warn(sv.not_reconstructable, "VCF record(s) that could not be spelled for some carrier");
        warn(unspelled_lost["sequence"], "event(s) left by the VCF only, at a record whose sequence call left out above its "
             "--max-seq-bp (call --max-seq-bp 0 writes them all)");
        warn(unspelled_lost["no_copy"], "event(s) left by the VCF only, at a copy-number record with no reference copy to paste");
        warn(unspelled_lost["other"], "event(s) left by the VCF only, at a record that could not be spelled");
        warn(sc.not_in_vcf, "haplotype(s) with no VCF column");
        warn(sc.unchecked + sc.cn_unchecked, "carrier check(s) skipped: the record has no single SVLEN, or no CN or copy length to compare");
        staged.commit();
    };
    for (std::size_t m = 0; m < 2; ++m) score_mode(specs[m], stages[m]);
    write_outputs();
    log.wrote(finals);
    log.done();
    return 0;
}

} // namespace panvar
