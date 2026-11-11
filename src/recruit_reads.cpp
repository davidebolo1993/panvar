#include "panvar/recruit_reads.hpp"

#include "panvar/flat_map.hpp"
#include "panvar/read_input.hpp"
#include "panvar/kmer.hpp"

#include <zlib.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <thread>

namespace panvar::recruit {

namespace {

constexpr std::size_t kK = 31, kS = 11;
constexpr std::uint8_t kLabelP = 0, kLabelS = 1, kLabelR = 2, kLabelD = 3, kLabelN = 4, kMultiBit = 8, kBackboneBit = 16;

// One locus, as reads see it: each k-mer's label (P/S/R, D for decoys), with bits for multi-copy P and
// backbone k-mers (single-copy, in every haplotype: their depth is the sample's depth at the locus); the
// graph nodes each private k-mer lies on.
struct LocusTable {
    std::string name;
    FlatMap<std::uint8_t> label;
    FlatMap<std::uint32_t> nodeset;   // private k-mer -> its node set
    std::vector<std::uint32_t> nodeset_first, nodeset_nodes;
    std::vector<std::string> nodes;
    std::uint64_t backbone = 0;       // backbone k-mers in the index
};

struct Row {
    std::string fragment;
    std::uint32_t locus = 0;
    Evidence ev;
    std::string nodes;           // graph nodes of the private k-mers in the segments, "." if none
    std::uint32_t backbone = 0;  // backbone k-mers in the segments
    std::string seg[2];          // each mate's locus segment in read bases, "." when it has none
    const char* mates = "";      // each mate's status: in, part (a segment covering part of it), out
    const char* fate = "";
    const char* reason = "";
};

struct Kmer {
    std::uint64_t code;
    std::uint32_t start;
    std::uint8_t lab;
    bool backbone;
};

// One mate's k-mers, labelled, with the whole-read counts.
std::vector<Kmer> label_read(const std::string& s, const LocusTable& t, Evidence& ev) {
    ev.bases += static_cast<std::uint32_t>(s.size());
    if (s.size() >= kK) ev.windows += static_cast<std::uint32_t>(s.size() - kK + 1);
    std::vector<Kmer> v;
    kmer::for_each_canonical_kmer(s.data(), s.size(), kK, [&](std::uint64_t c, std::size_t start, bool) {
        ++ev.kmers;
        const std::uint8_t* it = t.label.find(c);
        std::uint8_t lab = kLabelN;
        bool backbone = false;
        if (it != nullptr) {
            lab = *it & 7;
            backbone = *it & kBackboneBit;
            if (lab == kLabelP && (*it & kMultiBit)) ++ev.multi;
        }
        switch (lab) {
            case kLabelP: ++ev.n_p; break;
            case kLabelS: ++ev.n_s; break;
            case kLabelR: ++ev.n_r; break;
            case kLabelD: ++ev.n_d; break;
            default: ++ev.n_n;
        }
        v.push_back({c, static_cast<std::uint32_t>(start), lab, backbone});
    });
    return v;
}

// The read's locus segment: the stretch of k-mers with the highest score, a locus k-mer (P, S or R) +1, a
// decoy k-mer -2, and a k-mer in neither table free for the first k of a run (one differing base makes
// that many) and `gap_cost` after it, so a long run means the read has left the locus. Returns [first,
// last] k-mer indices (of equal scores, the longest); first > last when the read has no locus k-mer.
std::pair<std::size_t, std::size_t> locus_segment(const std::vector<Kmer>& v, double gap_cost) {
    double best = 0, cur = 0;
    std::size_t best_first = 1, best_last = 0, first = 0, run = 0;
    for (std::size_t i = 0; i < v.size(); ++i) {
        const bool continues = i > 0 && v[i].start == v[i - 1].start + 1;
        double x;
        if (v[i].lab == kLabelN) {
            run = continues && v[i - 1].lab == kLabelN ? run + 1 : 1;
            x = run <= kK ? 0.0 : -gap_cost;
        } else {
            run = 0;
            x = v[i].lab == kLabelD ? -2.0 : 1.0;
        }
        if (cur <= 0) { cur = 0; first = i; }
        cur += x;
        if (cur > best || (cur == best && cur > 0 && i - first > best_last - best_first)) { best = cur; best_first = first; best_last = i; }
    }
    return {best_first, best_last};
}

// The deciding evidence of a segment. Private and decoy k-mers count as independent sites only k bases
// apart (greedily, in read order), so one differing base, whose k-mers may be split by a shared one, is one
// site; codes in `seen` (the other mate's) are not counted again. Decoy sites outside the segment are kept
// apart when `boundary` (the segment reaches an end of the read, as where a read leaves the locus into a
// homologous flank), else they count with the segment's. Codes counted are added to `mine`.
void count_segment(const std::vector<Kmer>& v, std::pair<std::size_t, std::size_t> seg, bool boundary, Evidence& ev,
                   const std::vector<std::uint64_t>* seen, std::vector<std::uint64_t>* mine) {
    auto counted_before = [&](std::uint64_t c) { return seen != nullptr && std::binary_search(seen->begin(), seen->end(), c); };
    // A site is a group of k-mers starting within k bases of the group's first; a group with a k-mer the
    // other mate counted is that mate's site, whichever end of it this mate reaches first.
    struct Sites {
        bool open = false, seen = false;
        std::size_t first = 0;
        std::uint32_t* counter = nullptr;
        void add(std::size_t start, bool already, std::uint32_t* to) {
            if (!open || start >= first + kK) {
                close();
                open = true;
                seen = false;
                first = start;
                counter = to;
            }
            seen = seen || already;
        }
        void close() {
            if (open && !seen) ++*counter;
            open = false;
        }
    } p, d, out;
    bool in_n = false;
    for (std::size_t i = 0; i < v.size(); ++i) {
        const Kmer& k = v[i];
        const bool inside = i >= seg.first && i <= seg.second;
        if (!inside) {
            if (k.lab == kLabelD) {
                out.add(k.start, counted_before(k.code), boundary ? &ev.decoy_out : &ev.sites_d);
                if (mine) mine->push_back(k.code);
            }
            continue;
        }
        ++ev.seg_kmers;
        if (in_n && (k.lab != kLabelN || k.start != v[i - 1].start + 1)) in_n = false;
        if (k.lab == kLabelP) {
            p.add(k.start, counted_before(k.code), &ev.sites_p);
            if (mine) mine->push_back(k.code);
        } else if (k.lab == kLabelD) {
            d.add(k.start, counted_before(k.code), &ev.sites_d);
            if (mine) mine->push_back(k.code);
        } else if (k.lab == kLabelN && !in_n) {
            ++ev.runs_n;
            in_n = true;
        }
    }
    p.close();
    d.close();
    out.close();
}

// The most mismatch runs a locus read of this many callable k-mers should show: an upper Poisson bound on
// sequencing errors plus about one new allele per kilobase (normal approximation, 3.1 standard deviations).
double tolerated_runs(std::uint32_t callable, double e) {
    const double lambda = (e + 0.001) * static_cast<double>(callable);
    return lambda + 3.1 * std::sqrt(lambda) + 2.0;
}

bool uncallable(const Evidence& v) { return v.windows > 0 && v.kmers * 2 < v.windows; }
bool divergent(const Evidence& v, const ReadOptions& o) { return static_cast<double>(v.runs_n) > tolerated_runs(v.seg_kmers, o.error_rate); }

bool is_certain(const Evidence& v, const ReadOptions& o) {
    if (v.seg_kmers == 0 || uncallable(v) || divergent(v, o)) return false;
    switch (o.method) {
        case Method::Runs:
            return v.sites_p >= o.min_sites && v.sites_d <= o.decoy_share * v.sites_p &&
                   v.sites_p >= o.min_density * static_cast<double>(v.seg_kmers) / static_cast<double>(kK);
        case Method::Kmers: return v.n_d == 0 && v.n_p >= std::max<std::uint32_t>(3, o.min_sites);
        case Method::Fraction:
            return v.kmers > 0 && v.n_p >= 1 && static_cast<double>(v.n_p + v.n_s + v.n_r) / v.kmers >= o.min_fraction;
    }
    return false;
}

// A fragment's fate. The pair's pooled evidence decides; failing that, one mate certain on its own makes
// the pair certain too ("rescued"): the mates lie an insert apart, so the other one is in the locus or its
// flank. Certain fragments with a mate outside every segment are rescued, with a mate partly outside partial.
void decide(Row& r, const Evidence& one, const Evidence* two, const char* status[2], const ReadOptions& o) {
    const Evidence& v = r.ev;
    // pooled evidence decides a pair only with a mate wholly in the locus: two mates that both leave it
    // would pool islands (a pair across an edge has a mate inside unless the mates overlap)
    const bool pooled = is_certain(v, o) &&
                        (two == nullptr || std::string(status[0]) == "in" || std::string(status[1]) == "in");
    // one mate certain on its own carries the pair unless the other leans to a homolog (more decoy sites
    // than private ones) away from a locus edge; at an edge (the certain mate runs out of the locus) the
    // other mate may well be in a homologous flank
    auto carries = [&](const Evidence& mine, const Evidence& other, const char* my_status) {
        return is_certain(mine, o) &&
               (other.sites_d + other.decoy_out <= other.sites_p || std::string(my_status) == "part");
    };
    const bool alone = two != nullptr && (carries(one, *two, status[0]) || carries(*two, one, status[1]));
    if (pooled || alone) {
        r.fate = "certain";
        const bool out = two != nullptr && (std::string(status[0]) == "out" || std::string(status[1]) == "out");
        const bool part = std::string(status[0]) == "part" || (two && std::string(status[1]) == "part");
        r.reason = !pooled || out ? "rescued" : part ? "partial" : ".";
        return;
    }
    const std::uint32_t decoy = v.sites_d + v.decoy_out;
    if (v.sites_p == 0 && decoy >= 2) { r.fate = "unrecruited"; r.reason = "homolog"; return; }
    if (v.sites_p >= 1 && decoy >= 1) { r.fate = "ambiguous"; r.reason = "conflict"; return; }
    if (uncallable(v) && v.sites_p >= 1) { r.fate = "ambiguous"; r.reason = "uncallable"; return; }
    if (divergent(v, o) && v.sites_p >= 1) { r.fate = "ambiguous"; r.reason = "divergent"; return; }
    if (v.sites_p >= o.min_sites) { r.fate = "ambiguous"; r.reason = "sparse"; return; }
    if (v.sites_p >= 1) { r.fate = "ambiguous"; r.reason = "single_site"; return; }
    if (v.n_s >= 1) { r.fate = "ambiguous"; r.reason = "shared_only"; return; }
    r.fate = "unrecruited";
    r.reason = "no_evidence";
}

// Scores one fragment against one locus.
Row score_fragment(const Fragment& fr, const LocusTable& t, std::uint32_t li, const ReadOptions& o, double gap_cost,
                   double flank_share) {
    Row r;
    r.fragment = fr.name;
    r.locus = li;
    const bool pair = fr.paired;
    Evidence ev[2], dedup2;
    const char* status[2] = {"out", "out"};
    std::vector<std::uint64_t> seen;
    std::vector<std::uint32_t> sets;   // node sets of the private k-mers in the segments
    for (int m = 0; m < (pair ? 2 : 1); ++m) {
        const std::string& s = m == 0 ? fr.seq1 : fr.seq2;
        const std::vector<Kmer> v = label_read(s, t, ev[m]);
        auto seg = locus_segment(v, gap_cost);
        // a segment shorter than k k-mers (under 2k - 1 bases of locus) is no evidence of the locus
        if (seg.first <= seg.second && seg.second - seg.first + 1 < kK) seg = {1, 0};
        // The segment alone is judged when the read leaves the locus at its ends: what lies beyond each end
        // (if anything) holds few locus k-mers, under a quarter of what locus sequence read at the error
        // rate gives, and decoy k-mers lie beyond one end at most (a homologous flank on one side).
        // Otherwise (an error cluster or a new allele stopped the segment) the whole read is judged.
        bool boundary = false;
        if (seg.first <= seg.second) {
            auto side = [&](std::size_t a, std::size_t b, bool& decoy) {   // locus k-mer share in [a, b)
                std::size_t locus = 0;
                for (std::size_t i = a; i < b; ++i) {
                    if (v[i].lab == kLabelD) decoy = true;
                    else if (v[i].lab != kLabelN) ++locus;
                }
                return b > a ? static_cast<double>(locus) / static_cast<double>(b - a) : 0.0;
            };
            bool decoy_l = false, decoy_r = false;
            const double left = side(0, seg.first, decoy_l), right = side(seg.second + 1, v.size(), decoy_r);
            boundary = left < flank_share && right < flank_share && !(decoy_l && decoy_r);
            if (!boundary) seg = {0, v.size() - 1};
        }
        if (m == 0) {
            count_segment(v, seg, boundary, ev[0], nullptr, &seen);
            std::sort(seen.begin(), seen.end());
        } else {
            count_segment(v, seg, boundary, ev[1], nullptr, nullptr);
            dedup2 = ev[1];
            dedup2.sites_p = dedup2.sites_d = dedup2.decoy_out = dedup2.runs_n = dedup2.seg_kmers = 0;
            count_segment(v, seg, boundary, dedup2, &seen, nullptr);
        }
        for (std::size_t i = seg.first; i <= seg.second && i < v.size(); ++i) {
            if (v[i].backbone) ++r.backbone;
            if (v[i].lab == kLabelP)
                if (const std::uint32_t* ns = t.nodeset.find(v[i].code)) sets.push_back(*ns);
        }
        if (seg.first <= seg.second) {
            r.seg[m] = std::to_string(v[seg.first].start) + "-" + std::to_string(v[seg.second].start + kK);
            status[m] = ev[m].seg_kmers * 10 >= ev[m].kmers * 9 ? "in" : "part";
        } else {
            r.seg[m] = ".";
        }
    }
    std::sort(sets.begin(), sets.end());
    sets.erase(std::unique(sets.begin(), sets.end()), sets.end());
    std::vector<std::uint32_t> nodes;
    for (const std::uint32_t ns : sets)
        for (std::uint32_t i = t.nodeset_first[ns]; i < t.nodeset_first[ns + 1]; ++i) nodes.push_back(t.nodeset_nodes[i]);
    std::sort(nodes.begin(), nodes.end());
    nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
    for (const std::uint32_t nd : nodes) r.nodes += (r.nodes.empty() ? "" : ",") + t.nodes[nd];
    if (r.nodes.empty()) r.nodes = ".";
    r.ev = ev[0];
    if (pair) r.ev.add(dedup2);
    static const char* kMates[3][3] = {{"in,in", "in,part", "in,out"}, {"part,in", "part,part", "part,out"}, {"out,in", "out,part", "out,out"}};
    auto idx = [](const char* s) { return s[0] == 'i' ? 0 : s[0] == 'p' ? 1 : 2; };
    r.mates = pair ? kMates[idx(status[0])][idx(status[1])] : status[0];
    if (!pair) r.seg[1] = ".";
    decide(r, ev[0], pair ? &ev[1] : nullptr, status, o);
    return r;
}

}  // namespace

void Evidence::add(const Evidence& o) {
    sites_p += o.sites_p; sites_d += o.sites_d; decoy_out += o.decoy_out; runs_n += o.runs_n; seg_kmers += o.seg_kmers;
    n_p += o.n_p; n_s += o.n_s; n_r += o.n_r; n_d += o.n_d; n_n += o.n_n;
    kmers += o.kmers; windows += o.windows; bases += o.bases; multi += o.multi;
}

namespace {

// Recruited reads as gzip FASTQ, per locus and fate: <locus>.<fate>_1/_2.fq.gz for pairs, <locus>.<fate>.fq.gz
// for single reads. Files open when first written; past a cap they are closed and reopened to append (a
// gzip file of several members reads as one).
class ReadWriters {
public:
    explicit ReadWriters(std::string dir) : dir_(std::move(dir)) { std::filesystem::create_directories(dir_); }
    ~ReadWriters() { close_all(); }
    void write(const std::string& locus, const char* fate, const char* suffix, const std::string& header, const std::string& seq,
               const std::string& qual) {
        const std::string key = locus + "." + fate + suffix + ".fq.gz";
        auto it = open_.find(key);
        if (it == open_.end()) {
            if (open_.size() >= kMaxOpen) close_all();
            const bool again = started_.count(key) > 0;
            gzFile f = gzopen((std::filesystem::path(dir_) / key).string().c_str(), again ? "ab3" : "wb3");
            if (f == nullptr) throw std::runtime_error("recruit: cannot write " + dir_ + "/" + key);
            started_.insert(key);
            it = open_.emplace(key, f).first;
        }
        std::string rec = "@" + header + "\n" + seq + "\n+\n";
        rec += qual.size() == seq.size() ? qual : std::string(seq.size(), 'I');
        rec += '\n';
        if (gzwrite(it->second, rec.data(), static_cast<unsigned>(rec.size())) != static_cast<int>(rec.size()))
            throw std::runtime_error("recruit: cannot write " + dir_ + "/" + key);
    }
    void close_all() {
        for (auto& [key, f] : open_)
            if (gzclose(f) != Z_OK) throw std::runtime_error("recruit: cannot write " + dir_ + "/" + key);
        open_.clear();
    }

private:
    static constexpr std::size_t kMaxOpen = 200;
    std::string dir_;
    std::map<std::string, gzFile> open_;
    std::set<std::string> started_;
};

// The error rate from base qualities: the mean of 10^(-Q/10) over the bases that carry one, within [0.001,
// 0.08]; without qualities, or with qualities that say nothing (a mean above 30%, as from simulators that
// write Q0), 0.03 for long reads and 0.002 for short ones.
double estimated_error(const std::vector<Fragment>& batch, std::string& source) {
    double sum = 0;
    std::uint64_t n = 0, bases = 0, reads = 0;
    std::array<double, 94> p{};
    for (int q = 0; q < 94; ++q) p[q] = std::pow(10.0, -q / 10.0);
    for (const Fragment& f : batch)
        for (const std::string* q : {&f.qual1, &f.qual2}) {
            for (const char c : *q) {
                const int v = std::clamp(static_cast<int>(c) - 33, 0, 93);
                sum += p[v];
                ++n;
            }
        }
    for (const Fragment& f : batch) {
        bases += f.seq1.size() + f.seq2.size();
        reads += 1 + !f.seq2.empty();
    }
    const double raw = n ? sum / static_cast<double>(n) : 1.0;
    if (n >= 1000 && raw <= 0.3) {
        source = "estimated from " + std::to_string(n) + " base qualities";
        return std::clamp(raw, 0.001, 0.08);
    }
    const bool long_reads = reads > 0 && bases / reads > 1000;
    source = std::string(n >= 1000 ? "base qualities uninformative (all near Q0)" : "no base qualities") + ": the default for " +
             (long_reads ? "long" : "short") + " reads";
    return long_reads ? 0.03 : 0.002;
}

}  // namespace

ReadStats recruit_reads(const ReadOptions& options) {
    ReadOptions o = options;   // its error rate may be estimated below
    auto log = [&](const std::string& s) { if (o.log) o.log(s); };
    if (o.indexes.empty()) throw std::runtime_error("recruit: no locus index");
    std::vector<LocusTable> loci(o.indexes.size());
    std::vector<std::pair<std::uint64_t, std::uint32_t>> screen_items;
    for (std::size_t li = 0; li < o.indexes.size(); ++li) {
        const LocusIndex x = load(o.indexes[li]);
        LocusTable& t = loci[li];
        t.name = x.name;
        for (std::size_t lj = 0; lj < li; ++lj)
            if (loci[lj].name == x.name)
                throw std::runtime_error("recruit: two loci are named " + x.name + " (" + o.indexes[lj] + ", " + o.indexes[li] + "); rename one graph");
        t.label.reserve(x.codes.size() + x.decoy.size());
        std::size_t np = 0;
        for (std::size_t i = 0; i < x.codes.size(); ++i) np += x.cls[i] == kPrivate;
        t.nodeset.reserve(np);
        bool existing = false;
        for (std::size_t i = 0; i < x.codes.size(); ++i) {
            const bool backbone = x.backbone[i] != 0;
            t.backbone += backbone;
            t.label.slot(x.codes[i], existing) = static_cast<std::uint8_t>(
                x.cls[i] | (x.cls[i] == kPrivate && x.mult[i] > 1 ? kMultiBit : 0) | (backbone ? kBackboneBit : 0));
            if (x.cls[i] == kPrivate) t.nodeset.slot(x.codes[i], existing) = x.nodeset[i];
        }
        for (const std::uint64_t d : x.decoy) {
            std::uint8_t& v = t.label.slot(d, existing);
            if (!existing) v = kLabelD;
        }
        t.nodeset_first = x.nodeset_first;
        t.nodeset_nodes = x.nodeset_nodes;
        t.nodes = x.nodes;
        for (const std::uint64_t s : x.screen) screen_items.emplace_back(s, static_cast<std::uint32_t>(li));
        log(x.name + ": " + std::to_string(x.codes.size()) + " k-mers, " + std::to_string(x.decoy.size()) + " decoys, " +
            std::to_string(x.screen.size()) + " screen syncmers");
    }
    // the screen: syncmer -> an entry of `screen_loci`, the loci it belongs to (any number of loci)
    std::sort(screen_items.begin(), screen_items.end());
    std::vector<std::vector<std::uint32_t>> screen_loci;
    std::map<std::vector<std::uint32_t>, std::uint32_t> set_id;   // few distinct sets: loci rarely share syncmers
    FlatMap<std::uint32_t> screen;
    screen.reserve(screen_items.size());
    for (std::size_t i = 0; i < screen_items.size();) {
        std::size_t j = i;
        std::vector<std::uint32_t> ls;
        while (j < screen_items.size() && screen_items[j].first == screen_items[i].first) ls.push_back(screen_items[j++].second);
        ls.erase(std::unique(ls.begin(), ls.end()), ls.end());
        const auto [it, fresh] = set_id.emplace(ls, static_cast<std::uint32_t>(screen_loci.size()));
        if (fresh) screen_loci.push_back(std::move(ls));
        bool existing = false;
        screen.slot(screen_items[i].first, existing) = it->second;
        i = j;
    }
    screen_items = {};

    // whether a read hits some locus's screen (the input layer asks, for coordinate-sorted files)
    auto screened = [&](const std::string& seq) {
        bool hit = false;
        kmer::for_each_closed_syncmer(seq.data(), seq.size(), kK, kS, [&](std::uint64_t c, std::size_t, bool) {
            if (!hit && screen.find(c) != nullptr) hit = true;
        });
        return hit;
    };
    if (o.inputs.empty()) throw std::runtime_error("recruit: no reads given");
    InputOptions io;
    io.cram_reference = o.cram_reference;
    io.interleaved = o.interleaved;
    io.threads = std::max<std::size_t>(1, o.threads / 4);   // decompression: a quarter of the scorers keeps up
    io.qualities = !o.reads_out.empty() || o.error_rate < 0;
    io.screened = screened;
    io.log = o.log;
    std::size_t next_input = 0;
    std::unique_ptr<FragmentSource> source;
    ReadStats st;
    std::vector<std::pair<std::string, InputStats>> input_stats;
    auto close_source = [&] {
        if (!source) return;
        const InputStats& s = source->stats();
        st.fragments += s.unscreened + s.dropped;   // passed over before scoring, unscreened
        input_stats.emplace_back(source->describe(), s);
        log(source->describe() + ": " + std::to_string(s.records) + " records, " + std::to_string(s.fragments + s.unscreened + s.dropped) +
            " fragments (" + std::to_string(s.pairs + s.unscreened) + " pairs)" +
            (s.secondary + s.supplementary + s.qc_fail > 0 ? "; skipped " + std::to_string(s.secondary) + " secondary, " +
                 std::to_string(s.supplementary) + " supplementary, " + std::to_string(s.qc_fail) + " QC-failed" : std::string()) +
            (s.fetched > 0 ? "; " + std::to_string(s.fetched) + " mates fetched" : std::string()) +
            (s.unmatched > 0 ? "; " + std::to_string(s.unmatched) + " reads without their mate" : std::string()));
        source.reset();
    };
    constexpr std::size_t kBatch = 20000;
    // The next batch is read while the workers score the current one.
    auto read_batch = [&](std::vector<Fragment>& batch) {   // the inputs one after another
        batch.clear();
        std::vector<Fragment> part;
        while (batch.size() < kBatch) {
            if (!source) {
                if (next_input == o.inputs.size()) break;
                const auto& [a, b] = o.inputs[next_input++];
                source = b.empty() ? open_input(a, io) : open_fastq_pair(a, b);
            }
            if (source->next(part, kBatch - batch.size()) == 0) {
                close_source();
                continue;
            }
            for (auto& f : part) batch.push_back(std::move(f));
        }
    };
    std::vector<Fragment> batch, next;
    read_batch(batch);
    if (o.error_rate < 0) {
        o.error_rate = estimated_error(batch, st.error_source);
    } else {
        st.error_source = "given";
    }
    st.error_rate = o.error_rate;
    {
        char b[64];
        std::snprintf(b, sizeof b, "%.4f", o.error_rate);
        log(std::string("error rate ") + b + " (" + st.error_source + ")");
    }

    std::FILE* out = std::fopen(o.out.c_str(), "w");
    if (out == nullptr) throw std::runtime_error("recruit: cannot write " + o.out);
    std::fprintf(out, "fragment\tlocus\tfate\treason\tmates\tseg1\tseg2\tsites_p\tsites_d\tdecoy_out\truns_n\tseg_kmers\tn_p\tn_s\tn_r\tn_d\tn_n\tkmers\twindows\tbases\tmulti_p\tnodes\tother_loci\tcertain_loci\n");
    std::unique_ptr<ReadWriters> writers;
    if (!o.reads_out.empty()) writers = std::make_unique<ReadWriters>(o.reads_out);
    // a run of k-mers in neither table costs this much per k-mer past its first k: half the expected
    // score of a locus stretch read at the error rate (error-free k-mers +1, the rest at this cost)
    const double q = std::pow(std::max(0.0, 1.0 - o.error_rate - 0.001), static_cast<double>(kK));
    const double gap_cost = q >= 1.0 ? 1e9 : 0.5 * q / (1.0 - q);
    const double flank_share = 0.25 * q;   // locus k-mers beyond a segment's end below this: a flank

    // per locus: fragments by fate and reason, backbone k-mers in certain fragments, certain in several loci
    std::vector<std::map<std::string, std::uint64_t>> by_reason(loci.size());
    std::vector<std::uint64_t> backbone_seen(loci.size(), 0), multi_certain(loci.size(), 0);
    std::uint64_t read_bases = 0, reads_seen = 0;
    const std::size_t threads = std::max<std::size_t>(1, o.threads);
    std::vector<std::vector<Row>> rows;
    while (!batch.empty()) {
        rows.assign(batch.size(), {});
        std::atomic<std::size_t> cursor{0};
        std::vector<std::thread> pool;
        for (std::size_t t = 0; t < threads; ++t)
            pool.emplace_back([&] {
                std::vector<std::uint32_t> hit;
                for (std::size_t i; (i = cursor++) < batch.size();) {
                    const Fragment& fr = batch[i];
                    hit.clear();
                    for (const std::string* s : {&fr.seq1, &fr.seq2})
                        kmer::for_each_closed_syncmer(s->data(), s->size(), kK, kS, [&](std::uint64_t c, std::size_t, bool) {
                            if (const std::uint32_t* id = screen.find(c)) hit.insert(hit.end(), screen_loci[*id].begin(), screen_loci[*id].end());
                        });
                    std::sort(hit.begin(), hit.end());
                    hit.erase(std::unique(hit.begin(), hit.end()), hit.end());
                    for (const std::uint32_t li : hit) rows[i].push_back(score_fragment(fr, loci[li], li, o, gap_cost, flank_share));
                }
            });
        std::exception_ptr failed;
        try { read_batch(next); } catch (...) { failed = std::current_exception(); }
        for (auto& t : pool) t.join();
        if (failed) std::rethrow_exception(failed);
        for (std::size_t i = 0; i < batch.size(); ++i) {
            ++st.fragments;
            const Fragment& fr = batch[i];
            read_bases += fr.seq1.size() + fr.seq2.size();
            reads_seen += 1 + !fr.seq2.empty();
            if (rows[i].empty()) continue;
            ++st.screened;
            std::uint32_t certain_loci = 0;
            for (const Row& r : rows[i]) certain_loci += r.fate[0] == 'c';
            for (const Row& r : rows[i]) {
                const Evidence& v = r.ev;
                std::string others;
                for (const Row& x : rows[i])
                    if (x.locus != r.locus) others += (others.empty() ? "" : ",") + loci[x.locus].name;
                if (others.empty()) others = ".";
                std::fprintf(out, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%s\t%s\t%u\n", r.fragment.c_str(),
                             loci[r.locus].name.c_str(), r.fate, r.reason, r.mates, r.seg[0].c_str(), r.seg[1].c_str(), v.sites_p, v.sites_d,
                             v.decoy_out, v.runs_n, v.seg_kmers, v.n_p, v.n_s, v.n_r, v.n_d, v.n_n, v.kmers, v.windows, v.bases, v.multi,
                             r.nodes.c_str(), others.c_str(), certain_loci);
                ++by_reason[r.locus][std::string(r.fate) + (r.reason[0] == '.' ? "" : std::string("\t") + r.reason)];
                if (r.fate[0] == 'c') {
                    ++st.certain;
                    backbone_seen[r.locus] += r.backbone;
                    if (certain_loci > 1) ++multi_certain[r.locus];
                } else if (r.fate[0] == 'a') {
                    ++st.ambiguous;
                } else {
                    ++st.unrecruited;
                }
                if (writers && (r.fate[0] == 'c' || (r.fate[0] == 'a' && o.reads_ambiguous))) {
                    const std::string header = fr.name + " fate=" + r.fate + " reason=" + r.reason + " mates=" + r.mates +
                                               " sites_p=" + std::to_string(v.sites_p) + " other_loci=" + others;
                    const bool two = fr.paired && !fr.seq1.empty() && !fr.seq2.empty();
                    if (two) {
                        writers->write(loci[r.locus].name, r.fate, "_1", header, fr.seq1, fr.qual1);
                        writers->write(loci[r.locus].name, r.fate, "_2", header, fr.seq2, fr.qual2);
                    } else {
                        const bool first = !fr.seq1.empty();
                        writers->write(loci[r.locus].name, r.fate, "", header, first ? fr.seq1 : fr.seq2, first ? fr.qual1 : fr.qual2);
                    }
                }
            }
        }
        if (std::ferror(out) || std::fflush(out) != 0) throw std::runtime_error("recruit: cannot write " + o.out);
        batch.swap(next);
    }
    close_source();
    if (writers) writers->close_all();
    if (std::fclose(out) != 0) throw std::runtime_error("recruit: cannot write " + o.out);

    if (!o.summary.empty()) {
        std::FILE* sm = std::fopen(o.summary.c_str(), "w");
        if (sm == nullptr) throw std::runtime_error("recruit: cannot write " + o.summary);
        std::fprintf(sm, "section\tkey\tvalue\textra\n");
        std::fprintf(sm, "run\terror_rate\t%.5f\t%s\n", o.error_rate, st.error_source.c_str());
        std::fprintf(sm, "run\tfragments\t%llu\t\nrun\tscreened\t%llu\t\n", static_cast<unsigned long long>(st.fragments),
                     static_cast<unsigned long long>(st.screened));
        const double mean_read = reads_seen ? static_cast<double>(read_bases) / static_cast<double>(reads_seen) : 0;
        std::fprintf(sm, "run\tmean_read_length\t%.1f\t\n", mean_read);
        for (const auto& [name, s] : input_stats)
            std::fprintf(sm, "input\t%s\t%llu\trecords; skipped %llu secondary, %llu supplementary, %llu QC-failed; %llu duplicates, %llu unmapped kept; %llu fragments (%llu pairs), %llu unscreened pairs; %llu mates fetched, %llu without their mate\n",
                         name.c_str(), static_cast<unsigned long long>(s.records), static_cast<unsigned long long>(s.secondary),
                         static_cast<unsigned long long>(s.supplementary), static_cast<unsigned long long>(s.qc_fail),
                         static_cast<unsigned long long>(s.duplicates), static_cast<unsigned long long>(s.unmapped),
                         static_cast<unsigned long long>(s.fragments + s.unscreened + s.dropped), static_cast<unsigned long long>(s.pairs + s.unscreened),
                         static_cast<unsigned long long>(s.unscreened), static_cast<unsigned long long>(s.fetched), static_cast<unsigned long long>(s.unmatched));
        for (std::size_t li = 0; li < loci.size(); ++li) {
            for (const auto& [fate_reason, n] : by_reason[li]) {
                const auto tab = fate_reason.find('\t');
                const std::string fate = fate_reason.substr(0, tab), reason = tab == std::string::npos ? "." : fate_reason.substr(tab + 1);
                std::fprintf(sm, "locus\t%s\t%s\t%s %llu\n", loci[li].name.c_str(), fate.c_str(), reason.c_str(), static_cast<unsigned long long>(n));
            }
            // depth: backbone k-mers seen in certain fragments per backbone k-mer of the index, in k-mers
            // and, for reads of the mean length, in bases
            const double kdepth = loci[li].backbone ? static_cast<double>(backbone_seen[li]) / static_cast<double>(loci[li].backbone) : 0;
            const double bdepth = mean_read > kK ? kdepth * mean_read / (mean_read - kK + 1) : kdepth;
            std::fprintf(sm, "locus\t%s\tbackbone_depth\t%.2f k-mer depth, %.2f base depth over %llu backbone k-mers\n", loci[li].name.c_str(), kdepth, bdepth,
                         static_cast<unsigned long long>(loci[li].backbone));
            std::fprintf(sm, "locus\t%s\tcertain_in_other_loci\t%llu\n", loci[li].name.c_str(), static_cast<unsigned long long>(multi_certain[li]));
        }
        if (std::fclose(sm) != 0) throw std::runtime_error("recruit: cannot write " + o.summary);
    }
    return st;
}

}  // namespace panvar::recruit
