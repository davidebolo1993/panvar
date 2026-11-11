#include "panvar/recruit_scoring.hpp"

#include "panvar/kmer.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace panvar::recruit::scoring {

namespace {




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

bool is_certain(const Evidence& v, const ReadOptions& o, std::uint32_t min_sites) {
    if (v.seg_kmers == 0 || uncallable(v) || divergent(v, o)) return false;
    switch (o.method) {
        case Method::Runs:
            return v.sites_p >= min_sites && v.sites_d <= o.decoy_share * v.sites_p &&
                   v.sites_p >= o.min_density * static_cast<double>(v.seg_kmers) / static_cast<double>(kK);
        case Method::Kmers: return v.n_d == 0 && v.n_p >= std::max<std::uint32_t>(3, min_sites);
        case Method::Fraction:
            return v.kmers > 0 && v.n_p >= 1 && static_cast<double>(v.n_p + v.n_s + v.n_r) / v.kmers >= o.min_fraction;
    }
    return false;
}

// A fragment's fate. The pair's pooled evidence decides; failing that, one mate certain on its own makes
// the pair certain too ("rescued"): the mates lie an insert apart, so the other one is in the locus or its
// flank. Certain fragments with a mate outside every segment are rescued, with a mate partly outside partial.
void decide(Row& r, const Evidence& one, const Evidence* two, const char* status[2], const ReadOptions& o, std::uint32_t min_sites) {
    const Evidence& v = r.ev;
    // pooled evidence decides a pair only with a mate wholly in the locus: two mates that both leave it
    // would pool islands (a pair across an edge has a mate inside unless the mates overlap)
    const bool pooled = is_certain(v, o, min_sites) &&
                        (two == nullptr || std::string(status[0]) == "in" || std::string(status[1]) == "in");
    // one mate certain on its own carries the pair unless the other leans to a homolog (more decoy sites
    // than private ones) away from a locus edge; at an edge (the certain mate runs out of the locus) the
    // other mate may well be in a homologous flank
    auto carries = [&](const Evidence& mine, const Evidence& other, const char* my_status) {
        return is_certain(mine, o, min_sites) &&
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
    if (v.sites_p >= min_sites) { r.fate = "ambiguous"; r.reason = "sparse"; return; }
    if (v.sites_p >= 1) { r.fate = "ambiguous"; r.reason = "single_site"; return; }
    if (v.n_s >= 1) { r.fate = "ambiguous"; r.reason = "shared_only"; return; }
    r.fate = "unrecruited";
    r.reason = "no_evidence";
}

// Scores one fragment against one locus.
}  // namespace

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
    decide(r, ev[0], pair ? &ev[1] : nullptr, status, o, t.min_sites ? t.min_sites : o.min_sites);
    return r;
}


Costs costs_for(double error_rate) {
    // a run of k-mers in neither table costs this much per k-mer past its first k: half the expected
    // score of a locus stretch read at the error rate (error-free k-mers +1, the rest at this cost)
    Costs c;
    const double q = std::pow(std::max(0.0, 1.0 - error_rate - 0.001), static_cast<double>(kK));
    c.gap_cost = q >= 1.0 ? 1e9 : 0.5 * q / (1.0 - q);
    c.flank_share = 0.25 * q;   // locus k-mers beyond a segment's end below this: a flank
    return c;
}

LocusTable table_of(const LocusIndex& x) {
    LocusTable t;
    t.name = x.name;
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
    return t;
}

}  // namespace panvar::recruit::scoring
