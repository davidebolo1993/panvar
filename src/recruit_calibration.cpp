#include "panvar/recruit_calibration.hpp"

#include "panvar/graph_utils.hpp"
#include "panvar/kmer.hpp"
#include "panvar/read_input.hpp"
#include "panvar/recruit_scoring.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <cmath>
#include <cstring>
#include <map>
#include <string>
#include <thread>
#include <vector>

namespace panvar::recruit {

namespace {

using scoring::kK;
using scoring::kLabelN;
using scoring::kS;

struct Class {
    const char* name;
    double error, indel_share;            // per-base error rate; the share of errors that are indels
    std::uint64_t length, insert, insert_sd;   // insert 0: single reads
    double coverage;                      // of each held-out haplotype
};
constexpr Class kClasses[] = {
    {"short", 0.002, 0.05, 150, 400, 50, 3.0},
    {"hifi", 0.001, 0.5, 15000, 0, 0, 10.0},
    {"ont", 0.03, 0.5, 10000, 0, 0, 10.0},
};
constexpr std::size_t kMaxHeldOut = 8;
// off-target fragments per source (flanks, decoy windows), about: the strides widen to stay near them
constexpr std::uint64_t kOffShort = 100000, kOffLong = 4000;

// A small deterministic generator (splitmix64): one per simulated fragment, seeded by its number, so the
// reads do not depend on the thread count.
struct Rng {
    std::uint64_t s;
    std::uint64_t next() {
        std::uint64_t z = (s += 0x9E3779B97F4A7C15ULL);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }
    double unit() { return static_cast<double>(next() >> 11) * 0x1.0p-53; }
    double normal() {   // Box-Muller
        const double u = std::max(unit(), 1e-12), v = unit();
        return std::sqrt(-2.0 * std::log(u)) * std::cos(6.283185307179586 * v);
    }
};

// Errors at rate e (bases upper-cased: a substitution must change the base whatever its case).
std::string with_errors(const std::string& s, double e, double indel_share, Rng& rng) {
    static const char kBases[] = "ACGT";
    std::string out;
    out.reserve(s.size() + s.size() / 50 + 4);
    for (const char lc : s) {
        const char c = static_cast<char>(std::toupper(static_cast<unsigned char>(lc)));
        if (rng.unit() >= e) { out += c; continue; }
        const double k = rng.unit();
        if (k < indel_share / 2) continue;                                       // deletion
        if (k < indel_share) { out += kBases[rng.next() & 3]; out += c; continue; }   // insertion
        char b;
        do b = kBases[rng.next() & 3]; while (b == c);                          // substitution
        out += b;
    }
    return out;
}

// The same errors, drawn by the gaps between them: one draw per error instead of one per base. Its reads are
// not those of with_errors, which the first rules' counts were made with.
std::string with_errors_by_gaps(const std::string& s, double e, double indel_share, Rng& rng) {
    static const char kBases[] = "ACGT";
    std::string out;
    out.reserve(s.size() + s.size() / 50 + 4);
    const double lq = std::log(1.0 - e);
    auto upper = [](char c) { return static_cast<char>(std::toupper(static_cast<unsigned char>(c))); };
    for (std::size_t at = 0;;) {
        const double gap = std::floor(std::log(std::max(rng.unit(), 1e-300)) / lq);
        const std::size_t next = gap >= static_cast<double>(s.size() - at) ? s.size() : at + static_cast<std::size_t>(gap);
        for (std::size_t i = at; i < next; ++i) out += upper(s[i]);
        if (next >= s.size()) break;
        const char c = upper(s[next]);
        const double k = rng.unit();
        if (k >= indel_share) {   // substitution
            char b;
            do b = kBases[rng.next() & 3]; while (b == c);
            out += b;
        } else if (k >= indel_share / 2) {   // insertion; below that, deletion
            out += kBases[rng.next() & 3];
            out += c;
        }
        at = next + 1;
    }
    return out;
}

// A fragment to simulate: `length` bases of `src` from `at` (a pair: its two ends, `length` the insert).
struct Job {
    std::size_t src;     // which sequence
    std::uint64_t at, length;
    std::uint8_t kind;   // 0 locus, 1 flank, 2 decoy
};

struct Tally {
    std::uint64_t n = 0, screened = 0;
    std::array<std::uint64_t, kMaxCut + 1> certain{};
};

// One simulated fragment: screened, the largest cut it is certain at by the first rules, what the score decides on.
struct Result {
    bool screened = false;
    std::uint8_t through = 0;
    scoring::Judged judged;
};

// Simulates and scores the jobs. With `score_only`, what the score decides on and no more: errors drawn by their
// gaps, no count of the cuts the first rules would call the fragment certain at, and a fragment without a private
// k-mer not weighed.
std::vector<Result> run(const std::vector<Job>& jobs, const std::vector<const std::string*>& srcs, std::uint64_t seed, const Class& c,
                        const scoring::LocusTable& t, const std::vector<std::uint64_t>& screen, std::size_t threads, bool score_only = false) {
    ReadOptions o;
    o.error_rate = c.error;
    const scoring::Costs costs = scoring::costs_for(c.error);
    std::vector<Result> res(jobs.size());
    std::atomic<std::size_t> next{0};
    const auto errors = score_only ? with_errors_by_gaps : with_errors;
    auto work = [&] {
        for (std::size_t i; (i = next++) < jobs.size();) {
            const Job& j = jobs[i];
            const std::string* src = srcs[j.src];
            Rng rng{seed ^ (0xD1B54A32D192ED03ULL * (i + 1))};
            Fragment f;
            f.name = "f";
            if (c.insert > 0) {
                f.paired = true;
                const std::uint64_t r = std::min<std::uint64_t>(c.length, j.length);
                f.seq1 = errors(src->substr(j.at, r), c.error, c.indel_share, rng);
                f.seq2 = errors(reverse_complement(src->substr(j.at + j.length - r, r)), c.error, c.indel_share, rng);
            } else {
                f.seq1 = errors(src->substr(j.at, j.length), c.error, c.indel_share, rng);
            }
            // the screen, as recruit applies it: a closed syncmer of the index (and not held out)
            bool hit = false;
            for (const std::string* s : {&f.seq1, &f.seq2})
                kmer::for_each_closed_syncmer(s->data(), s->size(), kK, kS, [&](std::uint64_t code, std::size_t, bool) {
                    if (hit || !std::binary_search(screen.begin(), screen.end(), code)) return;
                    const std::uint8_t* l = t.label.find(code);
                    hit = l != nullptr && (*l & 7) != kLabelN;
                });
            if (!hit) continue;
            res[i].screened = true;
            std::uint32_t m = 0;
            scoring::score_fragment(f, t, 0, o, costs.gap_cost, costs.flank_share, score_only ? nullptr : &m, &res[i].judged, score_only);
            res[i].through = static_cast<std::uint8_t>(m);
        }
    };
    std::vector<std::thread> pool;
    for (std::size_t k = 1; k < std::max<std::size_t>(1, threads); ++k) pool.emplace_back(work);
    work();
    for (auto& th : pool) th.join();
    return res;
}

// Each kind's tally by the first rules.
void tally_of(const std::vector<Job>& jobs, const std::vector<Result>& res, Tally tally[3]) {
    for (std::size_t i = 0; i < jobs.size(); ++i) {
        Tally& x = tally[jobs[i].kind];
        ++x.n;
        x.screened += res[i].screened;
        for (std::uint32_t m = 1; m <= res[i].through; ++m) ++x.certain[m];
    }
}

// An off-target stretch as fragments are laid over it: its length, the window in it and whether it is one.
struct Extent {
    std::uint64_t n, core_start, core_end;
    bool decoy;
};

// Off-target fragments from flanks (tiled) and windows with context (tiled with reads overlapping the window: fragments
// overlapping a window share its sites, so false certain fragments cluster by site, and the bound treats them as
// independent), about `want` from each kind of source. With `every_overlap`, long reads are laid over a window so
// that every overlap with it, from a k-mer at one end of the read to the whole read, comes up as often as
// sequencing gives it: one start drawn in each stride from where the read's last k-mer is the window's first to
// where its first is the window's last. A read that touches a homolog only at its end is judged on a short
// stretch, where two sites score high; a fixed grid a set distance into the window never lays one.
std::vector<Job> off_jobs(const std::vector<Extent>& off, const Class& c, Rng& rng, std::uint8_t flank_kind, std::uint8_t decoy_kind,
                          bool every_overlap = false) {
    std::uint64_t flank_bases = 0, decoy_bases = 0;
    for (const Extent& r : off) (r.decoy ? decoy_bases : flank_bases) += r.decoy ? r.core_end - r.core_start + (c.insert ? 0 : c.length) : r.n;
    const std::uint64_t want = c.insert ? kOffShort : kOffLong;
    const std::uint64_t stride_flank = std::max<std::uint64_t>(c.insert ? 50 : c.length / 2, flank_bases / want);
    const std::uint64_t stride_decoy = std::max<std::uint64_t>(c.insert ? 10 : 500, decoy_bases / want);
    std::vector<Job> jobs;
    for (std::size_t ri = 0; ri < off.size(); ++ri) {
        const Extent& r = off[ri];
        const std::uint64_t n = r.n;
        const std::uint8_t kind = r.decoy ? decoy_kind : flank_kind;
        if (c.insert) {
            // decoy windows: fragments with at least a read length inside the window, in its context too
            const std::uint64_t max_l = c.insert + 3 * c.insert_sd;
            const std::uint64_t lo = r.decoy && r.core_start + c.length > max_l ? r.core_start + c.length - max_l : 0;
            const std::uint64_t hi = r.decoy ? (r.core_end > c.length ? r.core_end - c.length : 0) : n;
            for (std::uint64_t a = lo; a < hi; a += r.decoy ? stride_decoy : stride_flank) {
                const std::uint64_t L = static_cast<std::uint64_t>(std::clamp(static_cast<double>(c.insert) + c.insert_sd * rng.normal(),
                                                                              static_cast<double>(c.length), static_cast<double>(max_l)));
                if (a + L > n) continue;
                if (r.decoy && std::min(a + L, r.core_end) < std::max(a, r.core_start) + c.length) continue;
                jobs.push_back({ri, a, L, kind});
            }
        } else if (n <= c.length) {
            if (n >= 2 * kK) jobs.push_back({ri, 0, n, kind});
        } else if (r.decoy && every_overlap) {
            const std::uint64_t lo = std::min(n - c.length, r.core_start + kK > c.length ? r.core_start + kK - c.length : 0);
            const std::uint64_t hi = std::min(n - c.length, std::max(lo, r.core_end > kK ? r.core_end - kK : 0));
            for (std::uint64_t cell = lo; cell <= hi; cell += stride_decoy) {
                const std::uint64_t width = std::min<std::uint64_t>(stride_decoy, hi - cell + 1);
                jobs.push_back({ri, cell + std::min<std::uint64_t>(width - 1, static_cast<std::uint64_t>(rng.unit() * static_cast<double>(width))), c.length, kind});
            }
        } else if (r.decoy) {   // reads ending at least 1 kb into the window, starting at least 1 kb before its end
            const std::uint64_t lo = std::min(n - c.length, r.core_start + 1000 > c.length ? r.core_start + 1000 - c.length : 0);
            const std::uint64_t hi = std::min(n - c.length, std::max(lo, r.core_end > 1000 ? r.core_end - 1000 : 0));
            for (std::uint64_t a = lo; a <= hi; a += stride_decoy) jobs.push_back({ri, a, c.length, kind});
        } else {
            for (std::uint64_t a = 0; a + c.length <= n; a += stride_flank) jobs.push_back({ri, a, c.length, kind});
        }
    }
    return jobs;
}

// FNV-1a: the same seed on every platform (std::hash differs between standard libraries).
std::uint64_t fnv1a(const std::string& s) {
    std::uint64_t h = 0xCBF29CE484222325ULL;
    for (const char c : s) h = (h ^ static_cast<unsigned char>(c)) * 0x100000001B3ULL;
    return h;
}

// Sites per k k-mers, with half a site and one k-mer stretch of prior so that no rate is 0.
double site_rate(double sites, double kmers) { return (sites + 0.5) / (kmers / static_cast<double>(kK) + 1.0); }

ScoreHistogram histogram(const std::vector<const scoring::Judged*>& js, const scoring::Rates& r) {
    std::map<std::int32_t, std::uint64_t> h;
    for (const scoring::Judged* j : js) {
        double s = 0;
        if (scoring::deciding_score(*j, r, s)) ++h[scoring::score_bin(s)];
    }
    ScoreHistogram out;
    for (const auto& [b, n] : h) {
        out.bin.push_back(b);
        out.count.push_back(n);
    }
    return out;
}

}  // namespace

std::uint64_t ScoreHistogram::at_least(std::int32_t b) const {
    std::uint64_t n = 0;
    for (std::size_t i = std::lower_bound(bin.begin(), bin.end(), b) - bin.begin(); i < bin.size(); ++i) n += count[i];
    return n;
}

std::string calibration_sample(const std::string& path) {
    const std::size_t h = path.find('#');
    return h == std::string::npos ? path : path.substr(0, h);
}

DecoyFilter::DecoyFilter(const std::vector<std::uint64_t>& codes, const std::vector<std::uint64_t>& decoys) {
    std::uint64_t bits = 1024;
    while (bits < 16 * (codes.size() + decoys.size())) bits <<= 1;
    bits_.assign(bits / 64, 0);
    mask_ = bits - 1;
    for (const auto* v : {&codes, &decoys})
        for (const std::uint64_t d : *v) {
            const std::uint64_t h = kmer::mix64(d) & mask_;
            bits_[h >> 6] |= 1ULL << (h & 63);
        }
}

bool DecoyFilter::maybe(std::uint64_t code) const {
    const std::uint64_t h = kmer::mix64(code) & mask_;
    return (bits_[h >> 6] >> (h & 63)) & 1;
}

namespace {

// Sorted positions, counted over ranges.
struct Marks {
    std::vector<std::uint32_t> at;
    std::size_t count(std::uint64_t a, std::uint64_t b) const {   // in [a, b)
        return static_cast<std::size_t>(std::lower_bound(at.begin(), at.end(), b) - std::lower_bound(at.begin(), at.end(), a));
    }
};

// Sorted, disjoint stretches of positions.
struct Runs {
    std::vector<std::pair<std::uint64_t, std::uint64_t>> r;
    void add(std::uint64_t a, std::uint64_t b) {   // in increasing order of a
        if (!r.empty() && a <= r.back().second) r.back().second = std::max(r.back().second, b);
        else r.emplace_back(a, b);
    }
    std::uint64_t covered(std::uint64_t a, std::uint64_t b) const {   // positions of [a, b) in a stretch
        std::uint64_t n = 0;
        auto it = std::upper_bound(r.begin(), r.end(), a, [](std::uint64_t v, const std::pair<std::uint64_t, std::uint64_t>& x) { return v < x.second; });
        for (; it != r.end() && it->first < b; ++it) n += std::min(b, it->second) - std::max(a, it->first);
        return n;
    }
};

constexpr std::uint8_t kNoKmer = 7, kScreenBit = 8;

// A block of a held-out genome as reads of it are labelled: per k-mer start its label (kNoKmer where no k-mer
// starts, kScreenBit on a screen k-mer), and where the private k-mers, the screen k-mers and the decoy sites
// (decoy k-mers k apart) lie; the starts without a k-mer; the stretches that read as a copy of the locus.
struct BlockLabels {
    std::vector<std::uint8_t> lab;
    Marks priv, screen, dsite;
    Runs none, copy;
};

// A mate read off a block's labels, [at, at + n) of it, with errors at rate e: each error leaves the k-mers over
// it unknown (one that makes a known k-mer by chance is not looked for). True when a screen k-mer is read whole.
bool read_labels(const std::string& seq, const BlockLabels& L, std::uint64_t at, std::uint64_t n, bool reverse, double e, Rng& rng,
                 std::vector<scoring::Kmer>& out) {
    out.clear();
    if (n < kK) return false;
    std::vector<std::uint64_t> err;   // error positions: the gaps between them are geometric
    const double lq = std::log(1.0 - e);
    auto gap = [&] { return std::floor(std::log(std::max(rng.unit(), 1e-300)) / lq); };
    for (double p = gap(); p < static_cast<double>(n); p += 1.0 + gap()) err.push_back(static_cast<std::uint64_t>(p));
    bool screened = false;
    std::size_t ei = 0;
    out.reserve(n - kK + 1);
    kmer::for_each_canonical_kmer(seq.data() + at, n, kK, [&](std::uint64_t code, std::size_t start, bool) {
        while (ei < err.size() && err[ei] < start) ++ei;
        const bool hit = ei < err.size() && err[ei] < start + kK;
        const std::uint8_t l = L.lab[at + start];
        screened = screened || (!hit && (l & kScreenBit));
        out.push_back({code, static_cast<std::uint32_t>(start), static_cast<std::uint8_t>(hit ? kLabelN : (l & 7)), false});
    });
    if (reverse) {
        std::reverse(out.begin(), out.end());
        for (scoring::Kmer& k : out) k.start = static_cast<std::uint32_t>(n - kK - k.start);
    }
    return screened;
}

double median(std::vector<double> v) {
    std::sort(v.begin(), v.end());
    const std::size_t n = v.size();
    return n == 0 ? 0.0 : n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

}  // namespace

GenomeScan scan_genome(const LocusIndex& x, const FlatMap<std::uint32_t>& code_index, const DecoyFilter& filter,
                       const std::vector<std::uint64_t>& decoy_counts, const std::string& genome,
                       const std::vector<std::string>& blocks, const std::vector<OffStretch>& stretches) {
    GenomeScan out;
    out.stats.genome = genome;
    for (const OffStretch& s : stretches) out.stats.windows += s.window;
    const std::size_t nc = x.codes.size(), nd = x.decoy.size();
    // the k-mer starts holding a locus k-mer (its index) or a decoy (its index after them): few among the genome's,
    // which the filter passes over
    struct Hit {
        std::uint32_t pos, id;
    };
    std::vector<std::vector<Hit>> hits(blocks.size());
    std::vector<BlockLabels> labels(blocks.size());
    std::vector<std::uint32_t> in_locus(nc, 0), in_decoy(nd, 0);
    for (std::size_t b = 0; b < blocks.size(); ++b) {
        const std::string& seq = blocks[b];
        BlockLabels& L = labels[b];
        out.stats.bases += seq.size();
        const std::size_t n = seq.size() >= kK ? seq.size() - kK + 1 : 0;
        L.lab.assign(n, kLabelN);
        std::size_t next = 0;   // the starts before this one are done
        auto no_kmers = [&](std::size_t from, std::size_t to) {
            if (to <= from) return;
            std::fill(L.lab.begin() + static_cast<std::ptrdiff_t>(from), L.lab.begin() + static_cast<std::ptrdiff_t>(to), kNoKmer);
            L.none.add(from, to);
        };
        kmer::for_each_canonical_kmer(seq.data(), seq.size(), kK, [&](std::uint64_t c, std::size_t start, bool) {
            no_kmers(next, start);
            next = start + 1;
            if (!filter.maybe(c)) return;
            if (const std::uint32_t* i = code_index.find(c)) {
                ++in_locus[*i];
                hits[b].push_back({static_cast<std::uint32_t>(start), *i});
                return;
            }
            const auto it = std::lower_bound(x.decoy.begin(), x.decoy.end(), c);
            if (it == x.decoy.end() || *it != c) return;
            const std::size_t j = static_cast<std::size_t>(it - x.decoy.begin());
            ++in_decoy[j];
            hits[b].push_back({static_cast<std::uint32_t>(start), static_cast<std::uint32_t>(nc + j)});
        });
        no_kmers(next, n);
    }
    // held out, the genome's own homolog k-mers change label: a locus k-mer whose every occurrence outside the locus
    // is in it would be private, a decoy only it holds unknown
    std::vector<char> to_private(nc, 0), to_unknown(nd, 0);
    out.homolog.assign((nc + nd + 63) / 64, 0);
    auto hold = [&](std::size_t id) {
        out.homolog[id >> 6] |= 1ULL << (id & 63);
        ++out.stats.homolog_kmers;
    };
    for (std::size_t i = 0; i < nc; ++i) {
        if (in_locus[i] == 0 || x.cls[i] == kPrivate) continue;
        hold(i);
        if (x.outside[i] == in_locus[i]) { to_private[i] = 1; ++out.stats.relabelled; }
    }
    for (std::size_t j = 0; j < nd; ++j) {
        if (in_decoy[j] == 0) continue;
        hold(nc + j);
        if (decoy_counts[j] == in_decoy[j]) { to_unknown[j] = 1; ++out.stats.unlabelled; }
    }
    const std::size_t copy_need = static_cast<std::size_t>(std::ceil(kCopyDensity * static_cast<double>(kCopyWindow)));
    for (std::size_t b = 0; b < blocks.size(); ++b) {
        BlockLabels& L = labels[b];
        bool open = false;
        std::uint64_t first = 0;
        std::vector<std::uint32_t> like;   // the starts that read as the locus's own sequence
        for (const Hit& h : hits[b]) {
            std::uint8_t lab;
            if (h.id < nc) {
                lab = to_private[h.id] ? scoring::kLabelP : x.cls[h.id];
                if (x.cls[h.id] != kRepeat && kmer::is_closed_syncmer(x.codes[h.id], kK, kS)) lab |= kScreenBit;
            } else {
                lab = to_unknown[h.id - nc] ? kLabelN : scoring::kLabelD;
            }
            const std::uint8_t l = lab & 7;
            if (l == scoring::kLabelP) L.priv.at.push_back(h.pos);
            if (lab & kScreenBit) L.screen.at.push_back(h.pos);
            if (l == scoring::kLabelD && (!open || h.pos >= first + kK)) {
                L.dsite.at.push_back(h.pos);
                open = true;
                first = h.pos;
            }
            if (l == scoring::kLabelP || l == scoring::kLabelS) like.push_back(h.pos);
            L.lab[h.pos] = lab;
        }
        std::vector<Hit>().swap(hits[b]);
        // a stretch as dense in the locus's k-mers as a copy of it is one: its reads are the locus's, not off-target
        // (windows of kCopyWindow starts, each from one of the locus's k-mers)
        for (std::size_t i = 0, j = 0; i < like.size(); ++i) {
            while (j < like.size() && like[j] < like[i] + kCopyWindow) ++j;
            if (j - i >= copy_need && like[i] + kCopyWindow <= L.lab.size()) L.copy.add(like[i], like[i] + kCopyWindow);
        }
        out.stats.copy_bases += L.copy.covered(0, L.lab.size());
    }
    std::vector<Extent> extents;
    for (const OffStretch& s : stretches) extents.push_back({s.end - s.start, s.core_start - s.start, s.core_end - s.start, s.window});
    const ReadOptions defaults;
    for (const Class& c : kClasses) {
        GenomeClassScan cs;
        Rng rng{0x5EEDULL ^ fnv1a(x.name + c.name + "#" + genome)};
        const std::vector<Job> jobs = off_jobs(extents, c, rng, 0, 0, true);
        const std::uint64_t seed = rng.next();
        const scoring::Costs costs = scoring::costs_for(c.error);
        const double clean = std::pow(1.0 - c.error, static_cast<double>(kK));   // a k-mer read without an error
        std::vector<scoring::Kmer> mate[2];
        for (std::size_t i = 0; i < jobs.size(); ++i) {
            const Job& j = jobs[i];
            const OffStretch& st = stretches[j.src];
            const BlockLabels& L = labels[st.block];
            const std::uint64_t at = st.start + j.at;
            if (L.copy.covered(at, at + j.length) > 0) { ++cs.copies; continue; }
            ++cs.fragments;
            // the mates' bases, and their k-mer starts
            const std::uint64_t r = c.insert ? std::min<std::uint64_t>(c.length, j.length) : j.length;
            const int mates = c.insert ? 2 : 1;
            const std::uint64_t from[2] = {at, at + j.length - r};
            if (r < kK) continue;
            std::size_t priv = 0, hits = 0;
            for (int m = 0; m < mates; ++m) {
                priv += L.priv.count(from[m], from[m] + r - kK + 1);
                hits += L.screen.count(from[m], from[m] + r - kK + 1);
            }
            if (hits == 0) continue;   // not screened: no screen k-mer to read
            Rng frng{seed ^ (0xD1B54A32D192ED03ULL * (i + 1))};
            if (priv == 0) {   // it cannot be called: counted when an error-free screen k-mer is likely enough
                if (frng.unit() >= 1.0 - std::pow(1.0 - clean, static_cast<double>(hits))) continue;
                ++cs.screened;
                for (int m = 0; m < mates; ++m) {
                    const std::uint64_t a = from[m], b = from[m] + r - kK + 1;
                    cs.kmers += static_cast<double>(b - a - L.none.covered(a, b));
                    cs.sites_d += static_cast<double>(L.dsite.count(a, b));
                }
                continue;
            }
            bool screened = false;
            for (int m = 0; m < mates; ++m) screened = read_labels(blocks[st.block], L, from[m], r, m == 1, c.error, frng, mate[m]) || screened;
            if (!screened) continue;
            ++cs.screened;
            const scoring::Judged jd = scoring::judge_labelled(mate[0], r, mates == 2 ? &mate[1] : nullptr, mates == 2 ? r : 0, c.error,
                                                               r >= scoring::kLongRead ? defaults.decoy_share_long : defaults.decoy_share,
                                                               costs.gap_cost, costs.flank_share);
            cs.with_private += jd.pooled.p > 0;
            cs.sites_p += jd.pooled.p;
            cs.sites_d += jd.pooled.d;
            cs.kmers += jd.pooled.n;
            if (jd.pooled_ok || jd.mate_ok[0] || jd.mate_ok[1]) cs.judged.push_back(jd);
        }
        out.classes.push_back(std::move(cs));
    }
    return out;
}

void calibrate(LocusIndex& x, const Graph& g, const std::vector<OffTarget>& off, const std::vector<GenomeScan>& genomes,
               const std::vector<std::int32_t>& owner, std::size_t threads) {
    scoring::LocusTable t = scoring::table_of(x);
    // samples and their haplotypes, numbered by first path; the reference path's sample is in the background and is
    // not held out
    std::map<std::string, std::vector<std::size_t>> samples;
    std::map<std::string, std::int32_t> sample_id;
    for (std::size_t p = 0; p < g.paths.size(); ++p) {
        const std::string s = calibration_sample(g.paths[p].name);
        samples[s].push_back(p);
        sample_id.emplace(s, static_cast<std::int32_t>(sample_id.size()));
    }
    samples.erase(calibration_sample(x.reference_path));
    std::vector<std::string> names;
    for (const auto& [s, ps] : samples) names.push_back(s);
    std::vector<std::string> held;
    for (std::size_t i = 0; i < std::min(kMaxHeldOut, names.size()); ++i) held.push_back(names[i * names.size() / std::min(kMaxHeldOut, names.size())]);
    // the held-out samples' haplotypes
    std::vector<std::string> seqs(g.paths.size());
    std::vector<const std::string*> path_srcs(g.paths.size()), off_srcs;
    for (std::size_t p = 0; p < g.paths.size(); ++p) path_srcs[p] = &seqs[p];
    for (const std::string& s : held)
        for (const std::size_t p : samples[s]) seqs[p] = spell_path_steps_sequence(g, g.paths[p].steps);
    std::vector<Extent> off_extents;
    for (const OffTarget& r : off) {
        off_srcs.push_back(&r.seq);
        off_extents.push_back({r.seq.size(), r.core_start, r.core_end, r.decoy});
    }
    // labels changed for a while, then put back
    std::vector<std::pair<std::uint8_t*, std::uint8_t>> saved;
    auto relabel = [&](std::uint64_t code, std::uint8_t to) {
        bool existing = false;
        std::uint8_t& l = t.label.slot(code, existing);
        saved.emplace_back(&l, l);
        l = to;
    };
    auto restore = [&] {
        for (auto it = saved.rbegin(); it != saved.rend(); ++it) *it->first = it->second;
        saved.clear();
    };

    x.calibration.clear();
    for (std::size_t ci = 0; ci < std::size(kClasses); ++ci) {
        const Class& c = kClasses[ci];
        Calibration cal;
        cal.read_class = c.name;
        cal.error_rate = c.error;
        cal.read_length = c.length;
        cal.insert = c.insert;
        cal.held_out = held.size();
        Rng rng{0x5EEDULL ^ fnv1a(x.name + c.name)};
        Tally tally[3];
        std::vector<scoring::Judged> locus_j, off_j;   // screened fragments
        auto keep = [](const std::vector<Result>& res, std::vector<scoring::Judged>& to) {
            for (const Result& r : res) if (r.screened) to.push_back(r.judged);
        };
        // held-out samples: their own k-mers read as new (N), one sample at a time
        for (const std::string& s : held) {
            std::vector<Job> jobs;
            for (const std::size_t p : samples[s]) {
                const std::string& seq = seqs[p];
                const std::uint64_t span = c.insert ? c.insert + 3 * c.insert_sd : c.length;
                if (seq.size() < 2 * kK) continue;
                const std::uint64_t reads_per = c.insert ? 2 : 1;
                const std::uint64_t len_one = std::min<std::uint64_t>(c.length, seq.size());
                const std::uint64_t n = std::max<std::uint64_t>(3, static_cast<std::uint64_t>(c.coverage * static_cast<double>(seq.size()) / static_cast<double>(len_one * reads_per)));
                for (std::uint64_t k = 0; k < n; ++k) {
                    std::uint64_t L = c.insert ? static_cast<std::uint64_t>(std::clamp(static_cast<double>(c.insert) + c.insert_sd * rng.normal(), static_cast<double>(c.length), static_cast<double>(span))) : len_one;
                    L = std::min<std::uint64_t>(L, seq.size());
                    jobs.push_back({p, rng.next() % (seq.size() - L + 1), L, 0});
                }
            }
            const std::int32_t id = sample_id.at(s);
            for (std::size_t i = 0; i < x.codes.size(); ++i)
                if (owner[i] == id) relabel(x.codes[i], kLabelN);
            const auto res = run(jobs, path_srcs, rng.next(), c, t, x.screen, threads);
            restore();
            tally_of(jobs, res, tally);
            keep(res, locus_j);
        }
        // off target: the flanks of the coordinate reference and the decoy windows, read as the index labels them
        // The first rules count the long reads over the decoy windows as they were first laid; the score takes them
        // laid at every overlap, so those are simulated apart, and it screens its own number.
        std::uint64_t windows_screened = 0;
        {
            const std::vector<Job> jobs = off_jobs(off_extents, c, rng, 1, 2);
            const auto res = run(jobs, off_srcs, rng.next(), c, t, x.screen, threads);
            tally_of(jobs, res, tally);
            for (std::size_t i = 0; i < jobs.size(); ++i) {
                if (!res[i].screened || (c.insert == 0 && jobs[i].kind == 2)) continue;
                off_j.push_back(res[i].judged);
                windows_screened += jobs[i].kind == 2;
            }
        }
        if (c.insert == 0) {
            std::vector<Extent> windows;
            std::vector<const std::string*> window_srcs;
            for (std::size_t i = 0; i < off_extents.size(); ++i)
                if (off_extents[i].decoy) {
                    windows.push_back(off_extents[i]);
                    window_srcs.push_back(off_srcs[i]);
                }
            Rng wrng{0x5EEDULL ^ fnv1a(x.name + c.name + "#windows")};
            const std::vector<Job> jobs = off_jobs(windows, c, wrng, 1, 2, true);
            const auto res = run(jobs, window_srcs, wrng.next(), c, t, x.screen, threads, true);
            for (const Result& r : res) {
                if (!r.screened) continue;
                off_j.push_back(r.judged);
                ++windows_screened;
            }
        }
        cal.locus = tally[0].n;
        cal.flank = tally[1].n;
        cal.decoy = tally[2].n;
        cal.flank_screened = tally[1].screened;
        cal.decoy_screened = tally[2].screened;
        cal.off_screened = tally[1].screened + windows_screened;
        cal.locus_certain = tally[0].certain;
        cal.flank_certain = tally[1].certain;
        cal.decoy_certain = tally[2].certain;
        // the score: site rates of the pooled stretches, then each fragment's deciding score
        double sums[2][3] = {};   // locus, simulated off-target: private sites, decoy sites, k-mers
        for (const auto* v : {&locus_j, &off_j})
            for (const scoring::Judged& j : *v) {
                if (j.pooled.n == 0) continue;
                double* to = sums[v == &locus_j ? 0 : 1];
                to[0] += j.pooled.p;
                to[1] += j.pooled.d;
                to[2] += j.pooled.n;
            }
        // off target from each calibration genome held out of the background, scanned: its own site rates, to which the
        // private sites sequencing errors make (the simulated off-target fragments', which hold no other) add
        const double r_errors = sums[1][2] > 0 ? sums[1][0] / (sums[1][2] / static_cast<double>(kK) + 1.0) : 0.0;
        std::vector<double> r_outs, q_outs;
        for (const GenomeScan& gs : genomes) {
            const GenomeClassScan& cs = gs.classes.at(ci);
            GenomeTally gt;
            gt.genome = gs.stats.genome;
            gt.fragments = cs.fragments;
            gt.screened = cs.screened;
            gt.with_private = cs.with_private;
            gt.copies = cs.copies;
            if (cs.screened > 0) {
                gt.r_out = site_rate(cs.sites_p, cs.kmers) + r_errors;
                gt.q_out = site_rate(cs.sites_d, cs.kmers);
                r_outs.push_back(gt.r_out);
                q_outs.push_back(gt.q_out);
            }
            cal.genome_fragments += gt.fragments;
            cal.genome_screened += gt.screened;
            cal.genomes.push_back(gt);
        }
        if (!locus_j.empty() && !(off_j.empty() && r_outs.empty())) {
            cal.r_in = site_rate(sums[0][0], sums[0][2]);
            cal.q_in = site_rate(sums[0][1], sums[0][2]);
            // the median genome's rates: one genome holding much of the locus elsewhere does not set them
            cal.r_out = r_outs.empty() ? site_rate(sums[1][0], sums[1][2]) : median(r_outs);
            cal.q_out = q_outs.empty() ? site_rate(sums[1][1], sums[1][2]) : median(q_outs);
            const scoring::Rates r{cal.r_in, cal.r_out, cal.q_in, cal.q_out};
            std::vector<const scoring::Judged*> lp, op, gp;
            for (const auto& j : locus_j) lp.push_back(&j);
            for (const auto& j : off_j) op.push_back(&j);
            for (const GenomeScan& gs : genomes)
                for (const auto& j : gs.classes.at(ci).judged) gp.push_back(&j);
            cal.locus_score = histogram(lp, r);
            cal.off_score = histogram(op, r);
            cal.genome_score = histogram(gp, r);
        }
        x.calibration.push_back(cal);
    }
}

double upper_bound(std::uint64_t events, std::uint64_t n) {
    if (n == 0) return 1.0;
    // the mean at which `events` or fewer happen with probability 0.05 (exact Poisson, by bisection; above a
    // few hundred events, Wilson-Hilferty's chi-square approximation)
    const double k = static_cast<double>(events);
    double u;
    if (events > 300) {
        const double a = k + 1.0;
        u = a * std::pow(1.0 - 1.0 / (9.0 * a) + 1.645 / (3.0 * std::sqrt(a)), 3.0);
    } else {
        auto cdf = [&](double lambda) {   // P(X <= events)
            double term = std::exp(-lambda), sum = term;
            for (std::uint64_t i = 1; i <= events; ++i) { term *= lambda / static_cast<double>(i); sum += term; }
            return sum;
        };
        double lo = k, hi = k + 10.0 + 10.0 * std::sqrt(k + 1.0);
        for (int it = 0; it < 100; ++it) {
            const double mid = 0.5 * (lo + hi);
            (cdf(mid) > 0.05 ? lo : hi) = mid;
        }
        u = hi;
    }
    return std::min(1.0, u / static_cast<double>(n));
}

std::uint32_t lowest_cut(const std::string& read_class) { return read_class == "ont" ? 3 : 2; }

const Calibration* calibration_for(const std::vector<Calibration>& cs, const std::string& read_class, double mean_read_length,
                                   double error_rate) {
    const std::string want = !read_class.empty() ? read_class : mean_read_length < 1000 ? "short" : error_rate <= 0.01 ? "hifi" : "ont";
    for (const Calibration& c : cs)
        if (c.read_class == want) return &c;
    return nullptr;
}

bool has_score(const Calibration& c) { return c.r_in > 0 && c.r_out > 0 && c.q_in > 0 && c.q_out > 0; }

Tier score_tier(const Calibration& c, double rate) {
    Tier t;
    t.rate = rate;
    // the two kinds of off-target fragment, simulated and scanned: each one's share is of its own screened fragments
    const std::uint64_t screened[2] = {c.off_screened, c.genome_screened};
    const double w[2] = {screened[0] ? 1.0 / static_cast<double>(screened[0]) : 0.0, screened[1] ? 1.0 / static_cast<double>(screened[1]) : 0.0};
    t.off_screened = screened[0] + screened[1];
    t.locus_total = c.locus;
    std::map<std::int32_t, std::array<std::uint64_t, 2>> off;   // every off-target fragment that can be called, by score bin
    std::uint64_t callable = 0;
    for (int k = 0; k < 2; ++k) {
        const ScoreHistogram& h = k == 0 ? c.off_score : c.genome_score;
        for (std::size_t i = 0; i < h.bin.size(); ++i) {
            off[h.bin[i]][k] += h.count[i];
            callable += h.count[i];
        }
    }
    auto share_of = [&](const std::array<std::uint64_t, 2>& n) { return static_cast<double>(n[0]) * w[0] + static_cast<double>(n[1]) * w[1]; };
    // observed: from the top down, the bins whose off-target fragments stay within the rate
    t.bin = kEveryBin;
    {
        double share = 0;
        for (auto it = off.rbegin(); it != off.rend(); ++it) {
            share += share_of(it->second);
            if (share > rate * (1.0 + 1e-9)) {
                t.bin = it->first + 1;
                break;
            }
        }
    }
    // Too few fragments to show the rate (none beyond a threshold would still leave the bound above it): the scores
    // are taken to go on beyond those seen as an exponential tail. Its scale is fitted from what makes a higher
    // score: one more private site adds ln(r_in / r_out) and is as rare as a private site is among the off-target
    // fragments (the share that can be called, at its 95% bound). The threshold is the lowest at which a tail of
    // that scale lies above every score seen.
    t.observed_bin = t.bin;
    double least = 0, site_share = 0;
    std::array<std::uint64_t, 2> all{};
    for (const auto& [bin, n] : off)
        for (int k = 0; k < 2; ++k) all[k] += n[k];
    for (int k = 0; k < 2; ++k)
        if (screened[k]) {
            least += upper_bound(0, screened[k]);
            site_share += upper_bound(all[k], screened[k]);
        }
    if (rate > 0 && rate < least && callable > 0 && c.r_in > c.r_out && site_share < 1.0) {
        t.tail_scale = std::log(c.r_in / c.r_out) / std::log(1.0 / site_share);
        t.tail_site_share = site_share;
        double share = 0, s = 0;
        bool any = false;
        for (auto it = off.rbegin(); it != off.rend(); ++it) {
            share += share_of(it->second);
            const double at = static_cast<double>(it->first + 1) * kScoreBin + t.tail_scale * std::log(share / rate);
            if (!any || at > s) {
                s = at;
                t.tail_from = it->first;
                any = true;
            }
        }
        const std::int32_t bin = scoring::score_bin(s) + 1;
        if (t.bin == kEveryBin || bin > t.bin) {
            t.bin = bin;
            t.extrapolated = true;
        }
    }
    std::array<std::uint64_t, 2> reached{};
    for (auto it = off.lower_bound(t.bin); it != off.end(); ++it)
        for (int k = 0; k < 2; ++k) reached[k] += it->second[k];
    t.off = reached[0] + reached[1];
    t.share = share_of(reached);
    t.share_bound = t.off_screened == 0 ? 1.0 : 0.0;
    for (int k = 0; k < 2; ++k)
        if (screened[k]) t.share_bound += upper_bound(reached[k], screened[k]);
    t.share_bound = std::min(1.0, t.share_bound);
    t.locus = c.locus_score.at_least(t.bin);
    t.shown = t.share_bound <= rate;
    return t;
}

Cut calibrated_cut(const std::vector<Calibration>& cs, const std::string& read_class, double mean_read_length, double error_rate,
                   double rate) {
    const std::string want = !read_class.empty() ? read_class : mean_read_length < 1000 ? "short" : error_rate <= 0.01 ? "hifi" : "ont";
    Cut none;
    none.floor = lowest_cut(want);
    for (const Calibration& c : cs) {
        if (c.read_class != want) continue;
        const std::uint64_t screened = c.flank_screened + c.decoy_screened;
        Cut cut = none;
        cut.by = &c;
        if (screened == 0) return cut;   // nothing measured: no calibrated cut
        cut.min_sites = kMaxCut;
        cut.met = false;
        for (std::uint32_t m = cut.floor; m <= kMaxCut; ++m) {
            const std::uint64_t fc = c.flank_certain[m] + c.decoy_certain[m];
            if (static_cast<double>(fc) <= rate * static_cast<double>(screened)) {
                cut.min_sites = m;
                cut.met = true;
                break;
            }
        }
        cut.shown = upper_bound(c.flank_certain[cut.min_sites] + c.decoy_certain[cut.min_sites], screened) <= rate;
        return cut;
    }
    return none;
}

}  // namespace panvar::recruit
