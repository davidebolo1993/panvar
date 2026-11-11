#include "panvar/recruit_reads.hpp"

#include "panvar/flat_map.hpp"
#include "panvar/read_input.hpp"
#include "panvar/kmer.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <exception>
#include <memory>
#include <stdexcept>
#include <thread>

namespace panvar::recruit {

namespace {

constexpr std::size_t kK = 31, kS = 11;
constexpr std::uint8_t kLabelP = 0, kLabelS = 1, kLabelR = 2, kLabelD = 3, kLabelN = 4, kMultiBit = 8;

// One locus, as reads see it: each k-mer's label (P/S/R, D for decoys), with a bit for multi-copy P.
struct LocusTable {
    std::string name;
    FlatMap<std::uint8_t> label;
};

struct Row {
    std::string fragment;
    std::uint32_t locus = 0;
    Evidence ev;
    std::string seg[2];          // each mate's locus segment in read bases, "." when it has none
    const char* mates = "";      // each mate's status: in, part (a segment covering part of it), out
    const char* fate = "";
    const char* reason = "";
};

struct Kmer {
    std::uint64_t code;
    std::uint32_t start;
    std::uint8_t lab;
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
        if (it != nullptr) {
            lab = *it & 7;
            if (lab == kLabelP && (*it & kMultiBit)) ++ev.multi;
        }
        switch (lab) {
            case kLabelP: ++ev.n_p; break;
            case kLabelS: ++ev.n_s; break;
            case kLabelR: ++ev.n_r; break;
            case kLabelD: ++ev.n_d; break;
            default: ++ev.n_n;
        }
        v.push_back({c, static_cast<std::uint32_t>(start), lab});
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
        if (seg.first <= seg.second) {
            r.seg[m] = std::to_string(v[seg.first].start) + "-" + std::to_string(v[seg.second].start + kK);
            status[m] = ev[m].seg_kmers * 10 >= ev[m].kmers * 9 ? "in" : "part";
        } else {
            r.seg[m] = ".";
        }
    }
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

ReadStats recruit_reads(const ReadOptions& o) {
    auto log = [&](const std::string& s) { if (o.log) o.log(s); };
    if (o.indexes.empty()) throw std::runtime_error("recruit: no locus index");
    if (o.indexes.size() > 64) throw std::runtime_error("recruit: at most 64 loci in one run");
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
        bool existing = false;
        for (std::size_t i = 0; i < x.codes.size(); ++i)
            t.label.slot(x.codes[i], existing) = static_cast<std::uint8_t>(x.cls[i] | (x.cls[i] == kPrivate && x.mult[i] > 1 ? kMultiBit : 0));
        for (const std::uint64_t d : x.decoy) {
            std::uint8_t& v = t.label.slot(d, existing);
            if (!existing) v = kLabelD;
        }
        for (const std::uint64_t s : x.screen) screen_items.emplace_back(s, static_cast<std::uint32_t>(li));
        log(x.name + ": " + std::to_string(x.codes.size()) + " k-mers, " + std::to_string(x.decoy.size()) + " decoys, " +
            std::to_string(x.screen.size()) + " screen syncmers");
    }
    FlatMap<std::uint64_t> screen;   // syncmer -> loci it belongs to
    screen.reserve(screen_items.size());
    for (const auto& [code, li] : screen_items) {
        bool existing = false;
        screen.slot(code, existing) |= 1ULL << li;
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
    io.screened = screened;
    io.log = o.log;
    std::size_t next_input = 0;
    std::unique_ptr<FragmentSource> source;
    ReadStats st;
    auto close_source = [&] {
        if (!source) return;
        const InputStats& s = source->stats();
        st.fragments += s.unscreened + s.dropped;   // passed over before scoring, unscreened
        log(source->describe() + ": " + std::to_string(s.records) + " records, " + std::to_string(s.fragments + s.unscreened + s.dropped) +
            " fragments (" + std::to_string(s.pairs + s.unscreened) + " pairs)" +
            (s.secondary + s.supplementary + s.qc_fail > 0 ? "; skipped " + std::to_string(s.secondary) + " secondary, " +
                 std::to_string(s.supplementary) + " supplementary, " + std::to_string(s.qc_fail) + " QC-failed" : std::string()) +
            (s.fetched > 0 ? "; " + std::to_string(s.fetched) + " mates fetched" : std::string()) +
            (s.unmatched > 0 ? "; " + std::to_string(s.unmatched) + " reads without their mate" : std::string()));
        source.reset();
    };
    std::FILE* out = std::fopen(o.out.c_str(), "w");
    if (out == nullptr) throw std::runtime_error("recruit: cannot write " + o.out);
    std::fprintf(out, "fragment\tlocus\tfate\treason\tmates\tseg1\tseg2\tsites_p\tsites_d\tdecoy_out\truns_n\tseg_kmers\tn_p\tn_s\tn_r\tn_d\tn_n\tkmers\twindows\tbases\tmulti_p\n");
    // a run of k-mers in neither table costs this much per k-mer past its first k: half the expected
    // score of a locus stretch read at the error rate (error-free k-mers +1, the rest at this cost)
    const double q = std::pow(std::max(0.0, 1.0 - o.error_rate - 0.001), static_cast<double>(kK));
    const double gap_cost = q >= 1.0 ? 1e9 : 0.5 * q / (1.0 - q);
    const double flank_share = 0.25 * q;   // locus k-mers beyond a segment's end below this: a flank

    const std::size_t threads = std::max<std::size_t>(1, o.threads);
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
    std::vector<std::vector<Row>> rows;
    read_batch(batch);
    while (!batch.empty()) {
        rows.assign(batch.size(), {});
        std::atomic<std::size_t> cursor{0};
        std::vector<std::thread> pool;
        for (std::size_t t = 0; t < threads; ++t)
            pool.emplace_back([&] {
                for (std::size_t i; (i = cursor++) < batch.size();) {
                    const Fragment& fr = batch[i];
                    std::uint64_t hit = 0;
                    for (const std::string* s : {&fr.seq1, &fr.seq2})
                        kmer::for_each_closed_syncmer(s->data(), s->size(), kK, kS, [&](std::uint64_t c, std::size_t, bool) {
                            if (const std::uint64_t* m = screen.find(c)) hit |= *m;
                        });
                    for (std::uint32_t li = 0; hit != 0; ++li, hit >>= 1)
                        if (hit & 1) rows[i].push_back(score_fragment(fr, loci[li], li, o, gap_cost, flank_share));
                }
            });
        std::exception_ptr failed;
        try { read_batch(next); } catch (...) { failed = std::current_exception(); }
        for (auto& t : pool) t.join();
        if (failed) std::rethrow_exception(failed);
        for (std::size_t i = 0; i < batch.size(); ++i) {
            ++st.fragments;
            if (rows[i].empty()) continue;
            ++st.screened;
            for (const Row& r : rows[i]) {
                const Evidence& v = r.ev;
                std::fprintf(out, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\n", r.fragment.c_str(),
                             loci[r.locus].name.c_str(), r.fate, r.reason, r.mates, r.seg[0].c_str(), r.seg[1].c_str(), v.sites_p, v.sites_d,
                             v.decoy_out, v.runs_n, v.seg_kmers, v.n_p, v.n_s, v.n_r, v.n_d, v.n_n, v.kmers, v.windows, v.bases, v.multi);
                if (r.fate[0] == 'c') ++st.certain;
                else if (r.fate[0] == 'a') ++st.ambiguous;
                else ++st.unrecruited;
            }
        }
        batch.swap(next);
    }
    close_source();
    if (std::fclose(out) != 0) throw std::runtime_error("recruit: cannot write " + o.out);
    return st;
}

}  // namespace panvar::recruit
