#include "panvar/recruit_calibration.hpp"

#include "panvar/graph_utils.hpp"
#include "panvar/kmer.hpp"
#include "panvar/read_input.hpp"
#include "panvar/recruit_scoring.hpp"

#include <algorithm>
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

// A fragment to simulate: `length` bases of `src` from `at` (a pair: its two ends, `length` the insert).
struct Job {
    const std::string* src;
    std::uint64_t at, length;
    std::uint8_t kind;   // 0 locus, 1 flank, 2 decoy, 3 a calibration genome's
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

// Simulates and scores the jobs.
std::vector<Result> run(const std::vector<Job>& jobs, std::uint64_t seed, const Class& c, const scoring::LocusTable& t,
                        const std::vector<std::uint64_t>& screen, std::size_t threads) {
    ReadOptions o;
    o.error_rate = c.error;
    const scoring::Costs costs = scoring::costs_for(c.error);
    std::vector<Result> res(jobs.size());
    std::atomic<std::size_t> next{0};
    auto work = [&] {
        for (std::size_t i; (i = next++) < jobs.size();) {
            const Job& j = jobs[i];
            Rng rng{seed ^ (0xD1B54A32D192ED03ULL * (i + 1))};
            Fragment f;
            f.name = "f";
            if (c.insert > 0) {
                f.paired = true;
                const std::uint64_t r = std::min<std::uint64_t>(c.length, j.length);
                f.seq1 = with_errors(j.src->substr(j.at, r), c.error, c.indel_share, rng);
                f.seq2 = with_errors(reverse_complement(j.src->substr(j.at + j.length - r, r)), c.error, c.indel_share, rng);
            } else {
                f.seq1 = with_errors(j.src->substr(j.at, j.length), c.error, c.indel_share, rng);
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
            scoring::score_fragment(f, t, 0, o, costs.gap_cost, costs.flank_share, &m, &res[i].judged);
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

// Off-target fragments from flanks (tiled) and windows with context (tiled with reads overlapping the window: fragments
// overlapping a window share its sites, so false certain fragments cluster by site, and the bound treats them as
// independent), about `want` from each kind of source.
std::vector<Job> off_jobs(const std::vector<OffTarget>& off, const Class& c, Rng& rng, std::uint8_t flank_kind, std::uint8_t decoy_kind) {
    std::uint64_t flank_bases = 0, decoy_bases = 0;
    for (const OffTarget& r : off) (r.decoy ? decoy_bases : flank_bases) += r.decoy ? r.core_end - r.core_start + (c.insert ? 0 : c.length) : r.seq.size();
    const std::uint64_t want = c.insert ? kOffShort : kOffLong;
    const std::uint64_t stride_flank = std::max<std::uint64_t>(c.insert ? 50 : c.length / 2, flank_bases / want);
    const std::uint64_t stride_decoy = std::max<std::uint64_t>(c.insert ? 10 : 500, decoy_bases / want);
    std::vector<Job> jobs;
    for (const OffTarget& r : off) {
        const std::uint64_t n = r.seq.size();
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
                jobs.push_back({&r.seq, a, L, kind});
            }
        } else if (n <= c.length) {
            if (n >= 2 * kK) jobs.push_back({&r.seq, 0, n, kind});
        } else if (r.decoy) {   // reads ending at least 1 kb into the window, starting at least 1 kb before its end
            const std::uint64_t lo = std::min(n - c.length, r.core_start + 1000 > c.length ? r.core_start + 1000 - c.length : 0);
            const std::uint64_t hi = std::min(n - c.length, std::max(lo, r.core_end > 1000 ? r.core_end - 1000 : 0));
            for (std::uint64_t a = lo; a <= hi; a += stride_decoy) jobs.push_back({&r.seq, a, c.length, kind});
        } else {
            for (std::uint64_t a = 0; a + c.length <= n; a += stride_flank) jobs.push_back({&r.seq, a, c.length, kind});
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

std::string sample_of(const std::string& path) {
    const std::size_t h = path.find('#');
    return h == std::string::npos ? path : path.substr(0, h);
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

void calibrate(LocusIndex& x, const Graph& g, const std::vector<OffTarget>& off, const std::vector<GenomeOffTarget>& genomes,
               std::size_t threads) {
    scoring::LocusTable t = scoring::table_of(x);
    // samples and their haplotypes; the reference path's sample is in the background and is not held out
    std::map<std::string, std::vector<std::size_t>> samples;
    for (std::size_t p = 0; p < g.paths.size(); ++p) samples[sample_of(g.paths[p].name)].push_back(p);
    samples.erase(sample_of(x.reference_path));
    std::vector<std::string> names;
    for (const auto& [s, ps] : samples) names.push_back(s);
    std::vector<std::string> held;
    for (std::size_t i = 0; i < std::min(kMaxHeldOut, names.size()); ++i) held.push_back(names[i * names.size() / std::min(kMaxHeldOut, names.size())]);
    // the sample each k-mer occurs in (-1: several)
    std::vector<std::string> seqs(g.paths.size());
    std::vector<std::int32_t> owner(x.codes.size(), -2);
    std::map<std::string, std::int32_t> sample_id;
    for (std::size_t p = 0; p < g.paths.size(); ++p) {
        seqs[p] = spell_path_steps_sequence(g, g.paths[p].steps);
        const std::int32_t id = sample_id.emplace(sample_of(g.paths[p].name), static_cast<std::int32_t>(sample_id.size())).first->second;
        kmer::for_each_canonical_kmer(seqs[p].data(), seqs[p].size(), kK, [&](std::uint64_t c, std::size_t, bool) {
            const std::int64_t i = x.find(c);
            if (i < 0) return;
            std::int32_t& w = owner[static_cast<std::size_t>(i)];
            w = w == -2 || w == id ? id : -1;
        });
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
    for (const Class& c : kClasses) {
        Calibration cal;
        cal.read_class = c.name;
        cal.error_rate = c.error;
        cal.read_length = c.length;
        cal.insert = c.insert;
        cal.held_out = held.size();
        Rng rng{0x5EEDULL ^ fnv1a(x.name + c.name)};
        Tally tally[3];
        std::vector<scoring::Judged> locus_j, off_j, genome_j;   // screened fragments
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
                    jobs.push_back({&seq, rng.next() % (seq.size() - L + 1), L, 0});
                }
            }
            const std::int32_t id = sample_id.at(s);
            for (std::size_t i = 0; i < x.codes.size(); ++i)
                if (owner[i] == id) relabel(x.codes[i], kLabelN);
            const auto res = run(jobs, rng.next(), c, t, x.screen, threads);
            restore();
            tally_of(jobs, res, tally);
            keep(res, locus_j);
        }
        // off target: the flanks of the coordinate reference and the decoy windows, read as the index labels them
        {
            const std::vector<Job> jobs = off_jobs(off, c, rng, 1, 2);
            const auto res = run(jobs, rng.next(), c, t, x.screen, threads);
            tally_of(jobs, res, tally);
            keep(res, off_j);
        }
        // off target from each calibration genome held out of the background: the locus k-mers only it holds outside
        // the locus read as private, the decoy k-mers only it holds as unknown, as a new sample's paralogs would
        for (const GenomeOffTarget& go : genomes) {
            GenomeTally gt;
            gt.genome = x.calibration_genomes.at(go.genome).genome;
            Rng grng{0x5EEDULL ^ fnv1a(x.name + c.name + "#" + gt.genome)};
            const std::vector<Job> jobs = off_jobs(go.off, c, grng, 3, 3);
            for (const std::uint32_t i : go.relabel)
                relabel(x.codes[i], static_cast<std::uint8_t>(scoring::kLabelP | (x.mult[i] > 1 ? scoring::kMultiBit : 0) |
                                                              (x.backbone[i] ? scoring::kBackboneBit : 0)));
            for (const std::uint64_t d : go.unlabel) relabel(d, kLabelN);
            const auto res = run(jobs, grng.next(), c, t, x.screen, threads);
            restore();
            gt.fragments = jobs.size();
            for (const Result& r : res) {
                if (!r.screened) continue;
                ++gt.screened;
                gt.with_private += r.judged.pooled.p > 0;
                genome_j.push_back(r.judged);
            }
            cal.genome_fragments += gt.fragments;
            cal.genome_screened += gt.screened;
            cal.genomes.push_back(gt);
        }
        cal.locus = tally[0].n;
        cal.flank = tally[1].n;
        cal.decoy = tally[2].n;
        cal.flank_screened = tally[1].screened;
        cal.decoy_screened = tally[2].screened;
        cal.locus_certain = tally[0].certain;
        cal.flank_certain = tally[1].certain;
        cal.decoy_certain = tally[2].certain;
        // the score: site rates of the pooled stretches, then each fragment's deciding score
        double sums[2][3] = {};   // locus, off-target: private sites, decoy sites, k-mers
        for (const auto* v : {&locus_j, &off_j, &genome_j})
            for (const scoring::Judged& j : *v) {
                if (j.pooled.n == 0) continue;
                double* to = sums[v == &locus_j ? 0 : 1];
                to[0] += j.pooled.p;
                to[1] += j.pooled.d;
                to[2] += j.pooled.n;
            }
        if (!locus_j.empty() && !(off_j.empty() && genome_j.empty())) {
            cal.r_in = site_rate(sums[0][0], sums[0][2]);
            cal.q_in = site_rate(sums[0][1], sums[0][2]);
            cal.r_out = site_rate(sums[1][0], sums[1][2]);
            cal.q_out = site_rate(sums[1][1], sums[1][2]);
            const scoring::Rates r{cal.r_in, cal.r_out, cal.q_in, cal.q_out};
            auto ptrs = [](const std::vector<scoring::Judged>& v) {
                std::vector<const scoring::Judged*> p;
                for (const auto& j : v) p.push_back(&j);
                return p;
            };
            cal.locus_score = histogram(ptrs(locus_j), r);
            cal.off_score = histogram(ptrs(off_j), r);
            cal.genome_score = histogram(ptrs(genome_j), r);
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

double Tier::bound() const { return upper_bound(off, off_screened); }

Tier score_tier(const Calibration& c, double rate) {
    Tier t;
    t.rate = rate;
    t.off_screened = c.flank_screened + c.decoy_screened + c.genome_screened;
    t.locus_total = c.locus;
    std::map<std::int32_t, std::uint64_t> off;   // every off-target fragment that can be called, by score bin
    for (const ScoreHistogram* h : {&c.off_score, &c.genome_score})
        for (std::size_t i = 0; i < h->bin.size(); ++i) off[h->bin[i]] += h->count[i];
    // from the top down, the bins whose off-target fragments stay within the rate
    const double allowed = rate * static_cast<double>(t.off_screened);
    std::uint64_t reached = 0;
    t.bin = kEveryBin;
    for (auto it = off.rbegin(); it != off.rend(); ++it) {
        if (static_cast<double>(reached + it->second) > allowed) {
            t.bin = it->first + 1;
            break;
        }
        reached += it->second;
    }
    t.off = reached;
    t.locus = c.locus_score.at_least(t.bin);
    t.shown = t.bound() <= rate;
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
