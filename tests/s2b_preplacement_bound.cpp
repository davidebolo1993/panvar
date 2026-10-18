// S2-B (MEASUREMENT ONLY): an admissible PRE-PLACEMENT upper bound on the implemented fragment likelihood, its
// independent validation on exhaustive fixtures, and its practical usefulness on real union dumps.
// The derivation is tests/analysis/s2b/DERIVATION.md; every section number below refers to it.
//
//   s2b_preplacement_bound fixtures
//   s2b_preplacement_bound measure <union_evidence.bin> <reads_1.fq> <reads_2.fq> <threads> <out.tsv> [construction_only [max_haplotypes]]
//
// The bound reads only: the read bases, the haplotype bytes (through an exact substring index), the library, and
// sequence-independent recursion tables. It never calls gapped_emission, the prefilter or any read-window DP -- except
// the LAZY variant, which first computes EXACT emissions at every qualifying anchor of the early incumbent's haplotypes
// (its own cost, timed and reported apart: production's early evaluation covers only seeded starts) and reuses them,
// byte-verified, wherever another haplotype has the identical window. Dumped masses and scores are only the evaluator.
#include "panvar/fragment_likelihood.hpp"
#include "panvar/genotype_model.hpp"

#include <sys/resource.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {
using panvar::LibraryModel;
using panvar::ReadPair;
constexpr double kNegInf = -std::numeric_limits<double>::infinity();
using Clock = std::chrono::steady_clock;
double since(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }
double peak_rss_gib() {
    rusage u{};
    getrusage(RUSAGE_SELF, &u);
    return static_cast<double>(u.ru_maxrss) / 1073741824.0;   // bytes on macOS
}
double lae(double a, double b) { return panvar::log_add_exp(a, b); }
std::string revcomp(const std::string& s) {
    std::string out(s.rbegin(), s.rend());
    for (char& c : out) c = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : c == 'T' ? 'A' : 'N';
    return out;
}

// ---- section 2: the per-edit shape ceiling W_e(L, a) and the normaliser Z(L, a) --------------------------------------
// The reference recursion's shape mode (gap_recursion_logspace, emitting = false), with the M step's mismatch branch
// weighted `mismatch_weight` (epsilon for Z, epsilon/3 for W) and the terminal mass kept per edit count.
std::vector<double> shape_by_edits(std::size_t L, std::size_t available, std::size_t band, const LibraryModel& lib,
                                   double mismatch_weight) {
    const std::size_t delta = lib.max_gap_bp, width = 2 * delta + 1, budget = band + 1;
    const double log_continue = lib.gap_open > 0.0 ? std::log(1.0 - 2.0 * lib.gap_open) : 0.0;
    const double log_open = lib.gap_open > 0.0 ? std::log(lib.gap_open) : kNegInf;
    const double log_extend = lib.gap_extend > 0.0 ? std::log(lib.gap_extend) : kNegInf;
    const double log_close = lib.gap_extend < 1.0 ? std::log(1.0 - lib.gap_extend) : kNegInf;
    const double log_insert_emission = mismatch_weight == lib.substitution_error ? 0.0 : std::log(0.25);
    enum { M = 0, I = 1, D = 2, S = 3 };
    const auto at = [&](std::size_t d, std::size_t st, std::size_t e) { return (d * S + st) * budget + e; };
    std::vector<double> prev(width * S * budget, kNegInf), cur(width * S * budget, kNegInf);
    prev[at(delta, M, 0)] = 0.0;
    for (std::size_t i = 1; i <= L; ++i) {
        std::fill(cur.begin(), cur.end(), kNegInf);
        for (std::size_t d = 0; d < width; ++d) {
            const long k = static_cast<long>(d) - static_cast<long>(delta);
            const long consumed_before = static_cast<long>(i) - 1 + k;
            for (std::size_t e = 0; e < budget; ++e) {
                if (consumed_before >= 0 && static_cast<std::size_t>(consumed_before) < available) {
                    const double w[2] = {std::log(1.0 - lib.substitution_error), std::log(mismatch_weight)};
                    const std::size_t cost[2] = {0, 1};
                    for (int b = 0; b < 2; ++b) {
                        if (e < cost[b]) continue;
                        const std::size_t from = e - cost[b];
                        double sum = kNegInf;
                        sum = lae(sum, prev[at(d, M, from)] + log_continue);
                        sum = lae(sum, prev[at(d, I, from)] + log_close);
                        sum = lae(sum, prev[at(d, D, from)] + log_close);
                        cur[at(d, M, e)] = lae(cur[at(d, M, e)], sum + w[b]);
                    }
                }
                if (d + 1 < width && e >= 1) {
                    const double opened = prev[at(d + 1, M, e - 1)] + log_open;
                    const double extended = prev[at(d + 1, I, e - 1)] + log_extend;
                    cur[at(d, I, e)] = lae(cur[at(d, I, e)], lae(opened, extended) + log_insert_emission);
                }
            }
        }
        for (std::size_t d = 1; d < width; ++d) {
            const long k = static_cast<long>(d) - static_cast<long>(delta);
            const long consumed_before = static_cast<long>(i) + k - 1;
            if (consumed_before < 0 || static_cast<std::size_t>(consumed_before) >= available) continue;
            for (std::size_t e = 1; e < budget; ++e) {
                const double opened = cur[at(d - 1, M, e - 1)] + log_open;
                const double extended = cur[at(d - 1, D, e - 1)] + log_extend;
                cur[at(d, D, e)] = lae(cur[at(d, D, e)], lae(opened, extended));
            }
        }
        prev.swap(cur);
    }
    std::vector<double> by_edit(budget, kNegInf);
    const std::size_t floor_bp = L - std::min(L, delta);
    for (std::size_t d = 0; d < width; ++d) {
        const long k = static_cast<long>(d) - static_cast<long>(delta);
        const long consumed = static_cast<long>(L) + k;
        if (consumed < 0 || static_cast<std::size_t>(consumed) > available || static_cast<std::size_t>(consumed) < floor_bp) continue;
        for (std::size_t st = 0; st < S; ++st) {
            for (std::size_t e = 0; e < budget; ++e) by_edit[e] = lae(by_edit[e], prev[at(d, st, e)]);
        }
    }
    return by_edit;
}

// Per (L, a = min(available, L + delta)): log of sum_{e >= d} W_e / Z for d = 0..band (index band+1: -inf).
struct ShapeTables {
    const LibraryModel* lib = nullptr;
    std::size_t band = 0;
    std::map<std::pair<std::size_t, std::size_t>, std::vector<double>> tail;
    std::mutex guard;
    double max_z_mismatch = 0.0;   // |log Z (ours) - log_gap_shape_mass| over every table built
    const std::vector<double>& get(std::size_t L, std::size_t available) {
        const std::size_t a = std::min(available, L + lib->max_gap_bp);
        std::lock_guard<std::mutex> lock(guard);
        auto found = tail.find({L, a});
        if (found != tail.end()) return found->second;
        const std::vector<double> z_edits = shape_by_edits(L, a, band, *lib, lib->substitution_error);
        const std::vector<double> w_edits = shape_by_edits(L, a, band, *lib, lib->substitution_error / 3.0);
        double log_z = kNegInf;
        for (double x : z_edits) log_z = lae(log_z, x);
        const double reference = panvar::log_gap_shape_mass(L, a, band, *lib);
        max_z_mismatch = std::max(max_z_mismatch, std::abs(log_z - reference) / std::max(1.0, std::abs(reference)));
        std::vector<double> t(band + 2, kNegInf);
        for (std::size_t d = band + 1; d-- > 0;) t[d] = lae(t[d + 1], w_edits[d] - reference);
        return tail.emplace(std::make_pair(L, a), t).first->second;
    }
};

// ---- section 3: an exact substring index over the whole haplotype ----------------------------------------------------
int code_of(char c) { return c == 'A' ? 0 : c == 'C' ? 1 : c == 'G' ? 2 : c == 'T' ? 3 : -1; }
struct HaplotypeIndex {
    const std::string* h = nullptr;
    std::size_t q = 0;
    std::vector<std::uint64_t> entries;   // (q-gram code << 32) | position, every position, sorted
    void build(const std::string& seq, std::size_t qq) {
        h = &seq;
        q = qq;
        entries.clear();
        if (seq.size() < q) return;
        entries.reserve(seq.size());
        std::uint64_t code = 0;
        std::size_t valid = 0;
        const std::uint64_t mask = q >= 32 ? ~0ULL : ((1ULL << (2 * q)) - 1);
        for (std::size_t i = 0; i < seq.size(); ++i) {
            const int c = code_of(seq[i]);
            if (c < 0) { valid = 0; code = 0; continue; }
            code = ((code << 2) | static_cast<std::uint64_t>(c)) & mask;
            if (++valid >= q) entries.push_back((code << 32) | static_cast<std::uint64_t>(i + 1 - q));
        }
        std::sort(entries.begin(), entries.end());
    }
    // Every exact occurrence of piece = read[off, off+len): appended to `out`.
    void occurrences(const std::string& read, std::size_t off, std::size_t len, std::vector<std::uint32_t>& out) const {
        if (len < q || h == nullptr || h->size() < len) return;
        std::uint64_t code = 0;
        for (std::size_t i = 0; i < q; ++i) {
            const int c = code_of(read[off + i]);
            if (c < 0) return;   // a non-ACGT piece never occurs (and the scorer refuses non-ACGT reads)
            code = (code << 2) | static_cast<std::uint64_t>(c);
        }
        auto it = std::lower_bound(entries.begin(), entries.end(), code << 32);
        for (; it != entries.end() && (*it >> 32) == code; ++it) {
            const std::size_t x = static_cast<std::size_t>(*it & 0xffffffffULL);
            if (x + len <= h->size() && std::memcmp(h->data() + x, read.data() + off, len) == 0) {
                out.push_back(static_cast<std::uint32_t>(x));
            }
        }
    }
};

struct Variant {   // deliberate defects for the validator (section: independent validation)
    bool no_shift = false;          // K forced to 0: under-counts admissible starts/edits under indels
    bool drop_second_arm = false;   // omits one library orientation's contribution
    bool first_hit_only = false;    // one occurrence per piece, as a mapper reporting a single hit would
};

// The qualifying anchors of one read on one haplotype with their edit lower bounds (section 3).
// right_anchored: anchor = alignment end r (the reverse read); else anchor = start s.
struct Anchor {
    std::size_t at = 0;
    std::size_t d = 0;
};
void qualifying_anchors(const HaplotypeIndex& index, const std::string& read, std::size_t band, std::size_t delta,
                        bool right_anchored, const Variant& v, std::vector<Anchor>& out,
                        std::vector<std::uint32_t>& scratch, std::vector<std::uint64_t>& votes) {
    out.clear();
    const std::size_t L = read.size(), n = index.h->size(), p = band + 1;
    if (n < L) return;
    const std::size_t K = v.no_shift ? 0 : std::min(delta, band);
    votes.clear();
    for (std::size_t j = 0; j < p; ++j) {
        const std::size_t o = j * L / p, len = (j + 1) * L / p - o;
        if (len == 0) continue;
        scratch.clear();
        index.occurrences(read, o, len, scratch);
        if (v.first_hit_only && scratch.size() > 1) scratch.resize(1);
        for (const std::uint32_t x : scratch) {
            for (std::size_t kk = 0; kk <= 2 * K; ++kk) {
                const long k = static_cast<long>(kk) - static_cast<long>(K);
                // left-anchored: x = s + o + k ; right-anchored: x = r - (L - o) - k
                const long a = right_anchored ? static_cast<long>(x) + static_cast<long>(L - o) + k
                                              : static_cast<long>(x) - static_cast<long>(o) - k;
                if (a < 0) continue;
                const std::size_t anchor = static_cast<std::size_t>(a);
                if (right_anchored ? (anchor < L || anchor > n) : (anchor + L > n)) continue;
                const std::uint64_t absk = static_cast<std::uint64_t>(k < 0 ? -k : k);
                votes.push_back((static_cast<std::uint64_t>(anchor) << 24) | (static_cast<std::uint64_t>(j) << 8) | absk);
            }
        }
    }
    std::sort(votes.begin(), votes.end());
    std::size_t i = 0;
    std::vector<std::size_t> min_k(p);
    while (i < votes.size()) {
        const std::size_t anchor = static_cast<std::size_t>(votes[i] >> 24);
        std::fill(min_k.begin(), min_k.end(), std::numeric_limits<std::size_t>::max());
        for (; i < votes.size() && (votes[i] >> 24) == anchor; ++i) {
            const std::size_t j = static_cast<std::size_t>((votes[i] >> 8) & 0xffff), k = static_cast<std::size_t>(votes[i] & 0xff);
            min_k[j] = std::min(min_k[j], k);
        }
        std::size_t d = std::numeric_limits<std::size_t>::max();
        for (std::size_t t = 0; t <= K; ++t) {
            std::size_t c = 0;
            for (std::size_t j = 0; j < p; ++j) c += min_k[j] <= t ? 1 : 0;
            d = std::min(d, std::max(p - c, t));
        }
        if (d <= band) out.push_back({anchor, d});
    }
}

// Exact emissions on the early incumbent's haplotypes (the lazy variant): keyed by (read form, window hash), every hit
// verified byte for byte against the stored window, so a reuse is exact (never a hash-equality assumption).
std::uint64_t unit_key(std::uint64_t rid, std::uint64_t window_hash) { return rid * 0x9e3779b97f4a7c15ULL ^ window_hash; }
struct ExactStore {
    std::unordered_map<std::uint64_t, std::vector<std::pair<std::string, double>>> map;
    std::size_t entries = 0, bytes = 0;
    void put(std::uint64_t key, const char* w, std::size_t n, double value) {
        auto& bucket = map[key];
        for (const auto& [s, v] : bucket) if (s.size() == n && std::memcmp(s.data(), w, n) == 0) return;
        bucket.emplace_back(std::string(w, n), value);
        ++entries;
        bytes += n + sizeof(double) + 32;
    }
    const double* find(std::uint64_t key, const char* w, std::size_t n) const {
        const auto it = map.find(key);
        if (it == map.end()) return nullptr;
        for (const auto& [s, v] : it->second) if (s.size() == n && std::memcmp(s.data(), w, n) == 0) return &v;
        return nullptr;
    }
};

struct CellStats {
    std::uint64_t forward_anchors = 0, reverse_anchors = 0, reused = 0, bounded = 0, pair_terms = 0;
};

// Section 1 relaxation with section 2/3 (and 4) emission ceilings: U_fh for one fragment on one haplotype, in log,
// computed once per anchor for both variants: .first pure (section 3 only), .second lazy (exact where the window was
// evaluated on the incumbent). `exact_out`: when non-null, the exact emission at every qualifying anchor is computed
// and stored (the incumbent's haplotypes only). `units`: when non-null, every (read form, window) unit is recorded.
std::pair<double, double> cell_bound(const std::string& h, const HaplotypeIndex& index, const ReadPair& f, std::size_t fid,
                                     const LibraryModel& lib, const panvar::DiscreteInsertDistribution& insert,
                                     std::size_t band, ShapeTables& tables, const Variant& v, ExactStore* exact_out,
                                     const ExactStore* exact_in, CellStats& stats, std::unordered_set<std::uint64_t>* units) {
    static thread_local std::vector<Anchor> fa, ra;
    static thread_local std::vector<std::uint32_t> scratch;
    static thread_local std::vector<std::uint64_t> votes;
    static thread_local std::vector<std::pair<std::size_t, double>> fv, rv, fvl, rvl;
    const std::size_t delta = lib.max_gap_bp;
    const std::string rc1 = revcomp(f.first.bases), rc2 = revcomp(f.second.bases);
    double pure = 0.0, lazy = 0.0;
    for (int arm = 0; arm < 2; ++arm) {
        if (v.drop_second_arm && arm == 1) continue;
        const std::string& F = arm == 0 ? f.first.bases : f.second.bases;
        const std::string& R = arm == 0 ? rc2 : rc1;
        qualifying_anchors(index, F, band, delta, false, v, fa, scratch, votes);
        if (fa.empty()) continue;
        qualifying_anchors(index, R, band, delta, true, v, ra, scratch, votes);
        if (ra.empty()) continue;
        stats.forward_anchors += fa.size();
        stats.reverse_anchors += ra.size();
        // (pure value, lazy value) of one anchor
        const auto value = [&](const std::string& read, std::size_t at, std::size_t d, bool right, int side) {
            const std::size_t available = right ? at : h.size() - at;
            const std::size_t a = std::min(available, read.size() + delta);
            const double bound = tables.get(read.size(), available)[d];
            if (units == nullptr && exact_in == nullptr && exact_out == nullptr) return std::make_pair(bound, bound);
            const std::uint64_t rid = static_cast<std::uint64_t>(fid) * 4 + static_cast<std::uint64_t>(arm) * 2 + static_cast<std::uint64_t>(side);
            const char* w = h.data() + (right ? at - a : at);
            const std::uint64_t key = unit_key(rid, std::hash<std::string_view>()(std::string_view(w, a)));
            if (units != nullptr) units->insert(key);
            if (exact_out != nullptr) {
                const auto fit = right ? panvar::gapped_emission_ending_at(h, read, at, lib, band)
                                       : panvar::gapped_emission(h, read, at, lib, band);
                const double e = fit ? fit->log_emission : kNegInf;
                exact_out->put(key, w, a, e);
                ++stats.reused;
                return std::make_pair(bound, e);
            }
            if (exact_in != nullptr) {
                if (const double* e = exact_in->find(key, w, a)) { ++stats.reused; return std::make_pair(bound, *e); }
            }
            ++stats.bounded;
            return std::make_pair(bound, bound);
        };
        fv.clear(); rv.clear(); fvl.clear(); rvl.clear();
        for (const Anchor& x : fa) {
            const auto [b, l] = value(F, x.at, x.d, false, 0);
            if (b != kNegInf) fv.push_back({x.at, std::exp(b)});
            if (l != kNegInf) fvl.push_back({x.at, std::exp(l)});
        }
        for (const Anchor& x : ra) {
            const auto [b, l] = value(R, x.at, x.d, true, 1);
            if (b != kNegInf) rv.push_back({x.at, std::exp(b)});
            if (l != kNegInf) rvl.push_back({x.at, std::exp(l)});
        }
        // sum over s, r with r - s inside the library support (section 1 relaxation)
        const auto pair_sum = [&](const std::vector<std::pair<std::size_t, double>>& fw,
                                  const std::vector<std::pair<std::size_t, double>>& rw) {
            double total = 0.0;
            for (const auto& [s, fe] : fw) {
                const std::size_t lo = s + insert.support.min_bp, hi = s + insert.support.max_bp;
                auto it = std::lower_bound(rw.begin(), rw.end(), std::make_pair(lo, -1.0));
                double inner = 0.0;
                for (; it != rw.end() && it->first <= hi; ++it) {
                    if (it->first <= s) continue;
                    inner += std::exp(insert.log_probability(it->first - s)) * it->second;
                    ++stats.pair_terms;
                }
                total += 0.5 * fe * inner;
            }
            return total;
        };
        pure += pair_sum(fv, rv);
        lazy += pair_sum(fvl, rvl);
    }
    return {pure > 0.0 ? std::log(pure) : kNegInf, lazy > 0.0 ? std::log(lazy) : kNegInf};
}

// exposure_for_length, exactly as the scorer computes it (map order)
double exposure_of(std::size_t length, std::size_t fragments, const std::map<std::size_t, std::size_t>& by_floor,
                   const panvar::DiscreteInsertDistribution& insert) {
    if (fragments == 0) return panvar::haplotype_exposure(length, insert, 0);
    double total = 0.0;
    for (const auto& item : by_floor) {
        total += static_cast<double>(item.second) / static_cast<double>(fragments) *
                 panvar::haplotype_exposure(length, insert, item.first);
    }
    return total;
}

// UB of one candidate (section 0) and its absolute term sum A (section 5)
std::pair<double, double> candidate_bound(const std::vector<double>& ua, const std::vector<double>& ub,
                                          const std::vector<double>& background, double exposure_sum, double prior,
                                          const LibraryModel& lib) {
    double sum = 0.0, absolute = 0.0;
    for (std::size_t f = 0; f < background.size(); ++f) {
        const double term = panvar::log_fragment_contribution(lae(ua[f], ub[f]), background[f], lib);
        sum += term;
        absolute += std::abs(term);
    }
    const double exposure_term = -panvar::expected_signal_fragments(exposure_sum, lib);
    return {exposure_term + sum + prior, absolute + std::abs(exposure_term) + std::abs(prior)};
}
double allowance(std::size_t fragments, double absolute) {
    return 64.0 * static_cast<double>(fragments + 8) * std::ldexp(1.0, -53) * absolute;
}

std::vector<panvar::SequencingRead> load_fastq(const std::string& file) {
    std::vector<panvar::SequencingRead> out;
    std::ifstream in(file);
    if (!in) throw std::runtime_error("cannot read " + file);
    for (std::string head, bases, plus, quals; std::getline(in, head) && std::getline(in, bases) && std::getline(in, plus) &&
                                               std::getline(in, quals);) {
        panvar::SequencingRead r;
        const std::size_t sp = head.find_first_of(" \t");
        r.name = head.substr(1, sp == std::string::npos ? std::string::npos : sp - 1);
        for (char& c : bases) c = static_cast<char>(std::toupper(c));
        r.bases = bases;
        r.qualities = quals;
        out.push_back(std::move(r));
    }
    return out;
}

// ============================================================ fixtures ================================================
std::mt19937_64 g_rng(20261018);
std::string random_bases(std::size_t n) {
    static const char* b = "ACGT";
    std::string s(n, 'A');
    for (char& c : s) c = b[g_rng() & 3];
    return s;
}
std::string mutate(const std::string& s, std::size_t subs, std::size_t ins, std::size_t dels) {
    std::string out = s;
    for (std::size_t i = 0; i < subs && !out.empty(); ++i) {
        const std::size_t at = g_rng() % out.size();
        out[at] = "ACGT"[(code_of(out[at]) + 1 + g_rng() % 3) % 4];
    }
    for (std::size_t i = 0; i < ins; ++i) out.insert(g_rng() % (out.size() + 1), random_bases(1 + g_rng() % 2));
    for (std::size_t i = 0; i < dels && out.size() > 3; ++i) out.erase(g_rng() % (out.size() - 2), 1 + g_rng() % 2);
    return out;
}

struct FixtureResult {
    std::size_t cells = 0, violations = 0, production_over_oracle = 0, candidates = 0, candidate_violations = 0;
    std::size_t lazy_violations = 0, tie_exclusions = 0;
    double worst_margin = std::numeric_limits<double>::infinity();   // min over cells of log U - log M (finite M)
};

FixtureResult run_fixture(const std::string& name, const std::vector<std::string>& haps, const std::vector<ReadPair>& frags,
                          const LibraryModel& lib, std::size_t band, const Variant& v, bool report) {
    FixtureResult out;
    const panvar::DiscreteInsertDistribution insert = panvar::build_library_insert_distribution(lib);
    ShapeTables tables;
    tables.lib = &lib;
    tables.band = band;
    std::size_t q = 64;
    for (const ReadPair& f : frags) {
        for (const auto* r : {&f.first.bases, &f.second.bases}) q = std::min(q, r->size() / (band + 1));
    }
    q = std::max<std::size_t>(1, std::min<std::size_t>(q, 10));
    std::vector<HaplotypeIndex> index(haps.size());
    for (std::size_t h = 0; h < haps.size(); ++h) index[h].build(haps[h], q);
    // oracle masses (and the production scorer's, which must not exceed them)
    std::vector<std::vector<double>> oracle(haps.size(), std::vector<double>(frags.size())), bound = oracle, lazy = oracle;
    panvar::InternedHaplotypes ih;
    std::vector<panvar::MosaicPath> paths;
    for (std::size_t h = 0; h < haps.size(); ++h) paths.push_back(panvar::MosaicPath{{h}});
    ih = panvar::intern_haplotypes(paths, haps);
    std::vector<panvar::HaplotypeCandidate> pairs;
    for (std::size_t a = 0; a < ih.sequences.size(); ++a) {
        for (std::size_t b = a; b < ih.sequences.size(); ++b) pairs.push_back({a, b});
    }
    panvar::FragmentScoringOptions so;
    so.threads = 1;
    so.mismatch_band = band;
    so.max_physical_origins = std::numeric_limits<std::uint64_t>::max() / 4;
    so.minimum_piece_bp = 4;   // the fixtures use 60 bp reads; the oracle has no piece floor
    const panvar::FragmentScoringResult prod = panvar::score_diplotypes(ih, pairs, frags, lib, so);
    for (std::size_t h = 0; h < haps.size(); ++h) {
        for (std::size_t f = 0; f < frags.size(); ++f) {
            const auto placements = panvar::enumerate_placements_exhaustive(haps[h], frags[f], panvar::fragment_insert_support(frags[f], lib), lib, band);
            oracle[h][f] = panvar::log_placement_mass(placements, frags[f], insert, lib);
            const double pm = prod.log_mass[f * prod.unique_haplotypes + ih.index_of[h]];
            if (pm > oracle[h][f] + 1e-9 * std::max(1.0, std::abs(oracle[h][f]))) ++out.production_over_oracle;
        }
    }
    // the bound, pure and lazy (incumbent: haplotypes 0 and 1)
    ExactStore exact;
    for (std::size_t f = 0; f < frags.size(); ++f) {
        for (std::size_t h = 0; h < 2 && h < haps.size(); ++h) {
            CellStats st;
            (void)cell_bound(haps[h], index[h], frags[f], f, lib, insert, band, tables, v, &exact, nullptr, st, nullptr);
        }
    }
    for (std::size_t h = 0; h < haps.size(); ++h) {
        for (std::size_t f = 0; f < frags.size(); ++f) {
            CellStats st;
            const auto both = cell_bound(haps[h], index[h], frags[f], f, lib, insert, band, tables, v, nullptr, &exact, st, nullptr);
            bound[h][f] = both.first;
            lazy[h][f] = both.second;
            ++out.cells;
            const double m = oracle[h][f];
            const double tol = 1e-9 * std::max(1.0, std::abs(m));
            if (m != kNegInf && !(bound[h][f] >= m - tol)) ++out.violations;
            if (m != kNegInf && !(lazy[h][f] >= m - tol)) ++out.lazy_violations;
            if (m != kNegInf) out.worst_margin = std::min(out.worst_margin, bound[h][f] - m);
        }
    }
    // candidates: UB >= S for every pair, S from the oracle masses by the scorer's own combination
    std::map<std::size_t, std::size_t> by_floor;
    std::vector<double> background;
    for (const ReadPair& f : frags) {
        ++by_floor[std::max(f.first.bases.size(), f.second.bases.size())];
        background.push_back(panvar::log_background_probability(f));
    }
    std::vector<double> exposure(haps.size());
    for (std::size_t h = 0; h < haps.size(); ++h) exposure[h] = exposure_of(haps[h].size(), frags.size(), by_floor, insert);
    std::vector<std::pair<double, double>> exact_scores;   // (S, A) per pair (a <= b)
    std::vector<std::pair<std::size_t, std::size_t>> cand;
    for (std::size_t a = 0; a < haps.size(); ++a) {
        for (std::size_t b = a; b < haps.size(); ++b) {
            cand.push_back({a, b});
            const auto s = candidate_bound(oracle[a], oracle[b], background, exposure[a] + exposure[b], 0.0, lib);
            const auto u = candidate_bound(bound[a], bound[b], background, exposure[a] + exposure[b], 0.0, lib);
            const auto ul = candidate_bound(lazy[a], lazy[b], background, exposure[a] + exposure[b], 0.0, lib);
            exact_scores.push_back(s);
            ++out.candidates;
            if (!(u.first + allowance(frags.size(), u.second) >= s.first - allowance(frags.size(), s.second))) ++out.candidate_violations;
            if (!(ul.first + allowance(frags.size(), ul.second) >= s.first - allowance(frags.size(), s.second))) ++out.candidate_violations;
        }
    }
    // ties: no candidate whose exact score equals the incumbent's (pair 0: (0,1) when present) is ever excluded
    std::size_t inc = 0;
    for (std::size_t c = 0; c < cand.size(); ++c) if (cand[c] == std::make_pair<std::size_t, std::size_t>(0, 1)) inc = c;
    for (std::size_t c = 0; c < cand.size(); ++c) {
        if (exact_scores[c].first != exact_scores[inc].first) continue;
        const auto u = candidate_bound(lazy[cand[c].first], lazy[cand[c].second], background,
                                       exposure[cand[c].first] + exposure[cand[c].second], 0.0, lib);
        if (u.first + allowance(frags.size(), u.second) < exact_scores[inc].first - allowance(frags.size(), exact_scores[inc].second)) {
            ++out.tie_exclusions;
        }
    }
    if (report) {
        std::cout << name << "\tcells=" << out.cells << "\tviolations=" << out.violations << "\tlazy_violations="
                  << out.lazy_violations << "\tproduction_over_oracle=" << out.production_over_oracle << "\tcandidates="
                  << out.candidates << "\tcandidate_violations=" << out.candidate_violations << "\ttie_exclusions="
                  << out.tie_exclusions << "\tworst_cell_margin_nats=" << std::setprecision(4) << out.worst_margin
                  << "\tshape_Z_rel_err=" << tables.max_z_mismatch << '\n';
    }
    return out;
}

ReadPair fragment_from(const std::string& h, std::size_t start, std::size_t insert, std::size_t l1, std::size_t l2,
                       bool second_forward, std::size_t subs, std::size_t indels) {
    insert = std::min(insert, h.size() - start);
    const std::string templ = h.substr(start, insert);
    std::string m1 = templ.substr(0, std::min(l1, templ.size()));
    std::string m2 = revcomp(templ.substr(templ.size() - std::min(l2, templ.size())));
    m1 = mutate(m1, subs, indels / 2, indels - indels / 2);
    m2 = mutate(m2, subs / 2, 0, 0);
    ReadPair f;
    f.first.bases = second_forward ? m2 : m1;
    f.second.bases = second_forward ? m1 : m2;
    f.first.name = f.second.name = "f";
    f.first.qualities = std::string(f.first.bases.size(), 'I');
    f.second.qualities = std::string(f.second.bases.size(), 'I');
    return f;
}

int fixtures() {
    // Libraries: the production defaults, and a stressed gap geometry (frequent gaps, a narrow band)
    LibraryModel base;
    base.insert_mean_bp = 260;
    base.insert_sd_bp = 25;
    base.support_sigmas = 4;
    base.substitution_error = 0.01;
    base.outlier_weight = 0.002;
    base.fragment_start_intensity = 0.05;
    base.provenance = "s2b fixture library";
    LibraryModel gappy = base;
    gappy.gap_open = 0.02;
    gappy.gap_extend = 0.6;
    gappy.max_gap_bp = 2;
    struct Setup { std::string name; LibraryModel lib; std::size_t band; };
    const std::vector<Setup> setups = {{"defaults_band6_gap3", base, 6}, {"gappy_band3_gap2", gappy, 3}};
    // Haplotypes: H0 random; H1 SNPs; H2 a tandem duplication (repeated copies); H3 a distant copy (an alternative
    // origin far away); H4 a mosaic junction of H0 and H1; H5 short indels; H6 an empty block (a deleted segment
    // joined across); H7 sequence-identical to H1 (an alternative labelling); H8 short (edge effects).
    const std::string h0 = random_bases(1400);
    const std::string h1 = mutate(h0, 6, 0, 0);
    const std::string h2 = h0.substr(0, 700) + h0.substr(400, 300) + h0.substr(700);
    const std::string h3 = h0 + random_bases(300) + h0.substr(100, 250);
    const std::string h4 = h0.substr(0, 650) + h1.substr(650);
    const std::string h5 = mutate(h0, 0, 3, 3);
    const std::string h6 = h0.substr(0, 500) + h0.substr(560);
    const std::string h7 = h1;
    const std::string h8 = h0.substr(300, 330);
    const std::vector<std::string> haps = {h0, h1, h2, h3, h4, h5, h6, h7, h8};
    std::vector<ReadPair> frags;
    const std::size_t lengths[] = {60, 100, 150};
    for (std::size_t i = 0; i < 40; ++i) {
        const std::string& src = haps[i % haps.size()];
        if (src.size() < 200) continue;
        const std::size_t l1 = lengths[i % 3], l2 = lengths[(i + 1) % 3];
        const std::size_t insert = 200 + (g_rng() % 120);
        const std::size_t start = g_rng() % (src.size() - std::min(src.size() - 1, insert));
        frags.push_back(fragment_from(src, start, insert, l1, l2, (i % 2) == 1, i % 4, (i % 5 == 0) ? 2 : 0));
    }
    frags.push_back(fragment_from(random_bases(400), 0, 260, 100, 100, false, 0, 0));   // no placement anywhere
    frags.push_back(fragment_from(h0, 100, 250, 100, 100, false, 12, 0));               // beyond the band
    frags.push_back(fragment_from(h3, 1450, 240, 60, 60, true, 0, 0));                  // inside a distant copy
    // near ties: a fragment set symmetric between H1 and H7 is identical by construction
    int failures = 0;
    for (const Setup& s : setups) {
        const FixtureResult r = run_fixture(s.name, haps, frags, s.lib, s.band, Variant{}, true);
        if (r.violations || r.lazy_violations || r.candidate_violations || r.tie_exclusions || r.production_over_oracle) ++failures;
        // the validator must reject each deliberate defect on at least one cell
        const std::pair<std::string, Variant> defects[] = {
            {"defect_no_shift", Variant{true, false, false}},
            {"defect_drop_second_arm", Variant{false, true, false}},
            {"defect_first_hit_only", Variant{false, false, true}}};
        for (const auto& [dname, dv] : defects) {
            const FixtureResult d = run_fixture(s.name + "/" + dname, haps, frags, s.lib, s.band, dv, true);
            // pure, lazy (exact reuse) or candidate level: the pure bound's own slack can mask a defect that the exact
            // reuse exposes (observed for the no-shift defect), so every level counts
            if (d.violations + d.lazy_violations + d.candidate_violations == 0) {
                std::cout << "VALIDATOR FAILED TO REJECT " << s.name << "/" << dname << '\n';
                ++failures;
            }
        }
    }
    std::cout << (failures == 0 ? "s2b fixtures: PASS" : "s2b fixtures: FAIL") << '\n';
    return failures == 0 ? 0 : 1;
}

// ============================================================ measure ================================================
struct Dump {
    std::uint64_t H = 0, F = 0, C = 0, band = 0;
    LibraryModel lib;
    std::vector<std::string> haps;
    std::vector<double> log_mass, exposure, background, hap_prior;
    struct Candidate { std::uint64_t a, b; double prior, ll, ft, et, abs; std::uint64_t fwm; };
    std::vector<Candidate> cand;
    std::vector<std::uint64_t> winner_set;
    std::uint64_t top = 0, early_a = 0, early_b = 0;
};
Dump read_dump(const std::string& file) {
    std::ifstream in(file, std::ios::binary);
    if (!in) throw std::runtime_error("cannot read " + file);
    std::string all((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    std::size_t at = 0;
    const auto get = [&](auto& v) {
        std::memcpy(&v, all.data() + at, sizeof v);
        at += sizeof v;
    };
    if (all.compare(0, 5, "PVEV1") != 0) throw std::runtime_error("not a PVEV1 dump");
    at = 8;
    Dump d;
    get(d.H); get(d.F); get(d.C);
    get(d.lib.insert_mean_bp); get(d.lib.insert_sd_bp); get(d.lib.support_sigmas); get(d.lib.substitution_error);
    get(d.lib.outlier_weight); get(d.lib.fragment_start_intensity); get(d.lib.gap_open); get(d.lib.gap_extend);
    std::uint64_t gap = 0;
    get(gap);
    d.lib.max_gap_bp = gap;
    get(d.band);
    d.lib.provenance = "union_evidence dump";
    for (std::uint64_t h = 0; h < d.H; ++h) {
        std::uint64_t n = 0;
        get(n);
        d.haps.push_back(all.substr(at, n));
        at += n;
    }
    d.log_mass.resize(d.F * d.H);
    for (double& x : d.log_mass) get(x);
    d.exposure.resize(d.H);
    for (double& x : d.exposure) get(x);
    d.background.resize(d.F);
    for (double& x : d.background) get(x);
    d.hap_prior.resize(d.H);
    for (double& x : d.hap_prior) get(x);
    d.cand.resize(d.C);
    for (auto& c : d.cand) { get(c.a); get(c.b); get(c.prior); get(c.ll); get(c.ft); get(c.et); get(c.abs); get(c.fwm); }
    std::uint64_t w = 0;
    get(w);
    d.winner_set.resize(w);
    for (auto& x : d.winner_set) get(x);
    get(d.top); get(d.early_a); get(d.early_b);
    return d;
}

template <class Fn>
void parallel_for(std::size_t n, std::size_t threads, Fn fn) {
    std::atomic<std::size_t> next{0};
    std::vector<std::thread> pool;
    for (std::size_t t = 0; t < threads; ++t) {
        pool.emplace_back([&]() {
            for (std::size_t i; (i = next.fetch_add(1)) < n;) fn(i);
        });
    }
    for (auto& t : pool) t.join();
}

double quantile(std::vector<double> v, double q) {
    if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(v.begin(), v.end());
    return v[static_cast<std::size_t>(q * static_cast<double>(v.size() - 1))];
}

int measure(const std::string& dump_file, const std::string& r1, const std::string& r2, std::size_t threads,
            const std::string& out_file, bool construction_only, std::size_t max_haplotypes) {
    std::map<std::string, std::string> row;
    const auto put = [&](const std::string& k, auto v) {
        std::ostringstream o;
        o << std::setprecision(10) << v;
        row[k] = o.str();
    };
    const Clock::time_point t_load = Clock::now();
    const Dump d = read_dump(dump_file);
    const auto a1 = load_fastq(r1), a2 = load_fastq(r2);
    std::vector<ReadPair> frags;
    for (std::size_t i = 0; i < a1.size() && i < a2.size(); ++i) frags.push_back({a1[i], a2[i]});
    put("seconds_load", since(t_load));
    put("haplotypes_dump", d.H);
    put("fragments_dump", d.F);
    put("fragments_reads", frags.size());
    put("candidates", d.C);
    if (!construction_only && frags.size() != d.F) throw std::runtime_error("the reads do not match the dump's fragments");
    const std::size_t F = frags.size(), band = d.band;
    const bool have_early = d.early_a < d.H && d.early_b < d.H;
    // the haplotypes bounded: all of them, or (construction only) the incumbent's plus the first max_haplotypes others
    std::vector<std::size_t> hap_ids;
    if (have_early) { hap_ids.push_back(d.early_a); if (d.early_b != d.early_a) hap_ids.push_back(d.early_b); }
    for (std::size_t h = 0; h < d.H && (max_haplotypes == 0 || hap_ids.size() < max_haplotypes); ++h) {
        if (std::find(hap_ids.begin(), hap_ids.end(), h) == hap_ids.end()) hap_ids.push_back(h);
    }
    const std::size_t H = d.H;
    put("haplotypes_bounded", hap_ids.size());
    const panvar::DiscreteInsertDistribution insert = panvar::build_library_insert_distribution(d.lib);
    ShapeTables tables;
    tables.lib = &d.lib;
    tables.band = band;
    const Clock::time_point t_index = Clock::now();
    std::size_t q = 64;
    for (const ReadPair& f : frags) q = std::min({q, f.first.bases.size() / (band + 1), f.second.bases.size() / (band + 1)});
    q = std::max<std::size_t>(1, std::min<std::size_t>(q, 12));
    std::vector<HaplotypeIndex> index(H);
    parallel_for(hap_ids.size(), threads, [&](std::size_t i) { index[hap_ids[i]].build(d.haps[hap_ids[i]], q); });
    std::uint64_t index_bytes = 0;
    for (const auto& x : index) index_bytes += x.entries.capacity() * sizeof(std::uint64_t);
    put("seconds_index", since(t_index));
    put("index_bytes", index_bytes);
    put("q", q);
    // the lazy variant's exact store: exact emissions at EVERY qualifying anchor of the incumbent's haplotypes (its
    // own cost: production's early evaluation verifies only seeded starts, not the shifted neighbours)
    ExactStore exact;
    const Clock::time_point t_exact = Clock::now();
    CellStats exact_stats;
    if (have_early) {
        for (std::size_t f = 0; f < F; ++f) {
            for (const std::uint64_t h : {d.early_a, d.early_b}) {
                if (h == d.early_b && d.early_b == d.early_a) continue;
                (void)cell_bound(d.haps[h], index[h], frags[f], f, d.lib, insert, band, tables, Variant{}, &exact, nullptr, exact_stats, nullptr);
            }
        }
    }
    put("seconds_lazy_exact_store", since(t_exact));
    put("lazy_exact_store_entries", exact.entries);
    put("lazy_exact_store_bytes", exact.bytes);
    put("lazy_exact_store_anchor_evaluations", exact_stats.reused);
    put("peak_rss_gib_after_exact_store", peak_rss_gib());
    // every cell, both variants, one anchor pass; distinct (read form, window) units per worker
    const Clock::time_point t_cells = Clock::now();
    std::vector<std::vector<double>> pure(H), lazy(H);
    std::vector<CellStats> stats(hap_ids.size());
    std::mutex units_guard;
    std::unordered_set<std::uint64_t> units_all;
    {
        std::atomic<std::size_t> next{0};
        std::vector<std::thread> pool;
        for (std::size_t t = 0; t < threads; ++t) {
            pool.emplace_back([&]() {
                std::unordered_set<std::uint64_t> local;
                for (std::size_t i; (i = next.fetch_add(1)) < hap_ids.size();) {
                    const std::size_t h = hap_ids[i];
                    pure[h].assign(F, kNegInf);
                    lazy[h].assign(F, kNegInf);
                    for (std::size_t f = 0; f < F; ++f) {
                        const auto both = cell_bound(d.haps[h], index[h], frags[f], f, d.lib, insert, band, tables, Variant{},
                                                     nullptr, have_early ? &exact : nullptr, stats[i], &local);
                        pure[h][f] = both.first;
                        lazy[h][f] = both.second;
                    }
                    if (local.size() > (1u << 22)) {
                        std::lock_guard<std::mutex> lock(units_guard);
                        units_all.insert(local.begin(), local.end());
                        local.clear();
                    }
                }
                std::lock_guard<std::mutex> lock(units_guard);
                units_all.insert(local.begin(), local.end());
            });
        }
        for (auto& t : pool) t.join();
    }
    put("seconds_cells_pure_and_lazy", since(t_cells));
    put("shape_Z_rel_err", tables.max_z_mismatch);
    put("shape_tables", tables.tail.size());
    std::uint64_t anchors = 0, reused = 0, bounded = 0, pair_terms = 0, cells_nonzero = 0;
    for (const auto& st : stats) { anchors += st.forward_anchors + st.reverse_anchors; reused += st.reused; bounded += st.bounded; pair_terms += st.pair_terms; }
    for (const std::size_t h : hap_ids) for (std::size_t f = 0; f < F; ++f) cells_nonzero += pure[h][f] != kNegInf;
    put("qualifying_anchors", anchors);
    put("anchors_per_fragment_haplotype", static_cast<double>(anchors) / static_cast<double>(std::max<std::size_t>(1, F * hap_ids.size())));
    put("pair_terms", pair_terms);
    put("lazy_anchor_values_reused_exact", reused);
    put("lazy_anchor_values_bounded", bounded);
    put("cells", static_cast<std::uint64_t>(hap_ids.size()) * F);
    put("cells_nonzero_bound", cells_nonzero);
    put("distinct_window_units_all", units_all.size());
    put("peak_rss_gib_after_cells", peak_rss_gib());
    if (construction_only) {
        put("mode", "construction_only: no reference, no validation, no exclusion or completion claim");
    } else {
        std::uint64_t violations = 0, lazy_violations = 0, zero_bound_positive_mass = 0;
        std::vector<double> cell_slack, cell_slack_lazy;
        for (std::size_t h = 0; h < H; ++h) {
            for (std::size_t f = 0; f < F; ++f) {
                const double m = d.log_mass[f * H + h];
                if (m == kNegInf) continue;
                const double tol = 1e-9 * std::max(1.0, std::abs(m));
                if (!(pure[h][f] >= m - tol)) ++violations;
                if (!(lazy[h][f] >= m - tol)) ++lazy_violations;
                if (pure[h][f] == kNegInf) ++zero_bound_positive_mass;
                if (f % 7 == 0) { cell_slack.push_back(pure[h][f] - m); cell_slack_lazy.push_back(lazy[h][f] - m); }
            }
        }
        put("cell_violations_pure", violations);
        put("cell_violations_lazy", lazy_violations);
        put("cells_zero_bound_positive_mass", zero_bound_positive_mass);
        put("cell_slack_pure_median_sampled", quantile(cell_slack, 0.5));
        put("cell_slack_pure_p90_sampled", quantile(cell_slack, 0.9));
        put("cell_slack_lazy_median_sampled", quantile(cell_slack_lazy, 0.5));
        put("cell_slack_lazy_p90_sampled", quantile(cell_slack_lazy, 0.9));
        std::map<std::size_t, std::size_t> by_floor;
        std::vector<double> background(F);
        for (std::size_t f = 0; f < F; ++f) {
            ++by_floor[std::max(frags[f].first.bases.size(), frags[f].second.bases.size())];
            background[f] = panvar::log_background_probability(frags[f]);
        }
        std::vector<double> exposure(H);
        double exposure_diff = 0.0, background_diff = 0.0;
        for (std::size_t h = 0; h < H; ++h) {
            exposure[h] = exposure_of(d.haps[h].size(), F, by_floor, insert);
            exposure_diff = std::max(exposure_diff, std::abs(exposure[h] - d.exposure[h]));
        }
        for (std::size_t f = 0; f < F; ++f) background_diff = std::max(background_diff, std::abs(background[f] - d.background[f]));
        put("exposure_max_abs_diff_vs_scorer", exposure_diff);
        put("background_max_abs_diff_vs_scorer", background_diff);
        const Clock::time_point t_cand = Clock::now();
        std::vector<double> ub_pure(d.C), ub_lazy(d.C), a_pure(d.C), a_lazy(d.C);
        parallel_for(d.C, threads, [&](std::size_t c) {
            const auto& x = d.cand[c];
            const double prior = d.hap_prior[x.a] + d.hap_prior[x.b] + (x.a != x.b ? std::log(2.0) : 0.0);
            const auto p = candidate_bound(pure[x.a], pure[x.b], background, exposure[x.a] + exposure[x.b], prior, d.lib);
            const auto l = candidate_bound(lazy[x.a], lazy[x.b], background, exposure[x.a] + exposure[x.b], prior, d.lib);
            ub_pure[c] = p.first; a_pure[c] = p.second;
            ub_lazy[c] = l.first; a_lazy[c] = l.second;
        });
        put("seconds_candidate_bounds", since(t_cand));
        std::uint64_t cand_violations = 0, prior_mismatch = 0;
        std::vector<double> slack_pure, slack_lazy;
        const auto total = [&](std::size_t c) { return d.cand[c].ll + d.cand[c].prior; };
        for (std::size_t c = 0; c < d.C; ++c) {
            const auto& x = d.cand[c];
            const double prior = d.hap_prior[x.a] + d.hap_prior[x.b] + (x.a != x.b ? std::log(2.0) : 0.0);
            if (prior != x.prior) ++prior_mismatch;
            if (ub_pure[c] + allowance(F, a_pure[c]) + allowance(F, x.abs) < total(c) ||
                ub_lazy[c] + allowance(F, a_lazy[c]) + allowance(F, x.abs) < total(c)) ++cand_violations;
            slack_pure.push_back(ub_pure[c] - total(c));
            slack_lazy.push_back(ub_lazy[c] - total(c));
        }
        put("candidate_violations", cand_violations);
        put("candidate_prior_mismatch", prior_mismatch);
        put("candidate_slack_pure_median", quantile(slack_pure, 0.5));
        put("candidate_slack_pure_p10", quantile(slack_pure, 0.1));
        put("candidate_slack_lazy_median", quantile(slack_lazy, 0.5));
        put("candidate_slack_lazy_p10", quantile(slack_lazy, 0.1));
        std::size_t inc = d.C;
        const std::size_t win = d.top;
        for (std::size_t c = 0; c < d.C; ++c) if (d.cand[c].a == d.early_a && d.cand[c].b == d.early_b) inc = c;
        std::vector<double> gap_to_winner;
        for (std::size_t c = 0; c < d.C; ++c) gap_to_winner.push_back(total(win) - total(c));
        put("true_gap_to_winner_median", quantile(gap_to_winner, 0.5));
        put("true_gap_to_winner_p10", quantile(gap_to_winner, 0.1));
        put("winner_slack_pure", slack_pure[win]);
        put("winner_slack_lazy", slack_lazy[win]);
        std::set<std::uint64_t> top_set(d.winner_set.begin(), d.winner_set.end());
        top_set.insert(win);
        put("top_set_size", top_set.size());
        const auto exclude = [&](std::size_t ref, const std::vector<double>& ub, const std::vector<double>& a, const std::string& tag) {
            std::vector<char> out(d.C, 0);
            if (ref >= d.C) { put("excluded_" + tag, "NA"); return out; }
            const double floor_score = total(ref) - allowance(F, d.cand[ref].abs);
            std::uint64_t n = 0, top_lost = 0;
            for (std::size_t c = 0; c < d.C; ++c) {
                if (ub[c] + allowance(F, a[c]) < floor_score) { out[c] = 1; ++n; if (top_set.count(c)) ++top_lost; }
            }
            put("excluded_" + tag, n);
            put("excluded_fraction_" + tag, static_cast<double>(n) / static_cast<double>(d.C));
            put("top_set_lost_" + tag, top_lost);
            return out;
        };
        put("incumbent_found", inc < d.C ? "yes" : "no");
        if (inc < d.C) {
            put("incumbent_gap_to_winner", total(win) - total(inc));
            put("incumbent_rank_by_total", 1 + std::count_if(gap_to_winner.begin(), gap_to_winner.end(),
                                                               [&](double g) { return g < total(win) - total(inc); }));
        }
        (void)exclude(inc, ub_pure, a_pure, "vs_incumbent_pure");
        const auto ex_inc_lazy = exclude(inc, ub_lazy, a_lazy, "vs_incumbent_lazy");
        (void)exclude(win, ub_pure, a_pure, "vs_winner_pure_OPTIMISTIC");
        const auto ex_win_lazy = exclude(win, ub_lazy, a_lazy, "vs_winner_lazy_OPTIMISTIC");
        // evidence still required by the survivors: haplotypes, nonzero cells, distinct window units -- and of those
        // units, the ones the lazy exact store does not already hold (new exact work)
        const auto required = [&](const std::vector<char>& excluded, const std::string& tag) {
            std::vector<char> need(H, 0);
            for (std::size_t c = 0; c < d.C; ++c) if (!excluded[c]) { need[d.cand[c].a] = 1; need[d.cand[c].b] = 1; }
            std::uint64_t haps = 0, cells = 0;
            std::unordered_set<std::uint64_t> seen;
            for (std::size_t h = 0; h < H; ++h) {
                if (!need[h]) continue;
                ++haps;
                for (std::size_t f = 0; f < F; ++f) {
                    cells += pure[h][f] != kNegInf;
                    CellStats st;
                    (void)cell_bound(d.haps[h], index[h], frags[f], f, d.lib, insert, band, tables, Variant{}, nullptr, nullptr, st, &seen);
                }
            }
            put("required_haplotypes_" + tag, haps);
            put("required_cells_nonzero_" + tag, cells);
            put("required_window_units_" + tag, seen.size());
            put("required_window_units_fraction_" + tag, static_cast<double>(seen.size()) / static_cast<double>(std::max<std::size_t>(1, units_all.size())));
        };
        required(ex_inc_lazy, "lazy_vs_incumbent");
        required(ex_win_lazy, "lazy_vs_winner_OPTIMISTIC");
    }
    put("peak_rss_gib", peak_rss_gib());
    std::ofstream out(out_file);
    out << "key\tvalue\n";
    for (const auto& [k, v] : row) out << k << '\t' << v << '\n';
    for (const auto& [k, v] : row) std::cout << k << '\t' << v << '\n';
    return 0;
}
}  // namespace

int main(int argc, char** argv) {
    try {
        const std::string mode = argc > 1 ? argv[1] : "";
        if (mode == "fixtures") return fixtures();
        if (mode == "measure" && argc >= 7) {
            return measure(argv[2], argv[3], argv[4], std::stoul(argv[5]), argv[6], argc > 7 && std::string(argv[7]) == "construction_only",
                           argc > 8 ? std::stoul(argv[8]) : 0);
        }
        std::cerr << "usage: s2b_preplacement_bound fixtures | measure <dump> <r1> <r2> <threads> <out.tsv> [construction_only "
                     "[max_haplotypes]]\n";
        return 2;
    } catch (const std::exception& error) {
        std::cerr << "s2b_preplacement_bound: " << error.what() << '\n';
        return 1;
    }
}
