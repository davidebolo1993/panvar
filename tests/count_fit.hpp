// M6-A: a depth-normalised 31-mer count-fit comparator (DIAGNOSTIC; Ctyper-inspired; not a likelihood, never combined
// with the fragment likelihood or the copying prior). Shared by genotype_probe count-fit and tests/count_fit_gates.cpp.
//
//   loss(candidate) = sum_{k in D} ( observed(k) / depth - predicted(k) )^2
//
// D: the canonical 31-mers of the complete spellings of the active panel (no truth, no candidate additions).
// predicted(k): the summed multiplicity of k on the candidate's two complete haplotypes, junction-spanning 31-mers
// included, every occurrence anywhere in the locus counted. observed(k): the count of k over every 31-mer window of both
// mates of every read pair (each mate independently, both strands via canonical form; a window containing any base other
// than A/C/G/T is skipped). depth: one haploid 31-mer depth per read pool from stable-copy k-mers (see estimate_depth).
//
// Candidates substitute allele pairs at ONE target block inside a fixed two-homologue context. A haplotype is
// left + allele + right; its counts are base(left, right) + middle(allele), where middle counts exactly the 31-mers that
// overlap the allele or span a junction (the window left[-30:] + allele + right[:30]), so no candidate is materialised.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace countfit {

constexpr int K = 31;

inline int code(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1;
    }
}

// Every canonical 31-mer window of s (2-bit packed, 62 bits), in order; windows with a non-ACGT base are skipped.
template <class F>
void for_each_kmer(const std::string& s, F&& f) {
    const std::uint64_t mask = (std::uint64_t(1) << (2 * K)) - 1;
    std::uint64_t fw = 0, rc = 0;
    int valid = 0;
    for (std::size_t i = 0; i < s.size(); ++i) {
        const int c = code(s[i]);
        if (c < 0) { valid = 0; fw = rc = 0; continue; }
        fw = ((fw << 2) | static_cast<std::uint64_t>(c)) & mask;
        rc = (rc >> 2) | (static_cast<std::uint64_t>(3 - c) << (2 * (K - 1)));
        if (++valid >= K) f(std::min(fw, rc), i + 1 - K);
    }
}

inline std::uint64_t mix(std::uint64_t x) {   // splitmix64 finaliser: the per-k-mer weight of the linear profile hash
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

struct Dictionary {
    std::unordered_map<std::uint64_t, std::uint32_t> id;
    std::vector<std::uint64_t> kmer;
    std::uint32_t add(std::uint64_t k) {
        const auto [it, fresh] = id.emplace(k, static_cast<std::uint32_t>(kmer.size()));
        if (fresh) kmer.push_back(k);
        return it->second;
    }
    long find(std::uint64_t k) const {
        const auto it = id.find(k);
        return it == id.end() ? -1 : static_cast<long>(it->second);
    }
};

// A sparse count vector over dictionary ids (sorted), plus the number of 31-mers outside the dictionary.
struct Sparse {
    std::vector<std::pair<std::uint32_t, std::uint32_t>> v;
    std::uint64_t outside = 0;
    std::uint64_t hash = 0;   // sum of count * mix(kmer), mod 2^64: equal profiles <=> equal hash (collisions negligible)
};

inline Sparse count_string(const std::string& s, const Dictionary& d) {
    std::map<std::uint32_t, std::uint32_t> m;
    Sparse out;
    for_each_kmer(s, [&](std::uint64_t k, std::size_t) {
        const long i = d.find(k);
        if (i < 0) ++out.outside;
        else ++m[static_cast<std::uint32_t>(i)];
    });
    out.v.assign(m.begin(), m.end());
    for (const auto& [i, c] : out.v) out.hash += static_cast<std::uint64_t>(c) * mix(d.kmer[i]);
    return out;
}

// The 31-mers of left+allele+right that are NOT wholly inside left or wholly inside right.
inline std::string middle_window(const std::string& left, const std::string& allele, const std::string& right) {
    const std::size_t l = std::min<std::size_t>(left.size(), K - 1), r = std::min<std::size_t>(right.size(), K - 1);
    return left.substr(left.size() - l) + allele + right.substr(0, r);
}

struct Context {
    std::string left[2], right[2];   // per homologue: the spelling before / after the target block
};

struct Candidate {
    std::size_t a = 0, b = 0;   // allele index on homologue 1, on homologue 2
    double loss = 0.0;
    std::uint64_t profile = 0;  // the candidate's whole predicted profile hash
    std::uint64_t outside = 0;  // middle 31-mers outside the dictionary (both homologues)
};

// Exact loss for every ordered pair (a on h1, b on h2) over the dictionary. y = observed / depth (dense over D).
// loss(a,b) = sum_k (r - M1a - M2b)^2 with r = y - B, computed as
//   R0 + sum_{supp M1a} [(r - M1a)^2 - r^2] + sum_{supp M2b} [(r - M1a - M2b)^2 - (r - M1a)^2]
inline std::vector<Candidate> score_all(const std::vector<double>& y, const Sparse& base1, const Sparse& base2,
                                        const std::vector<Sparse>& mid1, const std::vector<Sparse>& mid2) {
    const std::size_t n = y.size();
    std::vector<double> r(y);
    for (const auto& [i, c] : base1.v) r[i] -= c;
    for (const auto& [i, c] : base2.v) r[i] -= c;
    double r0 = 0.0;
    for (double x : r) r0 += x * x;
    std::vector<double> ma(n, 0.0);
    std::vector<Candidate> out;
    out.reserve(mid1.size() * mid2.size());
    for (std::size_t a = 0; a < mid1.size(); ++a) {
        double ta = 0.0;
        for (const auto& [i, c] : mid1[a].v) {
            ma[i] = c;
            ta += (r[i] - c) * (r[i] - c) - r[i] * r[i];
        }
        for (std::size_t b = 0; b < mid2.size(); ++b) {
            double tb = 0.0;
            for (const auto& [i, c] : mid2[b].v) {
                const double before = r[i] - ma[i], after = before - c;
                tb += after * after - before * before;
            }
            Candidate cand;
            cand.a = a;
            cand.b = b;
            cand.loss = r0 + ta + tb;
            cand.profile = base1.hash + base2.hash + mid1[a].hash + mid2[b].hash;
            cand.outside = mid1[a].outside + mid2[b].outside;
            out.push_back(cand);
        }
        for (const auto& [i, c] : mid1[a].v) ma[i] = 0.0;
    }
    return out;
}

// Candidates grouped by IDENTICAL predicted count profile (M6-B: exact). The 64-bit linear hash is only a lookup key
// (masked by `hash_mask`; a test sets 0 to force every candidate into one bucket); within a bucket, candidates are split
// by equality of the complete predicted profile M1[a] + M2[b] (the base is shared by every candidate). A group is ONE
// count-model state: every member (sequence pair, allocation) stays attached, and a multi-member group never establishes a
// unique allocation or genotype. `collisions` counts buckets that held more than one distinct profile.
inline Sparse sum_sparse(const Sparse& x, const Sparse& y) {
    Sparse out;
    std::size_t i = 0, j = 0;
    while (i < x.v.size() || j < y.v.size()) {
        if (j == y.v.size() || (i < x.v.size() && x.v[i].first < y.v[j].first)) out.v.push_back(x.v[i++]);
        else if (i == x.v.size() || y.v[j].first < x.v[i].first) out.v.push_back(y.v[j++]);
        else { out.v.push_back({x.v[i].first, x.v[i].second + y.v[j].second}); ++i; ++j; }
    }
    return out;
}
inline std::vector<std::vector<std::size_t>> group_profiles(const std::vector<Candidate>& c, const std::vector<Sparse>& mid1, const std::vector<Sparse>& mid2,
                                                            std::uint64_t hash_mask = ~std::uint64_t(0), std::size_t* collisions = nullptr) {
    std::map<std::uint64_t, std::vector<std::size_t>> bucket;
    for (std::size_t i = 0; i < c.size(); ++i) bucket[c[i].profile & hash_mask].push_back(i);
    std::vector<std::vector<std::size_t>> groups;
    std::size_t multi = 0;
    for (auto& [h, members] : bucket) {
        std::vector<std::pair<Sparse, std::vector<std::size_t>>> distinct;
        for (std::size_t i : members) {
            const Sparse p = sum_sparse(mid1[c[i].a], mid2[c[i].b]);
            auto it = std::find_if(distinct.begin(), distinct.end(), [&](const auto& d) { return d.first.v == p.v; });
            if (it == distinct.end()) distinct.push_back({p, {i}});
            else it->second.push_back(i);
        }
        multi += distinct.size() > 1 ? 1 : 0;
        for (auto& d : distinct) groups.push_back(std::move(d.second));
    }
    if (collisions) *collisions = multi;
    std::sort(groups.begin(), groups.end(), [&](const auto& x, const auto& y) {
        if (c[x[0]].loss != c[y[0]].loss) return c[x[0]].loss < c[y[0]].loss;
        if (c[x[0]].profile != c[y[0]].profile) return c[x[0]].profile < c[y[0]].profile;   // the M6-A tie order
        return x[0] < y[0];
    });
    return groups;
}

// Independent scalar oracle: materialise both complete haplotypes, count, and sum over the whole dictionary.
inline double naive_loss(const std::vector<double>& y, const std::string& h1, const std::string& h2, const Dictionary& d) {
    std::vector<double> pred(y.size(), 0.0);
    for (const std::string* h : {&h1, &h2}) {
        for_each_kmer(*h, [&](std::uint64_t k, std::size_t) {
            const long i = d.find(k);
            if (i >= 0) pred[i] += 1.0;
        });
    }
    double s = 0.0;
    for (std::size_t i = 0; i < y.size(); ++i) s += (y[i] - pred[i]) * (y[i] - pred[i]);
    return s;
}

}  // namespace countfit
