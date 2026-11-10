#pragma once

// Shared k-mer primitives for read recruitment: rolling canonical k-mers, closed syncmers chosen by a
// hash of canonical s-mers, and a low-complexity score. Inline, so the per-base loops compile into
// their callers. k is at most 31 (a k-mer is 2 bits per base in one uint64_t).

#include <array>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <string>
#include <utility>
#include <vector>

namespace panvar::kmer {

constexpr std::size_t kMaxK = 31;

// A/C/G/T in either case -> 0..3; anything else -> -1, which breaks the window.
inline int base_code(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1;
    }
}

inline std::uint64_t mask_of(std::size_t k) { return k >= 32 ? ~0ULL : (1ULL << (2 * k)) - 1ULL; }

// The reverse complement of a k-mer's code.
inline std::uint64_t reverse_complement(std::uint64_t code, std::size_t k) {
    std::uint64_t out = 0;
    for (std::size_t i = 0; i < k; ++i) {
        out = (out << 2) | (3ULL - (code & 3ULL));
        code >>= 2;
    }
    return out;
}

inline std::uint64_t canonical(std::uint64_t code, std::size_t k) {
    const std::uint64_t rc = reverse_complement(code, k);
    return code < rc ? code : rc;
}

// A fixed 64-bit mixer (the splitmix64 finalizer): a hash with no seed, the same on every platform.
inline std::uint64_t mix64(std::uint64_t x) {
    x ^= x >> 30; x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27; x *= 0x94d049bb133111ebULL;
    x ^= x >> 31;
    return x;
}

// Calls f(code, start, forward) for every canonical k-mer of seq[0, n), in order: code is the smaller of
// the forward and reverse-complement codes, start its 0-based offset, forward whether code is the forward
// strand's. A base other than ACGT breaks the window.
template <class F>
void for_each_canonical_kmer(const char* seq, std::size_t n, std::size_t k, F&& f) {
    if (k == 0 || k > kMaxK || n < k) return;
    const std::uint64_t mask = mask_of(k);
    const std::size_t rc_shift = 2 * (k - 1);
    std::uint64_t fwd = 0, rev = 0;
    std::size_t filled = 0;
    for (std::size_t i = 0; i < n; ++i) {
        const int b = base_code(seq[i]);
        if (b < 0) { fwd = rev = 0; filled = 0; continue; }
        fwd = ((fwd << 2) | static_cast<std::uint64_t>(b)) & mask;
        rev = (rev >> 2) | (static_cast<std::uint64_t>(3 - b) << rc_shift);
        if (filled < k) ++filled;
        if (filled == k) f(fwd <= rev ? fwd : rev, i + 1 - k, fwd <= rev);
    }
}

// A k-mer is a closed syncmer when the smallest hash among its canonical s-mers sits at its first or
// last s-mer position (Edgar 2021, https://doi.org/10.7717/peerj.10805). Hashing canonical s-mers, not
// ranking their raw codes, keeps the density even across base composition and makes the choice the same
// on both strands. Reference implementation, one k-mer at a time; `code` may be either strand.
inline bool is_closed_syncmer(std::uint64_t code, std::size_t k, std::size_t s) {
    if (s == 0 || s > k || k > kMaxK) return false;
    if (s == k) return true;
    const std::uint64_t smask = mask_of(s);
    const std::size_t w = k - s + 1;
    std::uint64_t best = ~0ULL, first = 0, last = 0;
    for (std::size_t j = 0; j < w; ++j) {
        const std::uint64_t sub = (code >> (2 * (k - s - j))) & smask;
        const std::uint64_t h = mix64(canonical(sub, s));
        if (j == 0) first = h;
        if (j + 1 == w) last = h;
        if (h < best) best = h;
    }
    return first == best || last == best;
}

// Calls f(code, start, forward) for every k-mer of seq[0, n) that is a closed syncmer (as
// is_closed_syncmer), in order, in linear time: the s-mer hashes of the current window are kept in a
// monotone queue.
template <class F>
void for_each_closed_syncmer(const char* seq, std::size_t n, std::size_t k, std::size_t s, F&& f) {
    if (k == 0 || k > kMaxK || s == 0 || s > k || n < k) return;
    const std::uint64_t kmask = mask_of(k), smask = mask_of(s);
    const std::size_t krc = 2 * (k - 1), src = 2 * (s - 1), w = k - s + 1;
    std::uint64_t kf = 0, kr = 0, sf = 0, sr = 0;
    std::size_t filled = 0;
    std::vector<std::uint64_t> ring(w);            // s-mer hashes by start position modulo w
    std::deque<std::pair<std::uint64_t, std::size_t>> mins;   // (hash, s-mer start), hashes increasing
    for (std::size_t i = 0; i < n; ++i) {
        const int b = base_code(seq[i]);
        if (b < 0) { kf = kr = sf = sr = 0; filled = 0; mins.clear(); continue; }
        const auto ub = static_cast<std::uint64_t>(b), cb = static_cast<std::uint64_t>(3 - b);
        kf = ((kf << 2) | ub) & kmask; kr = (kr >> 2) | (cb << krc);
        sf = ((sf << 2) | ub) & smask; sr = (sr >> 2) | (cb << src);
        ++filled;
        if (filled >= s) {
            const std::size_t sstart = i + 1 - s;
            const std::uint64_t h = mix64(sf <= sr ? sf : sr);
            ring[sstart % w] = h;
            while (!mins.empty() && mins.back().first > h) mins.pop_back();
            mins.emplace_back(h, sstart);
        }
        if (filled >= k) {
            const std::size_t kstart = i + 1 - k;
            while (mins.front().second < kstart) mins.pop_front();
            const std::uint64_t best = mins.front().first;
            if (ring[kstart % w] == best || ring[(kstart + w - 1) % w] == best)
                f(kf <= kr ? kf : kr, kstart, kf <= kr);
        }
    }
}

// Low complexity of a k-mer as DUST scores it (Morgulis et al. 2006): over its k-2 overlapping triplets,
// the sum of c(c-1)/2 for each triplet's count c, divided by (k-3). 0 when no triplet repeats; a
// homopolymer 31-mer scores 14.5, a dinucleotide repeat about 6.5, random sequence about 0.2.
inline double dust_score(std::uint64_t code, std::size_t k) {
    if (k < 4) return 0.0;
    std::array<int, 64> count{};
    long pairs = 0;
    for (std::size_t j = 0; j + 3 <= k; ++j) {
        const int t = static_cast<int>((code >> (2 * (k - 3 - j))) & 63ULL);
        pairs += count[t]++;
    }
    return static_cast<double>(pairs) / static_cast<double>(k - 3);
}

// The ACGT string a code spells.
inline std::string decode(std::uint64_t code, std::size_t k) {
    std::string out(k, 'A');
    for (std::size_t i = 0; i < k; ++i) out[i] = "ACGT"[(code >> (2 * (k - 1 - i))) & 3ULL];
    return out;
}

} // namespace panvar::kmer
