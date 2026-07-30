#include "panvar/syncmer.hpp"

#include <algorithm>
#include <limits>

namespace panvar {

int base_bits(char c) {
    switch (c) {
        case 'A':
        case 'a':
            return 0;
        case 'C':
        case 'c':
            return 1;
        case 'G':
        case 'g':
            return 2;
        case 'T':
        case 't':
            return 3;
        default:
            return -1;
    }
}

std::size_t default_syncmer_s(std::size_t k) {
    return std::max<std::size_t>(1, std::min<std::size_t>(11, (k + 2) / 3));
}

std::vector<KmerOccurrence> collect_canonical_kmer_occurrences(
    const std::string& sequence,
    std::size_t k) {

    std::vector<KmerOccurrence> out;
    if (k == 0 || k > 31 || sequence.size() < k) {
        return out;
    }

    const std::size_t possible = sequence.size() - k + 1;
    out.reserve(possible);

    const std::uint64_t mask = (1ULL << (2 * k)) - 1ULL;
    const std::size_t rc_shift = 2 * (k - 1);
    std::uint64_t fwd = 0;
    std::uint64_t rev = 0;
    std::size_t filled = 0;

    for (std::size_t pos = 0; pos < sequence.size(); ++pos) {
        const int b = base_bits(sequence[pos]);
        if (b < 0) {
            fwd = 0;
            rev = 0;
            filled = 0;
            continue;
        }
        const std::uint64_t ub = static_cast<std::uint64_t>(b);
        const std::uint64_t cb = static_cast<std::uint64_t>(3 - b);
        fwd = ((fwd << 2) | ub) & mask;
        rev = (rev >> 2) | (cb << rc_shift);
        if (filled < k) {
            ++filled;
        }
        if (filled >= k) {
            out.push_back(KmerOccurrence{std::min(fwd, rev), pos + 1 - k});
        }
    }

    return out;
}

bool is_closed_syncmer(std::uint64_t code, std::size_t k, std::size_t s) {
    if (s == 0 || s > k) {
        return false;
    }
    if (s == k) {
        return true;
    }

    const std::uint64_t mask = (1ULL << (2 * s)) - 1ULL;
    const std::size_t last_offset = k - s;
    std::uint64_t best = std::numeric_limits<std::uint64_t>::max();
    bool best_at_end = false;

    for (std::size_t offset = 0; offset <= last_offset; ++offset) {
        const std::size_t shift = 2 * (k - s - offset);
        const std::uint64_t sub = (code >> shift) & mask;
        const bool at_end = offset == 0 || offset == last_offset;
        if (sub < best) {
            best = sub;
            best_at_end = at_end;
        } else if (sub == best && at_end) {
            best_at_end = true;
        }
    }

    return best_at_end;
}

std::vector<KmerOccurrence> collect_syncmers(
    const std::string& sequence,
    std::size_t k,
    std::size_t s) {

    std::vector<KmerOccurrence> out;
    if (k == 0 || k > 31 || sequence.size() < k) {
        return out;
    }
    if (s == 0) {
        s = default_syncmer_s(k);
    }

    const std::uint64_t mask = (1ULL << (2 * k)) - 1ULL;
    const std::size_t rc_shift = 2 * (k - 1);
    std::uint64_t fwd = 0;
    std::uint64_t rev = 0;
    std::size_t filled = 0;

    out.reserve(sequence.size() / 4 + 1);
    for (std::size_t pos = 0; pos < sequence.size(); ++pos) {
        const int b = base_bits(sequence[pos]);
        if (b < 0) {
            fwd = 0;
            rev = 0;
            filled = 0;
            continue;
        }
        const std::uint64_t ub = static_cast<std::uint64_t>(b);
        const std::uint64_t cb = static_cast<std::uint64_t>(3 - b);
        fwd = ((fwd << 2) | ub) & mask;
        rev = (rev >> 2) | (cb << rc_shift);
        if (filled < k) {
            ++filled;
        }
        if (filled >= k) {
            const std::uint64_t code = std::min(fwd, rev);
            if (is_closed_syncmer(code, k, s)) {
                out.push_back(KmerOccurrence{code, pos + 1 - k});
            }
        }
    }

    return out;
}

std::string decode_kmer(std::uint64_t code, std::size_t k) {
    std::string out(k, 'A');
    for (std::size_t i = 0; i < k; ++i) {
        const std::size_t shift = 2 * (k - 1 - i);
        const std::uint64_t b = (code >> shift) & 0x3ULL;
        switch (b) {
            case 0:
                out[i] = 'A';
                break;
            case 1:
                out[i] = 'C';
                break;
            case 2:
                out[i] = 'G';
                break;
            default:
                out[i] = 'T';
                break;
        }
    }
    return out;
}

} // namespace panvar
