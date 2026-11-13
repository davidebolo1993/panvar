#pragma once

// A static Elias-Fano sequence of strictly increasing 64-bit values below 2^universe_bits, with a
// membership/rank query: each value splits into high bits, stored in unary in a bit vector (about 2 bits
// per value), and low bits, stored packed (universe_bits - log2 n of them). Finding a value is one select
// on the high bits and a scan of the few values sharing them. Serialized into a flat buffer of 64-bit
// words that is read in place (memory-mapped).

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace panvar::ef {

constexpr std::size_t kSelectSample = 256;   // one sample per this many zeros of the high bit vector

inline unsigned popcount64(std::uint64_t x) { return static_cast<unsigned>(__builtin_popcountll(x)); }
inline unsigned ctz64(std::uint64_t x) { return static_cast<unsigned>(__builtin_ctzll(x)); }   // x != 0
inline std::size_t ceil_log2(std::uint64_t x) { return x <= 1 ? 0 : 64 - static_cast<std::size_t>(__builtin_clzll(x - 1)); }

// Builds the serialized words for `values` (strictly increasing, each < 2^universe_bits).
// Layout: [n, low_bits, max_high, high_words, sample_count, low_words] then high, samples, low.
inline std::vector<std::uint64_t> build(const std::vector<std::uint64_t>& values, std::size_t universe_bits) {
    if (universe_bits == 0 || universe_bits > 63) throw std::runtime_error("elias-fano: the universe must be 1 to 63 bits");
    const std::uint64_t n = values.size();
    const std::size_t lb = n == 0 ? 0 : (universe_bits > ceil_log2(n) ? universe_bits - ceil_log2(n) : 0);
    const std::uint64_t max_high = n == 0 ? 0 : values.back() >> lb;
    const std::uint64_t high_bits = n + max_high + 1;
    const std::uint64_t high_words = (high_bits + 63) / 64;
    const std::uint64_t low_words = (n * lb + 63) / 64;
    std::vector<std::uint64_t> high(high_words, 0), low(low_words, 0);
    const std::uint64_t lmask = lb == 0 ? 0 : (lb == 64 ? ~0ULL : (1ULL << lb) - 1);
    for (std::uint64_t i = 0; i < n; ++i) {
        if (i > 0 && values[i] <= values[i - 1]) throw std::runtime_error("elias-fano: values not strictly increasing");
        if ((values[i] >> universe_bits) != 0) throw std::runtime_error("elias-fano: a value outside the universe");
        const std::uint64_t h = values[i] >> lb, p = h + i;
        high[p / 64] |= 1ULL << (p % 64);
        if (lb > 0) {
            const std::uint64_t v = values[i] & lmask, bit = i * lb;
            low[bit / 64] |= v << (bit % 64);
            if (bit % 64 + lb > 64) low[bit / 64 + 1] |= v >> (64 - bit % 64);
        }
    }
    // positions of zero number 0, kSelectSample, 2*kSelectSample, ...
    std::vector<std::uint64_t> samples;
    std::uint64_t zeros = 0;
    for (std::uint64_t p = 0; p < high_bits; ++p) {
        if (!((high[p / 64] >> (p % 64)) & 1ULL)) {
            if (zeros % kSelectSample == 0) samples.push_back(p);
            ++zeros;
        }
    }
    std::vector<std::uint64_t> out{n, lb, max_high, high_words, samples.size(), low_words};
    out.insert(out.end(), high.begin(), high.end());
    out.insert(out.end(), samples.begin(), samples.end());
    out.insert(out.end(), low.begin(), low.end());
    return out;
}

// A read-only view over serialized words.
class View {
public:
    View() = default;
    explicit View(const std::uint64_t* w) {
        n_ = w[0]; lb_ = w[1]; max_high_ = w[2];
        const std::uint64_t hw = w[3], sc = w[4], lw = w[5];
        high_ = w + 6; high_words_ = hw;
        samples_ = high_ + hw;
        low_ = samples_ + sc;
        words_ = 6 + hw + sc + lw;
        lmask_ = lb_ == 0 ? 0 : (lb_ == 64 ? ~0ULL : (1ULL << lb_) - 1);
    }
    std::uint64_t size() const { return n_; }
    std::uint64_t words() const { return words_; }   // serialized length in 64-bit words

    std::uint64_t low(std::uint64_t i) const {
        if (lb_ == 0) return 0;
        const std::uint64_t bit = i * lb_;
        std::uint64_t v = low_[bit / 64] >> (bit % 64);
        if (bit % 64 + lb_ > 64) v |= low_[bit / 64 + 1] << (64 - bit % 64);
        return v & lmask_;
    }
    std::uint64_t at(std::uint64_t i) const {   // the i-th value (linear-time select1, for tests)
        std::uint64_t ones = 0;
        for (std::uint64_t wi = 0; wi < high_words_; ++wi) {
            const std::uint64_t c = static_cast<std::uint64_t>(popcount64(high_[wi]));
            if (ones + c > i) {
                std::uint64_t word = high_[wi];
                for (std::uint64_t k = ones; k < i; ++k) word &= word - 1;
                const std::uint64_t p = wi * 64 + static_cast<std::uint64_t>(ctz64(word));
                return ((p - i) << lb_) | low(i);
            }
            ones += c;
        }
        throw std::out_of_range("elias-fano: index");
    }
    // Calls f(index, value) for every value, in order.
    template <class F>
    void for_each_value(F&& f) const {
        std::uint64_t i = 0;
        for (std::uint64_t wi = 0; wi < high_words_ && i < n_; ++wi) {
            std::uint64_t word = high_[wi];
            while (word != 0 && i < n_) {
                const std::uint64_t p = wi * 64 + ctz64(word);
                f(i, ((p - i) << lb_) | low(i));
                ++i;
                word &= word - 1;
            }
        }
    }
    // Index of `value` in the sequence, or -1 when absent.
    std::int64_t find(std::uint64_t value) const {
        if (n_ == 0) return -1;
        const std::uint64_t h = value >> lb_;
        if (h > max_high_) return -1;
        const std::uint64_t begin = h == 0 ? 0 : select0(h - 1) + 1 - h;
        const std::uint64_t end = select0(h) - h;
        const std::uint64_t lo = value & lmask_;
        for (std::uint64_t i = begin; i < end; ++i) {
            const std::uint64_t v = low(i);
            if (v == lo) return static_cast<std::int64_t>(i);
            if (v > lo) break;
        }
        return -1;
    }

    // Sequential reading of the values in order, without select.
    class Cursor {
    public:
        explicit Cursor(const View& v) : v_(&v), word_(v.n_ > 0 ? v.high_[0] : 0) {}
        bool next(std::uint64_t& value) {
            if (i_ >= v_->n_) return false;
            while (word_ == 0) word_ = v_->high_[++wi_];
            const std::uint64_t p = wi_ * 64 + ctz64(word_);
            word_ &= word_ - 1;
            value = ((p - i_) << v_->lb_) | v_->low(i_);
            ++i_;
            return true;
        }
    private:
        const View* v_;
        std::uint64_t i_ = 0, wi_ = 0, word_ = 0;
    };

private:
    // Position of zero number j (0-based) in the high bit vector.
    std::uint64_t select0(std::uint64_t j) const {
        std::uint64_t p = samples_[j / kSelectSample];
        std::uint64_t left = j % kSelectSample;
        std::uint64_t wi = p / 64;
        std::uint64_t word = ~high_[wi] & (~0ULL << (p % 64));
        while (true) {
            const std::uint64_t c = static_cast<std::uint64_t>(popcount64(word));
            if (left < c) {
                for (std::uint64_t k = 0; k < left; ++k) word &= word - 1;
                return wi * 64 + static_cast<std::uint64_t>(ctz64(word));
            }
            left -= c;
            word = ~high_[++wi];
        }
    }

    std::uint64_t n_ = 0, lb_ = 0, max_high_ = 0, high_words_ = 0, words_ = 0, lmask_ = 0;
    const std::uint64_t* high_ = nullptr;
    const std::uint64_t* samples_ = nullptr;
    const std::uint64_t* low_ = nullptr;
};

} // namespace panvar::ef
