#pragma once

// An open-addressing table from k-mer codes to values: lookups probe linearly from a mixed hash. Codes
// of k <= 31 are below 2^62, so the all-ones key marks an empty slot.

#include <cstddef>
#include <cstdint>
#include <vector>

namespace panvar {

// Sized by reserve() to at least twice the keys it will hold.
template <typename V>
class FlatMap {
public:
    void reserve(std::size_t n) {
        std::size_t cap = 16;
        while (cap < 2 * n) cap <<= 1;
        keys_.assign(cap, kEmpty);
        vals_.assign(cap, V{});
        mask_ = cap - 1;
    }
    // The slot of `key`, claimed if absent; existing tells which.
    V& slot(std::uint64_t key, bool& existing) {
        std::size_t i = hash(key) & mask_;
        while (keys_[i] != kEmpty && keys_[i] != key) i = (i + 1) & mask_;
        existing = keys_[i] == key;
        keys_[i] = key;
        return vals_[i];
    }
    const V* find(std::uint64_t key) const {
        for (std::size_t i = hash(key) & mask_;; i = (i + 1) & mask_) {
            if (keys_[i] == key) return &vals_[i];
            if (keys_[i] == kEmpty) return nullptr;
        }
    }

private:
    static constexpr std::uint64_t kEmpty = ~0ULL;
    static std::size_t hash(std::uint64_t x) {
        x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33;
        return static_cast<std::size_t>(x);
    }
    std::vector<std::uint64_t> keys_;
    std::vector<V> vals_;
    std::size_t mask_ = 0;
};

}  // namespace panvar
