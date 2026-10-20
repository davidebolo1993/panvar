#include "panvar/junction_proposals.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace panvar {

namespace {

int base_code(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1;
    }
}

// Sorted (inventory index, count) entries of one window.
using Window = std::vector<std::pair<std::uint32_t, std::uint32_t>>;

}  // namespace

std::vector<std::pair<std::uint64_t, std::size_t>> junction_kmers(const std::string& sequence) {
    std::vector<std::pair<std::uint64_t, std::size_t>> out;
    const std::size_t k = kJunctionK;
    if (sequence.size() < k) return out;
    const std::uint64_t mask = (1ULL << (2 * k)) - 1ULL;
    const std::size_t rc_shift = 2 * (k - 1);
    std::uint64_t fwd = 0, rev = 0;
    std::size_t filled = 0;
    out.reserve(sequence.size() - k + 1);
    for (std::size_t pos = 0; pos < sequence.size(); ++pos) {
        const int b = base_code(sequence[pos]);
        if (b < 0) {
            fwd = rev = 0;
            filled = 0;
            continue;
        }
        fwd = ((fwd << 2) | static_cast<std::uint64_t>(b)) & mask;
        rev = (rev >> 2) | (static_cast<std::uint64_t>(3 - b) << rc_shift);
        if (filled < k) ++filled;
        if (filled >= k) out.push_back({std::min(fwd, rev), pos + 1 - k});
    }
    return out;
}

JunctionReadKmers JunctionReadKmers::from_fragments(const std::vector<ReadPair>& fragments) {
    JunctionReadKmers t;
    std::vector<std::uint64_t> all;
    for (const ReadPair& f : fragments) {
        for (const SequencingRead* r : {&f.first, &f.second}) {
            ++t.reads;
            for (const auto& kmer : junction_kmers(r->bases)) all.push_back(kmer.first);
        }
    }
    t.occurrences = all.size();
    std::sort(all.begin(), all.end());
    for (std::size_t i = 0; i < all.size();) {
        std::size_t j = i + 1;
        while (j < all.size() && all[j] == all[i]) ++j;
        t.codes.push_back(all[i]);
        t.counts.push_back(static_cast<std::uint32_t>(std::min<std::size_t>(j - i, 0xffffffffu)));
        i = j;
    }
    return t;
}

std::uint32_t JunctionReadKmers::count(std::uint64_t code) const {
    const auto it = std::lower_bound(codes.begin(), codes.end(), code);
    return it != codes.end() && *it == code ? counts[static_cast<std::size_t>(it - codes.begin())] : 0u;
}

bool junction_block_eligible(const std::vector<double>& pair_emissions, std::size_t limit) {
    if (pair_emissions.size() <= limit) return false;
    for (const double e : pair_emissions) {
        if (e != pair_emissions.front()) return false;
    }
    return true;
}

std::vector<JunctionPairScore> junction_rank_all(const std::vector<std::string>& alleles, const JunctionContext& context,
                                                 const JunctionReadKmers& reads, double rate, double background,
                                                 JunctionProposal* summary) {
    if (!(rate > 0.0) || !(background > 0.0) || !std::isfinite(rate) || !std::isfinite(background)) {
        throw std::invalid_argument("junction proposals need a finite positive rate and background");
    }
    const std::size_t n = alleles.size();
    const std::size_t flank = kJunctionK - 1;
    for (int h = 0; h < 2; ++h) {
        if (context.incumbent[h] >= n) throw std::invalid_argument("junction proposals: incumbent allele outside the block");
    }
    // windows (codes first), then the inventory
    std::vector<std::vector<std::uint64_t>> raw[2];
    std::vector<std::uint64_t> inventory;
    for (int h = 0; h < 2; ++h) {
        const std::string& pre = context.prefix[h];
        const std::string& suf = context.suffix[h];
        const std::string left = pre.substr(pre.size() > flank ? pre.size() - flank : 0);
        const std::string right = suf.substr(0, std::min(flank, suf.size()));
        raw[h].resize(n);
        for (std::size_t a = 0; a < n; ++a) {
            for (const auto& kmer : junction_kmers(left + alleles[a] + right)) raw[h][a].push_back(kmer.first);
            std::sort(raw[h][a].begin(), raw[h][a].end());
            inventory.insert(inventory.end(), raw[h][a].begin(), raw[h][a].end());
        }
    }
    std::sort(inventory.begin(), inventory.end());
    inventory.erase(std::unique(inventory.begin(), inventory.end()), inventory.end());
    const auto index_of = [&](std::uint64_t code) {
        return static_cast<std::uint32_t>(std::lower_bound(inventory.begin(), inventory.end(), code) - inventory.begin());
    };
    std::vector<Window> win[2];
    for (int h = 0; h < 2; ++h) {
        win[h].resize(n);
        for (std::size_t a = 0; a < n; ++a) {
            const auto& r = raw[h][a];
            for (std::size_t i = 0; i < r.size();) {
                std::size_t j = i + 1;
                while (j < r.size() && r[j] == r[i]) ++j;
                win[h][a].push_back({index_of(r[i]), static_cast<std::uint32_t>(j - i)});
                i = j;
            }
        }
    }
    // elsewhere: inventory K-mers in both incumbent haplotypes outside the window span of each
    std::vector<std::uint32_t> elsewhere(inventory.size(), 0);
    for (int h = 0; h < 2; ++h) {
        const std::string& pre = context.prefix[h];
        const std::string full = pre + alleles[context.incumbent[h]] + context.suffix[h];
        const std::size_t lo = pre.size() > flank ? pre.size() - flank : 0;
        const std::size_t hi = pre.size() + alleles[context.incumbent[h]].size();
        for (const auto& kmer : junction_kmers(full)) {
            if (kmer.second >= lo && kmer.second < hi) continue;
            const auto it = std::lower_bound(inventory.begin(), inventory.end(), kmer.first);
            if (it != inventory.end() && *it == kmer.first) ++elsewhere[static_cast<std::size_t>(it - inventory.begin())];
        }
    }
    std::vector<double> observed(inventory.size());
    std::size_t seen = 0, else_n = 0;
    for (std::size_t i = 0; i < inventory.size(); ++i) {
        observed[i] = reads.count(inventory[i]);
        seen += observed[i] > 0 ? 1 : 0;
        else_n += elsewhere[i] > 0 ? 1 : 0;
    }
    const auto f = [&](std::uint32_t m, double copies) {
        const double lambda = rate * copies + background;
        return observed[m] * std::log(lambda) - lambda;
    };
    // exact delta over the touched K-mers, merged in ascending inventory order
    const auto score = [&](const Window& w1, const Window& w2) {
        double s = 0.0;
        std::size_t i = 0, j = 0;
        while (i < w1.size() || j < w2.size()) {
            std::uint32_t m;
            double c = 0.0;
            if (j >= w2.size() || (i < w1.size() && w1[i].first < w2[j].first)) {
                m = w1[i].first;
                c = w1[i++].second;
            } else if (i >= w1.size() || w2[j].first < w1[i].first) {
                m = w2[j].first;
                c = w2[j++].second;
            } else {
                m = w1[i].first;
                c = static_cast<double>(w1[i++].second) + w2[j++].second;
            }
            s += f(m, elsewhere[m] + c) - f(m, elsewhere[m]);
        }
        return s;
    };
    std::vector<JunctionPairScore> rows;
    rows.reserve(n * (n + 1) / 2);
    for (std::size_t x = 0; x < n; ++x) {
        for (std::size_t y = x; y < n; ++y) {
            JunctionPairScore r;
            r.x = x;
            r.y = y;
            r.ab = score(win[0][x], win[1][y]);
            r.ba = x == y ? r.ab : score(win[0][y], win[1][x]);
            r.best = std::max(r.ab, r.ba);
            rows.push_back(r);
        }
    }
    // score order, chained tie classes, then CONTENT order inside each class
    std::stable_sort(rows.begin(), rows.end(), [](const JunctionPairScore& a, const JunctionPairScore& b) {
        if (a.best != b.best) return a.best > b.best;
        return std::tie(a.x, a.y) < std::tie(b.x, b.y);
    });
    std::size_t cls = 0;
    double head = 0.0;
    for (std::size_t i = 0; i < rows.size(); ++i) {
        if (i == 0 || head - rows[i].best > kJunctionTieNats) {
            ++cls;
            head = rows[i].best;
        }
        rows[i].tie_class = cls;
    }
    const auto content = [&](const JunctionPairScore& r) {
        const std::string& a = alleles[r.x];
        const std::string& b = alleles[r.y];
        return a <= b ? std::make_pair(&a, &b) : std::make_pair(&b, &a);
    };
    std::stable_sort(rows.begin(), rows.end(), [&](const JunctionPairScore& a, const JunctionPairScore& b) {
        if (a.tie_class != b.tie_class) return a.tie_class < b.tie_class;
        const auto ca = content(a), cb = content(b);
        if (*ca.first != *cb.first) return *ca.first < *cb.first;
        if (*ca.second != *cb.second) return *ca.second < *cb.second;
        return std::tie(a.x, a.y) < std::tie(b.x, b.y);   // identical content only (equal sequences under two labels)
    });
    if (summary != nullptr) {
        summary->inventory = inventory.size();
        summary->observed = seen;
        summary->elsewhere = else_n;
        summary->candidate_pairs = rows.size();
        summary->distinct_scores = cls;
    }
    return rows;
}

JunctionProposal junction_local_pairs(const std::vector<std::string>& alleles, const JunctionContext& context,
                                      const JunctionReadKmers& reads, double rate, double background, std::size_t limit) {
    JunctionProposal out;
    const std::vector<JunctionPairScore> rows = junction_rank_all(alleles, context, reads, rate, background, &out);
    const std::size_t take = std::min(limit, rows.size());
    for (std::size_t i = 0; i < take; ++i) {
        out.pairs.push_back({rows[i].x, rows[i].y});
        out.scores.push_back(rows[i].best);
    }
    if (take > 0) {
        const std::size_t c = rows[take - 1].tie_class;
        out.cutoff_tie_size = static_cast<std::size_t>(std::count_if(rows.begin(), rows.end(), [&](const auto& r) { return r.tie_class == c; }));
        out.cutoff_tie_crosses = take < rows.size() && rows[take].tie_class == c;
    }
    return out;
}

JunctionProposalCache::Key JunctionProposalCache::key(std::size_t block, const std::vector<std::size_t>& first,
                                                      const std::vector<std::size_t>& second) {
    return first <= second ? Key{block, {first, second}} : Key{block, {second, first}};
}

const JunctionProposal* JunctionProposalCache::find(const Key& key) {
    const auto it = entries_.find(key);
    if (it == entries_.end()) {
        ++misses_;
        return nullptr;
    }
    ++hits_;
    return &it->second;
}

const JunctionProposal& JunctionProposalCache::insert(const Key& key, JunctionProposal value) {
    if (capacity_ != 0 && entries_.size() >= capacity_ && !order_.empty()) {
        entries_.erase(order_.front());
        order_.erase(order_.begin());
    }
    order_.push_back(key);
    return entries_[key] = std::move(value);
}

}  // namespace panvar
