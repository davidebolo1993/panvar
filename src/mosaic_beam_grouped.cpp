// The grouped-prefix proposer. THIS IS WHAT PRODUCTION PROPOSES WITH. propose_mosaic_diplotypes is
// retained unchanged as the fixture oracle this is checked against, and is not called by the genotyper.
//
// Its state_limit counts allele-prefix PAIRS, not (prefix, template) states.
//
// The unit of retention is the allele-prefix PAIR, and each pair carries the whole distribution over
// current template pairs rather than one attribution of it. States sharing a prefix pair are different
// latent template histories spelling the same allele prefix, and the current template is what the next
// block's transition reads, so dropping all but one of them approximates that prefix's future.
//
// It is affordable because the joint factorises. Within a prefix pair every history spells the same
// alleles, so the diploid emission is a single scalar for the whole group -- and the emission is the
// only coupling between the two haplotypes. Hence
//
//     W_p(i, j) = E(p1, p2) * u_p1(i) * u_p2(j)
//
// with u_p(i) the mass of template paths spelling prefix p and ending on template i.
//
// Two structural facts keep the cost linear in T rather than A*T:
//   * Li-Stephens is low rank. Staying costs (1-r) + r/T and moving costs r/T, so
//         u'(i') = (1 - r) * u_p(i')  +  (r / T) * sum_i u_p(i)
//     and the rank-1 term is one scalar -- the prefix's running total -- added to every entry.
//   * Each template spells exactly ONE allele at a block. So a single scan over the T templates
//     produces every child of a prefix at once, scattered into its allele class, and all A children
//     together occupy exactly T entries. One scan, not A scans.
//
// Storage is split by LIFETIME and each part is predicted and gated separately, because pooling them
// hides exactly the defect an earlier version of this file had: it kept every layer's numerical arrays
// alive while claiming a rolling budget, so the real cost was O(B * 2LT) against a claimed O(2LT).
#include "panvar/mosaic_beam.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace panvar {
namespace {

constexpr double kNegInf = -std::numeric_limits<double>::infinity();
constexpr std::uint32_t kNoParent = std::numeric_limits<std::uint32_t>::max();
constexpr std::uint64_t kCountCap = std::numeric_limits<std::uint64_t>::max();

double log_add(double lhs, double rhs) {
    if (lhs == kNegInf) return rhs;
    if (rhs == kNegInf) return lhs;
    const double hi = std::max(lhs, rhs);
    return hi + std::log1p(std::exp(std::min(lhs, rhs) - hi));
}

void add_sat(std::uint64_t& total, bool& saturated, std::uint64_t amount, bool amount_saturated) {
    saturated = saturated || amount_saturated;
    if (amount > kCountCap - total) {
        total = kCountCap;
        saturated = true;
    } else {
        total += amount;
    }
}

std::uint64_t mul_sat(std::uint64_t lhs, std::uint64_t rhs, bool& saturated) {
    if (lhs == 0 || rhs == 0) return 0;
    if (lhs > kCountCap / rhs) {
        saturated = true;
        return kCountCap;
    }
    return lhs * rhs;
}

// ROLLING storage: only the parent layer's and the child layer's numerical arrays are ever live, and
// the parent's are released as soon as the child is complete and compacted. Vectors are compact --
// prefix p owns entries [offset[p], offset[p+1]), one per template spelling p's last allele, in
// template order -- so all the children of one parent occupy T entries between them.
struct PrefixVectors {
    std::vector<std::uint32_t> offset;   // size() + 1 entries
    std::vector<double> sum;             // log u_p(i), sum-product
    std::vector<double> best;            // log max-product over single histories
    std::vector<std::uint64_t> count;    // template paths ending at i
    std::vector<double> log_total;       // log sum_i u_p(i)
    std::vector<double> best_total;      // max_i best(i)
    std::vector<std::uint64_t> count_total;
    std::vector<char> count_saturated;
    std::vector<std::uint32_t> allele;   // last allele, needed to scatter back at the next block

    std::size_t size() const { return allele.size(); }
    std::uint64_t cells() const { return static_cast<std::uint64_t>(sum.size()); }
    void reset() {
        offset.clear();
        sum.clear();
        best.clear();
        count.clear();
        log_total.clear();
        best_total.clear();
        count_total.clear();
        count_saturated.clear();
        allele.clear();
        offset.push_back(0);
    }
};

// One rolling cell costs a sum, a max and a count.
constexpr std::uint64_t kCellBytes = sizeof(double) * 2 + sizeof(std::uint64_t);

// PERSISTENT storage, deliberately minimal: (parent, allele) per retained prefix is all that is needed
// to spell a candidate's allele path at the end. Compacted after every layer to the prefixes retained
// pairs actually reference, so it grows by at most 2L records per block.
struct AncestryLayer {
    std::vector<std::uint32_t> parent;   // index into the previous layer, kNoParent at block 0
    std::vector<std::uint32_t> allele;
    std::vector<char> matches_first;     // target's first path, tracked incrementally
    std::vector<char> matches_second;

    std::size_t size() const { return allele.size(); }
};

constexpr std::uint64_t kAncestryBytes = sizeof(std::uint32_t) * 2 + sizeof(char) * 2;

// The trace target is an UNORDERED diplotype, and pairs are canonicalised by prefix index, so the
// target's two prefixes can occupy either slot and the canonical order can flip between layers.
// Checking a single orientation reports the truth as lost at the first layer where it does.
bool pair_is_target(const AncestryLayer& layer, std::uint32_t first, std::uint32_t second) {
    return (layer.matches_first[first] != 0 && layer.matches_second[second] != 0) ||
           (layer.matches_second[first] != 0 && layer.matches_first[second] != 0);
}

// A pair-prefix: the unit this proposer retains and ranks.
struct PairState {
    std::uint32_t first = 0;    // index into the current layer
    std::uint32_t second = 0;   // first <= second, always
    double emission = 0.0;      // sum of diploid emissions along the pair, shared by every history
    double mass = kNegInf;      // sum-product, heterozygous factor of two included
    double best = kNegInf;      // max-product; the swapped history scores the same, so no factor here
    double best_transition = 0.0;
    std::uint64_t histories = 0;
    bool saturated = false;
};

// Sum-product mass first, then the best single history, then prefix identity. Total, so no tie is ever
// left to arrival order.
bool better_pair(const PairState& lhs, const PairState& rhs) {
    if (lhs.mass != rhs.mass) return lhs.mass > rhs.mass;
    if (lhs.best != rhs.best) return lhs.best > rhs.best;
    return std::tie(lhs.first, lhs.second) < std::tie(rhs.first, rhs.second);
}

bool worse_pair_on_top(const PairState& lhs, const PairState& rhs) { return better_pair(lhs, rhs); }

// Keeps the state_limit best pair-prefixes of a layer, streaming: the L * A^2 children a layer can
// generate are never all materialised.
class PairSelector {
public:
    explicit PairSelector(std::size_t limit) : limit_(limit), heap_(worse_pair_on_top) {}

    void offer(const PairState& pair) {
        ++offered_;
        if (heap_.size() < limit_) {
            heap_.push(pair);
        } else if (better_pair(pair, heap_.top())) {
            heap_.pop();
            heap_.push(pair);
        }
    }

    std::vector<PairState> take(std::uint64_t& pruned) {
        pruned += offered_ - heap_.size();
        std::vector<PairState> out;
        out.reserve(heap_.size());
        while (!heap_.empty()) {
            out.push_back(heap_.top());
            heap_.pop();
        }
        std::sort(out.begin(), out.end(), better_pair);
        return out;
    }

    std::uint64_t offered() const { return offered_; }
    std::size_t live() const { return heap_.size(); }

private:
    std::size_t limit_;
    std::uint64_t offered_ = 0;
    std::priority_queue<PairState, std::vector<PairState>, bool (*)(const PairState&, const PairState&)>
        heap_;
};

std::uint64_t checked_add(std::uint64_t total, std::uint64_t amount) {
    if (amount > std::numeric_limits<std::uint64_t>::max() - total) {
        throw std::length_error("grouped mosaic beam allocation bound overflows uint64");
    }
    return total + amount;
}

std::uint64_t checked_mul(std::uint64_t lhs, std::uint64_t rhs) {
    if (lhs != 0 && rhs > std::numeric_limits<std::uint64_t>::max() / lhs) {
        throw std::length_error("grouped mosaic beam allocation bound overflows uint64");
    }
    return lhs * rhs;
}

}   // namespace

// The work this proposer will do, from the exact loops below: block 0 builds one vector per allele over
// the templates (a single T scan) and forms its pairs; every later block scans T templates once per
// retained haplotype prefix (at most 2L of them) and forms A_b^2 children per retained pair.
std::uint64_t grouped_extension_upper_bound(std::size_t templates, std::size_t state_limit,
                                            const std::vector<std::size_t>& allele_counts) {
    if (templates == 0 || allele_counts.empty()) return 0;
    const std::uint64_t t = static_cast<std::uint64_t>(templates);
    const std::uint64_t limit = static_cast<std::uint64_t>(state_limit);
    const std::uint64_t first_alleles = static_cast<std::uint64_t>(allele_counts[0]);
    std::uint64_t total = checked_add(t, checked_mul(first_alleles, first_alleles));
    for (std::size_t b = 1; b < allele_counts.size(); ++b) {
        const std::uint64_t a = static_cast<std::uint64_t>(allele_counts[b]);
        total = checked_add(total, checked_mul(checked_mul(2u, limit), t));
        total = checked_add(total, checked_mul(limit, checked_mul(a, a)));
    }
    return total;
}

GroupedAllocationAccounts grouped_allocation_bounds(std::size_t templates, std::size_t blocks,
                                                    std::size_t state_limit,
                                                    const std::vector<std::size_t>& allele_counts) {
    GroupedAllocationAccounts bounds;
    if (templates == 0 || blocks == 0 || allele_counts.empty()) return bounds;
    const std::uint64_t t = static_cast<std::uint64_t>(templates);
    const std::uint64_t limit = static_cast<std::uint64_t>(state_limit);
    // A retained pair names at most two haplotype prefixes, so 2L bounds a layer. The LOCUS bounds it
    // too, and for a small one far more tightly: a prefix at block b is an allele path, so there are at
    // most A_0 * ... * A_b of them however large L is.
    //
    // It is the allele product and NOT the template count, which an earlier version of this bound got
    // wrong. The tempting argument -- a prefix exists only if some template spells it, so a layer holds
    // at most T prefixes -- is false precisely because this is a MOSAIC proposer: whenever r > 0 the
    // rank-1 jump term gives a child positive mass even when no single template spells its whole allele
    // path. That is the point of the search. The measured-allocation gate caught it at 14 ancestry
    // records against a predicted 12, which is exactly 2 + 4 + 8 for that fixture's allele products.
    const std::uint64_t pair_cap = checked_mul(2u, limit);
    std::uint64_t a_max = 0;
    for (const std::size_t a : allele_counts) a_max = std::max(a_max, static_cast<std::uint64_t>(a));
    // The coarsest a child layer can be: every retained parent extended by every allele.
    const std::uint64_t child_cap = checked_mul(pair_cap, a_max);

    // Retained (post-compaction) and pre-compaction child counts are DIFFERENT quantities and a layer's
    // children are the larger of the two. Bounding children by the retained count -- as this did -- is
    // wrong whenever the parent layer is at its width and the next block has several reachable alleles:
    // the child layer is built in full and only then selected and compacted, so it is
    //     children_b <= min(parents_{b-1} * A_b, A_0 * ... * A_b)
    // and child_index is literally a [parents][A_b] table.
    std::uint64_t max_prefixes = 0;
    std::uint64_t max_pairs = 0;
    std::uint64_t max_children = 0;
    std::uint64_t max_child_index = 0;
    std::uint64_t ancestry_total = 0;
    std::uint64_t product = 1;    // saturating at child_cap, which is above every bound drawn from it
    std::uint64_t parents = 0;
    for (std::size_t b = 0; b < allele_counts.size(); ++b) {
        const std::uint64_t a = static_cast<std::uint64_t>(allele_counts[b]);
        if (product < child_cap) {
            product = (a != 0 && product > child_cap / a) ? child_cap : product * a;
        }
        const std::uint64_t index_here = b == 0 ? 0 : std::min(child_cap, checked_mul(parents, a));
        const std::uint64_t children_here =
            b == 0 ? std::min(pair_cap, product) : std::min(index_here, product);
        max_child_index = std::max(max_child_index, index_here);
        max_children = std::max(max_children, children_here);
        // Retained after compaction: the width, or the allele paths at this depth, whichever binds.
        const std::uint64_t here = std::min(pair_cap, product);
        max_prefixes = std::max(max_prefixes, here);
        ancestry_total = checked_add(ancestry_total, here);
        const std::uint64_t pairs_here = here > 0 ? checked_mul(here, here + 1u) / 2u : 0;
        max_pairs = std::max(max_pairs, std::min(limit, pairs_here));
        parents = here;
    }
    const std::uint64_t prefixes_per_layer = max_prefixes;
    const std::uint64_t pairs_per_layer = max_pairs;
    // Parent and child layers are the only ones live; a prefix's compact vector is at most T cells.
    const std::uint64_t per_layer_cells = checked_mul(prefixes_per_layer, t);
    bounds.rolling_vector_cells = checked_mul(2u, per_layer_cells);
    bounds.rolling_vector_bytes = checked_mul(bounds.rolling_vector_cells, kCellBytes);
    // Ancestry persists: one compacted layer per block, each bounded the same structural way.
    bounds.ancestry_records = ancestry_total;
    bounds.ancestry_bytes = checked_mul(bounds.ancestry_records, kAncestryBytes);
    // The retained pairs of the previous layer, and the selector's heap for this one.
    bounds.pair_records = checked_mul(2u, pairs_per_layer);
    bounds.pair_bytes = checked_mul(bounds.pair_records, sizeof(PairState));

    // Everything else the proposer owns. The three dominant terms all scale as 2L * A_max, because a
    // layer's children are formed before they are selected and compacted: the child index, the remap
    // built over that child layer, and the child ancestry itself. A locus with 126 alleles in a block
    // makes these millions of bytes, not a rounding error.
    const std::uint64_t children = max_children;
    bounds.child_records = children;
    std::uint64_t index = checked_mul(max_child_index, sizeof(std::uint32_t));   // child_index
    index = checked_add(index, checked_mul(children, sizeof(std::uint32_t)));    // remap
    index = checked_add(index, checked_mul(children, kAncestryBytes));           // pre-compaction ancestry
    // The allele-class index: each template appears in exactly one class per block.
    index = checked_add(index, checked_mul(checked_mul(static_cast<std::uint64_t>(blocks), t),
                                           sizeof(std::uint32_t)));
    // Dense scratch, and one parent's scatter vectors, both O(T) cells of the same three arrays.
    index = checked_add(index, checked_mul(checked_mul(2u, t), kCellBytes));
    // The spelled candidate paths: two allele vectors of one entry per block, per retained pair.
    index = checked_add(index, checked_mul(checked_mul(checked_mul(pairs_per_layer, 2u),
                                                       static_cast<std::uint64_t>(blocks)),
                                           sizeof(std::size_t)));
    bounds.index_bytes = index;

    bounds.total_bytes = checked_add(
        checked_add(bounds.rolling_vector_bytes, bounds.ancestry_bytes),
        checked_add(bounds.pair_bytes, bounds.index_bytes));
    return bounds;
}

MosaicBeamResult propose_mosaic_diplotypes_grouped(const PanelAlleleMatrix& panel,
                                                   const PreparedProposalEmissions& emissions,
                                                   const std::vector<double>& switch_probability,
                                                   const MosaicBeamParameters& parameters) {
    // The SAME validation the state-based proposer runs, not a second, weaker set.
    validate_proposal_inputs(panel, emissions, switch_probability, parameters);
    const std::size_t templates = panel.template_count();
    const std::size_t blocks = panel.block_count();
    const double log_templates = std::log(static_cast<double>(templates));

    std::vector<std::size_t> allele_counts(blocks, 0);
    for (std::size_t b = 0; b < blocks; ++b) allele_counts[b] = emissions.blocks[b].allele_count;

    MosaicBeamResult result;
    result.state_limit_used = parameters.state_limit;

    // Transactional refusal: every budget is computed and checked BEFORE anything is allocated, so a
    // refused search allocates nothing and leaves no partial state. Each account is gated separately.
    result.extension_upper_bound =
        grouped_extension_upper_bound(templates, parameters.state_limit, allele_counts);
    const GroupedAllocationAccounts bounds =
        grouped_allocation_bounds(templates, blocks, parameters.state_limit, allele_counts);
    if (result.extension_upper_bound > parameters.max_extension_evaluations) {
        throw std::length_error("grouped mosaic beam work exceeds the declared pre-allocation budget");
    }
    if (bounds.rolling_vector_cells > parameters.max_rolling_vector_cells) {
        throw std::length_error("grouped mosaic beam rolling vectors exceed the declared budget");
    }
    if (bounds.ancestry_records > parameters.max_ancestry_records) {
        throw std::length_error("grouped mosaic beam ancestry exceeds the declared budget");
    }
    if (bounds.pair_records > parameters.max_pair_records) {
        throw std::length_error("grouped mosaic beam pair records exceed the declared budget");
    }
    // The fixed production ceiling, over the whole predicted peak, checked BEFORE construction: no
    // layer, no index and no class table has been built at this point, so a refusal allocates nothing
    // and can leave no partial output behind it.
    if (bounds.total_bytes > parameters.max_proposal_payload_bytes) {
        throw std::length_error("grouped mosaic beam predicted peak " +
                                std::to_string(bounds.total_bytes) +
                                " bytes exceeds the declared proposal payload ceiling of " +
                                std::to_string(parameters.max_proposal_payload_bytes) + " bytes");
    }
    result.allocation_predicted = bounds;

    // Which templates spell each allele, per block: the scatter target of the single template scan.
    std::vector<std::vector<std::vector<std::uint32_t>>> class_templates(blocks);
    for (std::size_t b = 0; b < blocks; ++b) {
        class_templates[b].assign(allele_counts[b], {});
        for (std::size_t h = 0; h < templates; ++h) {
            const std::size_t allele = panel.template_alleles[h][b];
            if (allele >= allele_counts[b]) {
                throw std::invalid_argument("grouped mosaic beam: template allele outside the emission table");
            }
            class_templates[b][allele].push_back(static_cast<std::uint32_t>(h));
        }
    }

    const bool traced = parameters.trace_target.has_value();
    const MosaicDiplotype empty_target;
    const MosaicDiplotype& target = traced ? *parameters.trace_target : empty_target;

    std::uint64_t peak_cells = 0;
    std::uint64_t peak_pairs = 0;
    std::uint64_t peak_index_bytes = 0;
    std::uint64_t peak_children = 0;

    // Compacts a freshly built layer to the prefixes retained pairs reference, remapping the pairs.
    // Indices move only WITHIN this layer; a layer's parent pointers name the previous layer, which was
    // compacted before this one was built, so they stay valid. `keep` is ascending, so the remap is
    // monotonic and each pair's first <= second survives it.
    const auto compact = [](const AncestryLayer& ancestry_in, const PrefixVectors& vectors_in,
                            std::vector<PairState>& pairs, AncestryLayer& ancestry_out,
                            PrefixVectors& vectors_out) {
        std::vector<std::uint32_t> keep;
        keep.reserve(pairs.size() * 2);
        for (const PairState& pair : pairs) {
            keep.push_back(pair.first);
            keep.push_back(pair.second);
        }
        std::sort(keep.begin(), keep.end());
        keep.erase(std::unique(keep.begin(), keep.end()), keep.end());

        std::vector<std::uint32_t> remap(ancestry_in.size(), kNoParent);
        ancestry_out = AncestryLayer{};
        vectors_out.reset();
        for (std::size_t n = 0; n < keep.size(); ++n) {
            const std::uint32_t old = keep[n];
            remap[old] = static_cast<std::uint32_t>(n);
            ancestry_out.parent.push_back(ancestry_in.parent[old]);
            ancestry_out.allele.push_back(ancestry_in.allele[old]);
            ancestry_out.matches_first.push_back(ancestry_in.matches_first[old]);
            ancestry_out.matches_second.push_back(ancestry_in.matches_second[old]);
            for (std::uint32_t e = vectors_in.offset[old]; e < vectors_in.offset[old + 1]; ++e) {
                vectors_out.sum.push_back(vectors_in.sum[e]);
                vectors_out.best.push_back(vectors_in.best[e]);
                vectors_out.count.push_back(vectors_in.count[e]);
            }
            vectors_out.offset.push_back(static_cast<std::uint32_t>(vectors_out.sum.size()));
            vectors_out.log_total.push_back(vectors_in.log_total[old]);
            vectors_out.best_total.push_back(vectors_in.best_total[old]);
            vectors_out.count_total.push_back(vectors_in.count_total[old]);
            vectors_out.count_saturated.push_back(vectors_in.count_saturated[old]);
            vectors_out.allele.push_back(vectors_in.allele[old]);
        }
        for (PairState& pair : pairs) {
            pair.first = remap[pair.first];
            pair.second = remap[pair.second];
        }
    };

    // ---- block 0: one haplotype prefix per allele, every template weighted 1/T ---------------------
    std::vector<AncestryLayer> ancestry;
    ancestry.reserve(blocks);
    PrefixVectors parent_vectors;
    PrefixVectors child_vectors;
    PrefixVectors rolled_vectors;
    AncestryLayer built_ancestry;
    PrefixVectors built_vectors;
    built_vectors.reset();
    {
        for (std::size_t a = 0; a < allele_counts[0]; ++a) {
            const auto& members = class_templates[0][a];
            built_ancestry.parent.push_back(kNoParent);
            built_ancestry.allele.push_back(static_cast<std::uint32_t>(a));
            built_ancestry.matches_first.push_back(traced && target.first.alleles[0] == a ? 1 : 0);
            built_ancestry.matches_second.push_back(traced && target.second.alleles[0] == a ? 1 : 0);
            double total = kNegInf;
            double best = kNegInf;
            std::uint64_t count = 0;
            bool saturated = false;
            for (std::size_t m = 0; m < members.size(); ++m) {
                built_vectors.sum.push_back(-log_templates);
                built_vectors.best.push_back(-log_templates);
                built_vectors.count.push_back(1);
                total = log_add(total, -log_templates);
                best = std::max(best, -log_templates);
                add_sat(count, saturated, 1, false);
            }
            built_vectors.offset.push_back(static_cast<std::uint32_t>(built_vectors.sum.size()));
            built_vectors.log_total.push_back(total);
            built_vectors.best_total.push_back(best);
            built_vectors.count_total.push_back(count);
            built_vectors.count_saturated.push_back(saturated ? 1 : 0);
            built_vectors.allele.push_back(static_cast<std::uint32_t>(a));
        }
        result.extensions_evaluated += templates;
        peak_cells = std::max(peak_cells, built_vectors.cells());
    }

    // Pairs at block 0: every unordered pair of block-0 alleles.
    std::vector<PairState> pairs;
    {
        const BlockProposalEmissions& emission = emissions.blocks[0];
        PairSelector selector(parameters.state_limit);
        MosaicBeamTraceLayer record;
        record.block = 0;
        for (std::size_t i = 0; i < built_vectors.size(); ++i) {
            for (std::size_t j = i; j < built_vectors.size(); ++j) {
                ++result.extensions_evaluated;
                if (built_vectors.log_total[i] == kNegInf || built_vectors.log_total[j] == kNegInf) continue;
                PairState pair;
                pair.first = static_cast<std::uint32_t>(i);
                pair.second = static_cast<std::uint32_t>(j);
                pair.emission = emission.log_score[built_vectors.allele[i] * emission.allele_count +
                                                   built_vectors.allele[j]];
                const bool heterozygous = i != j;
                // The heterozygous factor of two, in the SELECTION score from the first layer on: an
                // unordered pair of distinct prefixes stands for two ordered histories.
                pair.mass = pair.emission + built_vectors.log_total[i] + built_vectors.log_total[j] +
                            (heterozygous ? std::log(2.0) : 0.0);
                pair.best = pair.emission + built_vectors.best_total[i] + built_vectors.best_total[j];
                pair.best_transition = built_vectors.best_total[i] + built_vectors.best_total[j];
                bool saturated =
                    built_vectors.count_saturated[i] != 0 || built_vectors.count_saturated[j] != 0;
                pair.histories =
                    mul_sat(built_vectors.count_total[i], built_vectors.count_total[j], saturated);
                if (heterozygous) pair.histories = mul_sat(pair.histories, 2, saturated);
                pair.saturated = saturated;
                selector.offer(pair);
                if (traced && pair_is_target(built_ancestry, pair.first, pair.second)) {
                    record.target_prefixes_present = true;
                    record.target_states_offered = 1;
                    record.target_logsum_log_mass = pair.mass;
                    record.target_best_log_mass = pair.best;
                }
            }
        }
        peak_pairs = std::max<std::uint64_t>(peak_pairs, selector.live());
        pairs = selector.take(result.states_pruned);
        if (traced) {
            record.states_offered = selector.offered();
            record.states_retained = pairs.size();
            record.prefix_pairs_offered = static_cast<std::size_t>(selector.offered());
            record.prefix_pairs_retained = pairs.size();
            record.max_states_per_prefix = pairs.empty() ? 0 : 1;
            record.states_pruned = selector.offered() - pairs.size();
            if (record.target_states_offered > 0) {
                std::size_t rank = 1;
                for (const PairState& kept : pairs) {
                    if (kept.mass > record.target_logsum_log_mass) ++rank;
                    if (pair_is_target(built_ancestry, kept.first, kept.second)) {
                        record.target_states_retained = 1;
                        record.target_retained_logsum_log_mass = kept.mass;
                    }
                }
                record.target_group_rank_by_offered_mass = rank;
            }
            result.trace.push_back(record);
        }
        AncestryLayer compacted;
        compact(built_ancestry, built_vectors, pairs, compacted, parent_vectors);
        ancestry.push_back(std::move(compacted));
    }

    // ---- blocks 1..B-1 -----------------------------------------------------------------------------
    std::vector<double> scratch_sum(templates, kNegInf);
    std::vector<double> scratch_best(templates, kNegInf);
    std::vector<std::uint64_t> scratch_count(templates, 0);

    for (std::size_t block = 1; block < blocks && !pairs.empty(); ++block) {
        const double r = switch_probability[block - 1];
        const double log_same = std::log((1.0 - r) + r / static_cast<double>(templates));
        const double log_other = r > 0.0 ? std::log(r) - log_templates : kNegInf;
        const double log_stay = r < 1.0 ? std::log1p(-r) : kNegInf;
        const BlockProposalEmissions& emission = emissions.blocks[block];
        const AncestryLayer& parent_ancestry = ancestry.back();

        AncestryLayer child_ancestry;
        child_vectors.reset();
        // child_index[parent][allele] -> index in the child layer, or kNoParent when that allele has no
        // reachable template. Every parent is referenced: the layer was compacted after its selection.
        std::vector<std::vector<std::uint32_t>> child_index(
            parent_vectors.size(), std::vector<std::uint32_t>(allele_counts[block], kNoParent));

        for (std::size_t p = 0; p < parent_vectors.size(); ++p) {
            const std::uint32_t from = parent_vectors.offset[p];
            const std::uint32_t to = parent_vectors.offset[p + 1];
            const auto& parent_members = class_templates[block - 1][parent_vectors.allele[p]];

            // Dense scratch for this parent only, filled and cleared over its own class: O(class size).
            double best_first = kNegInf;
            double best_second = kNegInf;
            for (std::uint32_t e = from; e < to; ++e) {
                const std::uint32_t i = parent_members[e - from];
                scratch_sum[i] = parent_vectors.sum[e];
                scratch_best[i] = parent_vectors.best[e];
                scratch_count[i] = parent_vectors.count[e];
                if (parent_vectors.best[e] > best_first) {
                    best_second = best_first;
                    best_first = parent_vectors.best[e];
                } else if (parent_vectors.best[e] > best_second) {
                    best_second = parent_vectors.best[e];
                }
            }

            // One scan over every template, scattered into the child of its own allele class.
            std::vector<std::vector<double>> child_sum(allele_counts[block]);
            std::vector<std::vector<double>> child_best(allele_counts[block]);
            std::vector<std::vector<std::uint64_t>> child_count(allele_counts[block]);
            for (std::size_t i = 0; i < templates; ++i) {
                ++result.extensions_evaluated;
                // Low-rank Li-Stephens: staying keeps (1-r) of this template's own mass, and the rank-1
                // jump adds r/T of the prefix's total to every template alike.
                const double sum_here =
                    log_add(log_stay + scratch_sum[i], log_other + parent_vectors.log_total[p]);
                // The max-product cannot use that decomposition: it needs the best OTHER template,
                // which is the top of the parent's two best values, excluding i itself.
                const double other_best = (scratch_best[i] == best_first && scratch_best[i] != kNegInf)
                                              ? best_second
                                              : best_first;
                const double best_here =
                    std::max(scratch_best[i] == kNegInf ? kNegInf : scratch_best[i] + log_same,
                             other_best == kNegInf ? kNegInf : other_best + log_other);
                // Every retained path can reach any template when a switch is possible; with r = 0 only
                // the path already on this template continues.
                const std::uint64_t count_here =
                    log_other == kNegInf ? scratch_count[i] : parent_vectors.count_total[p];
                // Every template of the class is pushed, unreachable ones included. The compact vector
                // must stay positionally aligned with class_templates, because that is how the NEXT
                // block scatters back from a compact entry to its template id; skipping a zero here
                // would shift every later entry of the class and silently corrupt the next layer.
                const std::size_t allele = panel.template_alleles[i][block];
                child_sum[allele].push_back(sum_here);
                child_best[allele].push_back(best_here);
                child_count[allele].push_back(count_here);
            }

            for (std::uint32_t e = from; e < to; ++e) {
                const std::uint32_t i = parent_members[e - from];
                scratch_sum[i] = kNegInf;
                scratch_best[i] = kNegInf;
                scratch_count[i] = 0;
            }

            for (std::size_t a = 0; a < allele_counts[block]; ++a) {
                if (child_sum[a].empty()) continue;
                double total = kNegInf;
                for (const double value : child_sum[a]) total = log_add(total, value);
                // A child none of whose templates is reachable carries no mass and is not a prefix the
                // search can continue. It is simply never created, and that allele has no child here.
                if (total == kNegInf) continue;
                child_index[p][a] = static_cast<std::uint32_t>(child_ancestry.size());
                child_ancestry.parent.push_back(static_cast<std::uint32_t>(p));
                child_ancestry.allele.push_back(static_cast<std::uint32_t>(a));
                child_ancestry.matches_first.push_back(
                    traced && parent_ancestry.matches_first[p] != 0 && target.first.alleles[block] == a ? 1 : 0);
                child_ancestry.matches_second.push_back(
                    traced && parent_ancestry.matches_second[p] != 0 && target.second.alleles[block] == a ? 1 : 0);
                double best = kNegInf;
                std::uint64_t count = 0;
                bool saturated = parent_vectors.count_saturated[p] != 0;
                for (std::size_t k = 0; k < child_sum[a].size(); ++k) {
                    child_vectors.sum.push_back(child_sum[a][k]);
                    child_vectors.best.push_back(child_best[a][k]);
                    child_vectors.count.push_back(child_count[a][k]);
                    best = std::max(best, child_best[a][k]);
                    add_sat(count, saturated, child_count[a][k], false);
                }
                child_vectors.offset.push_back(static_cast<std::uint32_t>(child_vectors.sum.size()));
                child_vectors.log_total.push_back(total);
                child_vectors.best_total.push_back(best);
                child_vectors.count_total.push_back(count);
                child_vectors.count_saturated.push_back(saturated ? 1 : 0);
                child_vectors.allele.push_back(static_cast<std::uint32_t>(a));
            }
        }

        // Both layers are live at exactly this point, which is what the rolling account bounds.
        peak_cells = std::max(peak_cells, parent_vectors.cells() + child_vectors.cells());
        {
            const std::uint64_t kids = static_cast<std::uint64_t>(child_ancestry.size());
            peak_children = std::max(peak_children, kids);
            const std::uint64_t here =
                static_cast<std::uint64_t>(parent_vectors.size()) * allele_counts[block] *
                    sizeof(std::uint32_t) +                     // child_index
                kids * sizeof(std::uint32_t) +                  // remap
                kids * kAncestryBytes +                         // child ancestry before compaction
                static_cast<std::uint64_t>(blocks) * templates * sizeof(std::uint32_t) +
                2u * static_cast<std::uint64_t>(templates) * kCellBytes;
            peak_index_bytes = std::max(peak_index_bytes, here);
        }
        if (child_ancestry.size() >= kNoParent) {
            throw std::length_error("grouped mosaic beam prefix index exceeds uint32");
        }

        // Pair extension: two scalar lookups per ordered allele pair, no cell materialised.
        PairSelector selector(parameters.state_limit);
        MosaicBeamTraceLayer record;
        record.block = block;
        for (const PairState& pair : pairs) {
            const bool symmetric = pair.first == pair.second;
            for (std::size_t a = 0; a < allele_counts[block]; ++a) {
                const std::uint32_t c1 = child_index[pair.first][a];
                if (c1 == kNoParent) continue;
                // A symmetric parent's (a, b) and (b, a) are one unordered child, so it is formed once.
                for (std::size_t b = symmetric ? a : 0; b < allele_counts[block]; ++b) {
                    ++result.extensions_evaluated;
                    const std::uint32_t c2 = child_index[pair.second][b];
                    if (c2 == kNoParent) continue;
                    PairState child;
                    child.first = std::min(c1, c2);
                    child.second = std::max(c1, c2);
                    child.emission = pair.emission + emission.log_score[a * emission.allele_count + b];
                    const bool heterozygous = child.first != child.second;
                    child.mass = child.emission + child_vectors.log_total[child.first] +
                                 child_vectors.log_total[child.second] +
                                 (heterozygous ? std::log(2.0) : 0.0);
                    child.best = child.emission + child_vectors.best_total[child.first] +
                                 child_vectors.best_total[child.second];
                    child.best_transition =
                        child_vectors.best_total[child.first] + child_vectors.best_total[child.second];
                    bool saturated = child_vectors.count_saturated[child.first] != 0 ||
                                     child_vectors.count_saturated[child.second] != 0;
                    child.histories = mul_sat(child_vectors.count_total[child.first],
                                              child_vectors.count_total[child.second], saturated);
                    if (heterozygous) child.histories = mul_sat(child.histories, 2, saturated);
                    child.saturated = saturated;
                    if (child.mass == kNegInf) continue;
                    selector.offer(child);
                    if (traced && pair_is_target(child_ancestry, child.first, child.second)) {
                        record.target_prefixes_present = true;
                        record.target_states_offered = 1;
                        record.target_logsum_log_mass = child.mass;
                        record.target_best_log_mass = child.best;
                    }
                }
            }
        }

        peak_pairs = std::max<std::uint64_t>(peak_pairs, pairs.size() + selector.live());
        pairs = selector.take(result.states_pruned);
        if (traced) {
            record.states_offered = selector.offered();
            record.states_retained = pairs.size();
            record.prefix_pairs_offered = static_cast<std::size_t>(selector.offered());
            record.prefix_pairs_retained = pairs.size();
            record.max_states_per_prefix = pairs.empty() ? 0 : 1;
            record.states_pruned = selector.offered() - pairs.size();
            if (record.target_states_offered > 0) {
                std::size_t rank = 1;
                for (const PairState& kept : pairs) {
                    if (kept.mass > record.target_logsum_log_mass) ++rank;
                    if (pair_is_target(child_ancestry, kept.first, kept.second)) {
                        record.target_states_retained = 1;
                        record.target_retained_logsum_log_mass = kept.mass;
                    }
                }
                record.target_group_rank_by_offered_mass = rank;
            }
            result.trace.push_back(record);
        }

        // Compact the child layer to what the retained pairs reference, then ROLL: the child becomes
        // the parent and the previous parent's numerical arrays are released here, not at the end.
        AncestryLayer compacted;
        compact(child_ancestry, child_vectors, pairs, compacted, rolled_vectors);
        ancestry.push_back(std::move(compacted));
        parent_vectors = std::move(rolled_vectors);
        rolled_vectors = PrefixVectors{};
        child_vectors.reset();
    }

    // ---- candidates: every retained pair-prefix IS one diplotype --------------------------------
    const auto spell = [&](std::uint32_t index) {
        std::vector<std::size_t> alleles(blocks, 0);
        for (std::size_t b = ancestry.size(); b > 0; --b) {
            const AncestryLayer& layer = ancestry[b - 1];
            alleles[b - 1] = layer.allele[index];
            index = layer.parent[index];
        }
        return alleles;
    };

    std::vector<MosaicCandidate> candidates;
    candidates.reserve(pairs.size());
    for (const PairState& pair : pairs) {
        MosaicCandidate candidate;
        candidate.diplotype.first.alleles = spell(pair.first);
        candidate.diplotype.second.alleles = spell(pair.second);
        if (candidate.diplotype.second.alleles < candidate.diplotype.first.alleles) {
            std::swap(candidate.diplotype.first, candidate.diplotype.second);
        }
        candidate.retained_logsum_proposal_score = pair.mass;
        candidate.best_proposal_score = pair.best;
        candidate.proposal_emission_component = pair.emission;
        candidate.proposal_transition_component = pair.best_transition;
        candidate.equivalent_retained_histories = pair.histories;
        candidate.history_count_saturated = pair.saturated;
        candidate.every_allele_scored = emissions.every_allele_scored();
        candidates.push_back(std::move(candidate));
    }

    std::sort(candidates.begin(), candidates.end(), [](const MosaicCandidate& lhs, const MosaicCandidate& rhs) {
        if (lhs.retained_logsum_proposal_score != rhs.retained_logsum_proposal_score) {
            return lhs.retained_logsum_proposal_score > rhs.retained_logsum_proposal_score;
        }
        if (lhs.best_proposal_score != rhs.best_proposal_score) {
            return lhs.best_proposal_score > rhs.best_proposal_score;
        }
        if (lhs.diplotype.first.alleles != rhs.diplotype.first.alleles) {
            return lhs.diplotype.first.alleles < rhs.diplotype.first.alleles;
        }
        return lhs.diplotype.second.alleles < rhs.diplotype.second.alleles;
    });

    std::uint64_t ancestry_records = 0;
    for (const AncestryLayer& layer : ancestry) ancestry_records += layer.size();
    result.allocation_observed.rolling_vector_cells = peak_cells;
    result.allocation_observed.rolling_vector_bytes = peak_cells * kCellBytes;
    result.allocation_observed.ancestry_records = ancestry_records;
    result.allocation_observed.ancestry_bytes = ancestry_records * kAncestryBytes;
    result.allocation_observed.pair_records = peak_pairs;
    result.allocation_observed.pair_bytes = peak_pairs * sizeof(PairState);
    const std::uint64_t candidate_path_bytes =
        static_cast<std::uint64_t>(pairs.size()) * 2u * blocks * sizeof(std::size_t);
    result.allocation_observed.child_records = peak_children;
    result.allocation_observed.index_bytes = peak_index_bytes + candidate_path_bytes;
    result.allocation_observed.total_bytes = result.allocation_observed.rolling_vector_bytes +
                                             result.allocation_observed.ancestry_bytes +
                                             result.allocation_observed.pair_bytes +
                                             result.allocation_observed.index_bytes;

    result.final_states = pairs.size();
    result.unique_candidates_found = candidates.size();
    result.search_exhaustive = result.states_pruned == 0;
    result.unique_target_reached = candidates.size() >= parameters.unique_candidate_target;
    result.candidate_set_conditional =
        !result.search_exhaustive || candidates.size() > parameters.unique_candidate_target;
    if (candidates.size() > parameters.unique_candidate_target) {
        candidates.resize(parameters.unique_candidate_target);
    }
    for (std::size_t rank = 0; rank < candidates.size(); ++rank) {
        candidates[rank].candidate_id = rank;
        candidates[rank].unique_rank = rank + 1;
    }
    result.candidates = std::move(candidates);
    return result;
}

}   // namespace panvar
