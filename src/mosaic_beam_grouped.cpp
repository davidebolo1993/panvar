// The grouped-prefix proposer. Implemented ALONGSIDE propose_mosaic_diplotypes, which is untouched and
// is the oracle this is tested against; nothing here is wired into calling.
//
// The unit of retention is the allele-prefix PAIR, and each pair carries the whole distribution over
// current template pairs rather than one attribution of it. That is exact where prefix-first selection
// approximates: states sharing a prefix pair are different latent template histories spelling the same
// allele prefix, and the current template is what the next block's transition reads.
//
// It is affordable because the joint factorises. Within a prefix pair every history spells the same
// alleles, so the diploid emission is a single scalar for the whole group -- and the emission is the
// only coupling between the two haplotypes. Hence
//
//     W_p(i, j) = E(p1, p2) * u_p1(i) * u_p2(j)
//
// with u_p(i) the mass of template paths spelling prefix p and ending on template i. A group therefore
// needs two vectors of length T, each shared by every pair using that haplotype prefix, plus scalars.
//
// Two structural facts keep the cost linear in T rather than A*T:
//   * Li-Stephens is low rank. Staying costs (1-r) + r/T and moving costs r/T, so
//         u'(i') = (1 - r) * u_p(i')  +  (r / T) * sum_i u_p(i)
//     and the rank-1 term is one scalar -- the prefix's running total -- added to every entry.
//   * Each template spells exactly ONE allele at a block. So a single scan over the T templates
//     produces every child of a prefix at once, scattered into its allele class, and all A children
//     together occupy exactly T entries. One scan, not A scans.
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

// One layer's haplotype prefixes. Vectors are compact: prefix p owns entries
// [offset[p], offset[p+1]), one per template spelling p's last allele, in template order. All the
// children of one parent therefore occupy T entries between them, which is what makes the layer
// O(2 L T) rather than O(2 L A T).
struct HaploidLayer {
    std::vector<std::uint32_t> parent;   // index into the previous layer, kNoParent at block 0
    std::vector<std::uint32_t> allele;
    std::vector<char> matches_first;     // target's first path, tracked incrementally
    std::vector<char> matches_second;
    std::vector<std::uint32_t> offset;   // size() + 1 entries
    std::vector<double> sum;             // log u_p(i), sum-product
    std::vector<double> best;            // log max-product over single histories
    std::vector<std::uint64_t> count;    // template paths ending at i
    std::vector<double> log_total;       // log sum_i u_p(i)
    std::vector<double> best_total;      // max_i best(i)
    std::vector<std::uint64_t> count_total;
    std::vector<char> count_saturated;

    std::size_t size() const { return allele.size(); }
};

// The trace target is an UNORDERED diplotype, and pairs are canonicalised by prefix index, so the
// target's two prefixes can occupy either slot and the canonical order can flip between layers.
// Checking a single orientation reports the truth as lost at the first layer where it does.
bool pair_is_target(const HaploidLayer& layer, std::uint32_t first, std::uint32_t second) {
    return (layer.matches_first[first] != 0 && layer.matches_second[second] != 0) ||
           (layer.matches_second[first] != 0 && layer.matches_first[second] != 0);
}

// A pair-prefix: the unit this proposer retains and ranks.
struct PairState {
    std::uint32_t first = 0;    // index into the current HaploidLayer
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

private:
    std::size_t limit_;
    std::uint64_t offered_ = 0;
    std::priority_queue<PairState, std::vector<PairState>, bool (*)(const PairState&, const PairState&)>
        heap_;
};

void validate(const PanelAlleleMatrix& panel, const PreparedProposalEmissions& emissions,
              const std::vector<double>& switch_probability, const MosaicBeamParameters& parameters) {
    if (panel.template_count() == 0 || panel.block_count() == 0) {
        throw std::invalid_argument("grouped mosaic beam needs a non-empty panel");
    }
    for (const auto& row : panel.template_alleles) {
        if (row.size() != panel.block_count()) {
            throw std::invalid_argument("grouped mosaic beam: every template must span every block");
        }
    }
    if (emissions.blocks.size() != panel.block_count()) {
        throw std::invalid_argument("grouped mosaic beam: one emission table per block is required");
    }
    if (switch_probability.size() + 1 != panel.block_count()) {
        throw std::invalid_argument("grouped mosaic beam: one switch probability per block boundary");
    }
    for (const double r : switch_probability) {
        if (!(r >= 0.0 && r <= 1.0)) {
            throw std::invalid_argument("grouped mosaic beam: switch probability outside [0, 1]");
        }
    }
    if (parameters.state_limit == 0) {
        throw std::invalid_argument("grouped mosaic beam: a zero state limit retains nothing");
    }
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
    const auto checked_add = [](std::uint64_t total, std::uint64_t amount) {
        if (amount > std::numeric_limits<std::uint64_t>::max() - total) {
            throw std::length_error("grouped mosaic beam work bound overflows uint64");
        }
        return total + amount;
    };
    const auto checked_mul = [](std::uint64_t lhs, std::uint64_t rhs) {
        if (lhs != 0 && rhs > std::numeric_limits<std::uint64_t>::max() / lhs) {
            throw std::length_error("grouped mosaic beam work bound overflows uint64");
        }
        return lhs * rhs;
    };

    const std::uint64_t first_alleles = static_cast<std::uint64_t>(allele_counts[0]);
    std::uint64_t total = checked_add(t, checked_mul(first_alleles, first_alleles));
    for (std::size_t b = 1; b < allele_counts.size(); ++b) {
        const std::uint64_t a = static_cast<std::uint64_t>(allele_counts[b]);
        // One template scan per retained haplotype prefix, of which a layer holds at most two per pair.
        total = checked_add(total, checked_mul(checked_mul(2u, limit), t));
        // One child per ordered allele pair per retained pair.
        total = checked_add(total, checked_mul(limit, checked_mul(a, a)));
    }
    return total;
}

// Vector cells held at once: the parent layer and the child layer being built, each at most 2L prefixes
// whose compact vectors together span the T templates.
std::uint64_t grouped_prefix_vector_cells(std::size_t templates, std::size_t state_limit) {
    const std::uint64_t t = static_cast<std::uint64_t>(templates);
    const std::uint64_t limit = static_cast<std::uint64_t>(state_limit);
    if (limit != 0 && t > std::numeric_limits<std::uint64_t>::max() / (4u * limit)) {
        throw std::length_error("grouped mosaic beam vector bound overflows uint64");
    }
    return 4u * limit * t;
}

MosaicBeamResult propose_mosaic_diplotypes_grouped(const PanelAlleleMatrix& panel,
                                                   const PreparedProposalEmissions& emissions,
                                                   const std::vector<double>& switch_probability,
                                                   const MosaicBeamParameters& parameters) {
    validate(panel, emissions, switch_probability, parameters);
    const std::size_t templates = panel.template_count();
    const std::size_t blocks = panel.block_count();
    const double log_templates = std::log(static_cast<double>(templates));

    std::vector<std::size_t> allele_counts(blocks, 0);
    for (std::size_t b = 0; b < blocks; ++b) allele_counts[b] = emissions.blocks[b].allele_count;

    MosaicBeamResult result;
    result.state_limit_used = parameters.state_limit;

    // Transactional refusal: both budgets are computed and checked BEFORE anything is allocated, so a
    // refused search allocates nothing and leaves no partial state.
    result.extension_upper_bound =
        grouped_extension_upper_bound(templates, parameters.state_limit, allele_counts);
    const std::uint64_t vector_cells = grouped_prefix_vector_cells(templates, parameters.state_limit);
    if (result.extension_upper_bound > parameters.max_extension_evaluations) {
        throw std::length_error("grouped mosaic beam work exceeds the declared pre-allocation budget");
    }
    if (vector_cells > parameters.max_prefix_vector_cells) {
        throw std::length_error("grouped mosaic beam vectors exceed the declared pre-allocation budget");
    }

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

    // ---- block 0: one haplotype prefix per allele, every template weighted 1/T ---------------------
    std::vector<HaploidLayer> layers;
    layers.reserve(blocks);
    {
        HaploidLayer layer;
        layer.offset.push_back(0);
        for (std::size_t a = 0; a < allele_counts[0]; ++a) {
            const auto& members = class_templates[0][a];
            layer.parent.push_back(kNoParent);
            layer.allele.push_back(static_cast<std::uint32_t>(a));
            layer.matches_first.push_back(traced && target.first.alleles[0] == a ? 1 : 0);
            layer.matches_second.push_back(traced && target.second.alleles[0] == a ? 1 : 0);
            double total = kNegInf;
            double best = kNegInf;
            std::uint64_t count = 0;
            bool saturated = false;
            for (const std::uint32_t member : members) {
                (void)member;
                layer.sum.push_back(-log_templates);
                layer.best.push_back(-log_templates);
                layer.count.push_back(1);
                total = log_add(total, -log_templates);
                best = std::max(best, -log_templates);
                add_sat(count, saturated, 1, false);
            }
            layer.offset.push_back(static_cast<std::uint32_t>(layer.sum.size()));
            layer.log_total.push_back(total);
            layer.best_total.push_back(best);
            layer.count_total.push_back(count);
            layer.count_saturated.push_back(saturated ? 1 : 0);
        }
        result.extensions_evaluated += templates;
        layers.push_back(std::move(layer));
    }

    // Pairs at block 0: every unordered pair of block-0 alleles.
    std::vector<PairState> pairs;
    {
        const HaploidLayer& layer = layers.back();
        const BlockProposalEmissions& emission = emissions.blocks[0];
        PairSelector selector(parameters.state_limit);
        MosaicBeamTraceLayer record;
        record.block = 0;
        for (std::size_t i = 0; i < layer.size(); ++i) {
            for (std::size_t j = i; j < layer.size(); ++j) {
                ++result.extensions_evaluated;
                if (layer.log_total[i] == kNegInf || layer.log_total[j] == kNegInf) continue;
                PairState pair;
                pair.first = static_cast<std::uint32_t>(i);
                pair.second = static_cast<std::uint32_t>(j);
                pair.emission = emission.log_score[layer.allele[i] * emission.allele_count + layer.allele[j]];
                const bool heterozygous = i != j;
                // The heterozygous factor of two, in the SELECTION score from the first layer on: an
                // unordered pair of distinct prefixes stands for two ordered histories.
                pair.mass = pair.emission + layer.log_total[i] + layer.log_total[j] +
                            (heterozygous ? std::log(2.0) : 0.0);
                pair.best = pair.emission + layer.best_total[i] + layer.best_total[j];
                pair.best_transition = layer.best_total[i] + layer.best_total[j];
                bool saturated = layer.count_saturated[i] != 0 || layer.count_saturated[j] != 0;
                pair.histories = mul_sat(layer.count_total[i], layer.count_total[j], saturated);
                if (heterozygous) pair.histories = mul_sat(pair.histories, 2, saturated);
                pair.saturated = saturated;
                selector.offer(pair);
                if (traced && pair_is_target(layer, static_cast<std::uint32_t>(i),
                                             static_cast<std::uint32_t>(j))) {
                    record.target_prefixes_present = true;
                    record.target_states_offered = 1;
                    record.target_logsum_log_mass = pair.mass;
                    record.target_best_log_mass = pair.best;
                }
            }
        }
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
                    if (pair_is_target(layers.back(), kept.first, kept.second)) {
                        record.target_states_retained = 1;
                        record.target_retained_logsum_log_mass = kept.mass;
                    }
                }
                record.target_group_rank_by_offered_mass = rank;
            }
            result.trace.push_back(record);
        }
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
        const HaploidLayer& parent_layer = layers.back();

        // Which parent prefixes are still referenced, and where each one's children begin.
        std::vector<std::uint32_t> referenced;
        referenced.reserve(pairs.size() * 2);
        for (const PairState& pair : pairs) {
            referenced.push_back(pair.first);
            referenced.push_back(pair.second);
        }
        std::sort(referenced.begin(), referenced.end());
        referenced.erase(std::unique(referenced.begin(), referenced.end()), referenced.end());

        HaploidLayer child_layer;
        child_layer.offset.push_back(0);
        // child_index[parent slot][allele] -> index in child_layer, or kNoParent when that allele has
        // no template and therefore no child.
        std::vector<std::vector<std::uint32_t>> child_index(
            referenced.size(), std::vector<std::uint32_t>(allele_counts[block], kNoParent));
        std::vector<std::uint32_t> parent_slot(parent_layer.size(), kNoParent);
        for (std::size_t s = 0; s < referenced.size(); ++s) parent_slot[referenced[s]] = static_cast<std::uint32_t>(s);

        for (std::size_t s = 0; s < referenced.size(); ++s) {
            const std::uint32_t p = referenced[s];
            const std::uint32_t from = parent_layer.offset[p];
            const std::uint32_t to = parent_layer.offset[p + 1];
            const auto& parent_members = class_templates[block - 1][parent_layer.allele[p]];

            // Dense scratch for this parent only, filled and cleared over its own class: O(class size).
            double best_first = kNegInf;
            double best_second = kNegInf;
            for (std::uint32_t e = from; e < to; ++e) {
                const std::uint32_t i = parent_members[e - from];
                scratch_sum[i] = parent_layer.sum[e];
                scratch_best[i] = parent_layer.best[e];
                scratch_count[i] = parent_layer.count[e];
                if (parent_layer.best[e] > best_first) {
                    best_second = best_first;
                    best_first = parent_layer.best[e];
                } else if (parent_layer.best[e] > best_second) {
                    best_second = parent_layer.best[e];
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
                    log_add(log_stay + scratch_sum[i], log_other + parent_layer.log_total[p]);
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
                    log_other == kNegInf ? scratch_count[i] : parent_layer.count_total[p];
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
                child_index[s][a] = static_cast<std::uint32_t>(child_layer.size());
                child_layer.parent.push_back(p);
                child_layer.allele.push_back(static_cast<std::uint32_t>(a));
                child_layer.matches_first.push_back(
                    traced && parent_layer.matches_first[p] != 0 && target.first.alleles[block] == a ? 1 : 0);
                child_layer.matches_second.push_back(
                    traced && parent_layer.matches_second[p] != 0 && target.second.alleles[block] == a ? 1 : 0);
                double best = kNegInf;
                std::uint64_t count = 0;
                bool saturated = parent_layer.count_saturated[p] != 0;
                for (std::size_t k = 0; k < child_sum[a].size(); ++k) {
                    child_layer.sum.push_back(child_sum[a][k]);
                    child_layer.best.push_back(child_best[a][k]);
                    child_layer.count.push_back(child_count[a][k]);
                    best = std::max(best, child_best[a][k]);
                    add_sat(count, saturated, child_count[a][k], false);
                }
                child_layer.offset.push_back(static_cast<std::uint32_t>(child_layer.sum.size()));
                child_layer.log_total.push_back(total);
                child_layer.best_total.push_back(best);
                child_layer.count_total.push_back(count);
                child_layer.count_saturated.push_back(saturated ? 1 : 0);
            }
        }

        // Pair extension: two scalar lookups per ordered allele pair, no cell materialised.
        PairSelector selector(parameters.state_limit);
        MosaicBeamTraceLayer record;
        record.block = block;
        for (const PairState& pair : pairs) {
            const std::uint32_t slot1 = parent_slot[pair.first];
            const std::uint32_t slot2 = parent_slot[pair.second];
            const bool symmetric = pair.first == pair.second;
            for (std::size_t a = 0; a < allele_counts[block]; ++a) {
                const std::uint32_t c1 = child_index[slot1][a];
                if (c1 == kNoParent) continue;
                // A symmetric parent's (a, b) and (b, a) are one unordered child, so it is formed once.
                for (std::size_t b = symmetric ? a : 0; b < allele_counts[block]; ++b) {
                    ++result.extensions_evaluated;
                    const std::uint32_t c2 = child_index[slot2][b];
                    if (c2 == kNoParent) continue;
                    PairState child;
                    child.first = std::min(c1, c2);
                    child.second = std::max(c1, c2);
                    child.emission =
                        pair.emission + emission.log_score[a * emission.allele_count + b];
                    const bool heterozygous = child.first != child.second;
                    child.mass = child.emission + child_layer.log_total[child.first] +
                                 child_layer.log_total[child.second] +
                                 (heterozygous ? std::log(2.0) : 0.0);
                    child.best = child.emission + child_layer.best_total[child.first] +
                                 child_layer.best_total[child.second];
                    child.best_transition =
                        child_layer.best_total[child.first] + child_layer.best_total[child.second];
                    bool saturated = child_layer.count_saturated[child.first] != 0 ||
                                     child_layer.count_saturated[child.second] != 0;
                    child.histories = mul_sat(child_layer.count_total[child.first],
                                              child_layer.count_total[child.second], saturated);
                    if (heterozygous) child.histories = mul_sat(child.histories, 2, saturated);
                    child.saturated = saturated;
                    if (child.mass == kNegInf) continue;
                    selector.offer(child);
                    if (traced && pair_is_target(child_layer, child.first, child.second)) {
                        record.target_prefixes_present = true;
                        record.target_states_offered = 1;
                        record.target_logsum_log_mass = child.mass;
                        record.target_best_log_mass = child.best;
                    }
                }
            }
        }

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
                    if (pair_is_target(child_layer, kept.first, kept.second)) {
                        record.target_states_retained = 1;
                        record.target_retained_logsum_log_mass = kept.mass;
                    }
                }
                record.target_group_rank_by_offered_mass = rank;
            }
            result.trace.push_back(record);
        }
        layers.push_back(std::move(child_layer));
    }

    // ---- candidates: every retained pair-prefix IS one diplotype --------------------------------
    const auto spell = [&](std::uint32_t index) {
        std::vector<std::size_t> alleles(blocks, 0);
        for (std::size_t b = blocks; b > 0; --b) {
            const HaploidLayer& layer = layers[b - 1];
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
