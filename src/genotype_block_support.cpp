#include "panvar/genotype_block_support.hpp"

#include <set>
#include <stdexcept>

namespace panvar {

DataSupportResult block_data_support_v1(const std::vector<SupportGroup>& groups, std::size_t winner, std::size_t block_count) {
    if (winner >= groups.size()) throw std::invalid_argument("block support: winner outside the groups");
    for (const SupportGroup& g : groups) {
        if (g.rows.empty()) throw std::invalid_argument("block support: a group without members");
        for (const BlockGenotypeRow& row : g.rows) {
            if (row.size() != block_count) throw std::invalid_argument("block support: a genotype row with the wrong block count");
        }
    }
    const auto likelihood = [](const SupportGroup& g) { return g.fragment + g.exposure; };
    const auto above = [&](const SupportGroup& x, const SupportGroup& y) {   // x certifiably above y in L
        return likelihood(x) - x.bound > likelihood(y) + y.bound;
    };
    const SupportGroup& w = groups[winner];
    DataSupportResult out;
    out.blocks.resize(block_count);
    for (std::size_t b = 0; b < block_count; ++b) {
        BlockDataSupport& r = out.blocks[b];
        std::set<std::pair<std::size_t, std::size_t>> called;
        for (const BlockGenotypeRow& row : w.rows) called.insert(row[b]);
        if (called.size() != 1) {
            r.label = "AMBIGUOUS_PROJECTION";
            continue;
        }
        const auto g = *called.begin();
        long s = -1, d = -1;
        std::vector<std::size_t> mixed;
        for (std::size_t k = 0; k < groups.size(); ++k) {
            bool carries = false, lacks = false;
            for (const BlockGenotypeRow& row : groups[k].rows) (row[b] == g ? carries : lacks) = true;
            if (carries && !lacks) {
                if (s < 0 || likelihood(groups[k]) > likelihood(groups[static_cast<std::size_t>(s)])) s = static_cast<long>(k);
            } else if (lacks && !carries) {
                ++r.disagreeing;
                if (d < 0 || likelihood(groups[k]) > likelihood(groups[static_cast<std::size_t>(d)])) d = static_cast<long>(k);
            } else {
                mixed.push_back(k);
            }
        }
        r.mixed = mixed.size();
        if (d < 0) {
            r.label = "NO_EVALUATED_ALTERNATIVE";
            continue;
        }
        const SupportGroup& sg = groups[static_cast<std::size_t>(s)];   // s >= 0: w itself carries g
        const SupportGroup& dg = groups[static_cast<std::size_t>(d)];
        r.competitor = d;
        r.margin = likelihood(sg) - likelihood(dg);
        r.delta_fragment = w.fragment - dg.fragment;
        r.delta_exposure = w.exposure - dg.exposure;
        r.delta_likelihood = likelihood(w) - likelihood(dg);
        r.delta_prior = w.prior - dg.prior;
        bool mixed_below = true;
        for (const std::size_t k : mixed) mixed_below = mixed_below && above(sg, groups[k]);
        if (above(sg, dg) && mixed_below) {
            r.label = "DATA_SUPPORTED";
            r.certified = true;
        } else if (above(dg, sg)) {
            r.label = "PRIOR_DEPENDENT";
            r.certified = true;
        } else if (!mixed_below) {
            r.label = "AMBIGUOUS_COMPETITOR";
        } else {
            r.label = "DATA_TIED";
        }
    }
    // phase: classify the existing fragment contenders (the legacy definition), never relabel phase
    for (std::size_t k = 0; k < groups.size(); ++k) {
        if (k == winner) continue;
        const SupportGroup& c = groups[k];
        if (c.fragment + c.bound < w.fragment - w.bound) continue;   // certifiably below on the fragment term: not a contender
        bool same = false, differs = false;
        for (const BlockGenotypeRow& row : c.rows) {
            bool equal_to_some = false;
            for (const BlockGenotypeRow& wrow : w.rows) equal_to_some = equal_to_some || row == wrow;
            (equal_to_some && w.rows.size() == 1 ? same : differs) = true;
        }
        if (same && !differs) {
            ++out.phase.genotype_preserving;
        } else if (differs && !same) {
            ++out.phase.genotype_changing;
        } else {
            ++out.phase.mixed;
        }
    }
    return out;
}

std::string data_support_quality_reason(const std::string& label) {
    if (label == "DATA_SUPPORTED" || label == "NO_EVALUATED_ALTERNATIVE") return "";
    if (label == "DATA_TIED") return "data_tied";
    if (label == "PRIOR_DEPENDENT") return "prior_dependent";
    if (label == "AMBIGUOUS_PROJECTION") return "ambiguous_projection";
    if (label == "AMBIGUOUS_COMPETITOR") return "ambiguous_competitor";
    return "support_not_evaluated";
}

}  // namespace panvar
